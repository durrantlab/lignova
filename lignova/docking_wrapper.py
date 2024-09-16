r"Implementation for the wrapper for the docking program Schrodinger GLIDE."
from typing import Optional

import os

import MDAnalysis as mda
import pandas as pd
import pyarrow as pa
from loguru import logger

from lignova.docking.combind import Combind
from lignova.docking.contexts import GlideContext
from lignova.docking.contexts.combind import CombindContext
from lignova.docking.glide import Glide
from lignova.hdf5.parquet import ParquetParser
from lignova.io import write_text
from lignova.structure.ligand import DockedLigand, Ligand, PreparedLigand
from lignova.structure.protein import PreparedProtein, Protein
from lignova.structure.utils import (
    chery_pick_ligand,
    convert_cif2pdb,
    separate_protein_ligand,
    validate_ligands,
    validate_pdb,
    write_mda_universe,
)

# OUTLINE
# 1. Importing necessary modules
# 2. Reading the input files to get the ligand and receptor files
# 3. Running the preparation of the ligand and receptor files
# 4. Running the docking program
# 5. Running the post-processing of the docking results to get the top poses using combind score
# 6. Calculating the RMSD of the top poses with the reference ligand pose
# 7. Writing the output files with which proteins passed the docking
# and the RMSD values <= 2.5 Angstroms


# 2. Reading input files to get the ligand and receptor files


def get_pdb_ids_from_parquet(
    file_path: str, schema: Optional[pa.schema] = None
) -> list:
    r"""
    Get the pdb ids from the parquet file
    Parameters
    ----------
    file_path : str
        The path to the parquet file
    Returns
    -------
    pdb_ids : list
        The list of pdb ids
    """
    if not os.path.exists(file_path):
        logger.error(f"The file {file_path} does not exist")
        raise FileNotFoundError(f"The file {file_path} does not exist")
    if schema is None:
        schema = pa.schema(
            [
                ("Protein Cluster number", pa.int64()),
                ("PDB/Gene ID", pa.string()),
                ("Compound ID", pa.string()),
                ("Smiles", pa.string()),
                ("Ligand Cluster number", pa.int64()),
            ]
        )
    new_parquet = ParquetParser(file_path, schema)
    raw_prot_ids = new_parquet.convert_to_pandas()["PDB/Gene ID"].unique()
    pdb_ids = [
        raw_prot_ids[i]
        for i in range(len(raw_prot_ids))
        if any(char.isalpha() for char in raw_prot_ids[i])
    ]
    return pdb_ids


def extract_parquet_clusters(
    file_path: str, pdb_id: str, same_ligand_cluster: bool = True
) -> pd.DataFrame:
    r"""
    Read the parquet file containing the clustered protein and ligand information
    Parameters
    ----------
    file_path : str
        The path to the parquet file
    pdb_id : str
        The pdb id of the protein of interest
    same_ligand_cluster : bool (default=True)
        If true, only return the members in the same ligand cluster as the protein of interest
        not only the same protein cluster
    Returns
    -------
    members : pd.DataFrame
        The dataframe containing the protein and ligand information
    """
    if not os.path.exists(file_path):
        logger.error(f"The file {file_path} does not exist")
        raise FileNotFoundError(f"The file {file_path} does not exist")
    data = pd.read_parquet(file_path)
    prot_data = data[data["PDB/Gene ID"] == pdb_id]
    logger.info(f"The protein data is {prot_data}")
    # check if the protein is in the dataframe
    if prot_data.empty:
        logger.error(f"The protein {pdb_id} is not in the dataframe")
        raise ValueError(f"The protein {pdb_id} is not in the dataframe")
    prot_cluster_number = prot_data["Protein Cluster number"].values[0]
    lig_cluster_number = prot_data["Ligand Cluster number"].values[0]
    if same_ligand_cluster:
        members = data[
            (data["Ligand Cluster number"] == lig_cluster_number)
            & (data["Protein Cluster number"] == prot_cluster_number)
        ]
    else:
        members = data[data["Protein Cluster number"] == prot_cluster_number]
    return members


def parse_ligand_members(
    cluster_members: pd.DataFrame,
    pdb_id: str,
    find_pdb_ligand: bool = False,
    input_dir: str | None = None,
    water: bool = True,
) -> None:
    r"""Parse the cluster members information to write the ligand file
    Parameters
    ----------
    cluster_members : pd.DataFrame
        The dataframe containing the information about the ligand members in the cluster
        extracted from the parquet file
    pdb_id : str
        The pdb id of the protein of interest
    output_path : str
        The path to the output ligand file to be written
    find_pdb_ligand : bool (default=False)
        If true, we extract the crystallographic ligand from the pdb file
        if false, we extract the pubchem ligand from the smiles string
        and write it to the output file
    input_dir : str | None (default=None)
        The path to the directory containing the pdb files
    water : bool (default=True)
        If true, we remove the water molecules from the ligand file
    Returns
    -------
    ligand : pd.DataFrame | mda.Universe
    """
    # split the cluster members into pdb and pubchem ligands
    pdb_ligands = cluster_members[
        cluster_members["PDB/Gene ID"].apply(
            lambda x: any(char.isalpha() for char in x)
        )
        & (cluster_members["PDB/Gene ID"] == pdb_id)
    ]
    logger.info(f"The pdb ligands are {pdb_ligands}")
    pubchem_ligands = cluster_members[
        cluster_members["PDB/Gene ID"].apply(
            lambda x: all(char.isdigit() for char in x)
        )
    ].drop_duplicates()
    logger.info(f"The pubchem ligands are {pubchem_ligands}")
    if find_pdb_ligand:
        if input_dir is None:
            logger.error("The input directory is not provided")
            raise ValueError("The input directory is not provided")
        # loop through the pdb ligands and extract the ligand from the pdb file
        for lig_id in pdb_ligands["Compound ID"]:
            if os.path.exists(os.path.join(input_dir, f"{pdb_id.lower()}.pdb")):
                protein, ligand = chery_pick_ligand(
                    os.path.join(input_dir, f"{pdb_id.lower()}.pdb"),
                    lig_id,
                    remove_water=water,
                )
            else:
                logger.error(f"The file {pdb_id.lower()}.pdb does not exist")
                raise FileNotFoundError(f"The file {pdb_id.lower()}.pdb does not exist")
    else:
        ligand = (
            pubchem_ligands[["Compound ID", "Smiles"]]
            .drop_duplicates()
            .reset_index(drop=True)
        )
    return ligand


def get_pdb_coordinates(pdb_id: str, work_dir: str):
    """
    This function takes a list of PDB IDs and downloads the PDB files to a specified directory
    if they pass the validation test. (no mutation, has ligand, no covalent bond, x-ray structures,)
    Parameters
    ----------
    pdb_id : str|
        The PDB ids to be downloaded
    work_dir : str
        The working directory where the PDB file will be downloaded.
    Returns
    -------
    None.
    """
    current_dir = os.getcwd()
    protein = Protein()
    # check if the output directory exists
    if not os.path.exists(work_dir):
        logger.info("Output_Dir not found,Creating it in working directory")
        os.mkdir(os.path.join(current_dir, work_dir))
    if (
        not os.path.exists(os.path.join(work_dir, pdb_id.lower() + ".pdb"))
        and validate_pdb(pdb_id)
        and validate_ligands(pdb_id)
    ):
        logger.info(f"Downloading PDB file for {pdb_id}")
        file_ext = (
            "pdb" if protein.get_pdb_from_rcsb(pdb_id).startswith("HEADER") else "cif"
        )
        protein.load(
            pdb_id=pdb_id,
            write=True,
            write_path=os.path.join(work_dir, pdb_id.lower() + "." + file_ext),
        )
        if file_ext == "cif":
            logger.info(f"Converting {pdb_id} to pdb format")
            convert_cif2pdb(
                os.path.join(work_dir, pdb_id.lower() + ".cif"),
                os.path.join(work_dir, pdb_id.lower() + ".pdb"),
            )
    else:
        logger.warning(f"{pdb_id} failed validation test")


def prep_ligands(
    ligand_object: pd.DataFrame | mda.core.groups.AtomGroup,
    context: GlideContext | None,
    file_name: str,
) -> PreparedLigand:
    r"""
    Prepare the ligands for docking
    Parameters
    ----------
    ligand_file : pd.DataFrame | mda.Universe
        The ligand data to be prepared
    context : GlideContext | None
        The context object with information about the preparation
    file_name : str
        The name of the file to be written
    Returns
    -------
    PreparedLigand
    """
    glide = Glide()
    if not os.path.exists(context.write_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    if isinstance(ligand_object, pd.DataFrame):
        file_path = write_text(ligand_object, file_ext=".csv")
        file_name = file_name + "_lig_pubchem.csv"
    else:
        file_path = write_text(ligand_object, file_ext=".pdb")
        file_name = file_name + "_lig.pdb"
    temp_path = os.path.dirname(file_path)
    logger.debug(f"Temp location is {temp_path}")
    # Rename the temporary file
    new_path = os.path.join(temp_path, f"{file_name}")
    os.rename(file_path, new_path)
    logger.debug(f"Renamed temporary file to {new_path}")
    ligand = Ligand(new_path)
    try:
        glide.PrepLigand(ligand, context)
        prepped_lig = PreparedLigand(
            os.path.join(context.write_dir, ligand.file_id + "_prepared.mae")
        )
    except Exception as e:
        logger.error(f"Error in preparing ligand {ligand.file_id}")
        raise e
    os.remove(new_path)
    if os.path.exists(os.path.join(context.write_dir, ligand.file_id + ".mae")):
        os.remove(os.path.join(context.write_dir, ligand.file_id + ".mae"))
    return prepped_lig


def prep_proteins(pdb_file: str, context: GlideContext):
    """
    Prepare the protein for docking
    Parameters
    ----------
    pdb_file : str
        The path to the pdb file
    context : GlideContext
        The context object with information about the preparation
    Returns
    -------
    PreparedProtein
    """
    temp_prot = Protein(file_path=pdb_file)
    if not os.path.exists(context.write_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    protein, ligand = separate_protein_ligand(pdb_file, remove_water=False)
    write_mda_universe(
        protein, os.path.join(context.write_dir, f"{temp_prot.file_id}.pdb")
    )
    input_prot = Protein(
        file_path=os.path.join(context.write_dir, f"{temp_prot.file_id}.pdb")
    )
    glide = Glide()
    try:
        glide.PrepProtein(input_prot, context)
        prepped_prot = PreparedProtein(
            os.path.join(context.write_dir, input_prot.file_id + "_grid.zip")
        )
    except Exception as e:
        logger.error(f"Error in preparing protein {input_prot.file_id}")
        raise e
    os.remove(os.path.join(context.write_dir, f"{temp_prot.file_id}.pdb"))
    return prepped_prot


def dock_ligands(
    prepped_protein: PreparedProtein,
    prepped_ligand: PreparedLigand,
    context: GlideContext,
) -> None:
    r"""
    Dock the ligands to the protein
    Parameters
    ----------
    prepped_protein : PreparedProtein
        The prepared protein to dock the ligand to
    prepped_ligand : PreparedLigand
        The prepared ligand to be docked
    context : GlideContext
        The context object with information about the docking
    Returns
    -------
    DockedLigand
    """
    glide = Glide()
    if not os.path.exists(context.write_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    logger.info(f"Docking {prepped_ligand.file_id} to {prepped_protein.file_id}")
    try:
        glide.run(prepped_protein, prepped_ligand, context)
        docked_ligand = DockedLigand(
            os.path.join(
                context.write_dir, prepped_ligand.file_id + "_docking_pv.maegz"
            )
        )
    except Exception as e:
        logger.error(f"Error in docking ligand {prepped_ligand.file_id}")
        raise e
    return docked_ligand


# NOTE: The function run_combind is NOT READY FOR USE NOR WRITTEN correctly
def run_combind(docked_ligand: DockedLigand, context: CombindContext) -> DockedLigand:
    """
    Run the combind program to get the top poses
    Parameters
    ----------
    docked_ligand : DockedLigand
        The docked ligand to be scored
    context : GlideContext
        The context object with information about the scoring
    Returns
    -------
    DockedLigand
    """
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    if not os.path.exists(context.write_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    logger.info(f"Generating feautes for {docked_ligand.file_id}")
    try:
        combind.featurize(
            docking_filepaths=docked_ligand, file_name=docked_ligand.file_id
        )
        combind.select_pose(
            docked_ligand.file_id,
            os.path.join(context.write_dir, docked_ligand.file_id + "_features"),
        )
        combind.get_3d_top_pose(
            docked_ligand.file_path,
            os.path.join(context.write_dir, docked_ligand.file_id + ".csv"),
            dock_ligands.file_id,
        )
    except Exception as e:
        logger.error(f"Error in scoring ligand {docked_ligand.file_id}")
        raise e
    return top_poses


if __name__ == "__main__":
    PARQUET_FILENAME = "all_compounds_with_smiles_cluster.parquet"
    # PROOF OF CONCEPT FOR EACH FUNCTION
    pdbids = get_pdb_ids_from_parquet(PARQUET_FILENAME)
    """
    for pdbid in pdbids:
        logger.info(f"Getting the pdb coordinates for {pdbid}")
        get_pdb_coordinates(pdbid, "raw")
    logger.info(f"Example {extract_parquet_clusters(PARQUET_FILENAME, '5FTO')}")
    logger.info(f"Length of the pdb ids is {len(pdbids)}")
    """
    cluster_members = extract_parquet_clusters(PARQUET_FILENAME, "5FTO")
    ligand_info = parse_ligand_members(
        cluster_members, "5FTO", input_dir="raw", find_pdb_ligand=False
    )
    logger.info(f"The ligand information is\n {ligand_info}")
    prep_context = GlideContext.get_current()
    prep_context.write_dir = "./trial"
    prep_context.samplewater = True
    prep_context.set_current(prep_context)
    result = prep_ligands(ligand_info, prep_context, "5fto")
    logger.info(f"The prepared ligand is {result.file_path}")
    prep_prot = prep_proteins("raw/5fto.pdb", prep_context)
    logger.info(f"The prepared protein is {prep_prot.file_path}")
    final_lig = dock_ligands(prep_prot, result, prep_context)
    logger.info(f"The docked ligand is {final_lig.file_path}")
