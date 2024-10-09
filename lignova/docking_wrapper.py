r"Implementation for the wrapper for the docking program Schrodinger GLIDE."
from typing import Optional

import os
import random

import MDAnalysis as mda
import pandas as pd
import pyarrow as pa
from loguru import logger

from lignova.analysis.rmsd import RMSD
from lignova.docking.combind import Combind
from lignova.docking.contexts import GlideContext
from lignova.docking.contexts.combind import CombindContext
from lignova.docking.glide import Glide
from lignova.docking.utils import convert_to_pdb, manipulate_complexes
from lignova.hdf5.parquet import ParquetParser
from lignova.io import write_text
from lignova.structure.editing import (
    convert_cif2pdb,
    remove_residues,
    select_residues,
    write_mda_universe,
)
from lignova.structure.ligand import DockedLigand, Ligand, PreparedLigand
from lignova.structure.protein import PreparedProtein, Protein
from lignova.structure.utils import (
    chery_pick_ligand,
    separate_protein_ligand,
    validate_ligands,
    validate_pdb,
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
    combine: bool = False,
    datatype: str = "pdb",
) -> pd.DataFrame | mda.Universe:
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
        If true, we remove the water molecules from the ligand file extracted from the pdb file
    combine : bool (default=False)
        If true, we combine PDB and PubChem ligands into one csv file
    datatype : str (default="pdb")
        The type of data to be written to the output file when find_pdb_ligand is True
        if can be "pdb" or "panadas"
    Returns
    -------
    ligand : pd.DataFrame | mda.Universe
    """
    # split the cluster members into pdb and pubchem ligands
    all_pdb_ligands = cluster_members[
        cluster_members["PDB/Gene ID"].apply(
            lambda x: any(char.isalpha() for char in x)
        )
        & (cluster_members["PDB/Gene ID"] == pdb_id)
    ].drop_duplicates()
    logger.info(f"The pdb ligands are {all_pdb_ligands}")
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
        if datatype == "pdb":
            # loop through the pdb ligands and extract the ligand from the pdb file
            for lig_id in all_pdb_ligands["Compound ID"]:
                if os.path.exists(os.path.join(input_dir, f"{pdb_id.lower()}.pdb")):
                    _protein, ligand = chery_pick_ligand(
                        os.path.join(input_dir, f"{pdb_id.lower()}.pdb"),
                        lig_id,
                        remove_water=water,
                    )
                else:
                    logger.error(f"The file {pdb_id.lower()}.pdb does not exist")
                    raise FileNotFoundError(
                        f"The file {pdb_id.lower()}.pdb does not exist"
                    )
        elif datatype == "pandas":
            ligand = all_pdb_ligands[["Smiles", "Compound ID"]].rename(
                columns={"Smiles": "SMILES", "Compound ID": "s_m_title"}
            )
        else:
            logger.error(f"Invalid datatype {datatype}")
            raise ValueError(f"Invalid datatype {datatype}")
    else:
        ligand = (
            pubchem_ligands[["Smiles", "Compound ID"]]
            .drop_duplicates()
            .reset_index(drop=True)
        )
        ligand = ligand.rename(columns={"Smiles": "SMILES", "Compound ID": "s_m_title"})
    if combine:
        all_pdb_ligands_renamed = all_pdb_ligands[["Smiles", "Compound ID"]].rename(
            columns={"Smiles": "SMILES", "Compound ID": "s_m_title"}
        )
        ligand = pd.concat(
            [all_pdb_ligands_renamed, ligand],
            axis=0,
            ignore_index=True,
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
    elif os.path.exists(os.path.join(work_dir, pdb_id.lower() + ".pdb")):
        logger.info(f"{pdb_id} already exists in the directory")
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
        if not os.path.exists(
            os.path.join(context.write_dir, ligand.file_id + "_prepared.mae")
        ):
            glide.PrepLigand(ligand, context)
        prepped_lig = PreparedLigand(
            os.path.join(context.write_dir, ligand.file_id + "_prepared.mae")
        )
    except Exception as exc:
        logger.error(f"Error in preparing ligand {ligand.file_id}")
        raise exc
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
    if not os.path.exists(os.path.join(context.write_dir, temp_prot.file_id + ".pdb")):
        protein, _ligand = separate_protein_ligand(
            pdb_file, remove_water=not (context.samplewater), keep_het_chain="A"
        )
        write_mda_universe(
            protein, os.path.join(context.write_dir, f"{temp_prot.file_id}.pdb")
        )
    input_prot = Protein(
        file_path=os.path.join(context.write_dir, f"{temp_prot.file_id}.pdb")
    )
    glide = Glide()
    try:
        if not os.path.exists(
            os.path.join(context.write_dir, input_prot.file_id + "_grid.zip")
        ):
            glide.PrepProtein(input_prot, context)
        prepped_prot = PreparedProtein(
            os.path.join(context.write_dir, input_prot.file_id + "_grid.zip")
        )
    except Exception as exc:
        logger.error(f"Error in preparing protein {input_prot.file_id}")
        raise exc
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
        if not os.path.exists(
            os.path.join(
                context.write_dir,
                prepped_ligand.file_id.replace("_prepared", "_docking_pv.maegz"),
            )
        ):
            glide.run(prepped_protein, prepped_ligand, context)
        docked_ligand = DockedLigand(
            os.path.join(
                context.write_dir,
                prepped_ligand.file_id.replace("_prepared", "_docking_pv.maegz"),
            )
        )
    except Exception as exc:
        logger.error(f"Error in docking ligand {prepped_ligand.file_id}")
        raise exc
    return docked_ligand


def run_combind(docked_ligand: DockedLigand, context: CombindContext) -> str:
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
    str
    """
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    if not os.path.exists(context.work_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    logger.info(f"Generating features for {docked_ligand.file_id}")
    try:
        if not os.path.exists(
            os.path.join(context.work_dir, docked_ligand.file_id + "_features")
        ):
            combind.featurize(
                docking_filepaths=docked_ligand.file_path,
                file_name=docked_ligand.file_id,
            )

        logger.debug(
            f"Fetures generated for {docked_ligand.file_id}.Selecting top poses"
        )
        if not os.path.exists(
            os.path.join(context.work_dir, docked_ligand.file_id + "_poses.csv")
        ):
            combind.select_pose(
                docked_ligand.file_id + "_poses",
                os.path.join(context.work_dir, docked_ligand.file_id + "_features"),
            )
    except Exception as exc:
        logger.error(f"Error in scoring ligand {docked_ligand.file_id}")
        raise exc
    return os.path.join(context.work_dir, docked_ligand.file_id + "_poses.csv")


def get_top_combind_pose(
    combind_csv: str, glide_docking_file: DockedLigand, context: CombindContext
) -> DockedLigand:
    r"""Parse the combind docking results to get the top pose
    Parameters
    ----------
    combind_csv : str
        The path to the combind docking results CSV file
    glide_docking_file : DockedLigand
        The DockedLigand object containing the glide docking results
    context : CombindContext
        The context object with information about the combind docking
    Returns
    -------
    DockedLigand
    """
    if not os.path.exists(context.work_dir):
        logger.warning(f"The directory {context.work_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    if not os.path.exists(combind_csv):
        logger.error(f"The file {combind_csv} does not exist")
        raise FileNotFoundError(f"The file {combind_csv} does not exist")
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    try:
        if not os.path.exists(
            os.path.join(
                context.work_dir, glide_docking_file.file_id + "_top_poses.maegz"
            )
        ):
            combind.get_3d_top_pose(
                glide_docking_file.file_path,
                combind_csv,
                glide_docking_file.file_id + "_top_poses",
            )
        logger.debug(f"Top poses selected for {glide_docking_file.file_id}")
        top_poses = DockedLigand(
            os.path.join(
                context.work_dir, glide_docking_file.file_id + "_top_poses_pv.maegz"
            )
        )
    except Exception as exc:
        logger.error(f"Error in selecting top poses for {glide_docking_file.file_id}")
        raise exc
    return top_poses


def extract_pdb_top_poses(
    combind_result: DockedLigand, context: CombindContext, pdb_lig: str | None = None
):
    r"""Parse the combind docking results top poses to get each complex in a separate file
    Parameters
    ----------
    combind_result : DockedLigand
        The DockedLigand object containing the combind docking results
    context : CombindContext
        The context object with information about the combind docking
    pdb_lig : str | None (default=None)
        The id of the pdb ligand to be extracted from the complex
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
    glide_context = GlideContext.get_current()
    glide_context.write_dir = context.work_dir
    glide_context.set_current(glide_context)
    if not os.path.exists(context.work_dir):
        logger.warning(f"The directory {context.work_dir} does not exist.Creating it")
        os.makedirs(context.work_dir)
    if not os.path.exists(combind_result.file_path):
        logger.error(f"The file {combind_result.file_path} does not exist")
        raise FileNotFoundError(f"The file {combind_result.file_path} does not exist")
    try:
        if not os.path.exists(
            os.path.join(context.work_dir, combind_result.file_id + "_merge.maegz")
        ):
            manipulate_complexes(
                input_file=combind_result.file_path,
                context=glide_context,
                outfile_name=combind_result.file_id + "_merge.maegz",
                mode="merge",
            )
        if pdb_lig is not None:
            if not os.path.exists(
                os.path.join(context.work_dir, combind_result.file_id + "_merge.csv")
            ):
                combind.extract_data_csv(
                    os.path.join(
                        context.work_dir, combind_result.file_id + "_merge.maegz"
                    ),
                    combind_result.file_id + "_merge",
                    filter_data=False,
                )
            # read the csv file and find the index of the pdb ligand from the s_m_title column
            # after splitting the string by ":" and selecting the second element
            combind_data = pd.read_csv(
                os.path.join(context.work_dir, combind_result.file_id + "_merge.csv")
            )
            pdb_lig_index = (
                combind_data[
                    combind_data["s_m_title"].apply(lambda x: x.split(":")[1])
                    == pdb_lig.upper()
                ].index[0]
                + 1
            )
            logger.debug(f"The index of the pdb ligand is {pdb_lig_index}")
            if not os.path.exists(
                os.path.join(
                    context.work_dir, f"{combind_result.file_id }_{str(pdb_lig)}.pdb"
                )
            ):
                convert_to_pdb(
                    os.path.join(
                        context.work_dir, combind_result.file_id + "_merge.maegz"
                    ),
                    glide_context,
                    [pdb_lig_index],
                )
                os.rename(
                    os.path.join(
                        context.work_dir, combind_result.file_id + "_merge.pdb"
                    ),
                    os.path.join(
                        context.work_dir,
                        f"{combind_result.file_id }_{str(pdb_lig)}.pdb",
                    ),
                )
            if os.path.exists(
                os.path.join(context.work_dir, combind_result.file_id + "_merge.csv")
            ):
                os.remove(
                    os.path.join(
                        context.work_dir, combind_result.file_id + "_merge.csv"
                    )
                )
            if os.path.exists(
                os.path.join(context.work_dir, combind_result.file_id + "_merge.maegz")
            ):
                os.remove(
                    os.path.join(
                        context.work_dir, combind_result.file_id + "_merge.maegz"
                    )
                )
        final_docked_lig = DockedLigand(
            os.path.join(
                context.work_dir,
                f"{combind_result.file_id }_{str(pdb_lig)}.pdb",
            )
        )
        logger.debug(f"Top poses pdb complex for {combind_result.file_id} is extracted")
    except Exception as exc:
        logger.error(f"Error in selecting top poses for {combind_result.file_id}")
        raise exc
    return final_docked_lig


def calc_rmsd_spyrmsd(
    reference_file: DockedLigand,
    target_file: DockedLigand,
    context: GlideContext,
    ligand_name: str = None,
):
    """
    Calculate the symmetry corrected RMSD between the reference and target ligands
    Parameters
    ----------
    reference_file : DockedLigand
        The object containing the reference complex i.e protein and ligand
    target_file : DockedLigand
        The object containing the target complex i.e protein and ligand
    context : GlideContext
        The context object with information about the pre calculation prossessing
    ligand_name : str (default=None)
        The name of the ligand to be extracted from the complex
    Returns
    -------
    rmsd : float
    """
    if not os.path.exists(reference_file.file_path):
        logger.error(f"The file {reference_file.file_path} does not exist")
        raise FileNotFoundError(f"The file {reference_file.file_path} does not exist")
    if not os.path.exists(target_file.file_path):
        logger.error(f"The file {target_file.file_path} does not exist")
        raise FileNotFoundError(f"The file {target_file.file_path} does not exist")
    if not os.path.exists(context.write_dir):
        logger.warning(f"The directory {context.write_dir} does not exist.Creating it")
        os.makedirs(context.write_dir)
    # if file extension is not pdb, convert to pdb
    if not reference_file.file_ext == "pdb":
        logger.debug(f"Converting {reference_file.file_path} to pdb format")
        convert_to_pdb(reference_file.file_path, context)
        reference_file = DockedLigand(
            os.path.join(context.write_dir, reference_file.file_id + ".pdb")
        )
    if not target_file.file_ext == "pdb":
        logger.debug(f"Converting {target_file.file_path} to pdb format")
        convert_to_pdb(target_file.file_path, context)
        target_file = DockedLigand(
            os.path.join(context.write_dir, target_file.file_id + ".pdb")
        )
    # sepate the protein and ligand from the reference and target files
    _ref_prot, ref_lig = separate_protein_ligand(
        reference_file.file_path, remove_water=True
    )
    # get the residue names in ref_lig and exclude the ligand_name
    all_residues = set(list(ref_lig.resnames))
    ligand_residues = [residue for residue in all_residues if residue != ligand_name]
    logger.debug(f"The ligand residues are {ligand_residues}")
    write_mda_universe(
        select_residues(ref_lig, ligand_name),
        os.path.join(context.write_dir, reference_file.file_id + "_ref_lig.pdb"),
    )
    ref_lig_filepath = os.path.join(
        context.write_dir, reference_file.file_id + "_ref_lig.pdb"
    )
    logger.debug(f"Reference ligand is written to {ref_lig_filepath}")
    _tar_prot, tar_lig = separate_protein_ligand(
        target_file.file_path, remove_water=True
    )
    if len(ligand_residues) != 0:
        write_mda_universe(
            remove_residues(tar_lig, ligand_residues),
            os.path.join(context.write_dir, target_file.file_id + "_tar_lig.pdb"),
        )
    else:
        write_mda_universe(
            tar_lig,
            os.path.join(context.write_dir, target_file.file_id + "_tar_lig.pdb"),
        )
    tar_lig_filepath = os.path.join(
        context.write_dir, target_file.file_id + "_tar_lig.pdb"
    )
    logger.debug(f"Target ligand is written to {tar_lig_filepath}")
    new_ref_lig = DockedLigand(
        os.path.join(context.write_dir, reference_file.file_id + "_ref_lig.pdb")
    )
    new_tar_lig = DockedLigand(
        os.path.join(context.write_dir, target_file.file_id + "_tar_lig.pdb")
    )
    rmsd = RMSD(new_tar_lig, new_ref_lig, context)
    res = rmsd.symmetry_rmsd()
    logger.info(f"The RMSD between the reference and target ligands is {res}")
    if os.path.exists(
        os.path.join(context.write_dir, reference_file.file_id + "_ref_lig.pdb")
    ):
        os.remove(
            os.path.join(context.write_dir, reference_file.file_id + "_ref_lig.pdb")
        )
    if os.path.exists(
        os.path.join(context.write_dir, target_file.file_id + "_tar_lig.pdb")
    ):
        os.remove(os.path.join(context.write_dir, target_file.file_id + "_tar_lig.pdb"))
    if os.path.exists(reference_file.file_path):
        os.remove(reference_file.file_path)
    return res


if __name__ == "__main__":
    PARQUET_FILENAME = "all_compounds_with_smiles_cluster.parquet"
    RMSD_FILE_WATER = "rmsd_values_water.csv"
    MRSD_FILE_NO_WATER = "rmsd_values_no_water.csv"
    pdbids = get_pdb_ids_from_parquet(PARQUET_FILENAME)
    logger.info(f"The pdb ids are {pdbids}")
    logger.info(f"The number of pdb ids is {len(pdbids)}")
    failed = []
    rmsd_dict = {}
    for pdbid in pdbids:
        logger.info(f"Getting the pdb coordinates for {pdbid}")
        get_pdb_coordinates(pdbid, "raw")
        ligand_members = extract_parquet_clusters(PARQUET_FILENAME, pdbid.upper())
        ligand_info = parse_ligand_members(
            ligand_members,
            pdbid.upper(),
            input_dir="raw",
            find_pdb_ligand=False,
            combine=True,
        )
        pubchem_lig = ligand_info[
            ligand_info["s_m_title"].apply(lambda x: all(char.isdigit() for char in x))
        ]
        pdb_ligands = parse_ligand_members(
            ligand_members,
            pdbid.upper(),
            input_dir="raw",
            find_pdb_ligand=True,
            datatype="pandas",
        )
        if not (len(pubchem_lig) != 0 and len(pdb_ligands) >= 1):
            logger.error(f"The ligand information for {pdbid} is not complete")
            continue
        logger.info(f"The ligand information is\n {ligand_info}")
        prep_context = GlideContext.get_current()
        prep_context.write_dir = "./trial"
        prep_context.samplewater = False
        prep_context.set_current(prep_context)
        try:
            result = prep_ligands(ligand_info, prep_context, pdbid)
            ref_lig_obj = PreparedLigand(
                file_path=f"trial/{pdbid.lower()}_protein_prepared.mae"
            )
            logger.info(f"The prepared ligand is {result.file_path}")
            prep_prot = prep_proteins(f"raw/{pdbid.lower()}.pdb", prep_context)
            logger.info(f"The prepared protein is {prep_prot.file_path}")
            final_lig = dock_ligands(prep_prot, result, prep_context)
            logger.info(f"The docked ligand is {final_lig.file_path}")
            combind_context = CombindContext.get_current()
            combind_context.work_dir = "./trial"
            combind_context.set_current(combind_context)
            raw_combind_result = run_combind(final_lig, combind_context)
            final_combind_res = get_top_combind_pose(
                raw_combind_result, final_lig, combind_context
            )
            rmsd_target_file = extract_pdb_top_poses(
                final_combind_res, combind_context, pdb_ligands["s_m_title"].values[0]
            )
            RMSD_VAL = calc_rmsd_spyrmsd(
                ref_lig_obj,
                rmsd_target_file,
                prep_context,
                pdb_ligands["s_m_title"].values[0],
            )
            rmsd_dict[pdbid] = RMSD_VAL
            logger.info(f"The RMSD value is {RMSD_VAL}")
        except Exception as error:
            logger.error(f"Error in docking {pdbid}")
            failed.append(pdbid)
        pdb_ligand_name = pdb_ligands["s_m_title"].values[0]
        logger.debug(
            f"The pdb id is {pdbid} which had {pdb_ligand_name} pdb & {len(pubchem_lig)} pubchem ligands"
        )
    # save the failed pdb ids to a text file
    if len(failed) != 0 and not os.path.exists("trial/failed_pdb_ids.txt"):
        with open("trial/failed_pdb_ids.txt", "w", encoding="utf-8") as f:
            for item in failed:
                f.write(f"{item}\n")
    elif os.path.exists("trial/failed_pdb_ids.txt"):
        with open("trial/failed_pdb_ids.txt", "a", encoding="utf-8") as f:
            for item in failed:
                f.write(f"{item}\n")
    # save the rmsd values to a csv file
    rmsd_df = pd.DataFrame(rmsd_dict.items(), columns=["PDB_ID", "RMSD"])
    rmsd_df.to_csv("rmsd_values_no_water.csv", index=False)
