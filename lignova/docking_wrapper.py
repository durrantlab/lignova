r"Implementation for the wrapper for the docking program Schrodinger GLIDE."
from typing import Optional

import os

import pandas as pd
import pyarrow as pa
from loguru import logger

from lignova.docking.contexts import GlideContext
from lignova.hdf5.parquet import ParquetParser
from lignova.structure.protein import Protein
from lignova.structure.utils import (
    chery_pick_ligand,
    convert_cif2pdb,
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
    cluster_members : pd.DataFrame
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
        cluster_members = data[
            (data["Ligand Cluster number"] == lig_cluster_number)
            & (data["Protein Cluster number"] == prot_cluster_number)
        ]
    else:
        cluster_members = data[data["Protein Cluster number"] == prot_cluster_number]
    return cluster_members


def parse_ligand_members(
    cluster_members: pd.DataFrame,
    pdb_id: str,
    output_path: str,
    find_pdb_ligand: bool = False,
    input_dir: str | None = None,
    water: bool = False,
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
    water : bool (default=False)
        If true, we remove the water molecules from the ligand file
    Returns
    -------
    None
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
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    if find_pdb_ligand:
        if input_dir is None:
            logger.error("The input directory is not provided")
            raise ValueError("The input directory is not provided")
        # loop through the pdb ligands and extract the ligand from the pdb file
        for lig_id in pdb_ligands["Compound ID"]:
            if os.path.exists(os.path.join(input_dir, f"{pdb_id.lower()}.pdb")):
                ligand = chery_pick_ligand(
                    os.path.join(input_dir, f"{pdb_id.lower()}.pdb"),
                    lig_id,
                    remove_water=water,
                )
                write_mda_universe(
                    ligand[1], os.path.join(output_path, f"{pdb_id.lower()}_lig.pdb")
                )
            else:
                logger.error(f"The file {pdb_id.lower()}.pdb does not exist")
                raise FileNotFoundError(f"The file {pdb_id.lower()}.pdb does not exist")
    else:
        # save the Compound ID and Smiles columns to the output file
        pubchem_ligands[["Compound ID", "Smiles"]].to_csv(
            os.path.join(output_path, f"{pdb_id.lower()}_pubchem_lig.csv"),
            index=False,
            header=True,
        )


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
    ligand_file: str, output_dir: str, context: GlideContext | None
) -> None:
    r"""
    Prepare the ligands for docking
    Parameters
    ----------
    ligand_file : str
        The path to the ligand file
    output_dir : str
        The path to the output directory
    Returns
    -------
    None
    """
    pass


if __name__ == "__main__":
    PARQUET_FILENAME = "all_compounds_with_smiles_cluster.parquet"
    # PROOF OF CONCEPT FOR EACH FUNCTION
    pdbids = get_pdb_ids_from_parquet(PARQUET_FILENAME)
    for pdbid in pdbids:
        logger.info(f"Getting the pdb coordinates for {pdbid}")
        get_pdb_coordinates(pdbid, "raw")
    logger.info(f"The pdb ids are {pdbids}")
    logger.info(f"Example {extract_parquet_clusters(PARQUET_FILENAME, '5FTO')}")
    logger.info(f"Length of the pdb ids is {len(pdbids)}")
    cluster_members = extract_parquet_clusters(PARQUET_FILENAME, "5FTO")
    parse_ligand_members(
        cluster_members, "5FTO", "trial", input_dir="raw", find_pdb_ligand=False
    )
