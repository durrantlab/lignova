r"Implementation for the wrapper for the docking program Schrodinger GLIDE."
import os

import pandas as pd
import pyarrow as pa
from loguru import logger

from lignova.hdf5.parquet import ParquetParser

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


def get_pdb_ids_from_parquet(file_path: str) -> list:
    """
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


def parse_parquet_clusters(
    file_path: str, pdb_id: str, same_ligand_cluster: bool = True
) -> pd.DataFrame:
    """
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
    # check if the protein is in the dataframe
    if prot_data.empty:
        logger.error(f"The protein {pdb_id} is not in the dataframe")
        raise ValueError(f"The protein {pdb_id} is not in the dataframe")
    prot_cluster_number = data["Protein Cluster number"].values[0]
    lig_cluster_number = data["Ligand Cluster number"].values[0]
    if same_ligand_cluster:
        cluster_members = prot_data[
            (data["Ligand Cluster number"] == lig_cluster_number)
            & (data["Protein Cluster number"] == prot_cluster_number)
        ]
    else:
        cluster_members = data[data["Protein Cluster number"] == prot_cluster_number]
    return cluster_members


if __name__ == "__main__":
    PARQUET_FILENAME = "all_compounds_with_smiles_cluster.parquet"
    # PROOF OF CONCEPT FOR EACH FUNCTION
    pdbids = get_pdb_ids_from_parquet(PARQUET_FILENAME)
    logger.info(f"Example {parse_parquet_clusters(PARQUET_FILENAME, pdbids[0])}")
    logger.info(f"Length of the pdb ids is {len(pdbids)}")
