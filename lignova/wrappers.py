from typing import TextIO, Union

import csv
import os

from loguru import logger

from lignova.structure.base import Structure
from lignova.structure.ligand import Ligand
from lignova.structure.protein import Protein


def get_coordinates(pdb_ids: Union[str, TextIO, list], work_dir: str, limit: int = 50):
    """
    This function takes a list of PDB IDs and downloads the PDB files to a specified directory.
    Parameters
    ----------
    pdb_ids : Union[str,TextIO,list]
        A list of PDB IDs to be downloaded.
    work_dir : str
        The working directory where the PDB files will be downloaded.
    limit : int, optional
        The number of PDB IDs to be downloaded. The default is 50.
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
    # check if the input is a list or a file
    if isinstance(pdb_ids, list):
        logger.info("Input is a list")
        for pdb_id in pdb_ids:
            logger.info("Downloading PDB file for {}".format(pdb_id))
            protein._load_from_pdb_id(
                pdb_id, write=True, write_path=os.path.join(work_dir, pdb_id + ".pdb")
            )
    elif isinstance(pdb_ids, str):
        if os.path.exists(pdb_ids):
            logger.info("Input is a file")
            # read the file
            with open(pdb_ids, "r") as file:
                pdb_ids = file.readlines()
            # find the rows in the 1st columns with "Cluster:"
            pdb_ids = [row.split(",")[0] for row in pdb_ids if "Cluster:" in row]
            # remove the "Cluster:" from the pdb ids
            pdb_ids = [pdb_id.replace("Cluster: [", "") for pdb_id in pdb_ids]
            # remove the "\n" from the pdb ids
            pdb_ids = [pdb_id.replace("\n", "") for pdb_id in pdb_ids]
            # save the pdb ids to a file
            with open(os.path.join(work_dir, "pdb_ids.txt"), "w") as file:
                file.write("\n".join(pdb_ids))
            for pdb_id in pdb_ids:
                if "|" in pdb_id:
                    logger.debug(pdb_id)
                    pdb_id = pdb_id.split("_")[0]
                    logger.info(
                        "Downloading PDB file for {}".format(pdb_id.lstrip("'"))
                    )
                    protein._load_from_pdb_id(
                        pdb_id=pdb_id.lstrip("'"),
                        write=True,
                        write_path=os.path.join(work_dir, pdb_id.lstrip("'") + ".pdb"),
                    )
        else:
            # check if the input is a single pdb id
            logger.info("Input is a single PDB ID")
            pdb_id = pdb_ids
            logger.info("Downloading PDB file for {}".format(pdb_id))
            protein._load_from_pdb_id(
                pdb_id, write=True, write_path=os.path.join(work_dir, pdb_id + ".pdb")
            )


if __name__ == "__main__":
    csv_raw_file = "/home/mma121/PubChem_small/try_schrodinger/clusters.csv"
    csv_filter_file = "/home/mma121/PubChem_small/try_schrodinger/new_clusters_pdb_PostFilter_sep11.csv"

    get_coordinates(
        csv_raw_file, "/home/mma121/PubChem_small/representatives", limit=50
    )
