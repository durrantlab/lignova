r" Implemtnation for a wrapper to use lignova for validation"
from typing import TextIO, Union

import glob
import os

from loguru import logger

from lignova.docking import Glide
from lignova.docking.contexts import GlideContext
from lignova.structure.editing import write_mda_universe
from lignova.structure.ligand import Ligand
from lignova.structure.protein import Protein
from lignova.structure.utils import is_xray_structure, separate_protein_ligand


def clean_cluster_files(file_path: str, delim: list = ["[", "]"]):
    """
    This function takes a file containing a list of PDB IDs and returns a list of PDB IDs.
    Parameters
    ----------
    file_path : str
        The path to the file containing the PDB IDs.
    delim : str, optional
        The delimiter used in the file. The default is ":".
    Returns
    -------
    pdb_ids : list
        A list of PDB IDs.
    """
    with open(file_path, "r", encoding="utf-8") as file:
        pdb_ids = []
        for line in file:
            if line.startswith("Cluster: []"):
                continue
            if line.startswith("Cluster:"):
                cluster_id = (
                    line.split(":")[1]
                    .strip()
                    .replace(delim[0], "")
                    .replace(delim[1], "")
                    .replace("\n", "")
                )
                pdb_ids.extend(cluster_id.split(", "))
        pdb_ids = [x for x in pdb_ids if x not in (" ", "''")]
        return pdb_ids


def get_coordinates(
    pdb_ids: Union[str, TextIO, list], work_dir: str, limit: int = 1000
):
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
            logger.info(f"Downloading PDB file for {pdb_id}")
            file_ext = (
                "pdb"
                if protein.get_pdb_from_rcsb(pdb_id).startswith("HEADER")
                else "cif"
            )
            protein.load(
                pdb_id,
                write=True,
                write_path=os.path.join(work_dir, pdb_id + "." + file_ext),
            )
    elif isinstance(pdb_ids, str):
        if os.path.exists(pdb_ids):
            logger.info("Input is a file")
            pdb_ids = clean_cluster_files(
                pdb_ids, delim=["[", "]"]
            )  # NOTE:change this delim when dealing with filtered to ( and )
            logger.info(f"Found a total of {len(pdb_ids)} PDB IDs in the file")
            limit = len(pdb_ids) if len(pdb_ids) < limit else limit
            for i in range(limit):
                if "|" in pdb_ids[i]:
                    pdb_id = pdb_ids[i].split("_")[0]
                    tmp = pdb_id.strip("'")
                else:
                    pdb_id = pdb_ids[i]
                    tmp = pdb_id.lstrip("'").rstrip("'")
                if os.path.exists(os.path.join(work_dir, tmp.lower() + ".pdb")):
                    logger.info(f"PDB file for {tmp} already exists")
                    continue

                logger.info(f"Downloading PDB file for {tmp}")
                file_ext = (
                    "pdb"
                    if protein.get_pdb_from_rcsb(tmp).startswith("HEADER")
                    else "cif"
                )
                protein.load(
                    pdb_id=tmp,
                    write=True,
                    write_path=os.path.join(work_dir, tmp.lower() + "." + file_ext),
                )
        else:
            # check if the input is a single pdb id
            logger.info("Input is a single PDB ID")
            pdb_id = pdb_ids
            logger.info(f"Downloading PDB file for {pdb_id}")
            file_ext = (
                "pdb"
                if protein.get_pdb_from_rcsb(pdb_id).startswith("HEADER")
                else "cif"
            )
            protein.load(
                pdb_id=pdb_id,
                write=True,
                write_path=os.path.join(work_dir, pdb_id.lower() + "." + file_ext),
            )


def prep_structure(
    input_dir: str, output_dir: str, pdb_id: Union[str, list, None], limit: int = 500
):
    """
    This function takes a PDB ID and a ligand ID and returns a Structure object.
    Parameters
    ----------
    directory : str
        The directory where the PDB file is located.
    pdb_id : Optional[str]
        The PDB ID of the protein. The default is None.
    limit : int
        The number of PDB IDs to prep. The default is 50.
    Returns
    -------
    """
    glide = Glide()
    context = GlideContext.get_current()
    if not os.path.exists(input_dir):
        raise FileNotFoundError("Input directory not found")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if pdb_id is not None:
        if isinstance(pdb_id, str):
            pdb_id = [pdb_id]
        for ids in pdb_id:
            if not os.path.exists(
                os.path.join(input_dir, ids.lower() + ".pdb")
            ) and not os.path.exists(os.path.join(input_dir, ids.lower() + ".cif")):
                logger.warning(f"{ids} not found in the directory. Downloading it now")

                get_coordinates(ids, input_dir, limit=1)
    # NOTE:REMEMBER TO DEAL WITH THE CIF FILES
    pdb_files = glob.glob(os.path.join(input_dir, "*.pdb"))
    logger.info(f"Found a total of {len(pdb_files)} PDB files in the directory")
    limit = len(pdb_files) if len(pdb_files) < limit else limit
    pdb_files = pdb_files[:limit] if len(pdb_files) > limit else pdb_files
    for pdb_file in pdb_files:
        if "_lig" not in pdb_file:
            prot = Protein(pdb_file)
            prot.load(file_path=pdb_file)
            if os.path.exists(
                os.path.join(input_dir, pdb_file.replace(".pdb", "_lig.pdb"))
            ):
                continue
            if is_xray_structure(prot.file_path):
                logger.info(f"Prepping {pdb_file}")
                protein, ligand = separate_protein_ligand(prot.file_path)
                write_mda_universe(protein, os.path.join(input_dir, pdb_file))
                # NOTE:I can have this named with the lig name should i do it?
                write_mda_universe(
                    ligand,
                    os.path.join(input_dir, pdb_file.replace(".pdb", "_lig.pdb")),
                )
            raw_prot = Protein(pdb_file)
            raw_lig = Ligand(raw_prot.file_path.replace(".pdb", "_lig.pdb"))
            context.write_dir = output_dir
            prepared_lig = glide.PrepLigand(raw_lig, context)
            glide.convert_to_mae(raw_prot, context)
            prot_mae = Protein(
                file_path=os.path.join(
                    context.write_dir, raw_prot.file_name.replace(".pdb", ".mae")
                )
            )
            prepared_prot = glide.PrepProtein(prot_mae, context)


if __name__ == "__main__":
    RAW_FILE = "/home/mma121/PubChem_small/try_schrodinger/clusters.csv"
    FILTERED_FILE = (
        "/home/mma121/PubChem_small/try_schrodinger/new_clusters_pdb_postfilter.csv"
    )

    get_coordinates(RAW_FILE, "/home/mma121/PubChem_small/representatives", limit=20)
    prep_structure(
        "/home/mma121/PubChem_small/representatives",
        "./prepped",
        pdb_id=None,
        limit=20,
    )
