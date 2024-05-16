r" Implemtnation for a wrapper for MMseqs2 clustering"
from typing import TextIO, Union

import os
import time

import pandas as pd
from loguru import logger

from lignova.structure.utils import get_rcsb_data, validate_ligands, validate_pdb


def fasta_parser(fasta: Union[str, TextIO], delimiter: Union[str, None] = None) -> list:
    r"""Parse FASTA files to get the protein ids
    Parameters
    ----------
    fasta : Union[str, TextIO]
        Path to the FASTA file.
    delimiter : Union[str, None]
        Delimiter to split the protein id from the FASTA header. Default is None.
    Returns
    -------
    list
        List of protein ids.
    """
    # check if the fasta file exists and has .fasta extension
    if not os.path.exists(fasta):
        raise FileNotFoundError(f"FASTA file {fasta} not found.")
    if not fasta.endswith(".fasta"):
        raise ValueError(f"FASTA file {fasta} must have .fasta extension.")
    # read the fasta file
    with open(fasta, "r", encoding="utf-8") as file:
        lines = file.readlines()
    # get the protein ids by filtering the lines starting with ">" a
    # and splitting them by the delimiter
    if delimiter is not None:
        protein_ids = [
            line.split(delimiter)[0].strip(">")
            for line in lines
            if line.startswith(">")
        ]
    else:
        protein_ids = [line.strip(">") for line in lines if line.startswith(">")]
    return protein_ids


def pdb_validations(csvfilenames: str) -> dict:
    r"""Validate the PDB files for the proteins in the CSV file.
    Parameters
    ----------
    csvfilenames : str
        Path to the CSV file containing the protein ids.
    """
    start_time = time.time()
    # check if the CSV file exists
    if not os.path.exists(csvfilenames):
        raise FileNotFoundError(f"CSV file {csvfilenames} not found.")
    # read the CSV file using pandas
    protein_ids = pd.read_csv(csvfilenames)
    # get the PDB column
    pdb_ids = protein_ids["PDB"]
    # check if new_pdb_files exists and if so read it as a dictionary
    if os.path.exists("new_pdb_files"):
        new_file = pd.read_csv("new_pdb_files")
        new_file = new_file.to_dict()
    else:
        new_file = {}
    # iterate over the PDB ids and validate them while keeping track of the values of the From column of the CSV file
    for gene_id, pdb_id in zip(protein_ids["From"], pdb_ids):
        # split each pdb id by ; and skip the last one
        pdb_id = pdb_id.split(";")[:-1]
        logger.debug(f"Gene id: {gene_id} PDB id: {pdb_id}")
        # find the value of the From column
        if gene_id not in new_file:
            new_file[gene_id] = []
        else:
            continue
        for pdb in pdb_id:
            # validate the pdb id
            if len(get_rcsb_data(pdb)) == 0:
                logger.error(f"Invalid PDB id {pdb} for protein {gene_id}")
                continue
            if not validate_pdb(pdb) or not validate_ligands(pdb):
                logger.error(f"Invalid PDB id {pdb} for protein {gene_id}")
                continue
            logger.info(f"Valid PDB id {pdb} for gene_id {gene_id}")
            new_file[gene_id].append(pdb)
            # save this new file to a csv file every 1hr
            if time.time() - start_time > 60 * 60:
                df = pd.DataFrame(new_file.items(), columns=["Gene_id", "PDB"])
                df.to_csv("new_pdb_files", index=False)
                start_time = time.time()

    return new_file


if __name__ == "__main__":
    BIND_MOAD_FASTA = "/home/mma121/PubChem_small/try_schrodinger/MOAD_sequences.fasta"
    PUBCHEM_FASTA = "/home/mma121/PubChem_small/try_schrodinger/PUBCHEM_HDF5.fasta"
    """
    bind_moad_protein_ids = fasta_parser(BIND_MOAD_FASTA)
    PubChem_protein_ids = fasta_parser(PUBCHEM_FASTA)
    logger.info("Number of proteins in bind MOAD: {}", len(bind_moad_protein_ids))
    logger.info("Number of proteins in PubChem: {}", len(PubChem_protein_ids))
    # save the pubchem protein ids to a file
    with open("PubChem_protein_ids.txt", "w", encoding="utf-8") as f:
        for protein_id in PubChem_protein_ids:
            f.write(protein_id + "\n")
    """
    # validate the PDB files for the proteins in the CSV file
    new_file = pdb_validations(
        "/home/mma121/PubChem_small/try_schrodinger/gene-id_with_pdb.csv"
    )
    # save the new file which is a dictionary to a csv file
    df = pd.DataFrame(new_file.items(), columns=["Gene_id", "PDB"])
    df.to_csv("new_pdb_files", index=False)
