r" Implemtnation for a wrapper for MMseqs2 clustering"
from typing import TextIO, Union

import os

from loguru import logger


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


if __name__ == "__main__":
    BIND_MOAD_FASTA = "/home/mma121/PubChem_small/try_schrodinger/MOAD_sequences.fasta"
    PUBCHEM_FASTA = "/home/mma121/PubChem_small/try_schrodinger/PUBCHEM_HDF5.fasta"
    bind_moad_protein_ids = fasta_parser(BIND_MOAD_FASTA)
    PubChem_protein_ids = fasta_parser(PUBCHEM_FASTA)
    logger.info("Number of proteins in bind MOAD: {}", len(bind_moad_protein_ids))
    logger.info("Number of proteins in PubChem: {}", len(PubChem_protein_ids))
    # save the pubchem protein ids to a file
    with open("PubChem_protein_ids.txt", "w", encoding="utf-8") as f:
        for protein_id in PubChem_protein_ids:
            f.write(protein_id + "\n")
