r""" Utility functions for structure module. """
from typing import TextIO, Union

import os

from loguru import logger

from .editing import filter_hetatoms, get_mda_universe, remove_residues, select_chains


def is_xray_structure(pdb_file):
    """
    Check if the PDB file was generated from X-ray diffraction data.

    Parameters:
    -----------
    pdb_file : str
        Path to the PDB file.

    Returns:
    --------
    bool
        True if the PDB file was generated from X-ray data, False otherwise.
    """
    with open(pdb_file, "r") as file:
        lines = file.readlines()
    ext = os.path.splitext(pdb_file)[-1].lower()
    if ext == ".pdb":
        expdta_line = [line for line in lines if line.startswith("EXPDTA")]
        if expdta_line:
            return "X-RAY" in expdta_line[0]
        else:
            remark_200_line = [line for line in lines if line.startswith("REMARK 200")]
            return bool(remark_200_line)
    elif ext == ".cif":
        # find the _exptl.method line
        exptl_line = [line for line in lines if line.startswith("_exptl.method")]
        if exptl_line:
            return "X-RAY" in exptl_line[0]


def separate_protein_ligand(pdb: Union[str, TextIO]) -> tuple["Protein", "Ligand"]:
    r"""Separate protein and ligand from a PDB file.
    Parameters
    ----------
    pdb : str or file-like
        Path to the PDB file or file-like object.
    Returns
    -------
    Protein
        Universe object containing the protein.
    Ligand
        Universe object containing the ligand.
    """
    pdb_obj = get_mda_universe(pdb)
    selection = select_chains(pdb_obj)
    hetatm = filter_hetatoms(selection)
    ligand = remove_residues(hetatm, residues=["HOH"])
    if len(ligand.atoms) == 0 or len(hetatm.atoms) == 0:
        logger.warning("No ligand found in the selected chains.Trying another method.")
        hetatm = filter_hetatoms(pdb_obj)
        hetatm = remove_residues(hetatm, residues=["HOH"])
        for i in list(set(hetatm.resnames)):
            if len(i) != 3:
                hetatm = remove_residues(hetatm, residues=[i])
        selection_2 = select_chains(pdb_obj, chains=hetatm.segments.segids)
        ligand = filter_hetatoms(selection_2)
        return selection_2.atoms, ligand

    return selection.atoms, ligand
