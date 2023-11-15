r""" Utility functions for structure module. """
from typing import TextIO, Union

from .editing import filter_hetatoms, get_mda_universe, remove_residues, select_chains


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
    return selection.atoms, ligand
