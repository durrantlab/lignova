from typing import TextIO, Union

from collections.abc import Iterable

from .editing import filter_hetatoms, get_mda_universe, remove_residues, select_chains


def separate_protein_ligand(pdb: Union[str, TextIO]) -> tuple["Protein", "Ligand"]:
    # TODO: DONE? Function that takes in PDB, cleans (e.g., remove waters), and provides
    # Protein and Ligand objects
    pdb_obj = get_mda_universe(pdb)
    selection = select_chains(pdb_obj)
    # protein = remove_residues(selection, residues=["HOH"])
    hetatm = filter_hetatoms(selection)
    ligand = remove_residues(hetatm, residues=["HOH"])
    print(selection)
    return selection.atoms, ligand
