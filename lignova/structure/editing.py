r""" Implementation for editing protein structures using MDAnalysis."""
from typing import TextIO, Union

from collections.abc import Iterable

import MDAnalysis as mda
from loguru import logger


def get_mda_universe(pdb) -> mda.Universe:
    r"""Prepare MDAnalysis universe"""
    return mda.Universe(topology=pdb, format="PDB")


def select_chains(
    mda_univ: mda.Universe, chains: Union[str, Iterable[str], None] = None
) -> mda.Universe:
    r"""Select specific chains.

    Parameters
    ----------
    u
        MDAnalysis universe to process.
    chains
        Chains to keep.
    """
    n_chains = len(set(mda_univ.segments.segids))
    logger.info("There are {} chains in the structure", n_chains)
    if n_chains == 1:
        logger.info("Selecting the only chain available")
        return mda_univ

    # There are multiple chains in the structure
    if isinstance(chains, str):
        chains = [chains]
    if chains is None:
        chains = ["A"]
    selection = " or ".join([f"segid {c}" for c in chains])
    return mda_univ.select_atoms(selection)


def remove_residues(
    mda_univ: mda.Universe, residues: Union[str, Iterable[str]]
) -> mda.Universe:
    r"""Remove residues from structure.

    Parameters
    ----------
    u
        MDAnalysis universe to process.
    residues
        Names of residues to remove.
    """
    if isinstance(residues, str):
        residues = [residues]
    selection = " and ".join([f"not resname {r}" for r in residues])
    logger.info("MDAnalysis selection: {}", selection)
    return mda_univ.select_atoms(selection)


def select_residues(
    mda_univ: mda.Universe, residues: Union[str, Iterable[str]]
) -> mda.Universe:
    r"""Select residues from structure.

    Parameters
    ----------
    u
        MDAnalysis universe to process.
    residues
        Names of residues to select.
    """
    print(isinstance(residues, list))
    print(residues)
    if isinstance(residues, str):
        residues = list(residues)
    selection = " or ".join([f"record_type HETATM and resname {r}" for r in residues])
    logger.info("MDAnalysis selection: {}", selection)
    return mda_univ.select_atoms(selection)


def remove_hetatoms(mda_univ: mda.Universe) -> mda.Universe:
    r"""Remove hetero atoms.

    Parameters
    ----------
    u
        MDAnalysis universe to process.
    """
    return mda_univ.select_atoms("not record_type HETATM")


def filter_hetatoms(mda_univ: mda.Universe) -> mda.Universe:
    r"""Filter hetero atoms.

    Parameters
    ----------
    u
        MDAnalysis universe to process.
    """
    return mda_univ.select_atoms("record_type HETATM")


def find_common_atoms(mda_univ1: mda.Universe, mda_univ2: mda.Universe) -> Iterable:
    r"""Find common atoms between two MDAnalysis universes.

    Parameters
    ----------
    u1
        MDAnalysis universe to process.
    u2
        MDAnalysis universe to process.

    """
    # Get the atom names for each ligand
    ref_atom_names = [atom.name for atom in mda_univ1.atoms]
    dock_atom_names = [atom.name for atom in mda_univ2.atoms]

    # Find the common atoms
    common_atoms = list(set(ref_atom_names) & set(dock_atom_names))
    return common_atoms


def select_common_atoms(mda_univ: mda.Universe, common_atoms: Iterable) -> mda.Universe:
    r"""Select common atoms.

    Parameters
    ----------
    u
        MDAnalysis universe to process.
    common_atoms
        Names of atoms to select.
    """
    selection = " or ".join([f"name {atom}" for atom in common_atoms])
    return mda_univ.select_atoms(selection)


def write_mda_universe(mda_univ: mda.Universe, file_path: str) -> TextIO:
    r"""Write MDAnalysis universe to file.

    Parameters
    ----------
    u
        MDAnalysis universe to process.
    file_path
        File to write to.
    """
    return mda_univ.write(file_path)
