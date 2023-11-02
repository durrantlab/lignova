from typing import TextIO, Union

from collections.abc import Iterable

import MDAnalysis as mda
from loguru import logger


def get_mda_universe(pdb) -> mda.Universe:
    r"""Prepare MDAnalysis universe"""
    return mda.Universe(topology=pdb, format="PDB")


def select_chains(
    u: mda.Universe, chains: Union[str, Iterable[str], None] = None
) -> mda.Universe:
    r"""Select specific chains.

    Parameters
    ----------
    u
        MDAnalysis universe to process.
    chains
        Chains to keep.
    """
    n_chains = len(set(u.segments.segids))
    logger.info("There are {} chains in the structure", n_chains)
    if n_chains == 1:
        logger.info("Selecting the only chain available")
        return u

    # There are multiple chains in the structure
    if isinstance(chains, str):
        chains = [chains]
    if chains is None:
        chains = ["A"]
    selection = " and ".join([f"segid {c}" for c in chains])
    return u.select_atoms(selection)


def remove_residues(
    u: mda.Universe, residues: Union[str, Iterable[str]]
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
    return u.select_atoms(selection)


def remove_hetatoms(u: mda.Universe) -> mda.Universe:
    r"""Remove hetero atoms.

    Parameters
    ----------
    u
        MDAnalysis universe to process.
    """
    return u.select_atoms("not record_type HETATM")


def separate_protein_ligand(pdb: Union[str, TextIO]) -> tuple["Protein", "Ligand"]:
    # TODO: Function that takes in PDB, cleans (e.g., remove waters), and provides
    # Protein and Ligand objects.
    pass
