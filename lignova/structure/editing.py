r""" Implementation for editing protein structures using MDAnalysis."""

from typing import Literal

from collections.abc import Iterable

import gemmi
import MDAnalysis as mda
from loguru import logger
from MDAnalysis.core.groups import AtomGroup


def get_mda_universe(pdb: str) -> mda.Universe:
    r"""Prepare MDAnalysis universe
    Args:
        pdb : Path to PDB file.
    Returns:
        MDAnalysis universe for this pdb.
    """
    return mda.Universe(topology=pdb, format="PDB")


def select_water(
    pdb: mda.Universe | str,
    water_selection: Literal["surface", "interfacial"],
    ligand: str,
    water_distance: float = 3.6,
) -> mda.Universe:
    r"""Select specific water molecules from the structure.

    Args:
        pdb : MDAnalysis universe to process or path to PDB file.
        water_selection : Selection of water molecules to keep.
            surface: Water molecules on the second hydration shell.
                i.e 3.6 A from the ligand but not in direct contact with the protein.
            interfacial: Water molecules on the first hydration shell.
                i.e 3.6 A from the ligand and in direct contact with the protein.
        ligand : Ligand to calculate the distance from.
        water_distance : Distance to select water molecules. Default is 3.6 A.
    Returns:
        MDAnalysis universe with selected water molecules.
    """
    if isinstance(pdb, str):
        pdb = get_mda_universe(pdb)
    if len(pdb.select_atoms("resname HOH")) == 0:
        raise ValueError("No water molecules in the structure")
    electronegative = (
        "type O or type N or type F or type Cl or type Br or type I or type S"
    )
    water_electronegative = pdb.select_atoms("resname HOH and type O")
    logger.info(
        f"Total Number of electronegative atoms in the water: {len(water_electronegative)}"
    )
    surface_water = None
    interfacial_water = None
    if water_selection == "surface":
        surface_water = pdb.select_atoms(
            f"resname HOH and around {water_distance} ((resname {ligand} or protein) and ({electronegative}))"
        )
    elif water_selection == "interfacial":
        interfacial_water = pdb.select_atoms(
            f"resname HOH and around {water_distance} (resname {ligand} and {electronegative}) and (around {water_distance} protein and {electronegative})"
        )
    else:
        raise ValueError("Invalid water selection")
    return surface_water if water_selection == "surface" else interfacial_water


def select_chains(
    mda_univ: mda.Universe, chains: str | Iterable[str] | None = None
) -> mda.Universe:
    r"""Select specific chains.

    Args:
        u : MDAnalysis universe to process.
        chains : Chains to keep if None then select chain A. Default is None.
    """
    n_chains = len(set(mda_univ.segments.segids))
    logger.info(f"There are {n_chains} chains in the structure")
    if n_chains == 1:
        logger.warning("There is only one Chain available. Selecting it")
        return mda_univ

    # There are multiple chains in the structure
    if isinstance(chains, str):
        chains = [chains]
    if chains is None:
        chains = ["A"]
    selection = " or ".join([f"segid {c}" for c in chains])
    return mda_univ.select_atoms(selection)


def validate_chains(mda_univ: mda.Universe, chains: str | Iterable[str]) -> bool:
    r"""Validate chains.

    Args:
        u : MDAnalysis universe to process.
        chains : Chains to validate.
    Returns:
        True if chains are present in the structure.
    """
    if isinstance(chains, str):
        chains = [chains]
    n_chains = len(set(mda_univ.segments.segids))
    logger.info(f"There are {n_chains} chains in the structure")
    if set(chains).issubset(set(mda_univ.segments.segids)):
        return True
    return False


def remove_residues(mda_univ: mda.Universe, residues: str | Iterable[str]) -> AtomGroup:
    r"""Remove residues from structure.

    Args:
        u : MDAnalysis universe to process.
        residues : types of residues to remove.
    Returns:
        MDAnalysis universe for the structure without the specified residues.
    """
    if isinstance(residues, str):
        residues = [residues]
    selection = " and ".join([f"not resname {r}" for r in residues])
    return mda_univ.select_atoms(selection)


def merge_universes(mda_univs: Iterable[mda.Universe]) -> mda.Universe:
    r"""Merge multiple MDAnalysis universes.

    Args:
        mda_univs : list of MDAnalysis universes to merge.
    Returns:
        Merged MDAnalysis universe.

    """
    if isinstance(mda_univs, list) and len(mda_univs) > 1:
        logger.info(f"Merging {len(mda_univs)} MDAnalysis universes")
        merge_list: list[AtomGroup] = []
        for i in mda_univs:
            # check if the universe is empty
            if i.atoms is None or len(i.atoms) == 0:
                logger.warning("Empty universe found")
                continue
            # get the atoms from each universe and merge them
            merge_list.append(i.atoms)
        merged = mda.Merge(*merge_list)
    else:
        logger.warning("Only one MDAnalysis universe to merge")
        merged = mda_univs[0]
    if merged.atoms is None:
        raise ValueError("Merged universe has no atoms")
    return mda.Universe.empty(0) if len(merged.atoms) == 0 else merged


def select_residues(
    mda_univ: mda.Universe, residues: str | Iterable[str] | int | Iterable[int]
) -> mda.Universe:
    r"""Select residues from structure.

    Args:
        mda_univ : MDAnalysis universe to process.
        residues : types of residues or the residues id to select.
    Returns:
        MDAnalysis universe for the selected residues.
    """
    if isinstance(residues, str):
        residues = [residues]
    selection = ""
    # check if residues is a list of strings
    if all(residue.isdigit() and residue for residue in residues):
        residues = [int(residue) for residue in residues]
        selection = " or ".join([f"resid {r}" for r in residues])
    elif all(isinstance(residue, str) for residue in residues):
        selection = " or ".join([f"resname {r}" for r in residues])

    logger.debug(f"MDAnalysis selection: {selection}")
    return mda_univ.select_atoms(selection)


def remove_hetatoms(mda_univ: mda.Universe) -> mda.Universe:
    r"""Remove hetero atoms.

    Args:
        mda_univ : MDAnalysis universe to process.
    Returns:
        MDAnalysis universe with hetero atoms removed.
    """
    return mda_univ.select_atoms("not record_type HETATM")


def filter_hetatoms(
    mda_univ: mda.Universe | mda.AtomGroup,
    keep_het_chain: list[str] | str | None = None,
) -> mda.Universe:
    r"""Filter hetero atoms.

    Args:
        mda_univ : MDAnalysis universe to process.
        keep_het_chain : Chains to keep their HETATM and remove other HETATMs. Default is None.
    Returns:
        MDAnalysis universe with filtered hetero atoms.
    """
    if keep_het_chain is None:
        return mda_univ.select_atoms("record_type HETATM")
    if isinstance(keep_het_chain, str):
        keep_het_chain = [keep_het_chain]
    selection = " or ".join(
        [f"segid {c} and record_type HETATM" for c in keep_het_chain]
    )
    return mda_univ.select_atoms(selection)


def find_common_atoms(mda_univ1: AtomGroup, mda_univ2: AtomGroup) -> Iterable[str]:
    r"""Find common atoms between two MDAnalysis universes.

    Args:
        mda_univ1 : MDAnalysis universe to process.
        mda_univ : MDAnalysis universe to process.
    Returns:
        List of common atoms.

    """
    # Get the atom types for each ligand
    ref_atom_names = [atom.name for atom in mda_univ1.atoms] if mda_univ1.atoms else []
    dock_atom_names = [atom.name for atom in mda_univ2.atoms] if mda_univ2.atoms else []

    # Find the common atoms
    common_atoms = list(set(ref_atom_names) & set(dock_atom_names))
    return common_atoms


def select_common_atoms(
    mda_univ: mda.Universe, common_atoms: Iterable[str]
) -> mda.Universe:
    r"""Select common atoms.

    Args:
        mda_univ : MDAnalysis universe to process.
        common_atoms : types of atoms to select.
    Returns:
        MDAnalysis universe for the selected atoms in the structure.
    """
    selection = " or ".join([f"name {atom}" for atom in common_atoms])
    return mda_univ.select_atoms(selection)


def write_mda_universe(mda_univ: mda.Universe, file_path: str) -> None:
    r"""Write MDAnalysis universe to file.

    Args:
        mda_univ : MDAnalysis universe to process.
        file_path : File to write to.
    """
    with mda.Writer(file_path, n_atoms=mda_univ.atoms.n_atoms) as writer:
        writer.write(mda_univ.atoms)


def read_cif(file_path: str) -> "gemmi.Structure":
    r"""Read CIF file.

    Args:
        file_path : Path to CIF file.

    Returns:
        Structure object.
    """
    # pylint: disable=c-extension-no-member
    return gemmi.read_structure(file_path)


def convert_cif2pdb(file_path: str, write_path: str) -> None:
    r"""Convert CIT file to PDB file

    Args:
        file_path : Path to cif file.
        write_path : Path to write PDB file.
    """
    # pylint: disable=c-extension-no-member
    cif_structure = read_cif(file_path)
    cif_structure.setup_entities()
    cif_structure.shorten_chain_names()
    cif_structure.write_pdb(write_path)
