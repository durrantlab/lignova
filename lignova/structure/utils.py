r""" Utility functions for structure module. """
from typing import TextIO, Union

import ast
import os

import MDAnalysis as mda
import pandas as pd
import requests
from loguru import logger

from .editing import *


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
    with open(pdb_file, "r", encoding="utf-8") as file:
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


def separate_protein_ligand(
    pdb: Union[str, TextIO],
    reference: Union[str, TextIO] = None,
    remove_water: Union[bool, None] = True,
    keep_het_chain: Union[str, list, None] = None,
) -> tuple["Protein", "Ligand"]:
    r"""Separate protein and ligand from a PDB file.
    Parameters
    ----------
    pdb : str or file-like
        Path to the PDB file or file-like object.
    reference : str or file-like
        Path to the Reference file or file-like object.
    remove_water : bool
        Remove water molecules from the ligand structure. Default is True.
    keep_het_chain : str or list
        Chain(s) to keep their HETATM in the protein structure. Default is None. If None, all HETATM will be kept.
    Returns
    -------
    Protein
        Universe object containing the protein.
    Ligand
        Universe object containing the ligand.
    """
    pdb_obj = get_mda_universe(pdb)
    # check if the file has hetatoms in chain A or not
    if keep_het_chain is not None:
        selection = select_chains(pdb_obj, chains=keep_het_chain)
        # check if the selection has hetatoms
        while len(filter_hetatoms(selection)) == 0:
            logger.warning(
                "No HETATM found in the selected chains.Checking another chain."
            )
            # keep changing the chain by trying B,C,D,etc until hetatoms are found
            keep_het_chain = chr(ord(keep_het_chain) + 1)
            selection = select_chains(pdb_obj, chains=keep_het_chain)
        hetatm = filter_hetatoms(pdb_obj, keep_het_chain)
    else:
        logger.debug(
            f"No chain specified. Selecting all chains.i.e {pdb_obj.segments.segids}"
        )
        selection = select_chains(pdb_obj, chains=pdb_obj.segments.segids)
        hetatm = filter_hetatoms(pdb_obj)
    if reference is not None:
        reference_obj = get_mda_universe(reference)
        reference_chain = list(reference_obj.segments.segids)[0]
        logger.debug(f"The reference chain(s) : {reference_chain}")
        reference_ligand = set((reference_obj.residues.resids))
        # convert the set to a list of strings
        reference_ligand = [str(i) for i in reference_ligand]
        logger.debug(f"The reference resid(s) : {reference_ligand}")
        # check if the reference ligand more than one residue
        if len(reference_ligand) > 1:
            # get the resnames of the ligand
            reference_rename = list(set(reference_obj.residues.resnames))
            # remove resnames with less than 3 characters from the list in one line
            reference_resname = [i for i in reference_rename if len(i) == 3]
            # get the resids of the reference resname
            logger.debug(f"The reference resname(s) : {reference_resname}")
            # filter out that resname from the reference_obj
            with mda.Writer("tmp.pdb", multiframe=True) as writer:
                for ts in reference_obj.trajectory:
                    # Select atoms belonging to the specified residue name in the current frame
                    selection = " or ".join([f"resname {r}" for r in reference_resname])
                    selected_atoms = reference_obj.select_atoms(selection)
                    writer.write(selected_atoms)
            # rename tmp.pdb to reference.pdb
            os.rename("tmp.pdb", reference)
            reference_ligand = reference_resname
        if check_ligand(pdb, reference) is False:
            logger.warning("The ligand is not the same as the reference file.")
            # get the chains from the reference file
            """
            # select the chains from the pdb f
            reference = pd.read_csv(reference)
            print(os.path.basename(pdb).split("_")[0].upper())
            reference = reference[
                reference["PDB"] == os.path.basename(pdb).split("_")[0].upper()
            ]
            print(reference)
            reference_chain = reference["CHAIN"].values
            print(reference_chain)
            reference_ligand = list(set(reference["LIGAND"].values))
            print(reference_ligand[0])
            #ast.literal_eval(reference_ligand[0])
            if len(reference_chain) > 1:
                reference_chain = reference_chain[0]
            """
            selection = select_chains(pdb_obj, chains=reference_chain)
            ligand = select_residues(selection, residues=reference_ligand)
        ligand = select_residues(pdb_obj, residues=reference_ligand)
        return selection.atoms, ligand
    if remove_water:
        ligand = remove_residues(hetatm, residues=["HOH"])
    else:
        ligand = hetatm
        logger.warning("Water molecules are not removed from the ligand structure.")
        logger.debug(ligand.resnames.all())
        """
        for i in list(set(hetatm.resnames)):
            if len(i) != 3:
                hetatm = remove_residues(hetatm, residues=[i])
        selection_2 = select_chains(pdb_obj, chains=hetatm.segments.segids)
        return selection_2.atoms, ligand
        """
    protein = remove_hetatoms(pdb_obj)
    save_prot = merge_universes([protein, ligand])
    return save_prot, ligand


def check_ligand(pdb: Union[str, TextIO], reference: Union[str, TextIO]) -> bool:
    r"""Get ligand from a PDB file.
    Parameters
    ----------
    pdb : str or file-like
        Path to the PDB file or file-like object.
    reference : str or file-like
        Path to the csv file or file-like object.
    Returns
    -------
    bool
        True if the ligand is the same in the two PDB files, False otherwise.
    """
    # get the filename from the path
    pdb_obj = get_mda_universe(pdb)
    # get the chains id from the mdanalysis universe
    chains = list(set(pdb_obj.segments.segids))
    # get the ligand from the
    ligand = filter_hetatoms(pdb_obj)
    ligand_id = ligand.residues.resids
    # check the reference file extension
    if os.path.splitext(reference)[-1] == ".pdb":
        ref_obj = get_mda_universe(reference)
        reference_chain = list(set(ref_obj.segments.segids))
        reference_ligand = list(set(ref_obj.residues.resids))
    elif os.path.splitext(reference)[-1] == ".csv":
        reference = pd.read_csv(reference)
        # find the pdb in the reference file first column
        reference = pd.read_csv(reference)
        reference = reference[reference["PDB"] == os.path.basename(pdb).upper()]
        reference_ligand = reference["LIGAND"].values
        reference_chain = reference["CHAIN"].values
    else:
        logger.error("The reference file is not in the right format.")
    # check that all values in ligand_id are in reference_ligand and chains in reference_chain
    if all(i in reference_ligand for i in ligand_id) and all(
        i in reference_chain for i in chains
    ):
        return True
    else:
        return False


# TODO: FIX THIS FUNCTION AND TRY THE CHAIN EDITS
def check_pdb_mutation(pdb: str) -> bool:
    r"""Check if a PDB file has mutations or not
    Parameters
    ----------
    pdb : str
        Path to the PDB file.
    Returns
    -------
    bool
        True if the PDB file has mutations, False otherwise.
    """
    url = f"https://data.rcsb.org/rest/v1/core/entity_p/{pdb.upper()}"  # Assuming entity_id is 1
    response = requests.get(url)
    logger.debug(url)
    if response.status_code == 200:
        data = response.json()
        mutation_count = data.get("entity_poly", {}).get("rcsb_mutation_count", 0)
        logger.debug(f"Mutation count for PDB ID {pdb}: {mutation_count}")
        if mutation_count > 0:
            return True
        else:
            return False
    else:
        print(f"Error fetching data for PDB ID {pdb}: {response.status_code}")
        return False
