r""" Utility functions for structure module. """

from typing import TextIO

import os

import MDAnalysis as mda
import pandas as pd
import requests
from loguru import logger

from ..docking.contexts import ProteinContext
from .editing import *


def is_xray_structure(pdb: str | TextIO) -> bool:
    """
    Check if the PDB file was generated from X-ray diffraction data.

    Parameters:
    -----------
    pdb : str or file-like
        Path to the PDB file. or file-like object. or just pdb id
    Returns:
    --------
    bool
        True if the PDB was generated from X-ray data, False otherwise.
    """
    # check if the pdb is a file or a pdb id
    if os.path.isfile(pdb) and os.path.exists(pdb):
        with open(pdb, "r", encoding="utf-8") as file:
            lines = file.readlines()
        ext = os.path.splitext(pdb)[-1].lower()
        if ext == ".pdb":
            expdta_line = [line for line in lines if line.startswith("EXPDTA")]
            if expdta_line:
                return "X-RAY" in expdta_line[0]
            else:
                remark_200_line = [
                    line for line in lines if line.startswith("REMARK 200")
                ]
                return bool(remark_200_line)
        elif ext == ".cif":
            # find the _exptl.method line
            exptl_line = [line for line in lines if line.startswith("_exptl.method")]
            if exptl_line:
                return "X-RAY" in exptl_line[0]
    elif isinstance(pdb, str):
        raw_data = get_rcsb_data(pdb)
        if raw_data["exptl"][0]["method"] == "X-RAY DIFFRACTION":
            return True
        else:
            return False


def is_xray_structure(pdb: str | TextIO) -> bool:
    """
    Check if the PDB file was generated from X-ray diffraction data.

    Parameters:
    -----------
    pdb : str or file-like
        Path to the PDB file. or file-like object. or just pdb id
    Returns:
    --------
    bool
        True if the PDB was generated from X-ray data, False otherwise.
    """
    # check if the pdb is a file or a pdb id
    if os.path.isfile(pdb) and os.path.exists(pdb):
        with open(pdb, "r", encoding="utf-8") as file:
            lines = file.readlines()
        ext = os.path.splitext(pdb)[-1].lower()
        if ext == ".pdb":
            expdta_line = [line for line in lines if line.startswith("EXPDTA")]
            if expdta_line:
                return "X-RAY" in expdta_line[0]
            else:
                remark_200_line = [
                    line for line in lines if line.startswith("REMARK 200")
                ]
                return bool(remark_200_line)
        elif ext == ".cif":
            # find the _exptl.method line
            exptl_line = [line for line in lines if line.startswith("_exptl.method")]
            if exptl_line:
                return "X-RAY" in exptl_line[0]
    elif isinstance(pdb, str):
        raw_data = get_rcsb_data(pdb)
        if raw_data["exptl"][0]["method"] == "X-RAY DIFFRACTION":
            return True
        else:
            return False


def separate_protein_ligand(
    pdb: str | TextIO,
    reference: str | TextIO = None,
    remove_water: bool | None = True,
    keep_het_chain: str | list | None = None,
) -> tuple["Protein", "Ligand"]:
    r"""Separate protein and ligand from a PDB file.
    Parameters
    ----------
    pdb : str or file-like
        Path to the PDB file or file-like object.
    reference : str or file-like
        Path to the Reference file or file-like object.
    remove_water : bool
        Remove crystallographic waters from the protein structures. Default is True.
    keep_het_chain : str or list
        Chain(s) to keep their HETATM in the protein structure.
        Default is None. If None, all HETATM will be kept.
    reference : str or file-like
        Path to the Reference file or file-like object.
    remove_water : bool
        Remove crystallographic waters from the protein structures. Default is True.
    keep_het_chain : str or list
        Chain(s) to keep their HETATM in the protein structure.
        Default is None. If None, all HETATM will be kept.
    Returns
    -------
    Protein
        Universe object containing the protein.
    Ligand
        Universe object containing the ligand.
    """
    pdb_obj = get_mda_universe(pdb)
    water_object = select_residues(pdb_obj, residues=["HOH"])
    # check if the file has hetatoms in chain A or not
    if keep_het_chain is not None:
        selection = select_chains(pdb_obj, chains=keep_het_chain)
        # check if the hetatoms in the selection (the residue names)
        # are valid using the protein context impurities
        impurities = ProteinContext.get_current().impurities
        valid_hetatoms = [
            hetatom.resname
            for hetatom in filter_hetatoms(selection)
            if hetatom.resname not in impurities and len(hetatom.resname) == 3
        ]
        # check if the length of hetatoms line is < 4 using mda
        logger.debug((filter_hetatoms(selection).resnames.all()))
        logger.debug((filter_hetatoms(selection).atoms.resnames.all()))
        logger.debug(all(atom == "HOH" for atom in valid_hetatoms))
        while (
            len(filter_hetatoms(selection)) == 0
            or len(valid_hetatoms) == 0
            or all(atom == "HOH" for atom in list(set(valid_hetatoms))) is True
        ):
            logger.warning(
                "No HETATM found in the selected chains.Checking another chain."
            )
            # keep changing the chain by trying B,C,D,etc until hetatoms are found
            keep_het_chain = chr(ord(keep_het_chain) + 1)
            selection = select_chains(pdb_obj, chains=keep_het_chain)
            valid_hetatoms = [
                hetatom.resname
                for hetatom in filter_hetatoms(selection)
                if hetatom.resname not in impurities and len(hetatom.resname) == 3
            ]
            keep_het_chain = keep_het_chain
        hetatm = filter_hetatoms(pdb_obj, keep_het_chain)
    else:
        logger.debug(
            f"No chain specified. Selecting all chains.i.e {list(set(pdb_obj.segments.segids))}"
        )
        # make a list of all the chains in the pdb file using the pdb_obj object
        keep_het_chain = list(set(pdb_obj.segments.segids))
        # delete empty values from the list
        keep_het_chain = [i for i in keep_het_chain if i]
        logger.debug(f"Chains in the pdb file: {keep_het_chain}")
        selection = select_chains(pdb_obj, chains=keep_het_chain)
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
                for model in reference_obj.trajectory:
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
            selection = select_chains(pdb_obj, chains=reference_chain)
            ligand = select_residues(selection, residues=reference_ligand)
        ligand = select_residues(pdb_obj, residues=reference_ligand)
        return selection.atoms, ligand
    actual_ligand = remove_residues(hetatm, residues=["HOH"])
    if remove_water:
        # select the water molecules from the hetatm
        ligand = remove_residues(hetatm, residues=["HOH"])
    else:
        ligand = merge_universes(
            [hetatm, select_chains(water_object, chains=keep_het_chain)]
        )
        logger.warning(
            "Crystallographic Water molecules are not removed from the protein structure."
        )
        logger.debug(ligand.resnames.all())
    protein = remove_hetatoms(pdb_obj)
    save_prot = merge_universes([protein, ligand])
    return save_prot, actual_ligand


def check_ligand(pdb: str | TextIO, reference: str | TextIO) -> bool:
    r"""Compare the ligand in the PDB file with the ligand in the reference file.
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


def get_rcsb_data(pdb_id: str):
    r"""
    Get entry data for a given PDB ID using the RCSB API.
    Parameters:
    ----------
    pdb_id (str): The PDB ID to retrieve the entry data for.
    Returns:
    -------
    dict: The entry data for the given PDB ID.
    """

    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url, timeout=5)
    # check if the request was successful
    if response.status_code != 200:
        logger.error(f"Error fetching data for PDB ID {pdb_id}: {response.status_code}")
        return {}
    data = response.json()
    return data


def find_resolution(pdb_id: str, rcsb_data: dict | None = None) -> float:
    r"""Find the resolution of a PDB file using the RCSB API.
    Parameters
    ----------
    pdb_id : str
        The PDB ID to find the resolution for.
    rcsb_data : dict or None
        The data for the PDB ID from the RCSB API. If None, the data will be fetched.
    Returns
    -------
    float
        The resolution of the PDB file.
    """
    if rcsb_data is not None:
        data = rcsb_data
    else:
        data = get_rcsb_data(pdb_id)
    resolution = data["rcsb_entry_info"]["resolution_combined"]
    logger.debug(f"Resolution for {pdb_id} is {resolution}")
    return float(resolution[0])


def has_covalent_bonds(pdb: str, rcsb_data: dict | None = None) -> bool:
    r"""Check if the PDB file has covalent bonds or not.
    Parameters
    ----------
    pdb : str
        the PDB id needed to be checked.
    rcsb_data : dict or None
        The data for the PDB ID from the RCSB API. If None, the data will be fetched.
    Returns
    -------
    bool
        True if the PDB file has covalent bonds, False otherwise.
    """
    if rcsb_data is not None:
        data = rcsb_data
    else:
        data = get_rcsb_data(pdb)
    logger.debug(data["rcsb_entry_info"]["inter_mol_covalent_bond_count"])
    if data["rcsb_entry_info"]["inter_mol_covalent_bond_count"] > 0:
        return True
    else:
        return False


def has_ligands(pdb: str, rcsb_data: dict | None = None) -> bool:
    r"""Check if the PDB file has ligands or not.
    Parameters
    ----------
    pdb : str
        the PDB id needed to be checked.
    rcsb_data : dict or None
        The data for the PDB ID from the RCSB API. If None, the data will be fetched.
    Returns
    -------
    bool
        True if the PDB file has ligands which can be ions/additatives, False otherwise.
    """
    if rcsb_data is not None:
        data = rcsb_data
    else:
        data = get_rcsb_data(pdb)
    tmp = data["rcsb_entry_info"]["nonpolymer_entity_count"]
    logger.debug(f"ligand/non polymer count for {pdb} is {tmp}")
    if tmp > 0:
        return True
    else:
        return False


def get_entity_ids(pdb_id: str, rcsb_data: dict | None = None) -> dict[str, list[str]]:
    r"""Get the entity IDs for a given PDB ID using the RCSB API.
    Parameters
    ----------
    pdb_id : str
        The PDB ID to retrieve the entity IDs for.
    rcsb_data : dict or None
        The data for the PDB ID from the RCSB API. If None, the data will be fetched.
    Returns
    -------
    dict
        The entity IDs for the given PDB ID as a dictionary.
        The keys are the entity type and the values are the entity numbers.
    rcsb_data : dict or None
        The data for the PDB ID from the RCSB API. If None, the data will be fetched.
    """
    if rcsb_data is not None:
        data = rcsb_data
    else:
        data = get_rcsb_data(pdb_id)
    entity_ids = {}
    polymer_entity = data["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
    nonpolymer_entity = data["rcsb_entry_container_identifiers"][
        "non_polymer_entity_ids"
    ]
    entity_ids["polymer"] = polymer_entity
    entity_ids["nonpolymer"] = nonpolymer_entity
    logger.debug(f"Polymer entity IDs: {polymer_entity}")
    logger.debug(f"Non-polymer entity IDs: {nonpolymer_entity}")
    return entity_ids


def pdb_has_mutation(pdb_id: str) -> bool:
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
    polymer_ids = get_entity_ids(pdb_id)["polymer"]
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/"
    list_of_mutations = []
    if len(polymer_ids) == 1:
        response = requests.get(url + polymer_ids[0], timeout=5)
        data = response.json()
        logger.debug(
            f"Mutations in {pdb_id}: {data['entity_poly']['rcsb_mutation_count']}"
        )
        if data["entity_poly"]["rcsb_mutation_count"] > 0:
            return True
        else:
            return False
    else:
        for entity_id in polymer_ids:
            response = requests.get(url + entity_id, timeout=5)
            data = response.json()
            list_of_mutations.append(data["entity_poly"]["rcsb_mutation_count"])
    # check if the values of the list are 0
    if all(i == 0 for i in list_of_mutations):
        logger.debug(f"Mutations in {pdb_id}: {list_of_mutations}")
        return False
    else:
        return True


def get_nonpolymer_names(pdb_id: str, rcsb_data: dict | None = None) -> list:
    r"""Get the names of the non-polymer entities in a PDB file.
    Parameters
    ----------
    pdb_id : str
        The PDB ID to retrieve the non-polymer entity names for.
    rcsb_data : dict or None
        The data for the PDB ID from the RCSB API. If None, the data will be fetched.
    Returns
    -------
    list
        The names of the non-polymer entities in the PDB file.
    """
    if rcsb_data is not None:
        data = rcsb_data
    else:
        data = get_rcsb_data(pdb_id)
    if not has_ligands(pdb_id, data):
        logger.warning(f"No non-polymer entities found in {pdb_id}.")
        return []
    nonpolymer_ids = get_entity_ids(pdb_id, data)["nonpolymer"]
    nonpolymer_names = []
    url = f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/"
    for entity_id in nonpolymer_ids:
        response = requests.get(url + entity_id, timeout=10)
        data = response.json()
        nonpolymer_names.append(data["pdbx_entity_nonpoly"]["comp_id"])
    # exclude ligands with names less than 3 characters from the list
    # highly likely they are ions
    nonpolymer_names = [i for i in nonpolymer_names if len(i) == 3]
    logger.info(f"Non-polymer entity names: {nonpolymer_names}")
    return nonpolymer_names


def validate_ligands(
    pdb: str, impurities: list | None = ProteinContext.get_current().impurities
) -> bool:
    r"""Validate the ligands from pdb id using the impurities list.
    Parameters
    ----------
    pdb : str
        The PDB ID to validate.
    impurities : list or None
        List of impurities to check against. Default is impurities from the ProteinContext.
    Returns
    -------
    bool
        True if the ligands are valid, False otherwise.
    """
    ligands = get_nonpolymer_names(pdb)
    if not validate_pdb(pdb) or len(ligands) == 0:
        return False
    else:
        logger.debug(f"Ligands in {pdb}: {ligands}")
        logger.debug(all(i in impurities for i in ligands))
        # check if the ligands are in the impurities list
        if all(i in impurities for i in ligands):
            return False
        else:
            return True


def validate_pdb(pdb_id: str) -> bool:
    r"""Validate a PDB file using the RCSB API.
    Parameters
    ----------
    pdb_id : str
        The PDB ID to validate.
    Returns
    -------
    bool
        True if the PDB file is valid (i.e has ligand, no covalent bond and no mutation), False otherwise.
    """
    data = get_rcsb_data(pdb_id)
    if (
        has_ligands(pdb_id, data)
        and not has_covalent_bonds(pdb_id, data)
        and not pdb_has_mutation(pdb_id)
        and is_xray_structure(pdb_id)
        and find_resolution(pdb_id, data) <= 3.0
    ):
        logger.info(f"The PDB file {pdb_id} is valid.")
        return True
    else:
        logger.warning(
            f"The PDB file {pdb_id} is not valid. Check ./structure/utils.py functions for more details."
        )
        return False


def get_ligand_names(pdb: str | TextIO) -> list:
    r"""Get the names of the ligands in a PDB file.
    Parameters
    ----------
    pdb : str or file-like
        Path to the PDB file or file-like object.
    Returns
    -------
    list
        The names of the ligands in the PDB file.
    """
    ligand = separate_protein_ligand(pdb)[1]
    # get the residue name of the ligand
    ligand_resname = ligand.residues.resnames
    if len(ligand_resname) > 1:
        logger.warning("The ligand has more than one residue.")
        impurities = ProteinContext.get_current().impurities
        # delete ant values with less than 3 characters from the list
        ligand_resname = list(
            set([i for i in ligand_resname if len(i) == 3 and i not in impurities])
        )
    logger.debug(f"Ligand residue name: {ligand_resname}")
    return ligand_resname


def get_smiles(ligand_resname: str | TextIO) -> dict[str, str]:
    r"""Get the SMILES string of a ligand from a PDB file.
    Parameters
    ----------
    ligand_resname : str or file-like
        The residue name of the ligand.
    Returns
    -------
    str
        The SMILES string of the ligand.
    """
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_resname}"
    response = requests.get(url, timeout=5)
    data = response.json()
    # check if the data['rcsb_chem_comp_descriptor']['smiles'] is not found
    if "smiles" not in data["rcsb_chem_comp_descriptor"]:
        logger.error(
            f"SMILES not found for {ligand_resname},Checking pdbx_chem_comp_descriptor"
        )
        # get the pdbx_chem_comp_descriptor from the data
        if "pdbx_chem_comp_descriptor" not in data:
            logger.error(f"pdbx_chem_comp_descriptor not found for {ligand_resname}")
            return {"smiles": "", "stereo_smiles": ""}
        else:
            for descriptor in data["pdbx_chem_comp_descriptor"]:
                if descriptor["type"] == "SMILES_CANONICAL":
                    smiles = descriptor["descriptor"]
                    logger.debug(f"SMILES Canonical: {smiles}")
                    return {"smiles": smiles, "stereo_smiles": smiles}
            else:
                logger.error("SMILES Canonical not found for the ligand")
    smiles = data["rcsb_chem_comp_descriptor"]["smiles"]
    stereo_smiles = data["rcsb_chem_comp_descriptor"]["smilesstereo"]
    logger.debug(f"SMILES: {smiles}")
    logger.debug(f"Stereo SMILES: {stereo_smiles}")
    # make a dictionary of the two smiles
    smiles_dict = {"smiles": smiles, "stereo_smiles": stereo_smiles}
    return smiles_dict
