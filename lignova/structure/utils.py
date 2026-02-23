r""" Utility functions for structure module. """

from typing import Any, Literal, TextIO

import os
import time

import MDAnalysis as mda
import pandas as pd
import requests
from loguru import logger

from ..docking.contexts import ProteinContext
from .editing import (
    filter_hetatoms,
    get_mda_universe,
    merge_universes,
    remove_hetatoms,
    remove_residues,
    select_chains,
    select_residues,
    select_water,
    validate_chains,
)


def is_xray_structure(pdb: str | TextIO) -> bool:
    """
    Check if the PDB file was generated from X-ray diffraction data.

    Arg:
        pdb : Path to the PDB file. or file-like object. or just pdb id

    Returns:
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
            remark_200_line = [line for line in lines if line.startswith("REMARK 200")]
            return bool(remark_200_line)
        if ext == ".cif":
            # find the _exptl.method line
            exptl_line = [line for line in lines if line.startswith("_exptl.method")]
            if exptl_line:
                return "X-RAY" in exptl_line[0]
    if isinstance(pdb, str):
        raw_data = get_rcsb_data(pdb)
        return raw_data["exptl"][0]["method"] == "X-RAY DIFFRACTION"
    return False


def chery_pick_ligand(
    pdb: str | TextIO,
    ligand: str,
    remove_water: bool = True,
    water_selection: Literal["surface", "interfacial", "all"] | None = None,
) -> Literal["Protein", "Ligand"]:
    r"""Cherry pick a ligand from a PDB file.

    Args:
        pdb : str or file-like
            Path to the PDB file or file-like object.
        ligand : str
            The ligand to cherry pick.
        remove_water : bool
            Remove crystallographic waters from the protein structures. Default is True.
        water_selection : str
            The selection of water molecules to keep. Default is "none".

    Returns:
        Universe object containing the protein.

        Universe object containing the ligand.
    """
    if water_selection is None and remove_water is False:
        logger.warning("No water selection specified. Retaining all water molecules.")
        raise ValueError("Please specify a water selection.")
    pdb_obj = get_mda_universe(pdb)
    water_object = select_residues(pdb_obj, residues=["HOH"])
    # check if the ligand is in chain A and
    # if not change chains till it is found using filter_hetatoms
    selection = select_chains(pdb_obj, chains="A")
    ligand_obj = select_residues(selection, residues=ligand)
    chains = list(set(pdb_obj.segments.segids))
    chains.sort()
    index = 0
    save_prot = None
    while (
        len(ligand_obj) == 0
        or ligand not in ligand_obj.resnames
        and index < len(chains)
    ):
        logger.warning(f"No HETATM found in chain A. Checking chain.{chains[index]}")
        selection = select_chains(pdb_obj, chains=chains[index])
        ligand_obj = filter_hetatoms(selection)
        ligand_obj = remove_residues(ligand_obj, residues=["HOH"])
        water_object = select_residues(selection, residues=["HOH"])
        if index + 1 < len(chains):
            index += 1
        else:
            break
        # change the chain to the next chain in the pdb file
    if remove_water:
        save_prot = merge_universes([remove_hetatoms(pdb_obj), ligand_obj])
    else:
        if water_selection == "all":
            save_prot = merge_universes(
                [remove_hetatoms(pdb_obj), ligand_obj, water_object]
            )
        elif water_selection == "surface":
            surface_water = select_water(pdb_obj, ligand, water_selection)
            save_prot = merge_universes(
                [remove_hetatoms(pdb_obj), ligand_obj, surface_water]
            )
        elif water_selection == "interfacial":
            interfacial_water = select_water(pdb_obj, ligand, water_selection)
            save_prot = merge_universes(
                [remove_hetatoms(pdb_obj), ligand_obj, interfacial_water]
            )
    return save_prot, ligand_obj


def separate_protein_ligand(
    pdb: str | TextIO,
    remove_water: bool | None = True,
    keep_het_chain: str | list | None = None,
    water_selection: Literal["surface", "interfacial", "all"] | None = None,
    hetatm: Literal["valid_ligand", "no_hetam", "cofactors"] = "valid_ligand",
) -> tuple[mda.Universe, mda.Universe]:
    r"""Separate protein and ligand from a PDB file.

    Args:
        pdb : str or file-like
            Path to the PDB file or file-like object.
        remove_water : bool
            Remove crystallographic waters from the protein structures. Default is True.
        keep_het_chain : str or list
            Chain(s) to keep their HETATM in the protein structure.
            Default is None. If None, all HETATM will be kept.
        water_selection : str
            The selection of water molecules to keep if remove_water is False. Default is "none"
        hetatm : str
            The selection of hetatoms to keep in the protein structure. Default is "valid_ligand".
    Returns:
        Universe object containing the protein.

        Universe object containing the ligand.
    """
    save_prot = None
    if hetatm not in ["valid_ligand", "no_hetam", "cofactors"]:
        logger.warning(
            f"Invalid option for hetatm: {hetatm}. Defaulting to 'valid_ligand'."
        )
        hetatm = "valid_ligand"
    if water_selection is None and remove_water is False:
        logger.warning("No water selection specified. Retaining all water molecules.")
        raise ValueError("Please specify a water selection.")
    pdb_obj = get_mda_universe(pdb)
    logger.debug(f"Chains in the pdb file: {list(set(pdb_obj.segments.segids))}")
    water_object = select_residues(pdb_obj, residues=["HOH"])
    impurities = ProteinContext.get_current().impurities
    cofactor_obj = select_residues(
        pdb_obj, residues=ProteinContext.get_current().cofactors
    )
    metal_tmpobj = filter_hetatoms(pdb_obj)
    # Extract only the HETATMs with residue names of length less than 3
    metal_obj = select_residues(
        metal_tmpobj,
        residues=[res.resname for res in metal_tmpobj.residues if len(res.resname) < 3],
    )
    # check if the file has hetatoms in chain A or not
    if keep_het_chain is not None:
        if validate_chains(pdb_obj, keep_het_chain) is False:
            logger.warning(
                f"Chain {keep_het_chain} not found in the pdb file. Changing to the first chain."
            )
            # change the chain to the first chain in the pdb file
            keep_het_chain = list(set(pdb_obj.segments.segids))[0]
            logger.warning(f"Chain changed to {keep_het_chain}")
        logger.debug(f"Chain specified: {keep_het_chain}")
        selection = select_chains(pdb_obj, chains=keep_het_chain)
        # check if the hetatoms in the selection (the residue names)
        # are valid using the protein context impurities
        valid_hetatoms = [
            hetatom.resname
            for hetatom in filter_hetatoms(selection)
            if hetatom.resname not in impurities and len(hetatom.resname) == 3
        ]
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
            # Remove self-assignment
            # keep_het_chain = keep_het_chain
        hetatm_selection = filter_hetatoms(pdb_obj, keep_het_chain)
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
        hetatm_selection = filter_hetatoms(pdb_obj)
    # delete any hetatoms in hetatm that exist in the crystal additive list
    hetatm_selection = remove_residues(hetatm_selection, residues=impurities)
    actual_ligand = remove_residues(
        hetatm_selection,
        residues=["HOH"]
        + [res.resname for res in hetatm_selection.residues if len(res.resname) < 3],
    )
    ligand_name = actual_ligand.resnames.all()
    if remove_water:
        if hetatm == "valid_ligand":
            save_prot = merge_universes(
                [remove_hetatoms(pdb_obj), actual_ligand, cofactor_obj, metal_obj]
            )
        elif hetatm == "cofactors":
            save_prot = merge_universes(
                [remove_hetatoms(pdb_obj), cofactor_obj, metal_obj]
            )
        else:
            save_prot = remove_hetatoms(pdb_obj)
    else:
        logger.warning(
            "Crystallographic Water molecules are not removed from the structure."
        )
        if water_selection == "all":
            if hetatm == "valid_ligand":
                save_prot = merge_universes(
                    [
                        remove_hetatoms(pdb_obj),
                        actual_ligand,
                        water_object,
                        cofactor_obj,
                        metal_obj,
                    ]
                )
            elif hetatm == "cofactors":
                save_prot = merge_universes(
                    [remove_hetatoms(pdb_obj), water_object, cofactor_obj, metal_obj]
                )
            else:
                save_prot = merge_universes([remove_hetatoms(pdb_obj), water_object])
        elif water_selection == "surface":
            surface_water = select_water(
                pdb=pdb_obj, ligand=ligand_name, water_selection=water_selection
            )
            print(surface_water)
            if hetatm == "valid_ligand":
                save_prot = merge_universes(
                    [
                        remove_hetatoms(pdb_obj),
                        actual_ligand,
                        surface_water,
                        cofactor_obj,
                        metal_obj,
                    ]
                )
            elif hetatm == "cofactors":
                save_prot = merge_universes(
                    [remove_hetatoms(pdb_obj), surface_water, cofactor_obj, metal_obj]
                )
            else:
                save_prot = merge_universes([remove_hetatoms(pdb_obj), surface_water])
        elif water_selection == "interfacial":
            interfacial_water = select_water(
                pdb=pdb_obj, ligand=ligand_name, water_selection=water_selection
            )
            if hetatm == "valid_ligand":
                save_prot = merge_universes(
                    [
                        remove_hetatoms(pdb_obj),
                        actual_ligand,
                        interfacial_water,
                        cofactor_obj,
                        metal_obj,
                    ]
                )
            elif hetatm == "cofactors":
                save_prot = merge_universes(
                    [
                        remove_hetatoms(pdb_obj),
                        interfacial_water,
                        cofactor_obj,
                        metal_obj,
                    ]
                )
            else:
                save_prot = merge_universes(
                    [remove_hetatoms(pdb_obj), interfacial_water]
                )
    return save_prot, actual_ligand


def check_ligand(pdb: str | TextIO, reference: str | TextIO) -> bool:
    r"""Compare the ligand in the PDB file with the ligand in the reference file.

    Args:
        pdb : str or file-like
            Path to the PDB file or file-like object.
        reference : str or file-like
            Path to the csv file or file-like object.

    Returns:
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
    return False


def get_rcsb_data(pdb_id: str):
    r"""
    Get entry data for a given PDB ID using the RCSB API.

    Args:
        pdb_id (str): The PDB ID to retrieve the entry data for.

    Returns:
        The entry data for the given PDB ID.
    """

    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url, timeout=20)
    # check if the request was successful
    if response.status_code != 200:
        logger.error(f"Error fetching data for PDB ID {pdb_id}: {response.status_code}")
        return {}
    data = response.json()
    return data


def find_resolution(pdb_id: str, rcsb_data: dict | None = None) -> float:
    r"""Find the resolution of a PDB file using the RCSB API.

    Args:
        pdb_id : str
            The PDB ID to find the resolution for.
        rcsb_data : dict or None
            The data for the PDB ID from the RCSB API. If None, the data will be fetched.

    Returns:
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

    Args:
        pdb : str
            the PDB id needed to be checked.
        rcsb_data : dict or None
            The data for the PDB ID from the RCSB API. If None, the data will be fetched.

    Returns:
        True if the PDB file has covalent bonds, False otherwise.
    """
    if rcsb_data is not None:
        data = rcsb_data
    else:
        data = get_rcsb_data(pdb)
    # logger.debug(data["rcsb_entry_info"]["inter_mol_covalent_bond_count"])
    return data["rcsb_entry_info"]["inter_mol_covalent_bond_count"] > 0


def has_ligands(pdb: str, rcsb_data: dict | None = None) -> bool:
    r"""Check if the PDB file has ligands or not.

    Args:
        pdb : str
            the PDB id needed to be checked.
        rcsb_data : dict or None
            The data for the PDB ID from the RCSB API. If None, the data will be fetched.

    Returns:
        True if the PDB file has ligands which can be ions/additatives, False otherwise.
    """
    if rcsb_data is not None:
        data = rcsb_data
    else:
        data = get_rcsb_data(pdb)
    # tmp = data["rcsb_entry_info"]["nonpolymer_entity_count"]
    # logger.debug(f"ligand/non polymer count for {pdb} is {tmp}")
    return data["rcsb_entry_info"]["nonpolymer_entity_count"] > 0


def get_entity_ids(pdb_id: str, rcsb_data: dict | None = None) -> dict[str, list[str]]:
    r"""Get the entity IDs for a given PDB ID using the RCSB API.

    Args:
        pdb_id : str
            The PDB ID to retrieve the entity IDs for.
        rcsb_data : dict or None
            The data for the PDB ID from the RCSB API. If None, the data will be fetched.

    Returns:
        The entity IDs for the given PDB ID as a dictionary.
            The keys are the entity type and the values are the entity numbers.

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
    # logger.debug(f"Polymer entity IDs: {polymer_entity}")
    # logger.debug(f"Non-polymer entity IDs: {nonpolymer_entity}")
    return entity_ids


def pdb_has_mutation(pdb_id: str, rcsb_data: dict | None = None) -> bool:
    r"""Check if a PDB file has mutations or not

    Args:
        pdb : str
            Path to the PDB file.
        rcsb_data : dict or None
            The data for the PDB ID from the RCSB API. If None, the data will be fetched.

    Returns:
        True if the PDB file has mutations, False otherwise.
    """
    if rcsb_data is not None:
        polymer_ids = get_entity_ids(pdb_id, rcsb_data)["polymer"]
    polymer_ids = get_entity_ids(pdb_id)["polymer"]
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/"
    list_of_mutations = []
    if len(polymer_ids) == 1:
        response = requests.get(url + polymer_ids[0], timeout=20)
        data = response.json()
        logger.debug(
            f"Mutations in {pdb_id}: {data['entity_poly']['rcsb_mutation_count']}"
        )
        return data["entity_poly"]["rcsb_mutation_count"] > 0
    for entity_id in polymer_ids:
        response = requests.get(url + entity_id, timeout=20)
        data = response.json()
        list_of_mutations.append(data["entity_poly"]["rcsb_mutation_count"])
    # check if the values of the list are 0
    if all(i == 0 for i in list_of_mutations):
        # logger.debug(f"Mutations in {pdb_id}: {list_of_mutations}")
        return False
    return True


def get_nonpolymer_names(pdb_id: str, rcsb_data: dict | None = None) -> list:
    r"""Get the names of the non-polymer entities in a PDB file.

    Args:
        pdb_id : str
            The PDB ID to retrieve the non-polymer entity names for.
        rcsb_data : dict or None
            The data for the PDB ID from the RCSB API. If None, the data will be fetched.

    Returns:
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
        response = requests.get(url + entity_id, timeout=20)
        data = response.json()
        nonpolymer_names.append(data["pdbx_entity_nonpoly"]["comp_id"])
    # exclude ligands with names less than 3 characters from the list
    # highly likely they are ions
    nonpolymer_names = [i for i in nonpolymer_names if len(i) == 3]
    # logger.info(f"Non-polymer entity names: {nonpolymer_names}")
    return nonpolymer_names


def validate_ligands(
    pdb: str, impurities: list | None = ProteinContext.get_current().impurities
) -> bool:
    r"""Validate the ligands from pdb id using the impurities list.

    Args:
        pdb : str
            The PDB ID to validate.
        impurities : list or None
            List of impurities to check against. Default is impurities from the ProteinContext.

    Returns:
        True if the ligands are valid, False otherwise.
    """
    ligands = get_nonpolymer_names(pdb)
    if len(ligands) == 0:
        return False
    # logger.debug(f"Ligands in {pdb}: {ligands}")
    # logger.debug(all(i in impurities for i in ligands))
    return not all(i in impurities for i in ligands)


def validate_pdb(pdb_id: str) -> bool:
    r"""Validate a PDB file using the RCSB API.

    Args:
        pdb_id : str
            The PDB ID to validate.

    Returns:
        True if the PDB file is valid (i.e has ligand,
            no covalent bond and no mutation), False otherwise.
    """
    data = get_rcsb_data(pdb_id)
    if len(data) == 0:
        logger.error(f"Failed to fetch data for PDB ID {pdb_id}.")
        return False
    if (
        has_ligands(pdb_id, data)
        and not has_covalent_bonds(pdb_id, data)
        and not pdb_has_mutation(pdb_id, data)
        and is_xray_structure(pdb_id)
        and find_resolution(pdb_id, data) <= 3.0
    ):
        logger.info(f"The PDB file {pdb_id} is valid.")
        return True

    logger.warning(
        f"The PDB file {pdb_id} is not valid. "
        f"Check ./structure/utils.py functions for more details."
    )
    return False


def get_ligand_names(pdb: str | TextIO) -> list:
    r"""Get the names of the ligands in a PDB file.

    Args:
        pdb : str or file-like
            Path to the PDB file or file-like object.

    Returns:
        The names of the ligands in the PDB file.
    """
    _, ligand = separate_protein_ligand(pdb)
    # get the residue name of the ligand
    ligand_resname = ligand.residues.resnames
    if len(ligand_resname) > 1:
        logger.warning("The ligand has more than one residue.")
        impurities = ProteinContext.get_current().impurities
        # delete ant values with less than 3 characters from the list

        ligand_resname = {
            i for i in ligand_resname if len(i) == 3 and i not in impurities
        }
        ligand_resname = list(ligand_resname)
    # logger.debug(f"Ligand residue name: {ligand_resname}")
    return ligand_resname


def get_smiles(ligand_resname: str | TextIO) -> dict[str, str]:
    r"""Get the SMILES string of a ligand from a PDB file.

    Args:
        ligand_resname : str or file-like
            The residue name of the ligand.

    Returns:
        The SMILES string of the ligand.
    """
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_resname}"
    response = requests.get(url, timeout=20)
    data = response.json()
    # check if the data['rcsb_chem_comp_descriptor']['smiles'] is not found
    if "smiles" not in data["rcsb_chem_comp_descriptor"]:
        logger.error(
            f"SMILES not found for {ligand_resname},Checking pdbx_chem_comp_descriptor"
        )
        if "pdbx_chem_comp_descriptor" not in data:
            logger.error(f"pdbx_chem_comp_descriptor not found for {ligand_resname}")
            return {"smiles": "", "stereo_smiles": ""}

        for descriptor in data["pdbx_chem_comp_descriptor"]:
            if descriptor["type"] == "SMILES_CANONICAL":
                smiles = descriptor["descriptor"]
                logger.debug(f"SMILES Canonical: {smiles}")
                return {"smiles": smiles, "stereo_smiles": smiles}

        logger.error("SMILES Canonical not found for the ligand")
        return {"smiles": "", "stereo_smiles": ""}

    smiles = data["rcsb_chem_comp_descriptor"]["smiles"]
    stereo_smiles = data["rcsb_chem_comp_descriptor"]["smilesstereo"]
    # logger.debug(f"SMILES: {smiles}")
    # logger.debug(f"Stereo SMILES: {stereo_smiles}")
    return {"smiles": smiles, "stereo_smiles": stereo_smiles}


def map_genid_to_pdb(gene_ids: list[str]) -> list[dict]:
    r"""Map a list of gene IDs to PDB IDs using the UniProt ID Mapping API.

    Args:
        gene_ids : list of str
            The list of gene IDs to map.

    Returns:
        The mapping of each gene ID to the PDB ID and other attributes.
    """
    url = "https://rest.uniprot.org/idmapping/run"
    payload = {"from": "GeneID", "to": "UniProtKB", "ids": ",".join(gene_ids)}
    response = requests.post(url, data=payload, timeout=5)
    job_id = response.json().get("jobId")
    if not job_id:
        logger.error(f"Failed to retrieve job ID for gene IDs {gene_ids}")
        return []

    # Check the status of the job
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    status_response = requests.get(status_url, timeout=5)
    if status_response.status_code != 200:
        logger.error(
            f"Error checking job status for job ID {job_id}: {status_response.status_code}"
        )
        return []

    # Get the detailed results
    url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job_id}?format=json"
    response = requests.get(url, timeout=5)
    while response.status_code != 200:
        logger.error(f"Job ID {job_id} is not ready. Retrying in 5 seconds.")
        time.sleep(5)
        response = requests.get(url, timeout=5)

    logger.debug(f"Job ID {job_id} is ready.")
    results = response.json()
    if not results["results"]:
        logger.error(f"No results found for gene IDs {gene_ids}")
        return []
    logger.debug(f"Results found for {len(gene_ids)} gene IDs.")
    # Parse the results to extract all the attributes and save in a list of dictionaries
    uniprot_results = []
    for data in results["results"]:
        uniprot_result = {
            "Gene ID": data["from"],
            "UniprotID": data["to"]["primaryAccession"],
            "Organism": data["to"]["organism"]["scientificName"],
            "Protein Name": (
                data["to"]["proteinDescription"]["recommendedName"]["fullName"]["value"]
                if "recommendedName" in data["to"]["proteinDescription"]
                else data["to"]["proteinDescription"]["submissionNames"][0]["fullName"][
                    "value"
                ]
            ),
            "Gene Name": (
                data["to"]["genes"][0]["geneName"]["value"]
                if "genes" in data["to"] and "geneName" in data["to"]["genes"][0]
                else ""
            ),
            "PDB IDs": [
                ref["id"]
                for ref in data["to"]["uniProtKBCrossReferences"]
                if ref["database"] == "PDB"
            ],
            "AlphaFold IDs": [
                ref["id"]
                for ref in data["to"]["uniProtKBCrossReferences"]
                if ref["database"] == "AlphaFoldDB"
            ],
        }
        uniprot_results.append(uniprot_result)

    return uniprot_results
