import json
import os
import pickle

import rdkit
from rdkit import Chem
from tqdm import tqdm

from .utils import make_url, request_with_cache

# This list came from manually reviewing the 100 most common ligands in lig_id_counts.
ligs_to_ignore = [
    "MG",
    "SO4",
    "EDO",
    "GOL",
    "ZN",
    "CL",
    "CA",
    "CLA",
    "OHX",
    "NA",
    "K",
    "PO4",
    "MN",
    "ACT",
    "DMS",
    "UNX",
    "PEG",
    "PEB",
    "IOD",
    "HEM",
    "FMT",
    "SF4",
    "FE",
    "MPD",
    "CD",
    "CU",
    "UNL",
    "NI",
    "SR",
    "PG4",
    "BR",
    "HEC",
    "BCL",
    "PGE",
    "MES",
    "NO3",
    "CHL",
    "FE2",
    "FES",
    "LFA",
    "TRS",
    "CIT",
    "CO",
    "EPE",
    "ACY",
    "1PE",
    "IPA",
    "BME",
    "HG",
    "MLI",
    "C8E",
    "SCN",
    "FLC",
    "MRD",
    "CL7",
]

# Then I fed the above into chatgpt and asked it for additional ones. I reviewed these suggestions
# manually.
ligs_to_ignore += [
    "WAT",
    "NH4",
    "TAR",
    "CO3",
    "PO3",
    "AZI",
    "CAC",
    "IMD",
    "DTT",
    "TRI",
    "NO2",
    "ACN",
    "BCT",
    "BET",
    "CAC",
    "DMF",
    "EOH",
    "MOH",
    "MSE",
    "NCO",
    "NET",
    "PD",
    "TBU",
    "VO4",
]

ligs_to_ignore = set(ligs_to_ignore)


def _get_ligand_data_from_pdb(pdbid):
    """
    Get ligand data from PDB for a given pdbid.

    Args:
        pdbid (str): A pdbid.
    Returns:
        list: A list of lists, where each inner list is a ligand and its chain
            data.
    """

    pdbid = pdbid.upper()

    query = (
        '''{
  entry(entry_id:"'''
        + pdbid
        + """"){
    rcsb_id
    rcsb_entry_container_identifiers {
      entity_ids
    }
    polymer_entities {
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          auth_asym_id
          asym_id
          entity_id
        }
      }
      uniprots {
        rcsb_id
      }
    }
    nonpolymer_entities {
      nonpolymer_comp {
        chem_comp {
          id
        }
      }
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers {
          auth_asym_id
          asym_id
          entry_id
          entity_id
        }
        rcsb_target_neighbors {
          target_asym_id
          target_entity_id
          distance
          target_is_bound
        }
      }
    }
  }
}
"""
    )

    url = make_url("https://data.rcsb.org/graphql?query=", query)

    # Make the request
    results = json.loads(request_with_cache(url))

    entry = results["data"]["entry"]
    try:
        ligands = entry["nonpolymer_entities"]
    except:
        print("Should never get here")
        import pdb

        pdb.set_trace()
    if ligands is None:
        # There are no nonpolymer ligands.
        # In one example, ligand was a peptide (or at least a polymer).
        # In another, it was covalently bound to protein
        # In another, there truly was no ligand

        return None

    # You need to map the polymer asym_id_chains to auth_asym_id chains and entity_ids.
    polymers = entry["polymer_entities"]
    asym_id_to_auth_asym_and_entity_id = {}
    for polymer in polymers:
        for polymer_instance in polymer["polymer_entity_instances"]:
            identifiers = polymer_instance[
                "rcsb_polymer_entity_instance_container_identifiers"
            ]
            asym_id = identifiers["asym_id"]
            auth_asym_id_of_neighbor = identifiers["auth_asym_id"]
            entity_id_of_neighbor = identifiers["entity_id"]
            asym_id_to_auth_asym_and_entity_id[asym_id] = (
                auth_asym_id_of_neighbor,
                entity_id_of_neighbor,
            )

    data_to_return = []

    for ligand in ligands:
        ligand_id = ligand["nonpolymer_comp"]["chem_comp"]["id"]
        if ligand["nonpolymer_entity_instances"] is None:
            data_to_return.append(None)
            continue

        ligand_instances = ligand["nonpolymer_entity_instances"]

        # Remove instances that have no defined neighbors.
        ligand_instances = [
            x for x in ligand_instances if x["rcsb_target_neighbors"] is not None
        ]

        # Remove instances that are covalently bound to any neighbors.
        ligand_instances = [
            x
            for x in ligand_instances
            if all([y["target_is_bound"] == "N" for y in x["rcsb_target_neighbors"]])
        ]

        if len(ligand_instances) == 0:
            # There are no instances that are not covalently bound to any neighbors.
            data_to_return.append(None)
            continue

        for ligand_instance in ligand_instances:
            asym_ids_of_neighbors = [
                x["target_asym_id"] for x in ligand_instance["rcsb_target_neighbors"]
            ]

            for asym_id_of_neighbor in asym_ids_of_neighbors:
                if asym_id_of_neighbor not in asym_id_to_auth_asym_and_entity_id:
                    # print("This asym_id_of_neighbor is not in the polymers!", asym_id_of_neighbor, pdbid)
                    # In three cases, it was because the asym_id was a glycan, not a protein.
                    continue
                (
                    auth_asym_id_of_neighbor,
                    entity_id_of_neighbor,
                ) = asym_id_to_auth_asym_and_entity_id[asym_id_of_neighbor]
                data_to_return.append(
                    {
                        "pdbid": pdbid,
                        "ligand_id": ligand_id,
                        "auth_asym_id_of_neighbor": auth_asym_id_of_neighbor,
                        "asym_id_of_neighbor": asym_id_of_neighbor,
                        "entity_id_of_neighbor": entity_id_of_neighbor,
                    }
                )

    return data_to_return


def _check_if_cache():
    """
    Check if a cache exists for ligand data.

    Returns:
        str: The string "y" or "n" depending on if the cache should be used.
        ????: The data from the cache.
    """

    cache_filename = "ligand_data_resps.pk"
    use_cache = "n"
    uniq_ligs_to_keep, lig_data_resps = None, None
    if os.path.exists(cache_filename):
        # Ask if should use the cache.
        use_cache = input("Use the cache (good if debugging)? (y/n): ").lower()
        if use_cache == "y":
            with open(cache_filename, "rb") as f:
                uniq_ligs_to_keep, lig_data_resps = pickle.load(f)
        else:
            use_cache = "n"
    return use_cache, cache_filename, uniq_ligs_to_keep, lig_data_resps


def _get_lig_smiles_and_name(uniq_ligs_to_keep, use_cache="n"):
    """
    Get the smiles and name for a list of ligands.

    Args:
        uniq_ligs_to_keep (list): A list of unique ligands.
    Returns:
        dict: The ligand info. The keys are ligand ids and the values are
            dictionaries with keys "smiles" and "name".
    """

    if use_cache != "n" and os.path.exists("ligand_smiles_and_name.pk"):
        with open("ligand_smiles_and_name.pk", "rb") as f:
            return pickle.load(f)

    lig_info = {}
    fail_cnt = 0
    for lig_id in tqdm(uniq_ligs_to_keep):
        query = (
            '''{
            chem_comp(comp_id:"'''
            + lig_id
            + """"){
                chem_comp {
                    name
                }
                pdbx_reference_molecule {
                    name
                }
                rcsb_chem_comp_descriptor {
                    SMILES_stereo,
                    SMILES
                }
                pdbx_chem_comp_identifier {
                    identifier
                }
                rcsb_chem_comp_synonyms {
                    name
                }
            }
        }"""
        )

        url = make_url("https://data.rcsb.org/graphql?query=", query)
        resp = json.loads(request_with_cache(url))

        smiles = resp["data"]["chem_comp"]["rcsb_chem_comp_descriptor"]["SMILES_stereo"]
        if smiles is None:
            # Strange, but some entries don't have smiles. For example, https://www.rcsb.org/ligand/0DY
            # Let's just use "".
            fail_cnt += 1
            smiles = ""

        if resp["data"]["chem_comp"]["pdbx_chem_comp_identifier"] is None:
            name = resp["data"]["chem_comp"]["chem_comp"]["name"]
        else:
            name = resp["data"]["chem_comp"]["pdbx_chem_comp_identifier"][0][
                "identifier"
            ]

        lig_info[lig_id] = {"smiles": smiles, "name": name}

    with open("ligand_smiles_and_name.pk", "wb") as f:
        pickle.dump(lig_info, f)

    return lig_info


def _filter_ligands(lig_info):
    """
    Filter ligands.

    Args:
        lig_info (dict): The ligand info. The keys are ligand ids and the values are
            dictionaries with keys "smiles" and "name".
    Returns:
        ????: The ligand info (updated) and more ligands to ignore.
    """

    more_ligs_to_ignore = set([])

    for lig_id, lig_data in lig_info.items():
        smiles = lig_data["smiles"]

        # If the smiles string has no carbons in it, ignore this ligand.
        if "C" not in smiles and not "c" in smiles:
            more_ligs_to_ignore.add(lig_id)
            continue

        # If the ligand has fewer than 5 atoms, let's ignore it.
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # Can't process this molecule. Will have to skip. But still probably good, so store the info.
            continue

        if mol.GetNumAtoms() < 5:
            more_ligs_to_ignore.add(lig_id)
            continue

        # This is perhaps a bit controversial (good to discuss), but let's ignore large
        # ligands without rings. These tend to be extended aliphatic chains that I
        # don't think are very drug-like. Also, hard to dock.
        if mol.GetNumAtoms() > 30 and not mol.HasSubstructMatch(
            Chem.MolFromSmarts("[r]")
        ):
            more_ligs_to_ignore.add(lig_id)
            continue

    # Now remove the ligands we want to ignore.
    for lig_id in more_ligs_to_ignore:
        del lig_info[lig_id]

    return lig_info, more_ligs_to_ignore


def _get_ligs_per_chain(lig_data_resps):
    """
    Get ligands per chain.

    Args:
        lig_data_resps (???): The ligand data.
    Returns:
        ????: ???
    """

    pdbid_intid_to_pdbid_chain = {}
    pdbid_chain_to_ligs = {}
    lig_id_counts = {}

    for resp in tqdm(lig_data_resps):
        lig_id = resp["ligand_id"]
        pdbid = resp["pdbid"]
        entity_id = resp["entity_id_of_neighbor"]
        chain_id = resp["auth_asym_id_of_neighbor"]

        if lig_id in ligs_to_ignore:
            continue
        if lig_id not in lig_id_counts:
            lig_id_counts[lig_id] = 0
        lig_id_counts[lig_id] += 1

        pdbid_intid = pdbid + "_" + entity_id
        pdbid_chain = pdbid + "_" + chain_id

        if not pdbid_intid in pdbid_intid_to_pdbid_chain:
            pdbid_intid_to_pdbid_chain[pdbid_intid] = set([])
        pdbid_intid_to_pdbid_chain[pdbid_intid].add(pdbid_chain)

        if pdbid_chain not in pdbid_chain_to_ligs:
            pdbid_chain_to_ligs[pdbid_chain] = set([])
        pdbid_chain_to_ligs[pdbid_chain].add(lig_id)

    # Now map the pdbids (integer chain) to the associated ligands.
    pdbid_intid_to_ligs = {}
    for pdbid_intid in pdbid_intid_to_pdbid_chain:
        uniq_ligs = set([])
        for pdbid_chain in pdbid_intid_to_pdbid_chain[pdbid_intid]:
            if pdbid_chain in pdbid_chain_to_ligs:
                uniq_ligs.update(pdbid_chain_to_ligs[pdbid_chain])
        pdbid_intid_to_ligs[pdbid_intid] = uniq_ligs

    return pdbid_intid_to_pdbid_chain, pdbid_chain_to_ligs, pdbid_intid_to_ligs


def collect_all_lig_data(uniq_pdbids):
    """
    Collect all ligand data for a list of pdbids.

    Args:
        uniq_pdbids (list): A list of unique pdbids.
    Returns:
        ????: The ligand data.
    """

    # Let's get all the remaining ligands. This can take a while, so consider using a cache.

    use_cache, cache_filename, uniq_ligs_to_keep, lig_data_resps = _check_if_cache()

    if use_cache == "n":
        uniq_ligs_to_keep = set([])
        lig_data_resps = []
        failed_to_get_lig_data = set([])
        # for pdbid in tqdm(list(uniq_pdbids)[:10000]):
        for pdbid in tqdm(uniq_pdbids):
            resp = _get_ligand_data_from_pdb(pdbid)
            if resp is None:
                failed_to_get_lig_data.add(pdbid)
                continue
            for r in resp:
                if r is not None and r["ligand_id"] in ligs_to_ignore:
                    continue
                if r is not None:
                    uniq_ligs_to_keep.add(r["ligand_id"])
                lig_data_resps.append(r)

        if len(failed_to_get_lig_data) > 0:
            print("Failed to get ligand data for some pdbids.")
            print(
                100 * len(failed_to_get_lig_data) / len(uniq_pdbids),
                "% of pdbids failed to get ligand data.",
            )

        lig_data_resps_failed = len([x for x in lig_data_resps if x is None])
        if lig_data_resps_failed > 0:
            print("Failed to get ligand data for some ligands.")
            print(
                100 * lig_data_resps_failed / len(lig_data_resps),
                "ligands failed to get data.",
            )

        lig_data_resps = [x for x in lig_data_resps if x is not None]

        # Save the data to a cache.
        with open(cache_filename, "wb") as f:
            pickle.dump([uniq_ligs_to_keep, lig_data_resps], f)

    # Now get all the ligands smiles and names.
    lig_info = _get_lig_smiles_and_name(uniq_ligs_to_keep, use_cache)

    lig_info, more_ligs_to_ignore = _filter_ligands(lig_info)

    ligs_to_ignore.update(more_ligs_to_ignore)

    return _get_ligs_per_chain(lig_data_resps)
