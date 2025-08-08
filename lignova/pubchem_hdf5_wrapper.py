r""" Implementation of the PubChem wrapper."""

import time

import h5py
import requests
from loguru import logger

from lignova.hdf5.pubchem import PubChemAPI


def process_hdf5_file(hdf5_file_path: str) -> None:
    r"""Process the HDF5 file.

    Args:
        hdf5_file_path : Path to the HDF5 file.
    """

    # open hdf5 file
    hdf5_file = h5py.File(hdf5_file_path, "r")
    no_uniprot = []
    yes_uniprot = []
    no_targets = []
    no_cids = []
    no_smiles = []
    yes_smiles = []
    yes_cids = []
    yes_targets = []
    yes_seq = []
    total_cids = []

    # NOTE: run the following after the previous loop to check the current status of the file
    for aid in list(hdf5_file["/aids"].keys()):
        if "targets_gene_id" in hdf5_file["/aids/" + str(aid)].keys():
            yes_targets.append(aid)
        else:
            no_targets.append(aid)
            # print('aid '+str(aid)+' has no targets_gene_id')
        # print(hdf5_file['/aids/'+str(aid)+'/uniprot_pdb_alphafold'])
        if "protein_sequence" in hdf5_file["/aids/" + str(aid)].keys():
            # print(len(hdf5_file['/aids/'+str(aid)+'/protein_sequence'].asstr()[:]))
            yes_seq.append(aid)
        if "uniprot_pdb_alphafold" in hdf5_file["/aids/" + str(aid)].keys():
            yes_uniprot.append(aid)
            # print('aid '+str(aid)+' has no uniprot_pdb_alphafold')
        else:
            no_uniprot.append(aid)
        if "cids" in hdf5_file["/aids/" + str(aid)].keys():
            total_cids.extend((list(hdf5_file["/aids/" + str(aid) + "/cids"].keys())))
            for i in hdf5_file["/aids/" + str(aid) + "/cids"].keys():
                yes_cids.append(aid)
                if (
                    "smiles"
                    in hdf5_file["/aids/" + str(aid) + "/cids/" + str(i)].keys()
                ):
                    #    print('cid '+str(i)+' has no smiles')
                    yes_smiles.append(i)
                else:
                    no_smiles.append(i)
        else:
            no_cids.append(aid)
    aids = list(hdf5_file["/aids"].keys())
    hdf5_file.close()
    logger.info(f"aids with no targets_gene_id: {len(list(set(no_targets)))}")
    logger.info(f"aids with targets_gene_id: {len(list(set(yes_targets)))}")
    logger.info(f"aids with protein sequence: {len(list(set(yes_seq)))}")
    logger.info(f"aids with no uniprot_pdb_alphafold: {len(list(set(no_uniprot)))}")
    logger.info(f"aids with uniprot_pdb_alphafold: {len(list(set(yes_uniprot)))}")
    logger.info(f"aids: {len(list(set(aids)))}")
    logger.info(f"cids: {len(total_cids)}")
    logger.info(f"aids with cids: {len(list(set(yes_cids)))}")
    logger.info(f"aids with no cids: {len(list(no_cids))}")
    logger.info(f"cids with smiles: {len(list(set(yes_smiles)))}")
    logger.info(f"cids with no smiles: {len(list(set(no_smiles)))}")


def complete_hdf5(hdf5_file_path: str) -> None:
    r"""Complete the HDF5 file.

    Args:
        hdf5_file_path : Path to the HDF5 file.
    """
    try:
        # Open HDF5 file
        with h5py.File(hdf5_file_path, "r+") as hdf5_file:
            for aid in hdf5_file["/aids"]:
                for cids in hdf5_file[f"/aids/{aid}/cids"]:
                    if "smiles" not in hdf5_file[f"/aids/{aid}/cids/{cids}"]:
                        logger.info(f"Processing CID {cids} for AID {aid}.")
                        pubchem = PubChemAPI()
                        data = pubchem.get_cids_info(int(cids), ["SMILES", "ExactMass"])
                        smiles = data["SMILES"]
                        logger.info(f"Smiles: {smiles}")
                        mass = data["ExactMass"]
                        logger.info(f"Exact mass: {mass}")
                        # Save smiles to the HDF5 file
                        hdf5_file[f"/aids/{aid}/cids/{cids}"].create_dataset(
                            "smiles",
                            data=smiles,
                            dtype=h5py.special_dtype(vlen=str),
                            maxshape=(None,),
                            shape=(1,),
                        )
                        hdf5_file[f"/aids/{aid}/cids/{cids}"].create_dataset(
                            "exact_mass",
                            data=mass,
                            dtype="float64",
                            maxshape=(None,),
                            shape=(1,),
                        )
                    else:
                        logger.info(f"Skipping CID {cids} for AID {aid}.")

                    # add a timeout to avoid being blocked by the server
                    time.sleep(5)
        logger.info("HDF5 file completed.")
    except requests.RequestException as e:
        logger.error(f"Network error when accessing PubChem API: {e}")
    except (IOError, OSError) as e:
        logger.error(f"File operation or system error occurred: {e}")
    except KeyError as e:
        logger.error(f"Key error (possibly missing data from PubChemAPI): {e}")
    except ValueError as e:
        logger.error(f"Value error (possibly from PubChemAPI): {e}")
    except Exception as e:
        logger.error(f"Unexpected error occurred: {e}")
