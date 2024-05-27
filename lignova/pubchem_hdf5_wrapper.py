r""" Implementation of the PubChem wrapper."""
import time

import h5py
from hdf5.pubchem import PubChemAPI
from loguru import logger


def process_hdf5_file(hdf5_file_path: str):
    r"""Process the HDF5 file.
    Parameters:
    ----------
        hdf5_file_path (str): Path to the HDF5 file.
    """

    # open hdf5 file
    hdf5_file = h5py.File(hdf5_file_path, "r")
    no_seq = []
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
        else:
            no_seq.append(aid)
            # print('aid '+str(aid)+' has no protein sequence branch')
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
    a = list(hdf5_file["/aids"].keys())
    hdf5_file.close()
    logger.info(f"aids with no targets_gene_id: {len(list(set(no_targets)))}")
    logger.info(f"aids with targets_gene_id: {len(list(set(yes_targets)))}")
    logger.info(f"aids with no protein sequence: {len(list(set(no_seq)))}")
    logger.info(f"aids with protein sequence: {len(list(set(yes_seq)))}")
    logger.info(f"aids with no uniprot_pdb_alphafold: {len(list(set(no_uniprot)))}")
    logger.info(f"aids with uniprot_pdb_alphafold: {len(list(set(yes_uniprot)))}")
    logger.info(f"aids: {len(list(set(a)))}")
    logger.info(f"cids: {len(total_cids)}")
    logger.info(f"aids with cids: {len(list(set(yes_cids)))}")
    logger.info(f"aids with no cids: {len(list(no_cids))}")
    logger.info(f"cids with smiles: {len(list(set(yes_smiles)))}")
    logger.info(f"cids with no smiles: {len(list(set(no_smiles)))}")


def complete_hdf5(hdf5_file_path: str):
    r"""Complete the HDF5 file.
    Parameters:
    ----------
        hdf5_file_path (str): Path to the HDF5 file.
    """
    try:
        # Open HDF5 file
        with h5py.File(hdf5_file_path, "r+") as hdf5_file:
            for aid in hdf5_file["/aids"]:
                for cids in hdf5_file[f"/aids/{aid}/cids"]:
                    if "smiles" not in hdf5_file[f"/aids/{aid}/cids/{cids}"]:
                        logger.info(f"Processing CID {cids} for AID {aid}.")
                        pubchem = PubChemAPI()
                        data = pubchem.get_smiles(
                            int(cids), ["IsomericSMILES", "ExactMass"]
                        )
                        smiles = data["smiles"]
                        logger.info(f"Smiles: {smiles}")
                        mass = data["mass"]
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
                        # add a timeout to avoid being blocked by the server
                        time.sleep(5)

    except Exception as e:
        logger.error(f"Error occurred: {e}")
    else:
        logger.info("HDF5 file completed.")

    """
    # open hdf5 file
    hdf5_file = h5py.File(hdf5_file_path, "r+")
    for aid in list(hdf5_file["/aids"].keys()):
        for cids in list(hdf5_file["/aids/" + str(aid) + "/cids"].keys()):
            if (
                "smiles"
                not in hdf5_file["/aids/" + str(aid) + "/cids/" + str(cids)].keys()
            ):
                # get smiles from pubchem
                logger.info(f"Processing CID {cids} for AID {aid}.")
                data = get_smiles(int(cids), ["IsomericSMILES", "ExactMass"])
                smiles = data["smiles"]
                logger.info(f"Smiles: {smiles}")
                mass = data["mass"]
                logger.info(f"Exact mass: {mass}")
                # save smiles to the hdf5 file
                hdf5_file["/aids/" + str(aid) + "/cids/" + str(cids)].create_dataset(
                    "smiles",
                    data=smiles,
                    dtype=h5py.special_dtype(vlen=str),
                    maxshape=(None,),
                    shape=(1,),
                )
                hdf5_file["/aids/" + str(aid) + "/cids/" + str(cids)].create_dataset(
                    "exact_mass",
                    data=mass,
                    dtype="float64",
                    maxshape=(None,),
                    shape=(1,),
                )
    hdf5_file.close()
    logger.info("HDF5 file completed.")
"""


if __name__ == "__main__":
    hdf5_file_path = "../PubChem_data_edited.hdf5"
    process_hdf5_file(hdf5_file_path)
    complete_hdf5(hdf5_file_path)
    # read the hdf5 file and loop through the aids
    # for each aid, if there are no cids delete the aid
    """
    try:
        # Open HDF5 file
        with h5py.File(hdf5_file_path, "r+") as hdf5_file:
            for aid in hdf5_file["/aids"]:
                if "cids" not in hdf5_file[f"/aids/{aid}"].keys():
                    del hdf5_file[f"/aids/{aid}"]
    except Exception as e:
        logger.error(f"Error occurred: {e}")

    """
