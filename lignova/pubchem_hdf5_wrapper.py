r""" Implementation of the PubChem wrapper."""

import h5py
from loguru import logger

# open hdf5 file
hdf5_file = h5py.File("../PubChem_data_copy.hdf5", "r")
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
            if "smiles" in hdf5_file["/aids/" + str(aid) + "/cids/" + str(i)].keys():
                #    print('cid '+str(i)+' has no smiles')
                yes_smiles.append(i)
            else:
                no_smiles.append(i)
    else:
        no_cids.append(aid)
a = list(hdf5_file["/aids"].keys())

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
hdf5_file.close()
