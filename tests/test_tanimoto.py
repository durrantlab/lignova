import os

import pytest
from loguru import logger
from rdkit import Chem, DataStructs

from lignova.clustering.tanimoto import TanimotoClustering
from lignova.structure.utils import (
    get_smiles,
    separate_protein_ligand,
    write_mda_universe,
)

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_pubchem = {
    "write_dir": "./tmp/clustering",
    "pdb_file": "./files/6oav/6oav.pdb",
}


def prep_dirs():
    os.makedirs(context_pubchem["write_dir"])


if not os.path.exists(context_pubchem["write_dir"]):
    prep_dirs()

TanimotoClustering = TanimotoClustering()


def test_get_morgan_fingerprint():
    smiles = "C1=CC=CC=C1"
    fingerprint = TanimotoClustering.get_morgan_fingerprint(smiles)
    assert isinstance(fingerprint, DataStructs.cDataStructs.ExplicitBitVect)
    assert len(fingerprint) == 2048


def test_tanimoto_similarity_single_mol():
    smiles1 = "CC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5"
    smiles2 = "CC[C@@H](C(=O)N[C@@H](C1CCCCC1)C(=O)N2C[C@H]3CCCN3C[C@H]2C(=O)N[C@@H]4CCOC5=CC=CC=C45)NC"
    fingerprint1 = TanimotoClustering.get_morgan_fingerprint(smiles1)
    fingerprint2 = TanimotoClustering.get_morgan_fingerprint(smiles2)
    similarity = TanimotoClustering.tanimoto_similarity(fingerprint1, fingerprint2)
    assert isinstance(similarity, float)
    assert similarity == 0.17647058823529413


def test_tanimoto_similarity_multiple_mols():
    smiles1 = [
        "CC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5",
        "CC[C@@H](C(=O)N[C@@H](C1CCCCC1)C(=O)N2C[C@H]3CCCN3C[C@H]2C(=O)N[C@@H]4CCOC5=CC=CC=C45)NC",
        "CC(=O)O",
    ]
    smiles2 = "CC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5"
    fingerprint1 = [
        TanimotoClustering.get_morgan_fingerprint(smile) for smile in smiles1
    ]
    fingerprint2 = TanimotoClustering.get_morgan_fingerprint(smiles2)
    similarity = TanimotoClustering.tanimoto_similarity(fingerprint1, fingerprint2)
    assert similarity[0] == 1.0
    assert similarity[1] == 0.17647058823529413
    assert similarity[2] == 0.09433962264150944


def test_cal_distance():
    smiles1 = [
        "CC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5",
        "CC[C@@H](C(=O)N[C@@H](C1CCCCC1)C(=O)N2C[C@H]3CCCN3C[C@H]2C(=O)N[C@@H]4CCOC5=CC=CC=C45)NC",
        "CC(=O)O",
    ]
    smiles2 = "CC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5"
    fingerprint1 = [
        TanimotoClustering.get_morgan_fingerprint(smile) for smile in smiles1
    ]
    fingerprint2 = TanimotoClustering.get_morgan_fingerprint(smiles2)
    similarity = TanimotoClustering.tanimoto_similarity(fingerprint1, fingerprint2)
    logger.info(similarity)
    distance = TanimotoClustering.cal_distance(similarity)
    logger.info(distance)
    assert distance[0] == 0.0
    assert distance[1] == 0.8235294117647058
    assert distance[2] == 0.9056603773584906


def test_cluster_tanimoto():
    smiles = [
        "COC1=CC=C(C=C1)CNC(=O)CCCNS(=O)(=O)C2=CC=CC3=NON=C32",
        "CC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5",
        "CC(C)(C)OC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5",
        "C1CN(CCC1C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5)C(=O)O",
    ]
    # Get fingerprints for each molecule and tanimoto similarity
    fingerprints = [
        TanimotoClustering.get_morgan_fingerprint(smile) for smile in smiles
    ]
    similarity = []
    # loop over the fingerprints to get the similarity for each molecule with the previous molecules
    for i in range(1, len(fingerprints)):
        similarity.append(
            TanimotoClustering.tanimoto_similarity(fingerprints[:i], fingerprints[i])
        )
    logger.info(f"Similarity {similarity}")

    # Cluster the molecules using cluster_tanimoto
    clusters = TanimotoClustering.cluster_tanimoto(similarity, smiles, 0.5)
    assert len(clusters) == 2
    assert len(clusters[0]) == 3
    assert len(clusters[1]) == 1
    assert clusters[0][0] == smiles[3]
    assert clusters[0][1] == smiles[1]
    assert clusters[0][2] == smiles[2]
    assert clusters[1][0] == smiles[0]
