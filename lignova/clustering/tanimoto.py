r" Implementation of the tanimoto clustering algorithm."

from typing import List

import os

import rdkit


def calc_tanimoto_similarity(
    mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol
) -> float:
    r"""Calculate the Tanimoto similarity between two molecules.

    Parameters:
    ----------
        mol1: rdkit.Chem.rdchem.Mol
            First molecule.
        mol2: rdkit.Chem.rdchem.Mol
            Second molecule.

    Returns:
    ----------
        float: Tanimoto similarity.
    """
    fp1 = rdkit.Chem.RDKFingerprint(mol1)
    fp2 = rdkit.Chem.RDKFingerprint(mol2)
    return rdkit.DataStructs.FingerprintSimilarity(fp1, fp2)


def cluster_smiles(
    smiles_list: List[str], similarity_threshold: float
) -> List[List[str]]:
    r"""Cluster the SMILES based on Tanimoto similarity.

    Parameters:
    ----------
        smiles_list: List[str]
            List of SMILES strings.
        similarity_threshold: float
            Similarity threshold for clustering.

    Returns:
    ----------
        List[List[str]]: Clustered SMILES.
    """
    clusters = []
    for smile in smiles_list:
        mol = rdkit.Chem.MolFromSmiles(smile)
        if mol is None:
            continue
        added_to_cluster = False
        for cluster in clusters:
            for cluster_smile in cluster:
                cluster_mol = rdkit.Chem.MolFromSmiles(cluster_smile)
                if cluster_mol is None:
                    continue
                similarity = tanimoto_similarity(mol, cluster_mol)
                if similarity >= similarity_threshold:
                    cluster.append(smile)
                    added_to_cluster = True
                    break
            if added_to_cluster:
                break
        if not added_to_cluster:
            clusters.append([smile])
    return clusters
