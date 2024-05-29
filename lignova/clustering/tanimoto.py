r" Implementation of the tanimoto clustering algorithm."

from typing import Dict, List

import os
import subprocess

import rdkit
from loguru import logger
from rdkit import Chem, DataStructs
from rdkit.ML.Cluster import Butina


class TanimotoClustering:
    r"""Tanimoto clustering algorithm."""

    def __init__(self, similarity_threshold: float):
        r"""Initialize the Tanimoto clustering algorithm.
        parameters:
        ----------
            similarity_threshold: float
                Similarity threshold for clustering.
        """
        self.similarity_threshold = similarity_threshold

    def get_morgan_fingerprint(self, smiles: str, radius: int = 2):
        r"""Get the fingerprint of the molecule.
        parameters:
        ----------
            smiles: str
                Molecule.
            radius: int
        returns:
        ----------
            rdkit.DataStructs.cDataStructs.ExplicitBitVect: Fingerprint.
        """
        mol = rdkit.Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES: {smiles}")
            return None
        else:
            return rdkit.Chem.AllChem.GetMorganFingerprintAsBitVect(mol, radius)

    def tanimoto_similarity(
        self,
        mol1: DataStructs.cDataStructs.ExplicitBitVect
        | List[DataStructs.cDataStructs.ExplicitBitVect],
        mol2: DataStructs.cDataStructs.ExplicitBitVect,
    ):
        r"""Calculate the Tanimoto similarity between two fingerprints.
        parameters:
        ----------
            mol1: rdkit.DataStructs.cDataStructs.ExplicitBitVect | List[rdkit.DataStructs.cDataStructs.ExplicitBitVect]
                First molecule fingerprint or list of fingerprints.
            mol2: rdkit.DataStructs.cDataStructs.ExplicitBitVect | None
                Second molecule fingerprint or None.
        returns:
        ----------
            float: Tanimoto similarity.
        """
        if isinstance(mol1, list):
            return DataStructs.BulkTanimotoSimilarity(mol2, mol1)
        else:
            return DataStructs.TanimotoSimilarity(mol1, mol2)

    def cluster_tanimoto(self, tanimoto_score: Dict[str, float]) -> List[List[str]]:
        r"""Cluster the SMILES based on Tanimoto similarity.
        parameters:
        ----------
            tanimoto_score: Dict[str, float]
                Dictionary of SMILES or compound identifiers and their Tanimoto similarity scores.
        returns:
        ----------
            List[List[str]]: Clustered .
        """
        clusters = []
        for key, value in tanimoto_score.items():
            if not clusters:
                clusters.append([key])
            else:
                found = False
                for cluster in clusters:
                    for member in cluster:
                        if (
                            self.tanimoto_similarity(
                                self.get_morgan_fingerprint(key),
                                self.get_morgan_fingerprint(member),
                            )
                            > self.similarity_threshold
                        ):
                            cluster.append(key)
                            found = True
                            break
                    if found:
                        break
                if not found:
                    clusters.append([key])
