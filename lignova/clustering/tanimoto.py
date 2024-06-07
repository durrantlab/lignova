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

    def __init__(self, similarity_threshold: float = 0.2):
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

    def cal_distance(
        self, tanimoto_score: list[list[float]] | list[float]
    ) -> List[float]:
        r"""Calculate the distance between two fingerprints.
        parameters:
        ----------
            tanimoto_score: list[list[float]]
                List of Tanimoto similarity scores.
        returns:
        ----------
            float: Distance.
        """
        # knowing that the tanimoto_score is a list of list of similarity scores
        # we will calculate the distance between the fingerprints using the formula 1 - score
        distance = []
        if isinstance(tanimoto_score[0], float):
            for i in range(len(tanimoto_score)):
                distance.append(1 - tanimoto_score[i])
            return distance
        for i in range(len(tanimoto_score)):
            tmp = []
            for j in range(len(tanimoto_score[i])):
                tmp.append(1 - tanimoto_score[i][j])
            distance.extend(tmp)
        logger.debug(f"Distances: {distance}")
        return distance

    def cluster_tanimoto(
        self, tanimoto_score: list[list[float]], smiles: list[str]
    ) -> List[List[str]]:
        r"""Cluster the SMILES based on Tanimoto similarity using butina algorithm.
        parameters:
        ----------
            tanimoto_score: Dict[str, float]
                Dictionary of SMILES or compound identifiers and their Tanimoto similarity scores.
            smiles: List[str]
                List of SMILES or compound identifiers.
        returns:
        ----------
            List[List[str]]: Clustered .
        """
        clusters = []
        # calculate the distance between the fingerprints using cal_distance
        distance = self.cal_distance(tanimoto_score)
        clusters = Butina.ClusterData(
            distance, len(smiles), self.similarity_threshold, isDistData=True
        )
        logger.info(f"Number of clusters: {len(clusters)}")
        logger.debug(f"Clusters: {clusters}")
        # convert the clusters to smiles
        clusters = [[smiles[i] for i in cluster] for cluster in clusters]
        logger.info(f"Clusters: {clusters}")
        return clusters
