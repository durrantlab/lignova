r" Implementation of the tanimoto clustering algorithm."

from typing import List

from collections.abc import Iterable

from loguru import logger
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.ML.Cluster import Butina


class TanimotoClustering:
    r"""Tanimoto clustering algorithm."""

    def __init__(self):
        r"""Initialize the Tanimoto clustering algorithm."""

    # pylint: disable=c-extension-no-member,no-member
    def get_morgan_fingerprint(self, smiles: str, radius: int = 2) -> ExplicitBitVect:
        r"""Get the fingerprint of the molecule.

        Args:
            smiles: Molecule.
            radius: Radius of the mogan fingerprint. Default is 2.

        Returns:
            Morgan Fingerprint
        """
        morgan_fp = rdFingerprintGenerator.GetMorganGenerator(radius=radius)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES: {smiles}")
            return None

        return morgan_fp.GetFingerprint(mol)

    def tanimoto_similarity(
        self,
        mol1: ExplicitBitVect | Iterable[ExplicitBitVect],
        mol2: ExplicitBitVect,
    ) -> float | Iterable[float]:
        r"""Calculate the Tanimoto similarity between two fingerprints.

        Args:
            mol1 : First molecule fingerprint or List of fingerprints.
            mol2 : Second molecule fingerprint or None.

        Returns:
            Tanimoto similarity.
        """
        if isinstance(mol1, Iterable):
            return DataStructs.cDataStructs.BulkTanimotoSimilarity(mol2, mol1)

        return DataStructs.cDataStructs.TanimotoSimilarity(mol1, mol2)

    def cal_distance(
        self, tanimoto_score: Iterable[List[float]] | Iterable[float]
    ) -> List[float]:
        r"""Calculate the distance between two fingerprints.

        Args:
            tanimoto_score : List of Tanimoto similarity scores.

        Returns:
            Distance.
        """
        # knowing that the tanimoto_score is a List of List of similarity scores
        # we will calculate the distance between the fingerprints using the formula 1 - score
        distance = []
        if isinstance(tanimoto_score[0], float):
            distance = [1 - score for score in tanimoto_score]
        else:
            for item in tanimoto_score:
                distance.extend([1 - score for score in item])
        return distance

    def cluster_tanimoto(
        self,
        tanimoto_score: Iterable[float],
        smiles: Iterable[str],
        similarity_threshold: float,
    ) -> Iterable[List[str]]:
        r"""Cluster the SMILES based on Tanimoto similarity using butina algorithm.

        Args:
            tanimoto_score: Dictionary of SMILES or compound identifiers and
                their Tanimoto similarity scores.
            smiles: List of SMILES or compound identifiers.
            similarity_threshold: Similarity threshold.

        Returns:
            List of clusters, where each cluster is a List of SMILES or
                compoud identifiers depending on smiles input.
        """
        clusters = []
        # calculate the distance between the fingerprints using cal_distance
        distance = self.cal_distance(tanimoto_score)
        clusters = Butina.ClusterData(
            distance, len(smiles), 1 - similarity_threshold, isDistData=True
        )
        logger.info(f"Number of clusters: {len(clusters)}")
        # convert the clusters to smiles
        clusters = [[smiles[i] for i in cluster] for cluster in clusters]
        return clusters
