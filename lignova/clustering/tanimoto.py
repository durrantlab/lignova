r" Implementation of the tanimoto clustering algorithm."

from typing import List

from loguru import logger
from rdkit import Chem, DataStructs
from rdkit.ML.Cluster import Butina


class TanimotoClustering:
    r"""Tanimoto clustering algorithm."""

    def __init__(self):
        r"""Initialize the Tanimoto clustering algorithm."""

    # pylint: disable=c-extension-no-member,no-member
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
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES: {smiles}")
            return None

        return Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius)

    # pylint: disable=c-extension-no-member
    def tanimoto_similarity(
        self,
        mol1: (
            DataStructs.cDataStructs.ExplicitBitVect
            | List[DataStructs.cDataStructs.ExplicitBitVect]
        ),
        mol2: DataStructs.cDataStructs.ExplicitBitVect,
    ):
        r"""Calculate the Tanimoto similarity between two fingerprints.
        parameters:
        ----------
            mol1: rdkit.DataStructs.cDataStructs.ExplicitBitVect |
            List[rdkit.DataStructs.cDataStructs.ExplicitBitVect]
                First molecule fingerprint or list of fingerprints.
            mol2: rdkit.DataStructs.cDataStructs.ExplicitBitVect | None
                Second molecule fingerprint or None.
        returns:
        ----------
            float: Tanimoto similarity.
        """
        if isinstance(mol1, list):
            return DataStructs.cDataStructs.BulkTanimotoSimilarity(mol2, mol1)

        return DataStructs.cDataStructs.TanimotoSimilarity(mol1, mol2)

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
            distance = [1 - score for score in tanimoto_score]
        else:
            for item in tanimoto_score:
                distance.extend([1 - score for score in item])
        return distance

    def cluster_tanimoto(
        self,
        tanimoto_score: list[list[float]],
        smiles: list[str],
        similarity_threshold: float,
    ) -> List[List[str]]:
        r"""Cluster the SMILES based on Tanimoto similarity using butina algorithm.
        parameters:
        ----------
            tanimoto_score: Dict[str, float]
                Dictionary of SMILES or compound identifiers and their Tanimoto similarity scores.
            smiles: List[str]
                List of SMILES or compound identifiers.
            similarity_threshold: float
                Similarity threshold.
        returns:
        ----------
            List[List[str]]: Clustered .
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
