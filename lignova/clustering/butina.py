# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Implementation of the Butina clustering algorithm."""

from dataclasses import dataclass
from typing import ClassVar

from loguru import logger
from rdkit.ML.Cluster import Butina

from lignova.clustering import Clusterer, ClusterMethod, ClusterParams, ClusterResult

from .tanimoto import TanimotoSimilarities


@dataclass(frozen=True, slots=True)
class ButinaParams(ClusterParams):
    """Parameters for the Butina clustering algorithm."""

    method: ClassVar[ClusterMethod] = ClusterMethod.BUTINA
    """The clustering method used which is Butina."""

    similarity_cutoff: float = 0.7
    """The tanimoto similarity cutoff for clustering where higher values indicate more similar compounds. Must be between 0 and 1."""

    def __post_init__(self) -> None:
        if not (0.0 <= self.similarity_cutoff <= 1.0):
            raise ValueError(
                f"Similarity cutoff must be a float between 0 and 1, got {self.similarity_cutoff}"
            )


class ButinaClustering(Clusterer[ButinaParams]):
    """Butina clustering algorithm."""

    def __init__(self, params: ButinaParams):
        """Initialize the Butina clustering algorithm."""
        super().__init__(params)

    def cluster(self, sims: TanimotoSimilarities) -> ClusterResult:
        r"""Cluster the pass's compounds with RDKit Butina.

        Args:
            sims: Output of `compute_pairwise` over the compounds to cluster.

        Returns:
            A `ClusterResult` keyed by `sims.ids`.
        """
        ids = sims.ids
        n = sims.n

        if n == 0:
            return ClusterResult(labels={}, representatives={}, params=self.params)
        if n == 1:
            return ClusterResult(
                labels={ids[0]: 0}, representatives={0: ids[0]}, params=self.params
            )

        clusters = Butina.ClusterData(
            sims.condensed_distances,
            n,
            1.0 - self.params.similarity_cutoff,
            isDistData=True,
        )
        logger.info(
            f"Butina produced {len(clusters)} clusters from {n} compounds "
            f"at cutoff {self.params.similarity_cutoff}"
        )

        labels: dict[str, int] = {}
        representatives: dict[int, str] = {}
        for cid, members in enumerate(clusters):
            representatives[cid] = ids[members[0]]
            for idx in members:
                labels[ids[idx]] = cid

        return ClusterResult(
            labels=labels, representatives=representatives, params=self.params
        )
