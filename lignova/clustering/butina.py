# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Implementation of the Butina clustering algorithm."""

from dataclasses import dataclass
from typing import ClassVar

import numpy as np
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

    sparse_above: int | None = None
    """The largest cluster size above which the sparse (edge-based) path is used instead of the dense RDKit matrix. If None, the dense path is always used."""

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

        use_sparse = (not sims.is_dense) or (
            self.params.sparse_above is not None and n > self.params.sparse_above
        )
        if use_sparse:
            return self._cluster_sparse(sims)
        else:
            return self._cluster_dense(sims)

    def _cluster_dense(self, sims: TanimotoSimilarities) -> ClusterResult:
        """Cluster the compounds using the dense condensed distance array i.e, the RDKit Butina implementation.

        Args:
            sims: Output of `compute_pairwise` over the compounds to cluster.
            Returns:
                A `ClusterResult` keyed by `sims.ids`.
        """
        ids = sims.ids
        n = sims.n
        if sims.condensed_distances is None:
            raise ValueError(
                "dense Butina needs condensed_distances; got edges-only sims."
            )
        clusters = Butina.ClusterData(
            sims.condensed_distances,
            n,
            1.0 - self.params.similarity_cutoff,
            isDistData=True,
        )
        logger.info(
            "Butina produced {c} clusters from {n} compounds at cutoff {k}",
            c=len(clusters),
            n=n,
            k=self.params.similarity_cutoff,
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

    def _cluster_sparse(self, sims: TanimotoSimilarities) -> ClusterResult:
        """Cluster the compounds using the sparse edge list i.e, the Butina implementation on the edge graph.

        Args:
            sims: Output of `compute_pairwise` over the compounds to cluster.
            Returns:
                A `ClusterResult` keyed by `sims.ids`.
        """
        cutoff = self.params.similarity_cutoff
        if sims.min_sim > cutoff:
            raise ValueError(
                f"sparse Butina needs edge floor <= cutoff, but floor={sims.min_sim} "
                f"> cutoff={cutoff}: the edge graph is missing neighbors in "
                f"[{cutoff}, {sims.min_sim}). Rebuild edges with a lower --floor."
            )
        ids = sims.ids
        n = sims.n
        idx = {cid: k for k, cid in enumerate(ids)}
        dist_thresh = 1.0 - cutoff

        # Build the neighbor graph from edges only
        neighbors: list[list[int]] = [[] for _ in range(n)]
        for a, b, s in sims.edges:
            if (1.0 - s) <= dist_thresh:
                ia, ib = idx[a], idx[b]
                neighbors[ia].append(ib)
                neighbors[ib].append(ia)

        # RDKit sorts (neighbor_count, index) in a descending order and breaks ties by higher index.
        degrees = np.fromiter(
            (len(neigh) for neigh in neighbors), dtype=np.int64, count=n
        )
        order = np.lexsort(
            (-np.arange(n), -degrees)
        )  # primary: -degrees, tiebreak: -index

        # This is a Greedy assignment identical to RDKit Butina (non-reordering).
        seen = np.zeros(n, dtype=bool)
        clusters: list[list[int]] = []
        for i in order:
            if seen[i]:
                continue
            members = [int(i)]
            for m in neighbors[i]:
                if not seen[m]:
                    members.append(m)
            seen[members] = True  # NumPy index assignment to the whole cluster
            clusters.append(members)

        logger.info(
            "Butina produced {c} clusters from {n} compounds at cutoff {k} (sparse)",
            c=len(clusters),
            n=n,
            k=cutoff,
        )

        labels: dict[str, int] = {}
        representatives: dict[int, str] = {}
        for cid, members in enumerate(clusters):
            representatives[cid] = ids[members[0]]
            for m in members:
                labels[ids[m]] = cid

        return ClusterResult(
            labels=labels, representatives=representatives, params=self.params
        )
