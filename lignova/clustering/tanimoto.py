# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Compute pairwise Tanimoto similarities and distances for a set of fingerprints."""

from dataclasses import dataclass

from loguru import logger
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import ExplicitBitVect


@dataclass(frozen=True, slots=True)
class TanimotoSimilarities:
    """Dataclass holding the results of one similarity pass over a set of fingerprints."""

    ids: tuple[str, ...]
    """Frozen id order for the compounds set in the similarity pass. Indices in `condensed_distances` refer to this."""

    condensed_distances: list[float]
    """The lower-triangle of (1 - Tm) matrix where the compound at index i has distances to compounds at indices j in 0..i-1. Length is N(N-1)/2 where N is the number of compounds."""

    edges: list[tuple[str, str, float]]
    """List of tuples containing the compound ids and their Tanimoto similarity for every pair with Tm >= floor."""

    min_sim: float
    """The threshold applied to produce `edges`. Only pairs with Tm >= floor are included in `edges`."""

    def __post_init__(self) -> None:
        if len(self.condensed_distances) != self.n * (self.n - 1) // 2:
            raise ValueError(
                f"condensed_distances length {len(self.condensed_distances)} "
                f"does not match expected {self.n * (self.n - 1) // 2} for n={self.n}"
            )

        check_floor(self.min_sim)

    @property
    def n(self) -> int:
        """Number of compounds in the similarity pass."""
        return len(self.ids)


def check_floor(min_sim: float) -> None:
    """Check that the floor value is valid.

    Args:
        min_sim: The floor value to check.
    """
    if not 0.0 <= min_sim <= 1.0:
        raise ValueError(f"floor must be in [0, 1], got {min_sim}")


def compute_pairwise(
    items: dict[str, ExplicitBitVect], min_sim: float
) -> TanimotoSimilarities:
    """One pass: build the condensed distance array AND the >= min_sim edge list.

    Args:
        items: a dictionary with the  compound ids as the keys and their fingerprints as the values.
        min_sim: Minimum Tanimoto similarity for an edge to be kept

    Returns:
        A `TanimotoSimilarities` object containing the condensed distance array and the edge list.
    """
    check_floor(min_sim)

    ids = list(items)
    fps = [items[i] for i in ids]
    n = len(ids)

    condensed: list[float] = []
    edges: list[tuple[str, str, float]] = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        for j, s in enumerate(sims):
            condensed.append(1.0 - s)
            if s >= min_sim:
                edges.append((ids[i], ids[j], s))
    logger.info(
        "Computed pairwise similarities for {n} compounds producing {condensed} distances and {edges} edges at min similarity of {floor}",
        n=n,
        condensed=len(condensed),
        edges=len(edges),
        floor=min_sim,
    )
    return TanimotoSimilarities(
        ids=tuple(ids), condensed_distances=condensed, edges=edges, min_sim=min_sim
    )
