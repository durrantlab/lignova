# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Base class for clustering algorithms"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import StrEnum
from typing import ClassVar


class ClusterMethod(StrEnum):
    """An enum class to identify the cluster algorithm that produced a result."""

    BUTINA = "butina"
    LEADER = "leader"
    MMSEQS2 = "mmseqs2"


@dataclass(frozen=True, slots=True)
class ClusterParams(ABC):
    """Abstract base for a method's configuration."""

    method: ClassVar[ClusterMethod]


@dataclass(frozen=True, slots=True)
class ClusterResult:
    """The standardized output every clusterer returns."""

    labels: dict[str, int]
    """dictionary with the keys being the items id (e.g. cids) and the values is their cluster id."""

    representatives: dict[int, str]
    """ dictionary with the keys being the cluster id and the value is the id of the cluster representative id"""

    params: ClusterParams
    """The params object holding the metadata the clustering algorithm ran with."""

    def __post_init__(self) -> None:
        missing = set(self.representatives) - set(self.labels.values())
        """Ensure that cluster ids saved have a valid representative."""
        if missing:
            raise ValueError(
                f"Representatives with invalid cluster ids: {sorted(missing)}"
            )

    @property
    def method(self) -> ClusterMethod:
        """The algorithm generated this result."""
        return self.params.method

    @property
    def n_clusters(self) -> int:
        """Number of distinct clusters."""
        return len(set(self.labels.values()))

    def clusters(self) -> dict[int, list[str]]:
        """Converts the labels into a dictionary of clusters with the cluster id as the key and the cluster members as the value."""
        out: dict[int, list[str]] = {}
        for item_id, cid in self.labels.items():
            out.setdefault(cid, []).append(item_id)
        return out


class Clusterer(ABC):
    """Abstract base class for all clustering algorithms.

    Args:
        params: The ClusterParams object containing the configuration for the clustering algorithm.
    """

    def __init__(self, params: ClusterParams) -> None:
        self.params = params

    @property
    def method(self) -> ClusterMethod:
        """The algorithm this clustering method implements."""
        return self.params.method

    @abstractmethod
    def cluster(self, items: dict[str, object]) -> ClusterResult:
        """
        Abstract cluster method that must be implemented by subclasses.

        Args:
            items: A dictionary where the keys are stable identifiers (e.g., cids) and the values are the items to be clustered.
        """
        raise NotImplementedError
