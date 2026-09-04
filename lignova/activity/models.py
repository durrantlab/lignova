# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""The shared activity contract."""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class CompoundActivity:
    """Per-compound activity handed to cliff detection."""

    pActivity: float
    """Aggregated log-scale activity calculated from the median of the winning type's values."""

    winning_type: str
    """Winning assay type behind pActivity; enables the cross-type check."""

    passes_quality_gate: bool = True
    """False if the winning-type values disagreed by more than ~1 log unit. These cliffs are then
    flagged as low quality rather than excluded."""

    is_active: bool = True
    """False if the compound is considered inactive in its winning assay type against its target."""

    assay_ids: frozenset[str] = frozenset()
    """The assay ids (PubChem AID) behind the winning-type median. This enables the within-assay check for cliffs that share an assay id."""
