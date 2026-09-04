# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Test activity-cliff scoring: metric branches, severity labelling, and restriction."""

import numpy as np
import pytest

from lignova.activity.models import CompoundActivity
from lignova.clustering import (
    CliffParams,
    CliffSeverity,
    MorganFeaturizer,
    SaliUndefined,
    SeverityMetric,
    TanimotoSimilarities,
    compute_pairwise,
    find_activity_cliffs,
    label_cliff_severity,
    same_assay,
    same_assay_cliffs,
)

smiles1 = "CC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5"
smiles2 = "CC[C@@H](C(=O)N[C@@H](C1CCCCC1)C(=O)N2C[C@H]3CCCN3C[C@H]2C(=O)N[C@@H]4CCOC5=CC=CC=C45)NC"
ACETIC = "CC(=O)O"

FEATURIZER = MorganFeaturizer()


def _fingerprints(by_id: dict[str, str]) -> dict:
    """Featurize and assert nothing was skipped."""
    res = FEATURIZER.featurize(by_id)
    assert res.skipped == {}
    return res.items


def _sims(edges: list[tuple[str, str, float]]) -> TanimotoSimilarities:
    """Build a TanimotoSimilarities from an explicit edge list."""
    ids: list[str] = []
    seen: set[str] = set()
    for a, b, _ in edges:
        for x in (a, b):
            if x not in seen:
                seen.add(x)
                ids.append(x)
    n = len(ids)
    condensed = [0.0] * (n * (n - 1) // 2)
    return TanimotoSimilarities(
        ids=tuple(ids), condensed_distances=condensed, edges=list(edges), min_sim=0.0
    )


def _laddered_activity(
    n: int,
) -> tuple[TanimotoSimilarities, dict[str, CompoundActivity]]:
    """n independent pairs Ai-Bi (similarity 0.8) with increasing pActivity gaps."""
    edges = [(f"A{i}", f"B{i}", 0.8) for i in range(n)]
    activity: dict[str, CompoundActivity] = {}
    for i in range(n):
        activity[f"A{i}"] = CompoundActivity(pActivity=10.0 + i, winning_type="Ki")
        activity[f"B{i}"] = CompoundActivity(pActivity=0.0, winning_type="Ki")
    return _sims(edges), activity


def test_landscape_index():
    """Test that TS_SALI is defined for every pair, including Tm == 1 but SALI is undefined."""
    items = _fingerprints({"A": smiles1, "E": smiles1})
    sims = compute_pairwise(items, min_sim=0.0)
    activity = {
        "A": CompoundActivity(pActivity=9.0, winning_type="Ki"),
        "E": CompoundActivity(pActivity=4.0, winning_type="Ki"),
    }
    result = find_activity_cliffs(sims, activity, CliffParams(min_delta=2.0))
    (c,) = result.cliffs
    assert np.isclose(c.landscape_index, 25.0)
    assert result.n_undefined == 0
    result = find_activity_cliffs(
        sims, activity, CliffParams(min_delta=2.0, metric=SeverityMetric.SALI)
    )
    (c,) = result.cliffs
    assert c.landscape_index is None
    assert result.n_undefined == 1


def test_sali_value():
    """Test that SALI is calculated correctly when it is defined."""
    sims = _sims([("A", "B", 0.5)])
    activity = {
        "A": CompoundActivity(pActivity=9.0, winning_type="Ki"),
        "B": CompoundActivity(pActivity=6.0, winning_type="Ki"),
    }
    result = find_activity_cliffs(
        sims,
        activity,
        CliffParams(min_delta=2.0, metric=SeverityMetric.SALI, min_similarity=0.0),
    )
    (c,) = result.cliffs
    assert np.isclose(c.landscape_index, 6.0)


def test_severity_label():
    """Test that cliff severity is assigned based on percentiles."""
    sims, activity = _laddered_activity(100)
    result = label_cliff_severity(
        find_activity_cliffs(sims, activity, CliffParams(min_delta=2.0))
    )
    counts = result.severity_counts()
    assert counts[CliffSeverity.EXTREME] == 1
    assert counts[CliffSeverity.STRONG] == 4
    assert counts[CliffSeverity.MODERATE] == 95
    extreme = next(
        c for c in result.cliffs if c.severity_label is CliffSeverity.EXTREME
    )
    assert {extreme.id_a, extreme.id_b} == {"A99", "B99"}


def test_restrictions():
    """Test that severity labelling can be restricted to cliffs matching a condition."""
    sims = _sims([("A", "B", 0.8), ("A", "C", 0.8)])
    activity = {
        "A": CompoundActivity(9.0, "Ki", assay_ids=frozenset({"aid1"})),
        "B": CompoundActivity(4.0, "Ki", assay_ids=frozenset({"aid1"})),
        "C": CompoundActivity(4.0, "Ki", assay_ids=frozenset({"aid2"})),
    }
    result = find_activity_cliffs(sims, activity, CliffParams(min_delta=2.0))
    assert len(same_assay_cliffs(result.cliffs)) == 1
    result = label_cliff_severity(result, restrict=same_assay)
    labelled = [c for c in result.cliffs if c.severity_label is not None]
    assert len(labelled) == 1
    assert {labelled[0].id_a, labelled[0].id_b} == {"A", "B"}


def test_sali_undefined_treatment():
    """Test that the sali_undefined parameter controls how undefined SALI values are treated for severity labelling."""
    sims = _sims([("A", "B", 0.5), ("A", "E", 1.0)])
    activity = {
        "A": CompoundActivity(9.0, "Ki"),
        "B": CompoundActivity(6.0, "Ki"),
        "E": CompoundActivity(4.0, "Ki"),
    }
    result = label_cliff_severity(
        find_activity_cliffs(
            sims,
            activity,
            CliffParams(
                min_delta=2.0,
                metric=SeverityMetric.SALI,
                sali_undefined=SaliUndefined.NEXT_LARGEST,
                min_similarity=0.0,
            ),
        )
    )
    ae = next(c for c in result.cliffs if {c.id_a, c.id_b} == {"A", "E"})
    assert ae.severity_label is not None
    result = label_cliff_severity(
        find_activity_cliffs(
            sims,
            activity,
            CliffParams(
                min_delta=2.0,
                metric=SeverityMetric.SALI,
                sali_undefined=SaliUndefined.EXCLUDE,
                min_similarity=0.0,
            ),
        )
    )
    ae = next(c for c in result.cliffs if {c.id_a, c.id_b} == {"A", "E"})
    ab = next(c for c in result.cliffs if {c.id_a, c.id_b} == {"A", "B"})
    assert ae.severity_label is None
    assert ab.severity_label is not None
    result = label_cliff_severity(
        find_activity_cliffs(
            sims,
            activity,
            CliffParams(
                min_delta=2.0,
                metric=SeverityMetric.SALI,
                sali_undefined=SaliUndefined.MAX_SEVERITY,
            ),
        )
    )
    ae = next(c for c in result.cliffs if {c.id_a, c.id_b} == {"A", "E"})
    assert ae.severity_label is CliffSeverity.EXTREME
