# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Test Tanimoto similarity calculations and butina clustering."""

import numpy as np
import pytest
from rdkit.DataStructs.cDataStructs import ExplicitBitVect

from lignova.activity.models import CompoundActivity
from lignova.clustering import (
    ButinaClustering,
    ButinaParams,
    CliffParams,
    FeaturizeResult,
    MorganFeaturizer,
    TanimotoSimilarities,
    compute_pairwise,
    find_activity_cliffs,
    high_confidence_cliffs,
)

BENZENE = "C1=CC=CC=C1"
smiles1 = "CC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5"
smiles2 = "CC[C@@H](C(=O)N[C@@H](C1CCCCC1)C(=O)N2C[C@H]3CCCN3C[C@H]2C(=O)N[C@@H]4CCOC5=CC=CC=C45)NC"
ACETIC = "CC(=O)O"

CLUSTER_SMILES = {
    "s0": "COC1=CC=C(C=C1)CNC(=O)CCCNS(=O)(=O)C2=CC=CC3=NON=C32",
    "s1": "CC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5",
    "s2": "CC(C)(C)OC(=O)N1CCC(CC1)C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5",
    "s3": "C1CN(CCC1C2=NC3=CC=CC=C3N=C2OC4CN(C4)C5=NC6=CC=CC=C6C=C5)C(=O)O",
}

FEATURIZER = MorganFeaturizer()


def _fingerprints(smilesmiles2y_id: dict[str, str]) -> dict[str, ExplicitBitVect]:
    """Featurize and assert nothing was skipped."""
    res = FEATURIZER.featurize(smilesmiles2y_id)
    assert res.skipped == {}
    return res.items


def _sim_lookup(sims: TanimotoSimilarities) -> dict[frozenset, float]:
    return {frozenset((a, b)): s for a, b, s in sims.edges}


def test_featurize():
    """Test that the Morgan featurizer produces fingerprints of the expected length and type."""
    res = FEATURIZER.featurize({"benzene": BENZENE})
    assert isinstance(res, FeaturizeResult)
    fp = res.items["benzene"]
    assert isinstance(fp, ExplicitBitVect)
    assert len(fp) == 2048
    assert res.skipped == {}
    res_bad = FEATURIZER.featurize({"good": BENZENE, "bad": "not_a_smiles"})
    assert set(res_bad.items) == {"good"}
    assert "bad" in res_bad.skipped
    assert res_bad.n_compounds == 2


def test_pairwise_similarity_values():
    """Test that the pairwise similarity calculations."""
    items = _fingerprints({"A": smiles1, "B": smiles2, "C": ACETIC})
    sims = compute_pairwise(items, min_sim=0.0)
    assert isinstance(sims, TanimotoSimilarities)
    assert len(sims.condensed_distances) == 3
    lut = _sim_lookup(sims)
    assert np.isclose(lut[frozenset(("A", "B"))], 0.17647058823529413)
    assert np.isclose(lut[frozenset(("A", "C"))], 0.09433962264150944)
    items_same = _fingerprints({"x": smiles1, "y": smiles1})
    sims_1 = compute_pairwise(items_same, min_sim=0.0)
    ((_, _, s),) = sims_1.edges
    assert s == 1.0
    assert sims_1.condensed_distances[0] == 0.0


def test_condensed_distancesmiles1re_one_minus_similarity():
    """Test that the distances calculations"""
    items = _fingerprints({"A": smiles1, "B": smiles2, "C": ACETIC})
    sims = compute_pairwise(items, min_sim=0.0)
    lut = _sim_lookup(sims)
    assert np.isclose(1 - lut[frozenset(("A", "B"))], 0.8235294117647058)
    assert np.isclose(1 - lut[frozenset(("A", "C"))], 0.9056603773584906)


def test_min_sim_filters_edgesmiles2ut_not_condensed():
    """Test that the min_sim parameter filters edges but does not affect the distance matrix."""
    items = _fingerprints({"A": smiles1, "B": smiles2, "C": ACETIC})
    sims = compute_pairwise(items, min_sim=0.15)
    lut = _sim_lookup(sims)
    assert frozenset(("A", "B")) in lut
    assert frozenset(("A", "C")) not in lut
    assert len(sims.condensed_distances) == 3


def test_butina_clustering():
    """Test Butina clustering produces the expected partition and representatives."""
    items = _fingerprints(CLUSTER_SMILES)
    sims = compute_pairwise(items, min_sim=0.0)
    result = ButinaClustering(ButinaParams(similarity_cutoff=0.5)).cluster(sims)

    assert result.n_clusters == 2
    clusters = result.clusters()
    assert sorted(len(v) for v in clusters.values()) == [1, 3]

    singleton = next(v for v in clusters.values() if len(v) == 1)
    assert singleton == ["s0"]

    big_cid = next(cid for cid, v in clusters.items() if len(v) == 3)
    assert set(clusters[big_cid]) == {"s1", "s2", "s3"}
    assert result.representatives[big_cid] == "s3"


def test_butina_params_validation():
    """Test that ButinaParams raises a ValueError for invalid similarity_cutoff values."""
    with pytest.raises(ValueError):
        ButinaParams(similarity_cutoff=1.5)


def test_empty_and_single_compound():
    """Test that Butina clustering handles empty and single-compound cases."""
    empty = ButinaClustering(ButinaParams()).cluster(compute_pairwise({}, min_sim=0.5))
    assert empty.n_clusters == 0

    one = _fingerprints({"only": BENZENE})
    single = ButinaClustering(ButinaParams()).cluster(
        compute_pairwise(one, min_sim=0.5)
    )
    assert single.clusters() == {0: ["only"]}


def test_cliffs_read_edges_not_clusters():
    """Test that find_activity_cliffs reads the edges of the similarity graph"""
    items = _fingerprints({"A": smiles1, "B": smiles2, "C": ACETIC})
    sims = compute_pairwise(items, min_sim=0.15)
    activity = {
        "A": CompoundActivity(pActivity=9.0, winning_type="Ki"),
        "B": CompoundActivity(pActivity=5.0, winning_type="Ki"),
        "C": CompoundActivity(pActivity=6.0, winning_type="Ki"),
    }
    result = find_activity_cliffs(
        sims, activity, CliffParams(min_delta=2.0, min_similarity=0.15)
    )
    pairs = {frozenset((c.id_a, c.id_b)) for c in result.cliffs}
    assert frozenset(("A", "B")) in pairs


def test_cliff_flag_behavior():
    """Test that find_activity_cliffs flags low-quality and cross-type cliffs."""
    items = _fingerprints({"A": smiles1, "B": smiles2})
    sims = compute_pairwise(items, min_sim=0.15)
    activity = {
        "A": CompoundActivity(pActivity=9.0, winning_type="Ki"),
        "B": CompoundActivity(
            pActivity=5.0, winning_type="Ki", passes_quality_gate=False
        ),
    }
    result = find_activity_cliffs(
        sims, activity, CliffParams(min_delta=2.0, min_similarity=0.15)
    )

    pair = {frozenset((c.id_a, c.id_b)) for c in result.cliffs}
    assert frozenset(("A", "B")) in pair
    ab = next(c for c in result.cliffs if {c.id_a, c.id_b} == {"A", "B"})
    assert ab.is_low_quality is True

    clean = high_confidence_cliffs(result.cliffs)
    assert all(not c.is_low_quality for c in clean)
    assert not any({c.id_a, c.id_b} == {"A", "B"} for c in clean)
