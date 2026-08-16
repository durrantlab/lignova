# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Test PubChem API."""

import numpy as np
import pytest

from lignova.APIs.pubchem.client import PubChemAPI


@pytest.mark.asyncio
async def test_active_inactive_cids():
    r"""active_cids / inactive_cids match the dedicated /cids endpoint (AID 1000)."""
    async with PubChemAPI() as pubchem:
        assay = await pubchem.get_assay(1000)

    assert assay is not None
    active = assay.active_cids
    inactive = assay.inactive_cids

    assert len(active) == 36
    assert len(inactive) == 21
    assert 16749973 in active
    assert 16749973 not in inactive
    assert 730211 in inactive
    assert 730211 not in active


@pytest.mark.asyncio
async def test_get_cids_info():
    r"""Retrieve SMILES and ExactMass for a compound (aspirin, CID 2244)."""
    async with PubChemAPI() as pubchem:
        info = await pubchem._get_cids_info(2244, ["SMILES", "ExactMass"])

    assert info["SMILES"] == "CC(=O)OC1=CC=CC=C1C(=O)O"
    assert np.isclose(float(info["ExactMass"]), 180.04225873, rtol=0.01)


@pytest.mark.asyncio
async def test_binding_affinity():
    r"""Binding affinity for two CIDs in AID 1057958 ."""
    async with PubChemAPI() as pubchem:
        assay = await pubchem.get_assay(1057958)

    assert assay is not None
    affinity = assay.binding_affinity(cids=[135566761, 135566762])

    value_1, name_1, outcome_1 = affinity[135566761][0]
    value_2, name_2, outcome_2 = affinity[135566762][0]

    assert value_1 == 2.12
    assert value_2 == 0.88
    assert name_1 == "IC50"
    assert name_2 == "IC50"
    assert outcome_1 == "Active"
    assert outcome_2 == "Active"


@pytest.mark.asyncio
async def test_binding_affinity_outcome_filter():
    r"""outcome='active' keeps the Active IC50s; a bogus outcome yields nothing."""
    async with PubChemAPI() as pubchem:
        assay = await pubchem.get_assay(1057958)

    assert assay is not None
    active_only = assay.binding_affinity(cids=[135566761, 135566762], outcome="active")
    assert active_only[135566761][0][0] == 2.12

    # negative test for outcome filter
    inactive_only = assay.binding_affinity(
        cids=[135566761, 135566762], outcome="inactive"
    )
    assert inactive_only == {}


@pytest.mark.asyncio
async def test_pubmed_id():
    r"""PubMed ID is populated for a real assay and absent for a non-assay id."""
    async with PubChemAPI() as pubchem:
        assay = await pubchem.get_assay(1057958)
        non_assay = await pubchem.get_assay(2244)

    assert assay is not None
    assert assay.pubmed_id == 24188023
    assert non_assay is not None
    assert non_assay.pubmed_id is None


@pytest.mark.asyncio
async def test_enrichment_values():
    r"""Enriched record carries the exact expected CompoundProperties for AID 1000."""
    cid = 730195

    async with PubChemAPI() as pubchem:
        assay = await pubchem.get_assay(1000)
        assert assay is not None
        assert cid in assay.unique_cids
        enriched = await pubchem.enrich_cid_properties(assay)

    record = next(r for r in enriched.records if r.cid == cid)
    assert record.properties is not None
    p = record.properties

    assert p.smiles == "CC1=CC2=C(C=C1)N=C3CCCN3C2=O"
    assert p.inchikey == "WGMDMAZNZCEADV-UHFFFAOYSA-N"
    assert p.charge == 0
    assert np.isclose(p.xlogp, 1.3, rtol=0.01)
    assert np.isclose(p.exact_mass, 200.094963011, rtol=0.01)
    assert np.isclose(p.molecular_weight, 200.24, rtol=0.01)
    assert np.isclose(p.complexity, 324.0, rtol=0.01)
    assert np.isclose(p.tpsa, 32.7, rtol=0.01)
