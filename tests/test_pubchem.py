r"""Test PubChem API."""

import os

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from lignova.APIs.pubchem import PubChemAPI, PubChemBulk

# Column names exactly as they appear in bioactivities.tsv / the ingested parquet.
_DUMP_COLUMNS = [
    "AID",
    "CID",
    "Activity Outcome",
    "Activity Name",
    "Activity Value",
    "Activity Qualifier",
    "Activity Unit",
    "Protein Accession",
    "Gene ID",
    "PMID",
]

_TEST_AID = 999999


def _write_mini_parquet(path: str) -> None:
    """Write a small synthetic parquet using the real dump column names.

    Rows cover the interesting cases: a CID with both Active and Inactive rows
    (100), a blank-CID row (dropped by the model), a second CID (200), and an
    Unspecified row (300). All non-blank PMIDs agree (24188023), so the
    assay-level pubmed_id accessor returns that value.
    """
    rows = {
        "AID": [str(_TEST_AID)] * 5,
        "CID": ["100", "100", "200", "", "300"],
        "Activity Outcome": ["Active", "Inactive", "Active", "Active", "Unspecified"],
        "Activity Name": ["IC50", "IC50", "Ki", "IC50", "IC50"],
        "Activity Value": ["2.5", "40.0", "1.1", "5.0", "9.9"],
        "Activity Qualifier": ["=", "=", ">", "=", "="],
        "Activity Unit": ["uM", "uM", "nM", "uM", "uM"],
        "Protein Accession": ["P12345", "P12345", "", "P12345", "P12345"],
        "Gene ID": ["7157", "7157", "", "7157", "7157"],
        "PMID": ["24188023", "24188023", "", "", ""],
    }
    table = pa.table({c: pa.array(rows[c], type=pa.string()) for c in _DUMP_COLUMNS})
    pq.write_table(table, path)


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
        info = await pubchem._get_cids_info([2244], ["SMILES", "ExactMass"])
    props = info[2244]
    assert props["SMILES"] == "CC(=O)OC1=CC=CC=C1C(=O)O"
    assert np.isclose(float(props["ExactMass"]), 180.04225873, rtol=0.01)


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


@pytest.mark.asyncio
async def test_load_assay_pmid_per_row(tmp_path):
    r"""PMID rides per-record; blanks become None; assay-level agrees when unique."""
    pqpath = os.path.join(str(tmp_path), "mini.parquet")
    _write_mini_parquet(pqpath)

    assay = PubChemBulk().load_assay(_TEST_AID, pqpath)

    by_cid = {r.cid: r.pubmed_id for r in assay.records if r.cid is not None}
    assert by_cid[100] == 24188023
    assert by_cid[200] is None
    assert by_cid[300] is None
    assert assay.pubmed_id == 24188023
    assay = PubChemBulk().load_assay(111111, pqpath)
    assert assay.aid == 111111
    assert assay.records == []
    assert assay.pubmed_id is None


@pytest.mark.asyncio
async def test_load_assay_maps_columns(tmp_path):
    r"""load_assay maps dump columns onto AssayInfo records correctly."""
    pqpath = os.path.join(str(tmp_path), "mini.parquet")
    _write_mini_parquet(pqpath)

    assay = PubChemBulk().load_assay(_TEST_AID, pqpath)

    assert assay.aid == _TEST_AID
    assert len([r for r in assay.records if r.cid is not None]) == 4
    assert assay.pubmed_id == 24188023

    ki_row = next(r for r in assay.records if r.activity_name == "Ki")
    assert ki_row.activity_qualifier == ">"
    assert ki_row.activity_unit == "nM"


@pytest.mark.asyncio
async def test_load_assay_cid_accessors(tmp_path):
    r"""CID 100 is both Active and Inactive; accessors dedupe within each list."""
    pqpath = os.path.join(str(tmp_path), "mini.parquet")
    _write_mini_parquet(pqpath)

    assay = PubChemBulk().load_assay(_TEST_AID, pqpath)

    assert set(assay.active_cids) == {100, 200}
    assert set(assay.inactive_cids) == {100}
    assert 100 in assay.active_cids and 100 in assay.inactive_cids
    assert set(assay.unique_cids) == {100, 200, 300}

    aff = assay.binding_affinity(cids=[100])
    values = {(round(v, 2), name, outcome) for v, name, outcome in aff[100]}
    assert (2.5, "IC50", "Active") in values
    assert (40.0, "IC50", "Inactive") in values

    active = assay.binding_affinity(cids=[100], outcome="active")
    assert active[100] == [(2.5, "IC50", "Active")]


@pytest.mark.asyncio
async def test_load_assay_missing_file(tmp_path):
    r"""A missing parquet raises a clear FileNotFoundError, not a schema error."""
    with pytest.raises(FileNotFoundError):
        PubChemBulk().load_assay(
            _TEST_AID, os.path.join(str(tmp_path), "does_not_exist.parquet")
        )


@pytest.mark.asyncio
async def test_remote_signature():
    r"""remote_signature returns a Last-Modified|Content-Length signature."""
    async with PubChemBulk() as bulk:
        sig = await bulk.remote_signature()

    assert isinstance(sig, str)
    assert "|" in sig

    last_mod, length = sig.split("|")
    assert last_mod
    assert int(length) > 0


@pytest.mark.asyncio
async def test_load_assay_pmid_disagrees_is_none(tmp_path):
    r"""When records cite different PMIDs, assay-level pubmed_id is None; rows keep theirs."""
    pqpath = tmp_path / "disagree.parquet"
    rows = {
        c: v
        for c, v in {
            "AID": ["777777", "777777"],
            "CID": ["100", "200"],
            "Activity Outcome": ["Active", "Active"],
            "Activity Name": ["IC50", "IC50"],
            "Activity Value": ["2.5", "1.1"],
            "Activity Qualifier": ["=", "="],
            "Activity Unit": ["uM", "uM"],
            "Protein Accession": ["P1", "P1"],
            "Gene ID": ["7157", "7157"],
            "PMID": ["111", "222"],
        }.items()
    }
    pq.write_table(
        pa.table({c: pa.array(rows[c], type=pa.string()) for c in _DUMP_COLUMNS}),
        str(pqpath),
    )
    assay = PubChemBulk().load_assay(777777, pqpath)
    by_cid = {r.cid: r.pubmed_id for r in assay.records}
    assert by_cid[100] == 111 and by_cid[200] == 222
    assert assay.pubmed_id is None
