r"""Test PubChem API."""

import os

from lignova.APIs import PubChemAPI

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_pubchem = {
    "write_dir": "./tmp/pubchem",
}


def prep_dirs():
    os.makedirs(context_pubchem["write_dir"])


if not os.path.exists(context_pubchem["write_dir"]):
    prep_dirs()


def test_get_cids():
    r"""Retrieve CIDs from PubChem"""
    pubchem = PubChemAPI()
    cids_active = pubchem.get_cids("1000", True)
    assert len(cids_active) == 36
    assert 16749973 in cids_active
    assert 730211 not in cids_active
    cids_inactive = pubchem.get_cids("1000", False)
    assert len(cids_inactive) == 21
    assert 16749973 not in cids_inactive
    assert 730211 in cids_inactive


def test_get_cids_info():
    r"""Retrieve SMILES from PubChem"""
    pubchem = PubChemAPI()
    smiles_test = pubchem.get_cids_info(2244, ["SMILES", "ExactMass"])
    smiles_ref = "CC(=O)OC1=CC=CC=C1C(=O)O"
    mass_ref = "180.04225873"
    assert smiles_test["SMILES"] == smiles_ref
    assert smiles_test["ExactMass"] == mass_ref


def test_get_binding_affinity():
    r"""Retrieve binding affinity from PubChem"""
    pubchem = PubChemAPI()
    binding_affinity_test = pubchem.get_binding_affinity(
        1057958,
        [
            "135566761",
            "135566762",
        ],
    )
    value_1, type_1 = binding_affinity_test[135566761]
    value_2, type_2 = binding_affinity_test[135566762]

    assert value_1 == 2.12
    assert value_2 == 0.88
    assert type_1 == "IC50"
    assert type_2 == "IC50"


def test_get_pubmed_id():
    r"""Retrieve PubMed ID from PubChem"""
    pubchem = PubChemAPI()
    empty = pubchem.get_pubmed_id(2244)
    pubmed_id = pubchem.get_pubmed_id(1057924)
    assert pubmed_id == 24183742
    assert not empty
