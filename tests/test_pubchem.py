r"""Test PubChem API."""
import os

from lignova.hdf5.pubchem import PubChemAPI

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_pubchem = {
    "write_dir": "./tmp/pubchem",
}


def prep_dirs():
    os.makedirs(context_pubchem["write_dir"])


if not os.path.exists(context_pubchem["write_dir"]):
    prep_dirs()


def test_get_cids_info():
    r"""Retrieve SMILES from PubChem"""
    pubchem = PubChemAPI()
    smiles_test = pubchem.get_cids_info(2244, ["IsomericSMILES", "ExactMass"])
    smiles_ref = "CC(=O)OC1=CC=CC=C1C(=O)O"
    mass_ref = "180.04225873"
    assert smiles_test["IsomericSMILES"] == smiles_ref
    assert smiles_test["ExactMass"] == mass_ref


def test_get_binding_affinity():
    r"""Retrieve binding affinity from PubChem"""
    pubchem = PubChemAPI()
    binding_affinity_test = pubchem.get_binding_affinity(
        1057958,
        [
            135566761,
            135566762,
        ],
    )
    assert binding_affinity_test[135566761]["Activity Value [uM]"] == "2.12"
    assert binding_affinity_test[135566762]["Activity Value [uM]"] == "0.88"
    assert (
        binding_affinity_test[135566761]["Activity Name"]
        == "IC50"
        == binding_affinity_test[135566762]["Activity Name"]
    )
