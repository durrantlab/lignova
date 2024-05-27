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


def test_get_smiles():
    r"""Retrieve SMILES from PubChem"""
    pubchem = PubChemAPI()
    smiles_test = pubchem.get_smiles(2244, ["IsomericSMILES", "ExactMass"])
    smiles_ref = "CC(=O)OC1=CC=CC=C1C(=O)O"
    mass_ref = "180.04225873"
    assert smiles_test["smiles"] == smiles_ref
    assert smiles_test["mass"] == mass_ref
