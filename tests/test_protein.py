import os

from lignova.structure.protein import Protein

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
}


def prep_dirs():
    os.makedirs(context_protein_6Oav["write_dir"])


def test_get_pdb_6Oav():
    r"""Retrieve PDB from RCSB"""
    pdb_test = Protein.get_pdb_from_rcsb("6OAV")
    with open(context_protein_6Oav["file_path"], encoding="utf-8") as f:
        pdb_ref = f.read()
    assert pdb_test == pdb_ref
