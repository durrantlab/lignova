import os
import shutil

import pytest

from lignova.io import *
from lignova.structure.editing import *
from lignova.structure.protein import Protein
from lignova.structure.utils import (
    check_pdb_mutation,
    is_xray_structure,
    separate_protein_ligand,
)

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
}


def prep_dirs():
    os.makedirs(context_protein_6Oav["write_dir"])


if not os.path.exists(context_protein_6Oav["write_dir"]):
    prep_dirs()


def test_get_pdb_6oav():
    r"""Retrieve PDB from RCSB"""
    pdb_test = Protein.get_pdb_from_rcsb("6OAV")
    with open(context_protein_6Oav["file_path"], encoding="utf-8") as f:
        pdb_ref = f.read()
    assert pdb_test == pdb_ref


def test_load_from_pdb_id_6oav():
    r"""Load PDB from RCSB"""
    protein = Protein()
    # protein._load_from_pdb_id("6OAV",write=True,write_path=context_protein_6Oav["write_dir"]+'/6oav.pdb')
    protein._load_from_pdb_id(pdb_id="6OAV")
    with open(context_protein_6Oav["file_path"], encoding="utf-8") as f:
        pdb_ref = f.read()
    pdb_test = protein.pdb
    assert pdb_test == pdb_ref


def test_get_mda_universe():
    protein = Protein()
    protein._load_from_pdb_id(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    protein_p = get_mda_universe(protein._pdb_file_path)
    assert len(set(protein_p.segments.segids)) == 1


def test_select_chains():
    protein = Protein()
    protein._load_from_pdb_id(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    protein_p = get_mda_universe(protein._pdb_file_path)
    protein_p = select_chains(protein_p)
    assert set(protein_p.segments.segids) == set("A")


def test_is_xray_structure():
    protein = Protein()
    protein._load_from_pdb_id(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    assert is_xray_structure(protein._pdb_file_path)


def test_separate_protein_ligand():
    # rSeparate protein and ligand from PDB
    protein = Protein()
    protein._load_from_pdb_id(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    protein_p, ligand_p = separate_protein_ligand(protein._pdb_file_path)
    assert len(set(protein_p.segments.segids)) == 1
    assert ligand_p.resnames.all() == "M3A"


def test_select_residues():
    protein = Protein()
    protein._load_from_pdb_id(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    protein_p = get_mda_universe(protein._pdb_file_path)
    protein_p = select_residues(protein_p, residues=["M3A"])
    assert protein_p.resnames.all() == "M3A"


"""
def test_check_pdb_mutation():
    assert not check_pdb_mutation(context_protein_6Oav["id"])
    assert check_pdb_mutation("3c5e")
"""
