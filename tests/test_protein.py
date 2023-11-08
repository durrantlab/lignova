import os
import shutil
from lignova.structure.protein import Protein
import pytest
from lignova.structure.utils import separate_protein_ligand
from lignova.io import *
from lignova.structure.editing import *
# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
}


def prep_dirs():
    os.makedirs(context_protein_6Oav["write_dir"])


def test_get_pdb_6oav():
    r"""Retrieve PDB from RCSB"""
    pdb_test = Protein.get_pdb_from_rcsb("6OAV")
    with open(context_protein_6Oav["file_path"], encoding="utf-8") as f:
        pdb_ref = f.read()
    assert pdb_test == pdb_ref

def test_load_from_pdb_id_6oav():
    r"""Load PDB from RCSB"""
    protein = Protein()
    #protein._load_from_pdb_id("6OAV",write=True,write_path=context_protein_6Oav["write_dir"]+'/6oav.pdb')
    protein._load_from_pdb_id(pdb_id="6OAV")
    with open(context_protein_6Oav["file_path"], encoding="utf-8") as f:
        pdb_ref = f.read()
    pdb_test = protein.pdb
    assert pdb_test == pdb_ref


def test_get_mda_universe():
    protein=Protein()
    protein._load_from_pdb_id(pdb_id="6OAV",write=True,write_path=context_protein_6Oav["write_dir"]+'/6oav.pdb')
    protein_p=get_mda_universe(protein._pdb_file_path)
    assert len(set(protein_p.segments.segids)) == 1

def test_select_chains():
    protein=Protein()
    protein._load_from_pdb_id(pdb_id="6OAV",write=True,write_path=context_protein_6Oav["write_dir"]+'/6oav.pdb')
    protein_p=get_mda_universe(protein._pdb_file_path)
    protein_p=select_chains(protein_p)
    assert set(protein_p.segments.segids) == set("A")

def test_separate_protein_ligand():
    #rSeparate protein and ligand from PDB
    protein = Protein()
    protein._load_from_pdb_id(pdb_id="6OAV",write=True,write_path=context_protein_6Oav["write_dir"]+'/6oav.pdb')
    protein_p, ligand_p = separate_protein_ligand(protein._pdb_file_path)
    assert len(set(protein_p.segments.segids)) == 1
    assert ligand_p.resnames.all() == "M3A"

