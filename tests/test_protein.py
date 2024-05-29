import os
import shutil

import pytest

from lignova.io import *
from lignova.structure.editing import *
from lignova.structure.protein import Protein
from lignova.structure.utils import *

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


def test_get_rcsb_data():
    data = get_rcsb_data(context_protein_6Oav["id"])
    assert data["exptl"][0]["method"] == "X-RAY DIFFRACTION"


def test_find_resolution():
    data = get_rcsb_data(context_protein_6Oav["id"])
    assert find_resolution(context_protein_6Oav["id"], data) == 1.939


def test_has_covalent_bond():
    assert not has_covalent_bonds(context_protein_6Oav["id"])


def test_has_ligands():
    assert has_ligands(context_protein_6Oav["id"])


def test_get_entity_ids():
    entity = get_entity_ids(context_protein_6Oav["id"])
    assert entity["polymer"][0] == "1"
    assert entity["nonpolymer"][0] == "2"


def test_pdb_has_mutation():
    assert pdb_has_mutation(context_protein_6Oav["id"])
    assert pdb_has_mutation("3c5e")
    assert not pdb_has_mutation("4uxl")


def test_validate_pdb():
    assert not validate_pdb(context_protein_6Oav["id"])
    assert not validate_pdb("3c5e")
    assert validate_pdb("4uxl")


def test_validate_ligands():
    assert not validate_ligands(context_protein_6Oav["id"])
    assert validate_ligands("4uxl")


def test_get_smiles():
    protein = Protein()
    protein._load_from_pdb_id(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    smiles = get_smiles(protein._pdb_file_path)
    assert isinstance(smiles, dict)
    assert smiles["smiles"] == "c1ccc(cc1)NC(=O)n2c(nc(n2)Nc3ccc(cc3)C#N)N"
    assert smiles["stereo_smiles"] == "c1ccc(cc1)NC(=O)n2c(nc(n2)Nc3ccc(cc3)C#N)N"
