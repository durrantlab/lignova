r"""Tests for protein structure class and relevant functions."""

import os

from lignova.structure.editing import convert_cif2pdb, read_cif, select_water
from lignova.structure.protein import Protein
from lignova.structure.utils import (
    find_resolution,
    get_entity_ids,
    get_ligand_names,
    get_mda_universe,
    get_rcsb_data,
    get_smiles,
    has_covalent_bonds,
    has_ligands,
    is_xray_structure,
    pdb_has_mutation,
    select_chains,
    select_residues,
    separate_protein_ligand,
    validate_ligands,
    validate_pdb,
)

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
    "cif_file": "./files/6qsw.cif",
}


def prep_dirs():
    r"""Prepare directories."""
    os.makedirs(context_protein_6Oav["write_dir"])


if not os.path.exists(context_protein_6Oav["write_dir"]):
    prep_dirs()


def test_select_water():
    r"""test Select water atoms."""
    protein = Protein()
    protein.load(
        pdb_id="1GYY",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/1GYY.pdb",
    )
    protein_p = get_mda_universe(protein.file_path)

    surface_water = select_water(protein_p, water_selection="surface", ligand="FHC")
    bridging_water = select_water(
        protein_p, water_selection="interfacial", ligand="FHC"
    )
    assert bridging_water.atoms.issubset(surface_water.atoms)
    assert len(bridging_water.atoms) == 4
    # check that chain residd 2149 is in bridging water
    assert 2149 in [atom.resid for atom in bridging_water.atoms]
    assert 2107 in [atom.resid for atom in surface_water.atoms]
    assert 2008 in [atom.resid for atom in surface_water.atoms]
    assert 2150 in [atom.resid for atom in surface_water.atoms]


def test_get_pdb_6oav():
    r"""Retrieve PDB from RCSB"""
    pdb_test = Protein.get_pdb_from_rcsb("6OAV")
    with open(context_protein_6Oav["file_path"], encoding="utf-8") as f:
        pdb_ref = f.read()
    assert pdb_test == pdb_ref


def testload_6oav():
    r"""Load PDB from RCSB"""
    protein = Protein()
    # protein.load("6OAV",write=True,write_path=context_protein_6Oav["write_dir"]+'/6oav.pdb')
    protein.load(pdb_id="6OAV")
    with open(context_protein_6Oav["file_path"], encoding="utf-8") as f:
        pdb_ref = f.read()
    pdb_test = protein.pdb
    assert pdb_test == pdb_ref


def test_get_mda_universe():
    r"""Test get MDAnalysis universe"""
    protein = Protein()
    protein.load(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    protein_p = get_mda_universe(protein.file_path)
    assert len(set(protein_p.segments.segids)) == 1


def test_select_chains():
    r"""Test select chains"""
    protein = Protein()
    protein.load(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    protein_p = get_mda_universe(protein.file_path)
    protein_p = select_chains(protein_p)
    assert set(protein_p.segments.segids) == set("A")


def test_is_xray_structure():
    r"""Test if structure is X-ray"""
    protein = Protein()
    protein.load(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    assert is_xray_structure(protein.file_path)


def test_separate_protein_ligand():
    r"""Test separate protein and ligand"""
    protein = Protein()
    protein.load(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    protein_p, ligand_p = separate_protein_ligand(protein.file_path)
    assert len(set(protein_p.segments.segids)) == 1
    assert ligand_p.resnames.all() == "M3A"
    protein.load(
        pdb_id="4ZBG",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/4zbg.pdb",
    )
    protein_p, ligand_p = separate_protein_ligand(protein.file_path)
    assert len(set(protein_p.segments.segids)) == 1
    assert ligand_p.resnames.all() == "ACO"


def test_select_residues():
    r"""Test select residues"""
    protein = Protein()
    protein.load(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    protein_p = get_mda_universe(protein.file_path)
    protein_p = select_residues(protein_p, residues=["M3A"])
    assert protein_p.resnames.all() == "M3A"


def test_get_rcsb_data():
    r"""Test get RCSB data"""
    data = get_rcsb_data(context_protein_6Oav["id"])
    assert data["exptl"][0]["method"] == "X-RAY DIFFRACTION"


def test_find_resolution():
    r"""Test find resolution"""
    data = get_rcsb_data(context_protein_6Oav["id"])
    assert find_resolution(context_protein_6Oav["id"], data) == 1.939


def test_has_covalent_bond():
    r"""Test if structure has covalent bonds"""
    assert not has_covalent_bonds(context_protein_6Oav["id"])


def test_has_ligands():
    r"""Test if structure has ligands"""
    assert has_ligands(context_protein_6Oav["id"])


def test_get_entity_ids():
    r"""Test get entity ids"""
    entity = get_entity_ids(context_protein_6Oav["id"])
    assert entity["polymer"][0] == "1"
    assert entity["nonpolymer"][0] == "2"


def test_pdb_has_mutation():
    r"""Test if PDB has mutation"""
    assert pdb_has_mutation(context_protein_6Oav["id"])
    assert pdb_has_mutation("3c5e")
    assert not pdb_has_mutation("4uxl")


def test_validate_pdb():
    r"""Test validate PDB"""
    assert not validate_pdb(context_protein_6Oav["id"])
    assert not validate_pdb("3c5e")
    assert validate_pdb("7dbk")
    assert validate_pdb("4uxl")


def test_validate_ligands():
    r"""Test validate ligands"""
    assert validate_ligands(context_protein_6Oav["id"])
    assert validate_ligands("4uxl")
    assert not validate_ligands("7dbk")
    assert not validate_ligands("4zbg")


def test_get_ligand_names():
    r"""Test get ligand names"""
    protein = Protein()
    protein.load(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    ligands = get_ligand_names(protein.file_path)
    assert ligands == ["M3A"]


def test_get_smiles():
    r"""Test get SMILES"""
    protein = Protein()
    protein.load(
        pdb_id="6OAV",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/6oav.pdb",
    )
    ligand = get_ligand_names(protein.file_path)[0]
    smiles = get_smiles(ligand)
    assert isinstance(smiles, dict)
    assert smiles["smiles"] == "c1ccc(cc1)NC(=O)n2c(nc(n2)Nc3ccc(cc3)C#N)N"
    assert smiles["stereo_smiles"] == "c1ccc(cc1)NC(=O)n2c(nc(n2)Nc3ccc(cc3)C#N)N"


def test_read_cif():
    r"""Test read CIF"""
    protein = Protein(context_protein_6Oav["cif_file"])
    protein.load()
    data = read_cif(protein.file_path)
    assert protein.file_path == context_protein_6Oav["cif_file"]
    assert protein.file_id == "6qsw"
    assert protein.file_ext == "cif"
    assert data.resolution == 1.64


def test_convert_cif2pdb():
    r"""Test convert CIF to PDB"""
    protein = Protein(context_protein_6Oav["cif_file"])
    protein.load()
    convert_cif2pdb(
        protein.file_path, os.path.join(context_protein_6Oav["write_dir"], "6qsw.pdb")
    )
    assert os.path.exists(os.path.join(context_protein_6Oav["write_dir"], "6qsw.pdb"))


def test_separate_protein_ligand_water():
    r"""Test separate protein ligand with different water selection"""
    protein = Protein()
    protein.load(
        pdb_id="1GYY",
        write=True,
        write_path=context_protein_6Oav["write_dir"] + "/1GYY.pdb",
    )
    protein_surf, ligand_suf = separate_protein_ligand(
        protein.file_path, water_selection="surface", remove_water=False
    )
    protein_bridge, ligand_bridge = separate_protein_ligand(
        protein.file_path, water_selection="interfacial", remove_water=False
    )
    protein_all, ligand_all = separate_protein_ligand(
        protein.file_path, water_selection="all", remove_water=False
    )
    assert (
        ligand_suf.resnames.all()
        == "FHC"
        == ligand_bridge.resnames.all()
        == ligand_all.resnames.all()
    )
    assert (
        len(set(protein_surf.segments.segids))
        == 2
        == len(set(protein_bridge.segments.segids))
        == len(set(protein_all.segments.segids))
    )
    assert 2107 in [
        atom.resid for atom in protein_surf.atoms if atom.resname == "HOH"
    ] and 2107 in [atom.resid for atom in protein_all.atoms if atom.resname == "HOH"]
    assert 2008 in [
        atom.resid for atom in protein_surf.atoms if atom.resname == "HOH"
    ] and 2008 in [atom.resid for atom in protein_all.atoms if atom.resname == "HOH"]
    assert 2149 in [
        atom.resid for atom in protein_bridge.atoms if atom.resname == "HOH"
    ] and 2149 in [atom.resid for atom in protein_all.atoms if atom.resname == "HOH"]
