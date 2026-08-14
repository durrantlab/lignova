# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Test the docking analysis module for calculating RMSD."""

import os

import numpy as np
import pytest
from rdkit.Chem import rdmolfiles

from lignova.analysis import mdaRMSD, obabelRMSD, spyrmsdRMSD
from lignova.analysis.utils import (
    mae_convert,
    obabel_convert,
)
from lignova.structure.editing import write_mda_universe
from lignova.structure.ligand import DockedLigand, Ligand
from lignova.structure.utils import separate_protein_ligand

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
    "docked_ligand_filepath": "./files/6oav/6oav_m3a_top_pose_pv.maegz",
}


def prep_dirs():
    """Prepare directories for writing files."""
    os.makedirs(context_protein_6Oav["write_dir"])


if not os.path.exists(context_protein_6Oav["write_dir"]):
    prep_dirs()


@pytest.fixture
def docked_pdb():
    """Convert the .mae/.maegz docked ligand to PDB for tests that need it."""
    src = context_protein_6Oav["docked_ligand_filepath"]
    base = os.path.splitext(os.path.basename(src))[0]
    dst_dir = context_protein_6Oav["write_dir"]
    ligand_pdb = os.path.join(dst_dir, f"{base}.pdb")

    if not os.path.exists(ligand_pdb):
        mae_convert(src, ligand_pdb)

    return ligand_pdb


@pytest.fixture
def docked_pdb_with_protein():
    """Convert .maegz and also write the receptor to a separate file."""
    src = context_protein_6Oav["docked_ligand_filepath"]
    base = os.path.splitext(os.path.basename(src))[0]
    dst_dir = context_protein_6Oav["write_dir"]
    ligand_pdb = os.path.join(dst_dir, f"{base}_with_prot.pdb")
    protein_pdb = os.path.join(dst_dir, f"{base}_with_prot_protein.pdb")

    if not os.path.exists(ligand_pdb):
        mae_convert(src, ligand_pdb, protein=True)

    return ligand_pdb, protein_pdb


@pytest.fixture
def reference_pdb():
    """Prepare the reference PDB file for RMSD calculations."""
    src = context_protein_6Oav["file_path"]
    dst = os.path.join(context_protein_6Oav["write_dir"], "reference_lig.pdb")

    if not os.path.exists(dst):
        _, ligand = separate_protein_ligand(src)
        write_mda_universe(ligand, dst)

    return dst


@pytest.fixture
def docked_sdf(docked_pdb):
    """Convert docked PDB to SDF using RDKit."""
    dst = docked_pdb.replace(".pdb", ".sdf")
    if not os.path.exists(dst):
        mol = rdmolfiles.MolFromPDBFile(docked_pdb, sanitize=False, removeHs=True)
        if mol is None:
            raise ValueError(f"RDKit failed to load {docked_pdb}")
        rdmolfiles.MolToMolFile(mol, dst)
    return dst


@pytest.fixture
def reference_sdf(reference_pdb):
    """Convert reference PDB to SDF using RDKit (not obabel) for spyrmsd."""
    dst = os.path.join(context_protein_6Oav["write_dir"], "reference_lig.sdf")

    if not os.path.exists(dst):
        mol = rdmolfiles.MolFromPDBFile(reference_pdb, sanitize=False, removeHs=True)
        if mol is None:
            raise ValueError(f"RDKit failed to load {reference_pdb}")
        rdmolfiles.MolToMolFile(mol, dst)

    return dst


def test_obabel_convert():
    """Test the conversion of file formats using Open Babel."""
    output_filename = os.path.join(context_protein_6Oav["write_dir"], "6oav.sdf")
    obabel_convert(context_protein_6Oav["file_path"], output_filename)
    assert os.path.exists(output_filename)


def test_mae_convert_pdb(docked_pdb):
    """Test that mae_convert produces a valid ligand PDB with HETATM records."""
    assert os.path.exists(docked_pdb)
    with open(docked_pdb) as f:
        content = f.read()
    assert "HETATM" in content
    assert "ATOM  " not in content
    assert "M3A" in content


def test_mae_convert_protein_output(docked_pdb_with_protein):
    """Test that protein=True writes a separate receptor PDB with ATOM records."""
    ligand_pdb, protein_pdb = docked_pdb_with_protein
    assert os.path.exists(ligand_pdb)
    assert os.path.exists(protein_pdb)
    with open(protein_pdb) as f:
        content = f.read()
    assert "ATOM  " in content
    assert "CG  LYS A 677  " in content
    assert "M3A" not in content


def test_rmsd_mda(docked_pdb, reference_pdb):
    """Test the RMSD calculation using MDAnalysis."""
    ligand = DockedLigand(docked_pdb)
    reference = Ligand(reference_pdb)
    rmsd_calc = mdaRMSD(
        target=ligand,
        reference=reference,
    )
    value = rmsd_calc.calculate()
    assert isinstance(value, list)
    assert isinstance(value[0], float)
    assert np.isclose(value[0], 1.1236, atol=1e-4)


def test_rmsd_obabel(docked_sdf, reference_sdf):
    """Test the RMSD calculation using Open Babel."""
    ligand = DockedLigand(docked_sdf)
    reference = Ligand(reference_sdf)
    rmsd_calc = obabelRMSD(target=ligand, reference=reference)
    values = rmsd_calc.calculate()
    assert isinstance(values, list)
    assert len(values) == 1
    assert np.isclose(values[0], 0.58534, atol=1e-4)


def test_spyrmsd_calculation(docked_sdf, reference_sdf):
    """Test the RMSD calculation using spyrmsd."""
    ligand = DockedLigand(docked_sdf)
    reference = Ligand(reference_sdf)
    rmsd_calc = spyrmsdRMSD(target=ligand, reference=reference)
    values_api = rmsd_calc.calculate(backend="api")
    values_cli = rmsd_calc.calculate(backend="cli")
    values_api_mcs = rmsd_calc.calculate(backend="api", mcs=True)
    assert values_api == values_api_mcs
    assert np.allclose(values_api, values_cli, atol=1e-4)
    assert isinstance(values_api, list)
    assert len(values_api) == 1
    assert all(isinstance(v, float) for v in values_api)
    assert np.isclose(values_api[0], 0.58534, atol=1e-4)


def test_spyrmsd_no_symmetry(docked_sdf, reference_sdf):
    """Test spyrmsd without symmetry correction."""
    ligand = DockedLigand(docked_sdf)
    reference = Ligand(reference_sdf)
    rmsd_calc = spyrmsdRMSD(target=ligand, reference=reference)
    values_sym_api = rmsd_calc.calculate(symmetry=True, backend="api")
    values_sym_cli = rmsd_calc.calculate(symmetry=True, backend="cli")
    values_sym_api_mcs = rmsd_calc.calculate(symmetry=True, backend="api", mcs=True)
    assert values_sym_api == values_sym_api_mcs
    assert np.allclose(values_sym_api, values_sym_cli, atol=1e-4)
    assert np.isclose(values_sym_api[0], 0.58534, atol=1e-4)
    values_nosym = rmsd_calc.calculate(symmetry=False, backend="api")
    values_nosym_cli = rmsd_calc.calculate(symmetry=False, backend="cli")
    assert np.allclose(values_nosym, values_nosym_cli, atol=1e-4)
    assert np.isclose(values_nosym[0], 1.1236, atol=1e-4)
    assert len(values_sym_api) == len(values_nosym)
    assert values_nosym[0] > values_sym_api[0]
