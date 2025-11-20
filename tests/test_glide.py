"""
This module contains tests for the Glide docking process.
"""

import csv
import os

import numpy as np
import pytest
from loguru import logger

from lignova.structure.ligand import Ligand, PreparedLigand
from lignova.structure.protein import PreparedProtein, Protein
from lignova.structure.utils import separate_protein_ligand

# Skip this entire module as it had Schrödinger-dependent code
if os.getenv("SCHRODINGER") is None:
    pytest.skip(
        "Schrödinger is not installed or the $SCHRODINGER environment variable is not set. Skipping Glide tests.",
        allow_module_level=True,
    )

from lignova.docking import Glide
from lignova.docking.contexts import GlideContext
from lignova.structure.editing import write_mda_universe

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))


context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
    "ligand_file_path": "./files/6oav/6OAV_A_M3A_lig.mae",
    "ligand_prepared_path": "./files/6oav/6OAV_A_M3A_lig_prepared.mae",
    "protein_prepared_path": "./files/6oav/6OAV_A_protein_prepared.mae",
    "docking_results_path": "./files/6oav/6OAV_A_M3A_lig_docking.csv",
    "grid_log_path": "./files/6oav/6OAV_chA_grid.log",
}


def prep_dirs():
    """Prepare directories for writing files."""
    os.makedirs(context_protein_6Oav["write_dir"])


if not os.path.exists(context_protein_6Oav["write_dir"]):
    prep_dirs()


# to generate the protein and ligand files for 6oav
prot_raw = Protein(context_protein_6Oav["file_path"])
prot_raw.load(file_path=context_protein_6Oav["file_path"])
protein, lig = separate_protein_ligand(prot_raw.file_path)
write_mda_universe(protein, context_protein_6Oav["write_dir"] + "/6oav_chA.pdb")
write_mda_universe(lig, context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig.pdb")
protein_obj = Protein(file_path=context_protein_6Oav["write_dir"] + "/6oav_chA.pdb")
lig_object = Ligand(
    file_path=os.path.join(context_protein_6Oav["write_dir"], "6oav_A_M3A_lig.pdb")
)

glide = Glide()
context = GlideContext.get_current()


def filter_lines(lines):
    """Filter out empty lines and strip whitespace."""
    # If no file path to ignore specified, ignore trailing empty or '\n' lines
    filtered_lines = [line.strip() for line in lines if line.strip()]
    return filtered_lines


def extract_m_bond_section(file_content):
    """Extract the m_bond section from the file content."""
    start_marker = "m_bond["
    end_marker = "}"
    start_index = file_content.find(start_marker)
    end_index = file_content.find(end_marker, start_index) + len(end_marker)
    return file_content[start_index:end_index]


def parse_mae_atom_block(lines):
    """
    Extract the atom data rows from an MAE file's m_atom block.
    Returns a list of token lists, one per atom.
    """
    in_block = False
    in_data = False
    atoms = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("m_atom["):
            in_block = True
            continue
        if in_block and stripped == ":::":
            if not in_data:
                # first ':::' marks start of data
                in_data = True
                continue
            else:
                # second ':::' marks end of data
                break
        if in_data:
            if not stripped:
                continue
            tokens = stripped.split()
            if len(tokens) < 12:
                continue
            atoms.append(tokens)
    return atoms


def test_convert_protein_to_mae():
    """Test conversion of protein to MAE format."""
    glide.convert_to_mae(protein_obj, context)
    prot_test_mae = Protein(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_chA.mae"
    )
    prot_test_mae_content = prot_test_mae.pdb
    assert prot_test_mae_content is not None


def test_convert_ligand_to_mae():
    """Test conversion of ligand to MAE format."""
    glide.convert_to_mae(lig_object, context)
    lig_ref_mae = Ligand(context_protein_6Oav["ligand_file_path"])
    lig_test_mae = Ligand(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig.mae"
    )

    # Extract the m_bond section from both reference and test files
    lig_ref_content = lig_ref_mae.ligand_text
    lig_test_content = lig_test_mae.ligand_text
    lig_ref_m_bond = extract_m_bond_section(lig_ref_content)
    lig_test_m_bond = extract_m_bond_section(lig_test_content)

    # Compare the extracted m_bond sections
    assert lig_ref_m_bond == lig_test_m_bond


def test_prep_ligand():
    """Test preparation of ligand."""
    prepared = glide.PrepLigand(lig_object, context)
    prep_info = prepared.get_info(context)
    # read the ref and test files and compare the lines
    with open(
        context_protein_6Oav["ligand_prepared_path"], "r", encoding="utf-8"
    ) as file:
        lig_ref = file.readlines()
        lig_ref_atom_block = parse_mae_atom_block(lig_ref)
    with open(
        context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig_prepared.mae",
        "r",
        encoding="utf-8",
    ) as file:
        lig_test = file.readlines()
        lig_test_atom_block = parse_mae_atom_block(lig_test)
    assert lig_ref_atom_block == lig_test_atom_block
    assert prepared.file_name == "6oav_A_M3A_lig_prepared.mae"
    assert prep_info["pH"] == "7.0"
    assert prep_info["pHt"] == "2.0"
    assert prep_info["stereoisomers"] == "32"
    assert prep_info["forcefield"] == "14"


def test_prepprotein():
    """Test preparation of protein."""
    prot_mae = Protein(context_protein_6Oav["write_dir"] + "/6oav_chA.mae")
    prepared = glide.PrepProtein(prot_mae, context)
    prep_info = prepared.get_info(context)
    # read the grid log file and get the OUTERBOX value and the grid center
    outerbox_ref = None
    outerbox_test = None
    with open(context_protein_6Oav["grid_log_path"], "r", encoding="utf-8") as file:
        lines = file.readlines()
        for line in lines:
            if "OUTERBOX" in line:
                outerbox_ref = line.split()[1]
            if "GRID_CENTER" in line:
                grid_center_ref = [line.split()[1], line.split()[2], line.split()[3]]
                grid_center_ref = [
                    float(item.strip("[], ")) for item in grid_center_ref
                ]
                logger.info("The grid center reference numbers are", grid_center_ref)
    log_file = context_protein_6Oav["write_dir"] + "/6oav_chA_grid.log"
    with open(log_file, "r", encoding="utf-8") as file:
        lines = file.readlines()
        for line in lines:
            if "OUTERBOX" in line:
                outerbox_test = line.split()[1]
            if "GRID_CENTER" in line:
                temp = [line.split()[1], line.split()[2], line.split()[3]]
                grid_center_test = [float(item.strip("[], ")) for item in temp]
                logger.info("The grid center test numbers", grid_center_test)

    assert np.isclose(float(outerbox_ref), float(outerbox_test), rtol=0.1)
    assert np.all(
        [
            np.isclose(float(grid_center_ref[i]), float(grid_center_test[i]), rtol=0.1)
            for i in range(3)
        ]
    )
    assert prep_info["epik_pH"] == "7.0"
    assert prep_info["epik_pHt"] == "2.0"
    assert prep_info["forcefield"] == "OPLS_2005"
    assert prep_info["grid_innerbox"] == "10"
    assert prep_info["rmsd"] == "0.3"
    assert prep_info["propka_pH"] == "7.0"
    assert prep_info["fillsidechains"]
    assert prep_info["disulfides"]
    assert prep_info["water_distance"]
    assert prep_info["sample_water"]
    assert prep_info["minimize_adj_h"]
    assert prep_info["hydrogen_addition"]


def test_docking():
    """Test the docking process."""
    prep_prot = PreparedProtein(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_chA_grid.zip"
    )
    prep_lig = PreparedLigand(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig_prepared.mae"
    )
    glide.run(prep_prot, prep_lig, context)
    with open(
        context_protein_6Oav["docking_results_path"], "r", encoding="utf-8"
    ) as file:
        docking_results_ref_reader = csv.DictReader(file)
        docking_results_ref = list(docking_results_ref_reader)
    test_file_path = context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig_docking.csv"
    with open(test_file_path, "r", encoding="utf-8") as file:
        docking_results_test_reader = csv.DictReader(file)
        docking_results_test = list(docking_results_test_reader)
    assert np.isclose(
        float(docking_results_ref[0]["r_i_glide_gscore"]),
        float(docking_results_test[0]["r_i_glide_gscore"]),
        rtol=0.1,
    )
    assert len(docking_results_ref) == len(docking_results_test)
    assert context.docking_protocol == "SP"
    assert context.forcefield == "OPLS_2005"
    assert context.n_enhanced_sampling == "4"
    assert context.postdock_nposes == "100"
    assert context.docking_method == "confgen"
