import csv
import glob
import os
import shutil

import pytest

from lignova.docking import Glide
from lignova.docking.contexts import GlideContext
from lignova.docking.docking import Docking
from lignova.io import *
from lignova.structure.editing import *
from lignova.structure.ligand import Ligand, PreparedLigand
from lignova.structure.protein import PreparedProtein, Protein
from lignova.structure.utils import separate_protein_ligand

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
    os.makedirs(context_protein_6Oav["write_dir"])


if not os.path.exists(context_protein_6Oav["write_dir"]):
    prep_dirs()

# to generate the protein and ligand files for 6oav
prot_raw = Protein(context_protein_6Oav["file_path"])
prot_raw.load(file_path=context_protein_6Oav["file_path"])
protein, lig = separate_protein_ligand(prot_raw._pdb_file_path)
write_mda_universe(protein, context_protein_6Oav["write_dir"] + "/6oav_chA.pdb")
write_mda_universe(lig, context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig.pdb")
protein_obj = Protein(file_path=context_protein_6Oav["write_dir"] + "/6oav_chA.pdb")
lig_object = Ligand(
    file_path=os.path.join(context_protein_6Oav["write_dir"], "6oav_A_M3A_lig.pdb")
)

glide = Glide()
context = GlideContext.get_current()


def filter_lines(lines):
    # If no file path to ignore specified, ignore trailing empty or '\n' lines
    filtered_lines = [line.strip() for line in lines if line.strip()]
    return filtered_lines


def test_convert_protein_to_mae():
    glide.convert_to_mae(protein_obj, context)
    prot_test_mae = Protein(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_chA.mae"
    )
    prot_test_mae.load(file_path=context_protein_6Oav["write_dir"] + "/6oav_chA.mae")
    prot_test = prot_test_mae.pdb
    assert prot_test


def test_convert_ligand_to_mae():
    glide.convert_to_mae(lig_object, context)
    lig_ref_mae = Ligand(context_protein_6Oav["ligand_file_path"])
    lig_test_mae = Ligand(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig_prepared.mae"
    )
    lig_ref = filter_lines(lig_ref_mae.ligand_text)
    lig_test = filter_lines(lig_ref_mae.ligand_text)
    assert lig_ref == lig_test


def test_prep_Ligand():
    prepared = glide.PrepLigand(lig_object, context)
    prep_info = prepared.get_info(context)
    # read the ref and test files and compare the lines
    with open(context_protein_6Oav["ligand_prepared_path"], "r") as file:
        lig_ref = file.readlines()
    with open(
        context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig_prepared.mae", "r"
    ) as file:
        lig_test = file.readlines()
    assert lig_ref[113] == lig_test[113]
    assert prepared.file_name == "6oav_A_M3A_lig_prepared.mae"
    assert prep_info["pH"] == "7.0"
    assert prep_info["pHt"] == "2.0"
    assert prep_info["stereoisomers"] == "32"
    assert prep_info["forcefield"] == "14"


def test_prep_Protein():
    prot_mae = Protein(
        file_path=os.path.join(context_protein_6Oav["write_dir"], "6oav_chA.mae")
    )
    prepared = glide.PrepProtein(prot_mae, context)
    prepared.load(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_chA_protein_prepared.mae"
    )
    prep_info = prepared.get_info(context)
    # read the grid log file and get the OUTERBOX value and the grid center
    with open(context_protein_6Oav["grid_log_path"], "r") as file:
        lines = file.readlines()
        for line in lines:
            if "OUTERBOX" in line:
                outerbox_ref = line.split()[1]
            if "GRID_CENTER" in line:
                grid_center_ref = [line.split()[1], line.split()[2], line.split()[3]]
    with open(context_protein_6Oav["write_dir"] + "/6oav_chA_grid.log", "r") as file:
        lines = file.readlines()
        for line in lines:
            if "OUTERBOX" in line:
                outerbox_test = line.split()[1]
            if "GRID_CENTER" in line:
                grid_center_test = [line.split()[1], line.split()[2], line.split()[3]]
    assert outerbox_ref == outerbox_test
    assert grid_center_ref == grid_center_test
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
    prep_prot = Protein(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_chA_grid.zip"
    )
    prep_lig = Ligand(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig_prepared.mae"
    )
    glide.run(prep_prot, prep_lig, context)
    with open(context_protein_6Oav["docking_results_path"], "r") as file:
        docking_results_ref_reader = csv.DictReader(file)
        docking_results_ref = list(docking_results_ref_reader)
    test_file_path = context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig_docking.csv"
    with open(test_file_path, "r") as file:
        docking_results_test_reader = csv.DictReader(file)
        docking_results_test = list(docking_results_test_reader)
    assert (
        docking_results_ref[0]["r_i_glide_gscore"]
        == docking_results_test[0]["r_i_glide_gscore"]
    )
    assert len(docking_results_ref) == len(docking_results_test)
    assert context.docking_protocol == "SP"
    assert context.forcefield == "OPLS_2005"
    assert context.n_enhanced_sampling == "4"
    assert context.postdock_nposes == "100"
