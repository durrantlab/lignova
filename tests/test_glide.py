import glob
import os
import shutil

import pytest

from lignova.docking import Glide
from lignova.docking.contexts import GlideContext
from lignova.docking.docking import Docking
from lignova.structure.ligand import Ligand, PreparedLigand
from lignova.structure.protein import PreparedProtein, Protein
from lignova.structure.utils import separate_protein_ligand
from lignova.io import *
from lignova.structure.editing import *


# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))



def prep_dirs():
    os.makedirs(context_protein_6Oav["write_dir"])

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
}

prot_object = Protein(context_protein_6Oav["file_path"])

# ligand_file_mae = "./files/6OAV_A_M3A_lig.mae"
# lig_object_mae = Ligand(ligand_file_mae)
# prepared_ligand = PreparedLigand(file_path="./files/6OAV_A_M3A_lig_prepared.mae")
# prepared_protein = PreparedProtein(file_path="./files/6OAV_A_grid.zip")
# glide = Glide(prepared_ligand, prepared_protein)

glide = Glide()
context = GlideContext.get_current()


def test_convert_protein_to_mae():
    glide.convert_to_mae(prot_object, context)
    assert os.path.exists("6OAV_A.mae")

"""
def test_convert_Ligand_to_mae():
    glide.convert_to_mae(lig_object,context)
    assert os.path.exists("6OAV_A_M3A_lig.mae")


def test_prep_Protein():
    prepared = glide.PrepProtein(prot_object)
    assert os.path.exists("6OAV_A_grid.zip")
    assert os.path.exists("6OAV_A_protein_prepared.mae")
    assert os.path.exists("6OAV_A_grid.log")
    assert prepared.epik
    assert not prepared.water
    assert prepared.pH == 7.0 / 2.0
    assert prepared.propka
    assert prepared.forcefield == "OPLS_2005"
    assert prepared.RMSD == 0.3
    assert prepared.grid_center == "ligand"
    assert prepared.grid_innerbox == 10
    assert prepared.propka_pH == 7.0


def test_prep_Ligand():
    prepared = glide.PrepLigand(lig_object)
    assert prepared.epik
    assert prepared.pH == 7.0 / 2.0
    assert prepared.stereo == "All"
    assert prepared.forcefield == "OPLS_2005"
    assert prepared.lig_name == "M3A"
    assert prepared.stereo_num == 32
    assert os.path.exists("6OAV_A_M3A_lig_prepared.mae")
    assert os.path.exists("6OAV_A_M3A_lig_prepared.log")


def test_prep_Ligand_mae():
    prepared = glide.PrepLigand(lig_object_mae)
    assert prepared.epik
    assert prepared.pH == 7.0 / 2.0
    assert prepared.stereo == "All"
    assert prepared.forcefield == "OPLS_2005"
    assert prepared.lig_name == "M3A"
    assert prepared.stereo_num == 32
    assert os.path.exists("6OAV_A_M3A_lig_prepared.mae")
    assert os.path.exists("6OAV_A_M3A_lig_prepared.log")


def test_docking():
    glide.run()
    assert os.path.exists("M3A_docking_pv.maegz")
    assert os.path.exists("M3A_docking.in")
    assert os.path.exists("M3A_docking.csv")
    assert os.path.exists("M3A_docking_skip.csv")
    assert os.path.exists("M3A_docking.log")
"""
