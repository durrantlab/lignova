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
}


def prep_dirs():
    os.makedirs(context_protein_6Oav["write_dir"])


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

""" WIP
def filter_lines(lines):
    # Filter out lines starting with "./"
    filtered = [line for line in lines if not (line.startswith('  ./') or line.startswith('/'))]
    # Remove trailing empty lines or lines starting with "/n"
    while filtered and (not filtered[-1].strip() or filtered[-1].startswith("/n")):
        filtered.pop()
    return filtered
"""


def filter_lines(lines):
    # If no file path to ignore specified, ignore trailing empty or '\n' lines
    filtered_lines = [line.strip() for line in lines if line.strip()]
    print(len(filtered_lines))
    return filtered_lines


"""
def test_convert_protein_to_mae():
    glide.convert_to_mae(protein_obj, context)
    prot_test_mae=Protein(file_path=context_protein_6Oav["write_dir"] + "/6oav_chA.mae")
    prot_test_mae.load(file_path=context_protein_6Oav["write_dir"] + "/6oav_chA.mae")
    prot_test = (prot_test_mae.pdb)
    assert prot_test

def test_convert_ligand_to_mae():
    glide.convert_to_mae(lig_object, context)
    lig_ref_mae=Ligand(context_protein_6Oav["ligand_file_path"])
    lig_test_mae=Ligand(file_path=context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig_prepared.mae")
    lig_ref = filter_lines(lig_ref_mae.ligand_text)
    lig_test = filter_lines(lig_ref_mae.ligand_text)
    assert lig_ref == lig_test


def test_prep_Ligand():
    prepared = glide.PrepLigand(lig_object,context)
    lig_ref=Ligand(context_protein_6Oav["ligand_prepared_path"])
    lig_ref=filter_lines(lig_ref.ligand_text)
    lig_test=filter_lines(prepared.ligand_text)
    #assert lig_ref == lig_test
    assert prepared.file_name == "6oav_A_M3A_lig_prepared.mae"
    assert context.lig_ph == '7.0'
    assert context.lig_pht == '2.0'
    assert context.lig_stereoisomers == '32'
    assert context.lig_forcefield ==  "14"
    
def test_prep_Protein():
    prot_mae=Protein(file_path=os.path.join(context_protein_6Oav["write_dir"],"6oav_chA.mae"))
    prepared=glide.PrepProtein(prot_mae,context)
    prepared.load(file_path=context_protein_6Oav['write_dir']+"/6oav_chA_protein_prepared.mae")
    prot_ref=Protein(context_protein_6Oav["protein_prepared_path"])
    prot_ref.load(file_path=context_protein_6Oav["protein_prepared_path"])
    prot_ref=filter_lines(prot_ref.pdb)
    prot_test=filter_lines(prepared.pdb)
    #assert prot_ref == prot_test
    assert context.epik_ph == '7.0'
    assert context.epik_pht == '2.0'
    assert context.forcefield == 'OPLS_2005'
    assert context.grid_innerbox == '10'
    assert context.prot_rmsd == '0.3'
    assert context.propka_ph == '7.0'
"""


def test_docking():
    prep_prot = Protein(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_chA_grid.zip"
    )
    prep_lig = Ligand(
        file_path=context_protein_6Oav["write_dir"] + "/6oav_A_M3A_lig_prepared.mae"
    )
    glide.run(prep_prot, prep_lig, context)
    assert context.docking_protocol == "SP"
    assert context.forcefield == "OPLS_2005"
    assert context.n_enhanced_sampling == "4"
    assert context.postdock_nposes == "100"


"""
def test_docking():
    glide.run()
    assert os.path.exists("M3A_docking_pv.maegz")
    assert os.path.exists("M3A_docking.in")
    assert os.path.exists("M3A_docking.csv")
    assert os.path.exists("M3A_docking_skip.csv")
    assert os.path.exists("M3A_docking.log")
"""
