import os

import pytest

from lignova.analysis.rmsd import RMSD
from lignova.docking import Glide
from lignova.docking.contexts import GlideContext
from lignova.structure.ligand import DockedLigand
from lignova.structure.protein import Protein

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
    "docked_ligand_filepath": "/home/mma121/PubChem_small/try_schrodinger/6oav_validation/6oav_m3a_combind_sorted.maegz",
}


def prep_dirs():
    os.makedirs(context_protein_6Oav["write_dir"])


if not os.path.exists(context_protein_6Oav["write_dir"]):
    prep_dirs()


# convert pdb to mae using the convert_to_mae function in glide.py
protein = Protein(context_protein_6Oav["file_path"])
protein.load(file_path=context_protein_6Oav["file_path"])
glide = Glide()
context = GlideContext.get_current()
glide.convert_to_mae(protein, context)


def test_rmsd():
    ligand = DockedLigand(context_protein_6Oav["docked_ligand_filepath"])
    reference = Protein(os.path.join(context_protein_6Oav["write_dir"], "6oav.mae"))
    rmsd = RMSD(ligand, reference, context)
    output_filename = os.path.join(
        context_protein_6Oav["write_dir"], "6oav_m3a_rmsd.csv"
    )
    rmsd.rmsd_schrodinger(output_filename)
    # assert os.path.exists(output_filename)
    # os.remove(output_filename)
