import os

import pytest
from loguru import logger

from lignova.analysis.rmsd import RMSD
from lignova.analysis.utils import (
    interconvert_mae_sdf,
    obabel_convert,
    obabel_result_parser,
)
from lignova.docking import Glide
from lignova.docking.contexts import GlideContext
from lignova.docking.utils import get_complexes
from lignova.structure.ligand import DockedLigand
from lignova.structure.protein import Protein

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
    "docked_ligand_filepath": "/home/mma121/PubChem_small/try_schrodinger/6oav_validation/6oav_m3a_top_pose_pv.maegz",
}
# "/home/mma121/PubChem_small/try_schrodinger/6oav_validation/6oav_m3a_combind_sorted.maegz",


def prep_dirs():
    os.makedirs(context_protein_6Oav["write_dir"])


# get the complex for the docking file
# get_complexes(context_protein_6Oav["docked_ligand_filepath"])

if not os.path.exists(context_protein_6Oav["write_dir"]):
    prep_dirs()

"""
# convert pdb to mae using the convert_to_mae function in glide.py
protein = Protein(context_protein_6Oav["file_path"])
protein.load(file_path=context_protein_6Oav["file_path"])
glide = Glide()
"""
context = GlideContext.get_current()

"""
def test_rmsd():
    ligand = DockedLigand(context_protein_6Oav["docked_ligand_filepath"])
    reference = Protein(os.path.join(context_protein_6Oav["write_dir"], "6oav.mae"))
    rmsd = RMSD(ligand, reference, context)
    output_filename = os.path.join(
        context_protein_6Oav["write_dir"], "6oav_m3a_rmsd"
    )
    rmsd.rmsd_schrodinger(output_filename)
"""


def test_interconvert_mae_sdf():
    # Call the interconvert_mae_sdf function
    interconvert_mae_sdf(
        test_file=context_protein_6Oav["docked_ligand_filepath"],
        output_filename=os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv.sdf"
        ),
        context=context,
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv.sdf")
    )


def test_rmsd_mda():
    ligand = DockedLigand(context_protein_6Oav["docked_ligand_filepath"])
    reference = Protein(context_protein_6Oav["file_path"])
    # Create an instance of the RMSD class
    rmsd = RMSD(ligand, reference, context)
    logger.debug(context.write_dir)
    # Call the rmsd_mda function
    value = rmsd.rmsd_mda()
    assert isinstance(value[0], float)
    assert round(value[0], 4) == 1.1226


def test_obabel_convert():
    obabel_convert(
        test_file=os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_complexes.pdb"
        ),
        output_filename=os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_obabel.sdf"
        ),
    )
    assert os.path.exists(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_obabel.sdf"
        )
    )


def test_rmsd_obabel():
    ligand = DockedLigand(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv.sdf")
    )

    reference = DockedLigand(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv.sdf")
    )
    # Create an instance of the RMSD class
    rmsd = RMSD(ligand, reference, context)
    # Call the rmsd_obabel function
    output_filename = os.path.join(
        context_protein_6Oav["write_dir"], "6oav_m3a_rmsd_obabel"
    )
    values = rmsd.rmsd_obabel(output_filename=output_filename)
    assert os.path.exists(output_filename + ".txt")


def test_rmsd_obabel2():
    ligand = DockedLigand(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_obabel.sdf"
        )
    )

    reference = DockedLigand(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_obabel.sdf"
        )
    )
    # Create an instance of the RMSD class
    rmsd = RMSD(ligand, reference, context)
    # Call the rmsd_obabel function
    output_filename = os.path.join(
        context_protein_6Oav["write_dir"], "6oav_m3a_rmsd_obabel2"
    )
    values = rmsd.rmsd_obabel(output_filename=output_filename)
    print(values)
    assert os.path.exists(output_filename + ".txt")


def test_obabel_parser():
    ligand = DockedLigand(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_obabel.sdf"
        )
    )
    reference = DockedLigand(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_obabel.sdf"
        )
    )
    output_filename = os.path.join(
        context_protein_6Oav["write_dir"], "6oav_m3a_rmsd_obabel"
    )
    rmsd = RMSD(ligand, reference, context)
    res = rmsd.rmsd_obabel(output_filename=output_filename)
    values = obabel_result_parser(res)
    print(values)
    # assert isinstance(values, float)
    assert values == 0.0
