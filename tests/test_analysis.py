r"""Test the docking analysis module for calculating RMSD."""

import os

import pytest
from loguru import logger

# Skip this entire module as it had Schrödinger-dependent code
if os.getenv("SCHRODINGER") is None:
    pytest.skip(
        "Schrödinger is not installed or the $SCHRODINGER environment variable is not set. Skipping Glide tests.",
        allow_module_level=True,
    )

from lignova.analysis.rmsd import RMSD
from lignova.analysis.utils import (
    interconvert_mae_sdf,
    obabel_convert,
    obabel_result_parser,
)
from lignova.docking.contexts import GlideContext
from lignova.docking.utils import manipulate_complexes
from lignova.structure.editing import write_mda_universe
from lignova.structure.ligand import DockedLigand
from lignova.structure.protein import Protein
from lignova.structure.utils import separate_protein_ligand

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
    "docked_ligand_filepath": "/home/mma121/PubChem_small/try_schrodinger/6oav_validation/6oav_m3a_top_pose_pv.maegz",
    "complexes_filepath": "/home/mma121/PubChem_small/try_schrodinger/6oav_validation/6oav_m3a_top_pose_merge.maegz",
}


def prep_dirs():
    r"""Prepare directories for writing files."""
    os.makedirs(context_protein_6Oav["write_dir"])


if not os.path.exists(context_protein_6Oav["write_dir"]):
    prep_dirs()


if not os.path.exists(context_protein_6Oav["write_dir"]):
    prep_dirs()


context = GlideContext.get_current()
context.write_dir = context_protein_6Oav["write_dir"]
context.set_current(context)


# get the complex for the docking file
@pytest.mark.skipif(
    condition=os.getenv("SCHRODINGER") is None, reason="Schrödinger not installed"
)
def test_manipulate_complexes():
    r"""Test the manipulation of complexes as it is used in the rmsd calculation functions."""
    manipulate_complexes(
        context_protein_6Oav["docked_ligand_filepath"],
        context=context,
        mode="merge",
        outfile_name="6oav_m3a_top_pose_merge.maegz",
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_merge.maegz")
    )


@pytest.mark.skipif(
    condition=os.getenv("SCHRODINGER") is None, reason="Schrödinger not installed"
)
def test_interconvert_mae_sdf():
    r"""Test the interconversion of MAE to SDF."""
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


def test_obabel_convert():
    r"""Test the conversion of file formats using Open Babel."""
    output_filename = os.path.join(context_protein_6Oav["write_dir"], "6oav.sdf")
    obabel_convert(context_protein_6Oav["file_path"], output_filename)
    assert os.path.exists(os.path.join(context_protein_6Oav["write_dir"], "6oav.sdf"))


def test_rmsd_mda():
    r"""Test the RMSD calculation using MDAnalysis."""
    ligand = DockedLigand(context_protein_6Oav["docked_ligand_filepath"])
    reference = Protein(context_protein_6Oav["file_path"])
    # Create an instance of the RMSD class
    rmsd = RMSD(ligand, reference, context)
    logger.debug(context.write_dir)
    # Call the rmsd_mda function
    value = rmsd.rmsd_mda()
    assert isinstance(value[0], float)
    assert round(value[0], 4) == 1.1226


def test_rmsd_obabel():
    r"""Test the RMSD calculation using Open Babel."""
    ligand = DockedLigand(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_complexes.pdb"
        )
    )

    reference = DockedLigand(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_complexes.pdb"
        )
    )
    # Create an instance of the RMSD class
    rmsd = RMSD(ligand, reference, context)
    # Call the rmsd_obabel function
    output_filename = os.path.join(
        context_protein_6Oav["write_dir"], "6oav_m3a_rmsd_obabel"
    )
    _ = rmsd.rmsd_obabel(output_filename=output_filename, save=True)
    assert os.path.exists(output_filename + ".txt")


def test_obabel_parser():
    r"""Test the Open Babel result parser."""
    ligand = DockedLigand(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_complexes.pdb"
        )
    )
    reference = DockedLigand(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_complexes.pdb"
        )
    )
    output_filename = os.path.join(
        context_protein_6Oav["write_dir"], "6oav_m3a_rmsd_obabel.txt"
    )
    rmsd = RMSD(ligand, reference, context)
    res = rmsd.rmsd_obabel(output_filename=output_filename)
    values = obabel_result_parser(res)
    print(values)
    assert isinstance(values, dict)
    assert list(values.values()) == [[0.0]]


@pytest.mark.skipif(
    condition=os.getenv("SCHRODINGER") is None, reason="Schrödinger not installed"
)
def test_manipulate_complexes_lig_sep():
    r"""Test the manipulation of complexes to split ligand."""
    manipulate_complexes(
        context_protein_6Oav["complexes_filepath"],
        context=context,
        mode="split_ligand",
        outfile_name="6oav_m3a_top_pose_split_lig.maegz",
    )
    assert os.path.exists(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_split_lig.maegz"
        )
    )


"""
def test_spyrmsd_calculation():
    r""Test the RMSD calculation using spyrmsd.""
    ligand = DockedLigand(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_complexes.pdb"
        )
    )

    reference = DockedLigand(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv_complexes.pdb"
        )
    )
    # Create an instance of the RMSD class
    rmsd = RMSD(ligand, reference, context)
    # Call the rmsd_obabel function
    output_filename = os.path.join(
        context_protein_6Oav["write_dir"], "6oav_m3a_rmsd_spyrmsd"
    )
    _ = rmsd.symmetry_rmsd(output_filename=output_filename, save=True)
    assert os.path.exists(output_filename + ".txt")
"""
