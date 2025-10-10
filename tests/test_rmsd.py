r"""Test for the RMSD calculation class."""

import os

import numpy as np
from loguru import logger

from lignova.analysis.rmsdclass.mda import mdaRMSD
from lignova.analysis.rmsdclass.obabel import obabelRMSD
from lignova.analysis.rmsdclass.spyrmsd import spyrmsdRMSD
from lignova.analysis.utils import obabel_convert
from lignova.docking.contexts import GlideContext
from lignova.structure.ligand import DockedLigand
from lignova.structure.protein import Protein

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav_B_M3A.sdf",
    "target_path": "./files/6oav/6oav_m3a_top_pose_pv.sdf",
    "write_dir": "./tmp/6oav",
}


def test_prep_dirs():
    r"""Prepare directories."""
    os.makedirs(context_protein_6Oav["write_dir"], exist_ok=True)


context = GlideContext.get_current()


def test_rmsd_spyrmsd():
    r"""Test RMSD calculation using SpyRMSD."""
    logger.debug(context_protein_6Oav["file_path"])
    reference = Protein(file_path=context_protein_6Oav["file_path"])
    ligand = DockedLigand(file_path=context_protein_6Oav["target_path"])
    spy_obj = spyrmsdRMSD(
        reference=reference,
        target=ligand,
        context=context,
    )
    rmsd_val = spy_obj.calculate()
    assert isinstance(rmsd_val, float)
    assert np.isclose(rmsd_val, 0.58534)


def test_rmsd_obabel():
    r"""Test RMSD calculation using OpenBabel."""
    logger.debug(context_protein_6Oav["file_path"])
    reference = Protein(file_path=context_protein_6Oav["file_path"])
    ligand = DockedLigand(file_path=context_protein_6Oav["target_path"])
    obabel_obj = obabelRMSD(
        reference=reference,
        target=ligand,
        context=context,
    )
    rmsd_val = obabel_obj.calculate()
    assert isinstance(rmsd_val, float)
    assert np.isclose(rmsd_val, 0.58534)


def test_rmsd_mda():
    r"""Test RMSD calculation using MDAnalysis."""
    logger.debug(context_protein_6Oav["file_path"])
    # Convert files to pdb format for OpenBabel
    obabel_convert(
        test_file=context_protein_6Oav["file_path"],
        output_filename=os.path.join(
            context_protein_6Oav["write_dir"], "6oav_B_M3A.pdb"
        ),
    )
    obabel_convert(
        test_file=context_protein_6Oav["target_path"],
        output_filename=os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv.pdb"
        ),
    )
    reference = Protein(
        file_path=os.path.join(context_protein_6Oav["write_dir"], "6oav_B_M3A.pdb")
    )
    ligand = DockedLigand(
        file_path=os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv.pdb"
        )
    )
    mda_obj = mdaRMSD(
        reference=reference,
        target=ligand,
        context=context,
    )
    rmsd_val = mda_obj.calculate()
    assert isinstance(rmsd_val, list)
    assert np.isclose(rmsd_val[0], 1.0, atol=0.1)
