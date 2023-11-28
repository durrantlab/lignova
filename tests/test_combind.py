import csv
import glob
import os
import shutil

import pytest

from lignova.docking.combind import Combind
from lignova.docking.contexts.combind import CombindContext

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
    "docking_results_path": "/home/mma121/PubChem_small/try_schrodinger/6oav_validation/glide-dock_SP_1/glide-dock_SP_1_pv.maegz",
}


def prep_dirs():
    os.makedirs(context_protein_6Oav["write_dir"])


if not os.path.exists(context_protein_6Oav["write_dir"]):
    prep_dirs()


def test_combindcontext():
    combind = CombindContext.get_current()
    assert combind.command == "/home/mma121/PubChem_small/combind"
    assert combind.work_dir == context_protein_6Oav["write_dir"]
    assert combind.schrodinger == os.environ["SCHRODINGER"]
    assert combind.schrodinger_env == "./tmp/6oav/schrodinger.ve"


def test_featurize():
    context = CombindContext.get_current()
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    combind.featurize(
        docking_filepath=context_protein_6Oav["docking_results_path"],
        file_name="6oav_m3a",
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_features")
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_features.log")
    )


def test_select_pose():
    context = CombindContext.get_current()
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    combind.select_pose(
        file_name="6oav_m3a",
        features_dir=os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_features"
        ),
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a.csv")
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_pose_selection.log")
    )


def test_get_3d_top_pose():
    context = CombindContext.get_current()
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    combind.get_3d_top_pose(
        docking_filepath=context_protein_6Oav["docking_results_path"],
        combind_csv=os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a.csv"),
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_pv.maegz")
    )


def test_compute_combind_score():
    context = CombindContext.get_current()
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    combind.compute_combind_score(
        features_dir=os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_features"
        ),
        filename="6oav_m3a",
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_screen.npy")
    )


def test_apply_combind_score():
    context = CombindContext.get_current()
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    combind.apply_combind_score(
        docking_filepath=context_protein_6Oav["docking_results_path"],
        combind_score_file=os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_screen.npy"
        ),
        output_filename="6oav_m3a",
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_combind_score.maegz")
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_combind_score.log")
    )
