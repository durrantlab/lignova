import csv
import glob
import os
import shutil

import numpy as np
import pytest

from lignova.docking.combind import Combind
from lignova.docking.contexts.combind import CombindContext

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))
#    "docking_results_path": "/home/mma121/PubChem_small/try_schrodinger/lignova/tests/tmp/6oav/6oav_A_M3A_lig_docking_pv.maegz",

context_protein_6Oav = {
    "id": "6OAV",
    "file_path": "./files/6oav/6oav.pdb",
    "write_dir": "./tmp/6oav",
    "docking_results_path": "/home/mma121/PubChem_small/try_schrodinger/6oav_validation/glide-dock_SP_1/glide-dock_SP_1_pv.maegz",
    "ref_rmsd_file": "./files/6oav/rmsd1.npy",
    "ref_csv_file": "./files/6oav/6oav_m3a.csv",
    "ref_screen_file": "./files/6oav/6oav_m3a_screen.npy",
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
    assert os.path.basename(combind.schrodinger_env) == "schrodinger.ve"


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
    ref_rmsd = np.load(context_protein_6Oav["ref_rmsd_file"])
    rmsd = np.load(
        os.path.join(
            context_protein_6Oav["write_dir"], "6oav_m3a_features", "rmsd1.npy"
        )
    )
    assert np.allclose(ref_rmsd, rmsd)
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
    with open(context_protein_6Oav["ref_csv_file"], "r") as file:
        reader = csv.DictReader(file)
        ref_csv = list(reader)
    with open(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a.csv"), "r"
    ) as file:
        reader = csv.DictReader(file)
        test_csv = list(reader)
    assert ref_csv[0]["COMBIND_RMSD"] == test_csv[0]["COMBIND_RMSD"]
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a.csv")
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
        extract_filename="6oav_m3a_top_pose",
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_top_pose_pv.maegz")
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
    ref_screen_file = np.load(context_protein_6Oav["ref_screen_file"])
    test_screen_file = np.load(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_screen.npy")
    )
    assert np.allclose(ref_screen_file, test_screen_file)
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_combind_sorted.maegz")
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_combind_sort.log")
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_combind.maegz")
    )


def test_extract_data_csv():
    context = CombindContext.get_current()
    combind = Combind(
        command=context.command,
        work_dir=context.work_dir,
        schrodinger=context.schrodinger,
        schrodinger_env=context.schrodinger_env,
    )
    combind.extract_data_csv(
        docking_file=context_protein_6Oav["docking_results_path"],
        filename="6oav_m3a_filter",
        filter=False,
    )
    assert os.path.exists(
        os.path.join(context_protein_6Oav["write_dir"], "6oav_m3a_filter.csv")
    )
