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
    assert combind.command == "/home/mma121/PubChem_small/combind/"
    assert combind.work_dir == context_protein_6Oav["write_dir"]
    assert combind.schrodinger == os.environ["SCHRODINGER"]
    assert combind.schrodinger_env == "schrodinger.ve"
