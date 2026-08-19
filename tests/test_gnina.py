# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Test the YAML configuration handler for GNINA & proper file handling."""

import os
import shutil
from copy import deepcopy
from typing import Any, LiteralString

import pytest
import yaml
from loguru import logger

from lignova.docking.gnina import GNINA
from lignova.structure.ligand import PreparedLigand
from lignova.structure.protein import PreparedProtein
from lignova.yaml.docking_config import GninaConfig

os.chdir(os.path.dirname(os.path.realpath(__file__)))

PQR_FILE = os.path.join(os.path.dirname(__file__), "files", "receptor.pqr")


def read_yaml(file_path: str) -> dict[str, Any]:
    """Helper to read a YAML file into a dictionary."""
    with open(file_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _make_file(path: str, content: str = "dummy\n") -> str:
    """Create a minimal file and return its path."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    return path


def _make_gnina(tmp_path: str) -> GNINA:
    """Create a GNINA instance with a valid autobox ligand."""
    lig = _make_file(os.path.join(tmp_path, "box.sdf"), "$$$$\n")
    return GNINA(autobox=True, box_ligand=lig)


def _copy_pqr(tmp_path: str) -> str:
    """Copy the test PQR fixture into tmp_path."""
    dest = os.path.join(tmp_path, os.path.basename(PQR_FILE))
    shutil.copy2(PQR_FILE, dest)
    return dest


def make_cfg_dict() -> dict[str, Any]:
    """Start from defaults to keep tests focused; tweak per-test."""
    tmp = GninaConfig("does_not_exist_gnina.yaml", data_dict=None)
    d = tmp.declare_defaults()
    if os.path.exists("does_not_exist_gnina.yaml"):
        os.remove("does_not_exist_gnina.yaml")
    return d


def test_valid_init(tmp_path: str):
    """Defaults are created when no file or data_dict is given."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")

    cfg = GninaConfig(cfg_path)
    assert os.path.exists(cfg_path)

    config_file = read_yaml(cfg_path)
    assert config_file["gnina"]["cnn"]["cnn_scoring"] == "rescore"
    assert config_file["gnina"]["scoring"]["num_mc_saved"] == 50

    cli_args = cfg.to_cli()
    cli_str = " ".join(cli_args)

    assert "--cnn_scoring rescore" in cli_str
    assert "--num_mc_saved 50" in cli_str
    assert "--num_mc_steps" not in cli_str


def test_valid_data_population(tmp_path: str):
    """Missing sections/keys are populated with defaults."""
    base = make_cfg_dict()
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")

    # Drop a key + drop a whole section
    base2 = deepcopy(base)
    base2["gnina"]["scoring"].pop("num_mc_saved")
    base2["gnina"].pop("cnn")

    cfg = GninaConfig(str(cfg_path), data_dict=base2)
    cli_str = " ".join(cfg.to_cli())

    # num_mc_saved default restored
    assert "--num_mc_saved 50" in cli_str
    # cnn section restored; default cnn_scoring appears
    assert "--cnn_scoring rescore" in cli_str


def test_invalid_unknown_key_type(tmp_path: str):
    """Unknown (non-typed) keys must be string or None."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    # Inject unknown key under scoring; your validator treats unknown keys as string-ish
    d["gnina"]["scoring"]["my_custom_param"] = 123  # not str/None
    with pytest.raises(TypeError, match="must be a string or None"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_boolean_param_type(tmp_path: str):
    """Boolean parameters must be actual booleans."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["scoring"]["score_only"] = "false"
    with pytest.raises(ValueError, match="must be a boolean"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_numeric_param_type(tmp_path: str):
    """Numeric parameters must be int/float if provided."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["misc"]["exhaustiveness"] = "8"
    with pytest.raises(ValueError, match="must be numeric"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_int_only_param_type(tmp_path: str):
    """Int-only parameters must be integers (not floats)."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["scoring"]["num_mc_saved"] = 50.0  # should be int
    with pytest.raises(TypeError, match="must be an integer"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_file_path_input_must_exist(tmp_path: str):
    """Input-like file paths must exist if provided."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["input"]["autobox_ligand"] = os.path.join(tmp_path, "nope.sdf")
    with pytest.raises(FileNotFoundError, match="does not exist"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_file_path_output_need_not_exist(tmp_path: str):
    """Output-like file paths are allowed to not exist."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    # should not raise even though doesn't exist
    d["gnina"]["output"]["out"] = os.path.join(tmp_path, "will_be_created.sdf")
    _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_scoring_enum(tmp_path: str):
    """scoring must be allowed or 'default'/None."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["scoring"]["scoring"] = "totally_fake"
    with pytest.raises(ValueError, match="Invalid scoring"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_cnn_scoring_enum(tmp_path: str):
    """cnn_scoring must be one of allowed values."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["cnn"]["cnn_scoring"] = "banana"
    with pytest.raises(ValueError, match="Invalid cnn_scoring"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_approximation_enum(tmp_path: str):
    """approximation must be one of linear|spline|exact."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["scoring"]["approximation"] = "cubic"
    with pytest.raises(ValueError, match="approximation must be one of"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_box_sizes(tmp_path: str):
    """size_x/size_y/size_z must be > 0 if provided."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["docking_region"]["center_x"] = 0.0
    d["gnina"]["docking_region"]["center_y"] = 0.0
    d["gnina"]["docking_region"]["center_z"] = 0.0
    d["gnina"]["docking_region"]["size_x"] = -1.0
    d["gnina"]["docking_region"]["size_y"] = 10.0
    d["gnina"]["docking_region"]["size_z"] = 10.0

    with pytest.raises(ValueError, match="size_x must be > 0"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_exhaustiveness(tmp_path: str):
    """exhaustiveness must be >= 1."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["misc"]["exhaustiveness"] = 0
    with pytest.raises(ValueError, match="exhaustiveness must be >= 1"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_num_modes(tmp_path: str):
    """num_modes must be >= 1."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["misc"]["num_modes"] = 0
    with pytest.raises(ValueError, match="num_modes must be >= 1"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_cpu_zero(tmp_path: str):
    """cpu must be non-zero if provided."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["misc"]["cpu"] = 0
    with pytest.raises(ValueError, match="cpu must be non-zero"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_min_rmsd_filter(tmp_path: str):
    """min_rmsd_filter must be >= 0."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["misc"]["min_rmsd_filter"] = -0.1
    with pytest.raises(ValueError, match="min_rmsd_filter must be >= 0"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_pose_sort_order(tmp_path: str):
    """pose_sort_order must be"CNNscore" OR "CNNaffinity" OR "Energy"}"""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["output"]["pose_sort_order"] = "not_good"
    with pytest.raises(ValueError, match=r"Invalid pose_sort_order"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_user_grid_lambda(tmp_path: str):
    """user_grid_lambda must be >= -1."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["scoring"]["user_grid_lambda"] = -2
    with pytest.raises(ValueError, match="user_grid_lambda must be >= -1"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_covalent_bond_order(tmp_path: str):
    """covalent_bond_order must be >= 1."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["covalent"]["covalent_bond_order"] = 0
    with pytest.raises(ValueError, match="covalent_bond_order must be >= 1"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_autobox_conflicts_with_explicit_box(tmp_path: str):
    """Cannot set autobox_ligand and any center/size simultaneously."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    # Make a real dummy ligand file for existence check
    lig = os.path.join(tmp_path, "lig.sdf")
    with open(lig, "w", encoding="utf-8") as f:
        f.write("$$$$\n")

    d["gnina"]["input"]["autobox_ligand"] = lig
    d["gnina"]["docking_region"]["center_x"] = 0.0

    with pytest.raises(
        ValueError,
        match="choose either \\(center_\\* \\+ size_\\*\\) OR autobox_ligand",
    ):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_incomplete_explicit_box_rejected(tmp_path: str):
    """If explicit box is used, require all center_* and size_*."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["docking_region"]["center_x"] = 0.0
    d["gnina"]["docking_region"]["center_y"] = 0.0
    # center_z missing
    d["gnina"]["docking_region"]["size_x"] = 10.0
    d["gnina"]["docking_region"]["size_y"] = 10.0
    d["gnina"]["docking_region"]["size_z"] = 10.0

    with pytest.raises(ValueError, match="must set all center_\\* and all size_\\*"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_no_lig_incompatible_with_autobox_ligand(tmp_path: str):
    """no_lig=True cannot be used with autobox_ligand."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    lig = os.path.join(tmp_path, "lig.sdf")
    with open(lig, "w", encoding="utf-8") as f:
        f.write("$$$$\n")

    d["gnina"]["docking_region"]["no_lig"] = True
    d["gnina"]["input"]["autobox_ligand"] = lig

    with pytest.raises(
        ValueError, match="no_lig=True is incompatible with autobox_ligand"
    ):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_no_lig_incompatible_with_flexdist_ligand(tmp_path: str):
    """no_lig=True cannot be used with flexdist_ligand."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    # flexdist_ligand is treated as a file path in your validator, so create it to avoid FileNotFoundError
    flex_lig = os.path.join(tmp_path, "lig_for_flexdist.sdf")
    with open(flex_lig, "w", encoding="utf-8") as f:
        f.write("$$$$\n")

    d["gnina"]["docking_region"]["no_lig"] = True
    d["gnina"]["flexibility"]["flexdist_ligand"] = flex_lig

    with pytest.raises(
        ValueError, match="no_lig=True is incompatible with flexdist_ligand"
    ):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_score_only_and_randomize_only_conflict(tmp_path: str):
    """score_only and randomize_only cannot both be True."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["scoring"]["score_only"] = True
    d["gnina"]["scoring"]["randomize_only"] = True

    with pytest.raises(
        ValueError, match="score_only and randomize_only cannot both be True"
    ):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_local_only_and_score_only_conflict(tmp_path: str):
    """local_only and score_only cannot both be True."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["scoring"]["local_only"] = True
    d["gnina"]["scoring"]["score_only"] = True

    with pytest.raises(
        ValueError, match="local_only and score_only cannot both be True"
    ):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_valid_vector_expansion(tmp_path: str):
    """center/size vectors expand into center_x/y/z and size_x/y/z."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["docking_region"]["center"] = [0.0, 1.0, 2.0]
    d["gnina"]["docking_region"]["size"] = [10.0, 11.0, 12.0]

    cfg = GninaConfig(str(cfg_path), data_dict=d)
    region = cfg.data_dict["gnina"]["docking_region"]

    assert region["center_x"] == 0.0
    assert region["center_y"] == 1.0
    assert region["center_z"] == 2.0
    assert region["size_x"] == 10.0
    assert region["size_y"] == 11.0
    assert region["size_z"] == 12.0
    assert "center" not in region
    assert "size" not in region


def test_invalid_vector_expansion_type(tmp_path: str):
    """center must be a 3-item sequence, not a scalar/string."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["docking_region"]["center"] = "0,0,0"

    with pytest.raises(TypeError, match=r"docking_region\.center.*must be a list"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_vector_expansion_len(tmp_path: str):
    """size must have exactly 3 entries."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["docking_region"]["size"] = [10.0, 10.0]  # invalid length

    with pytest.raises(ValueError, match=r"docking_region\.size.* exactly 3 values"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_invalid_vector_expansion_population(tmp_path: str):
    """Vector elements must be numeric (or None) after expansion."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["docking_region"]["center"] = [0.0, "bad", 2.0]

    with pytest.raises(TypeError, match="'center_y' must be numeric or None"):
        _ = GninaConfig(str(cfg_path), data_dict=d)


def test_cpu_less_than_exhaustiveness_warns(tmp_path: str):
    """cpu < exhaustiveness should log a warning (not raise)."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")
    d = make_cfg_dict()

    d["gnina"]["misc"]["cpu"] = 2
    d["gnina"]["misc"]["exhaustiveness"] = 8

    # capture loguru output
    logs: list[str] = []
    sink_id = logger.add(lambda msg: logs.append(msg), level="WARNING")
    try:
        _ = GninaConfig(str(cfg_path), data_dict=d)
    finally:
        logger.remove(sink_id)

    assert any("cpu is less than exhaustiveness" in str(m) for m in logs)


def test_invalid_protein_prep(tmp_path: str):
    """run() must reject a receptor path that does not exist."""
    g = _make_gnina(tmp_path)
    lig = _make_file(os.path.join(tmp_path, "box.sdf"), "$$$$\n")
    pqr = os.path.join(tmp_path, "receptor.pqr")
    cfg_path = os.path.join(tmp_path, "gnina.yaml")
    with pytest.raises(FileNotFoundError, match="Receptor file .*does not exist"):
        g.run(PreparedProtein(pqr), PreparedLigand(lig), cfg_path)


def test_run_command(tmp_path: str):
    """run() with a .pdbqt target returns a single gnina command string."""
    box_lig = _make_file(os.path.join(tmp_path, "box.sdf"), "$$$$\n")
    g = GNINA(autobox=True, box_ligand=box_lig)

    rec = _make_file(os.path.join(tmp_path, "receptor.pdbqt"))
    lig = _make_file(os.path.join(tmp_path, "mol.sdf"), "$$$$\n")
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")

    cmd = g.run(target=rec, ligand=lig, context=cfg_path)

    assert isinstance(cmd, str)
    assert cmd.startswith("gnina ")
    assert "--autobox_ligand" in cmd
    assert box_lig in cmd


def test_multiple_ligands_run(tmp_path: str):
    """run() with multiple ligands returns one command per ligand."""
    box_lig = _make_file(os.path.join(tmp_path, "box.sdf"), "$$$$\n")
    g = GNINA(autobox=True, box_ligand=box_lig)

    rec = _make_file(os.path.join(tmp_path, "receptor.pdbqt"))
    ligs = [
        _make_file(os.path.join(tmp_path, f"mol{i}.sdf"), "$$$$\n") for i in range(3)
    ]
    cfg_path: LiteralString = os.path.join(tmp_path, "gnina.yaml")

    cmd = g.run(target=rec, ligand=ligs, context=cfg_path)

    assert isinstance(cmd, list)
    assert len(cmd) == 3
    assert all(cmd.startswith("gnina ") for cmd in cmd)
