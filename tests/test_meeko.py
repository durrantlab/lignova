# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""testing for the yaml class that writes Meeko configuration files."""

import copy
import os
from typing import Any

import pytest

from lignova.yaml.meeko_config import MeekoConfig

os.chdir(os.path.dirname(os.path.realpath(__file__)))


@pytest.fixture
def defaults() -> dict[str, Any]:
    """Provide a fresh copy of the default Meeko configuration."""
    return copy.deepcopy(MeekoConfig.declare_defaults(MeekoConfig))


@pytest.fixture
def pqr_file() -> str:
    """Path to PQR file for testing."""
    return os.path.join(os.path.dirname(__file__), "files", "receptor.pqr")


@pytest.fixture
def pdb_file(tmp_path) -> str:
    """Create a tmp PDB file on disk."""
    p = tmp_path / "receptor.pdb"
    p.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000\n")
    return str(p)


@pytest.fixture
def ligand_file(tmp_path) -> str:
    """Create a tmp SDF file for box_enveloping."""
    p = tmp_path / "ligand.sdf"
    p.write_text("$$$$\n")
    return str(p)


def build(
    tmp_path, defaults: dict[str, Any], **sections: dict[str, Any]
) -> MeekoConfig:
    """Build a MeekoConfig from the defaults with the given section overrides.
    Args:
        tmp_path: Pytest temporary directory.
        defaults : Default configuration to start from.
        sections : Per-section overrides, e.g. box={"box_size": [20, 20, 20]}.
    """
    config = copy.deepcopy(defaults)
    for section, params in sections.items():
        config["meeko"][section].update(params)
    return MeekoConfig(str(tmp_path / "meeko_config.yaml"), data_dict=config)


def test_default_behavior(tmp_path):
    """Test a default config file is created when none exists."""
    target = tmp_path / "meeko_config.yaml"
    assert not target.exists()
    cfg = MeekoConfig(str(target))
    assert target.exists()
    assert set(cfg.data_dict["meeko"]) == {
        "input_output",
        "receptor_perception",
        "box",
        "reactive",
    }
    cfg.update_config({"meeko": {"input_output": {"output_basename": "mine"}}})
    io = cfg.data_dict["meeko"]["input_output"]
    assert io["output_basename"] == "mine"


def test_input_rules(tmp_path, defaults, pqr_file, pdb_file):
    """Test that only one reader is allowed."""
    with pytest.raises(ValueError, match="Only one input may be given"):
        build(
            tmp_path,
            defaults,
            input_output={"read_pqr": pqr_file, "read_pdb": pdb_file},
        )
    with pytest.raises(FileNotFoundError, match="does not exist"):
        build(tmp_path, defaults, input_output={"read_pqr": str(tmp_path / "nope.pqr")})

    cfg = build(tmp_path, defaults)
    assert cfg.data_dict["meeko"]["input_output"]["read_pqr"] is None


def test_charge_model(tmp_path, defaults, pdb_file, pqr_file):
    """Test charge model rules"""
    with pytest.raises(ValueError, match="Invalid charge model"):
        build(tmp_path, defaults, receptor_perception={"charge_model": "mulliken"})

    with pytest.raises(ValueError, match="requires 'read_pqr'"):
        defaults["meeko"]["receptor_perception"]["charge_model"] = "read"
        build(tmp_path, defaults, input_output={"read_pdb": pdb_file})

    defaults = MeekoConfig.declare_defaults(MeekoConfig)
    cfg = build(tmp_path, defaults, input_output={"read_pqr": pqr_file})
    assert cfg.data_dict["meeko"]["receptor_perception"]["charge_model"] == "gasteiger"


def test_mk_config(tmp_path, defaults):
    """Test mk_config rules"""
    bad = str(tmp_path / "conf.yaml")
    open(bad, "w").close()
    with pytest.raises(ValueError, match="must be a .json file"):
        build(tmp_path, defaults, receptor_perception={"mk_config": bad})

    with pytest.raises(FileNotFoundError, match="does not exist"):
        build(
            tmp_path,
            defaults,
            receptor_perception={"mk_config": str(tmp_path / "nope.json")},
        )


def test_box_size_rules(tmp_path, defaults, ligand_file):
    """Test a box_size rules."""
    with pytest.raises(ValueError, match="exactly 3 values"):
        build(tmp_path, defaults, box={"box_size": [20.0, 20.0]})
    with pytest.raises(ValueError, match="must be floats or ints"):
        build(tmp_path, defaults, box={"box_size": [20.0, "20", 20.0]})
    with pytest.raises(ValueError, match="must be positive"):
        build(tmp_path, defaults, box={"box_size": [20.0, 0.0, 20.0]})
    with pytest.raises(TypeError, match="must be a list or None"):
        build(tmp_path, defaults, box={"box_size": 20.0})
    with pytest.raises(ValueError, match="cannot be combined"):
        build(
            tmp_path,
            defaults,
            box={"box_enveloping": ligand_file, "box_size": [20.0, 20.0, 20.0]},
        )


def test_box_center_off_reactive(tmp_path, defaults):
    with pytest.raises(ValueError, match="cannot be combined with 'box_center'"):
        build(
            tmp_path,
            defaults,
            box={
                "box_center": [0.0, 0.0, 0.0],
                "box_center_off_reactive_res": True,
            },
            reactive={"reactive_flexres": "A:42"},
        )
    residues = ",".join(f"A:{n}" for n in range(9))
    with pytest.raises(ValueError, match="At most 8 reactive"):
        build(tmp_path, defaults, reactive={"reactive_flexres": residues})
    residues = ",".join(f"A:{n}" for n in range(8))
    cfg = build(tmp_path, defaults, reactive={"reactive_flexres": residues})
    assert cfg.data_dict["meeko"]["reactive"]["reactive_flexres"] == residues


def test_to_cli(tmp_path, defaults, pqr_file):
    """Test the configuration renders the expected mk_prepare_receptor flags."""
    cfg = build(
        tmp_path,
        defaults,
        input_output={"read_pqr": pqr_file, "output_basename": "rec"},
        box={"box_size": [20.0, 20.0, 20.0], "box_center": [1.0, 2.0, 3.0]},
    )
    args = cfg.to_cli()
    cli = " ".join(args)
    assert f"--read_pqr {pqr_file}" in cli
    assert "--output_basename rec" in cli
    assert "--charge_model gasteiger" in cli
    assert "--write_pdbqt" in args
    assert "--box_size" in args
    assert args[args.index("--box_size") + 1 : args.index("--box_size") + 4] == [
        "20.0",
        "20.0",
        "20.0",
    ]
    assert "--allow_bad_res" in cli
    assert "--set_template" not in cli
