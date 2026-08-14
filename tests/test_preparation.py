"""Test the preparation module for protein-ligand systems pre-docking."""

import os
import re
from copy import deepcopy

import pytest
import yaml

from lignova.preparation.pdb2pqr import PDB2PQR
from lignova.yaml.protonation_config import ProtonationConfig

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

pdbfile = "./files/6oav/6oav.pdb"
tmp_path = "./tmp"


def read_yaml(file_path: str) -> str:
    """Helper to read a YAML file into a dictionary."""
    with open(file_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def make_cfg_dict():
    """Start from defaults to keep tests focused; tweak per-test."""
    tmp = ProtonationConfig(
        "does_not_exist.yaml", data_dict=None
    )  # will write defaults in CWD
    d = tmp.declare_defaults()
    # clean up the accidental file if created
    if os.path.exists("does_not_exist.yaml"):
        os.remove("does_not_exist.yaml")
    return d


def test_valid_init(tmp_path):
    """Test that defaults are created when no file or data_dict is given."""
    cfg_path = tmp_path / "protonation.yaml"
    # No file and no data_dict -> class should create defaults and validate
    pc = ProtonationConfig(str(cfg_path))
    assert cfg_path.exists(), "Expected defaults to be written to disk"
    on_disk = read_yaml(cfg_path)
    assert on_disk["pdb2pqr"]["general"]["ff"] == "PARSE"


def test__valid_data_population(tmp_path):
    """Test that missing sections/keys are populated with defaults."""
    cfg_path = tmp_path / "protonation.yaml"
    base = make_cfg_dict()
    base["pdb2pqr"].pop("propka")
    base["pdb2pqr"]["general"].pop("include-header")

    pc = ProtonationConfig(str(cfg_path), data_dict=base)
    merged = read_yaml(cfg_path)
    assert "propka" in merged["pdb2pqr"]
    assert "include-header" in merged["pdb2pqr"]["general"]
    assert isinstance(merged["pdb2pqr"]["propka"]["window"], list)


def test_valid_ff(tmp_path):
    """Test that force field parameter is valid."""
    cfg_path = tmp_path / "protonation.yaml"
    d = make_cfg_dict()
    d["pdb2pqr"]["general"]["ff"] = "BADFF"
    with pytest.raises(ValueError, match="Invalid force field"):
        ProtonationConfig(str(cfg_path), data_dict=d)


def test_valid_titration_method(tmp_path):
    """Test that titration method parameter is valid."""
    cfg_path = tmp_path / "protonation.yaml"
    d = make_cfg_dict()
    d["pdb2pqr"]["pka"]["titration-state-method"] = "something-else"
    with pytest.raises(ValueError, match="Invalid titration method"):
        ProtonationConfig(str(cfg_path), data_dict=d)


def test_valid_log_level(tmp_path):
    """Test that PropKa log level parameter is valid."""
    cfg_path = tmp_path / "protonation.yaml"
    d = make_cfg_dict()
    d["pdb2pqr"]["propka"]["log-level"] = "TRACE"
    with pytest.raises(ValueError, match="Invalid log level"):
        ProtonationConfig(str(cfg_path), data_dict=d)


def test_valid_mutator(tmp_path):
    """Test that PropKa mutator parameter is valid."""
    cfg_path = tmp_path / "protonation.yaml"
    d = make_cfg_dict()
    d["pdb2pqr"]["propka"]["mutator"] = "weird"
    with pytest.raises(ValueError, match="Invalid mutator"):
        ProtonationConfig(str(cfg_path), data_dict=d)


def test_valid_reference(tmp_path):
    """Test that PropKa reference parameter is valid."""
    cfg_path = tmp_path / "protonation.yaml"
    d = make_cfg_dict()
    d["pdb2pqr"]["propka"]["reference"] = "basic"
    with pytest.raises(ValueError, match="Invalid PropKa reference"):
        ProtonationConfig(str(cfg_path), data_dict=d)


def test__valid_boolean(tmp_path):
    """Test that boolean parameters are indeed booleans."""
    cfg_path = tmp_path / "protonation.yaml"
    d = make_cfg_dict()
    # choose a boolean param and make it a string
    d["pdb2pqr"]["general"]["whitespace"] = "false"
    with pytest.raises(ValueError, match="must be a boolean"):
        ProtonationConfig(str(cfg_path), data_dict=d)


def test_valid_ph_checks(tmp_path):
    """Test that pH values are within acceptable ranges."""
    cfg_path = tmp_path / "protonation.yaml"
    d = make_cfg_dict()

    d_bad1 = deepcopy(d)
    d_bad1["pdb2pqr"]["pka"]["with-ph"] = 20
    with pytest.raises(ValueError, match="between 0 and 14"):
        ProtonationConfig(str(cfg_path), data_dict=d_bad1)

    d_bad2 = deepcopy(d)
    d_bad2["pdb2pqr"]["propka"]["pH"] = -1
    with pytest.raises(ValueError, match="PropKa pH must be between 0 and 14"):
        ProtonationConfig(str(cfg_path), data_dict=d_bad2)


def test_valid_parse_only_neutral_flags(tmp_path):
    """test that neutraln and neutralc can only be True if ff is 'PARSE'."""
    cfg_path = tmp_path / "protonation.yaml"
    d = make_cfg_dict()
    d["pdb2pqr"]["general"]["ff"] = "AMBER"
    d["pdb2pqr"]["general"]["neutraln"] = True
    d["pdb2pqr"]["general"]["neutralc"] = False
    with pytest.raises(
        ValueError, match="neutraln and neutralc can only be True if ff is 'PARSE'"
    ):
        ProtonationConfig(str(cfg_path), data_dict=d)


def test_to_cli(tmp_path):
    """Test conversion of config to command-line arguments."""
    cfg_path = tmp_path / "protonation.yaml"
    d = make_cfg_dict()
    expect_pattern = r"--\w+ [^\s]+"
    d["pdb2pqr"]["general"]["clean"] = True
    d["pdb2pqr"]["propka"]["chain"] = "A"
    pc = ProtonationConfig(str(cfg_path), data_dict=d)
    cli_args = pc.to_cli()
    cli_str = " ".join(cli_args)
    matches = re.findall(expect_pattern, cli_str)
    assert matches
    assert "--ff AMBER" not in cli_str
    assert "--clean" in cli_str
    assert "--noopt True" not in cli_str
    assert "--chain A" in cli_str


def test_pdb2pqr_run(tmp_path):
    cfg_path = tmp_path / "protonation.yaml"
    d = make_cfg_dict()
    pc = ProtonationConfig(str(cfg_path), data_dict=d)
    output_pqr = os.path.join(tmp_path, "tmp.pqr")
    prep = PDB2PQR(pdb_file=pdbfile, outfile=output_pqr, config_obj=pc)
    prep.run()
    reference = "./files/6oav/6oav_reference.pqr"
    # read the output and reference files to compare
    with open(reference, "r") as ref_f:
        ref_data = ref_f.readlines()
    assert os.path.exists(output_pqr), "Expected output PQR file to be created."
    assert os.path.getsize(output_pqr) > 0, "Expected output PQR file to be non-empty."
    with open(output_pqr, "r") as out_f:
        out_data = out_f.readlines()
    assert len(out_data) == len(
        ref_data
    ), "Output PQR file line count does not match reference."
    assert out_data == ref_data, "Output PQR file content does not match reference."
