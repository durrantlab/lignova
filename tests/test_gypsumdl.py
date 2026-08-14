"""Test the YAML configuration handler for Gypsum-DL."""

import os
import re
from copy import deepcopy
from typing import Any, LiteralString

import numpy as np
import pytest
import yaml

from lignova.preparation.gypsumdl import Gypsum
from lignova.yaml.ligprep_config import GypsumDLConfig

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))
tmp_path = "./tmp"


def read_yaml(file_path: str) -> dict[str, Any]:
    """Helper to read a YAML file into a dictionary."""
    with open(file_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def make_cfg_dict() -> dict[str, Any]:
    """Start from defaults to keep tests focused; tweak per-test."""
    tmp = GypsumDLConfig(
        "does_not_exist_gypsum.yaml", data_dict=None
    )  # will write defaults in CWD
    d = tmp.declare_defaults()
    # clean up the accidental file if created
    if os.path.exists("does_not_exist_gypsum.yaml"):
        os.remove("does_not_exist_gypsum.yaml")
    return d


def test_valid_init(tmp_path: str):
    """Defaults are created when no file or data_dict is given."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")
    # No file and no data_dict -> class should create defaults and validate
    cfg = GypsumDLConfig(cfg_path)
    assert os.path.exists(cfg_path)
    config_file = read_yaml(cfg_path)
    assert config_file["gypsum_dl"]["job_specs"]["job_manager"] == "multiprocessing"

    # convert to_cli to string for easy checking
    cli_args = cfg.to_cli()
    cli_str = " ".join(cli_args)
    # Defaults should be reflected in CLI
    assert "--max_variants_per_compound 5" in cli_str
    assert "--num_processors 4" in cli_str


def test_valid_data_population(tmp_path: str):
    """Missing sections/keys are populated with defaults."""
    base = make_cfg_dict()
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")
    # Drop an entire section and a key inside another section
    base["gypsum_dl"].pop("format")
    base["gypsum_dl"]["ph"].pop("pka_precision")

    cfg = GypsumDLConfig(str(cfg_path), data_dict=base)
    cli_str = " ".join(cfg.to_cli())

    # Defaults should be brought back and appear in CLI
    # use_durrant_lab_filters default: True -> boolean flag present
    assert "--use_durrant_lab_filters" in cli_str
    # pka_precision default: 1 -> numeric argument present
    assert "--pka_precision 1" in cli_str


def test_invalid_keys(tmp_path: str):
    """Unknown (non-typed) keys must be string or None."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    # Inject a custom key that is not in the typed sets
    d["gypsum_dl"]["job_specs"]["custom_param"] = 123  # not str/None

    with pytest.raises(TypeError, match="must be a string or None"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d)


def test_valid_bool_param(tmp_path: str):
    """Boolean parameters must be actual booleans."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    # choose a boolean param and make it a string
    d["gypsum_dl"]["format"]["separate_output_files"] = "false"
    with pytest.raises(ValueError, match="must be a boolean"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d)


def test_valid_int_param(tmp_path: str):
    """Integer parameters must be integers."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    d["gypsum_dl"]["job_specs"]["num_processors"] = "4"
    with pytest.raises(ValueError, match="must be an integer"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d)


def test_valid_float_param(tmp_path: str):
    """Float parameters must be floats/ints."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    d["gypsum_dl"]["ph"]["min_ph"] = "six"
    with pytest.raises(ValueError, match="must be a float"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d)


def test_invalid_job_manager(tmp_path: str):
    """job_manager must be one of the allowed values."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    d["gypsum_dl"]["job_specs"]["job_manager"] = "slurm"
    with pytest.raises(ValueError, match="Invalid job_manager"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d)


def test_invalid_num_processors(tmp_path: str):
    """num_processors cannot be zero."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    d["gypsum_dl"]["job_specs"]["num_processors"] = 0
    with pytest.raises(ValueError, match="num_processors must be non-zero"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d)


def test_invalid_thoroughness(tmp_path: str):
    """thoroughness must be >= 1."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    d["gypsum_dl"]["job_specs"]["thoroughness"] = 0
    with pytest.raises(ValueError, match="thoroughness must be >= 1"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d)


def test_invalid_max_compound(tmp_path: str):
    """max_variants_per_compound must be >= 1."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    d["gypsum_dl"]["job_specs"]["max_variants_per_compound"] = 0
    with pytest.raises(ValueError, match="max_variants_per_compound must be >= 1"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d)


def test_ph_bounds_and_order(tmp_path: str):
    """pH values must be within 0–14 and min_ph <= max_ph."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    # min_ph below 0
    d_bad1 = deepcopy(d)
    d_bad1["gypsum_dl"]["ph"]["min_ph"] = -1.0
    with pytest.raises(ValueError, match="min_ph must be between 0 and 14"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d_bad1)

    # max_ph above 14
    d_bad2 = deepcopy(d)
    d_bad2["gypsum_dl"]["ph"]["max_ph"] = 20.0
    with pytest.raises(ValueError, match="max_ph must be between 0 and 14"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d_bad2)

    # min_ph > max_ph
    d_bad3 = deepcopy(d)
    d_bad3["gypsum_dl"]["ph"]["min_ph"] = 9.0
    d_bad3["gypsum_dl"]["ph"]["max_ph"] = 8.0
    with pytest.raises(ValueError, match="min_ph must be <= max_ph"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d_bad3)


def test_pka_precision_positive(tmp_path: str):
    """pka_precision must be > 0."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    d["gypsum_dl"]["ph"]["pka_precision"] = 0.0
    with pytest.raises(ValueError, match="pka_precision must be > 0"):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d)


def test_mpijobs_processor_type(tmp_path: str):
    """For mpi, tasks_per_processor must be an int (not string, etc.)."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    job_specs = d["gypsum_dl"]["job_specs"]
    job_specs["job_manager"] = "mpi"
    job_specs["tasks_per_processor"] = "4"

    with pytest.raises(
        TypeError,
        match="tasks_per_processor must be an integer for mpi job.",
    ):
        _ = GypsumDLConfig(str(cfg_path), data_dict=d)


def test_valid_mpijobs_processor(tmp_path: str):
    """For mpi, a valid tasks_per_processor (int > 1) should pass and appear in CLI."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()

    job_specs = d["gypsum_dl"]["job_specs"]
    job_specs["job_manager"] = "mpi"
    job_specs["tasks_per_processor"] = 8

    cfg = GypsumDLConfig(str(cfg_path), data_dict=d)
    cli_str = " ".join(cfg.to_cli())

    assert "--job_manager mpi" in cli_str
    assert "--tasks_per_processor 8" in cli_str


def test_to_cli(tmp_path: str):
    """Test conversion of config to CLI arguments."""
    cfg_path: LiteralString = os.path.join(tmp_path, "gypsum_dl.yaml")

    d = make_cfg_dict()
    expect_pattern = r"--\w+ [^\s]+"
    # Modify a few parameters away from defaults
    d["gypsum_dl"]["job_specs"]["num_processors"] = 64
    d["gypsum_dl"]["ph"]["min_ph"] = 5.5
    d["gypsum_dl"]["format"]["separate_output_files"] = False
    d["gypsum_dl"]["ph"]["pka_precision"] = 9.0

    cfg = GypsumDLConfig(str(cfg_path), data_dict=d)
    cli_args = cfg.to_cli()
    cli_str = " ".join(cli_args)
    matches = re.findall(expect_pattern, cli_str)
    assert matches
    assert "--num_processors 64" in cli_str
    assert "--max_ph 8" in cli_str
    assert "--min_ph 5.5" in cli_str
    assert "5.5" in cli_str
    assert "64" in cli_str
    # Boolean False should not emit a flag
    assert "--separate_output_files" not in cli_str


def test_gypsum(tmp_path: str):
    """Test the Gypsum class initialization and run method."""
    # Create a dummy .smi file
    smiles_file = os.path.join(tmp_path, "test.smi")
    with open(smiles_file, "w", encoding="utf-8") as f:
        # Use a space or tab between SMILES and ID, *not* "\s"
        f.write(
            "C1CCN(CC1)CCN2C3=CC=CC=C3N=C2NCC4=NC5=CC=CC=C5N4CC6=CC=CC=C6\t1043331\n"
            + "CCCC1=CC(=O)NC(=N1)N2C(=O)C3=C(N2)CCCC3\t135406797\n"
        )
    outfolder = os.path.join(tmp_path, "output")
    # make the config use serial job manager for testing
    d = make_cfg_dict()
    d["gypsum_dl"]["job_specs"]["job_manager"] = "serial"
    d["gypsum_dl"]["job_specs"]["num_processors"] = 1
    config = GypsumDLConfig(
        file_path=os.path.join(tmp_path, "gypsum_config.yaml"), data_dict=d
    )
    gypsum = Gypsum(
        smiles_file=str(smiles_file),
        outfolder=str(outfolder),
        config_obj=config,
    )
    gypsum.run()
    assert gypsum.smiles_file == str(smiles_file)
    assert gypsum.outfolder == str(outfolder)
    test_output_file = os.path.join(outfolder, "gypsum_dl_success.sdf")
    assert os.path.exists(test_output_file)
    with open(test_output_file, "r", encoding="utf-8") as out_f:
        out_data = out_f.read()
    assert "1043331" in out_data
    assert "135406797" in out_data
    assert np.isclose(len(out_data), 50000, rtol=0.05)


# test the gypsum with multiprocessing job manager
def test_gypsum_multi(tmp_path: str):
    """Test the Gypsum class initialization and run method with multiprocessing."""
    # Create a dummy .smi file
    smiles_file = os.path.join(tmp_path, "test_multi.smi")
    with open(smiles_file, "w", encoding="utf-8") as f:
        # Use a space or tab between SMILES and ID, *not* "\s"
        f.write(
            "C1CCN(CC1)CCN2C3=CC=CC=C3N=C2NCC4=NC5=CC=CC=C5N4CC6=CC=CC=C6\t1043331\n"
            + "CCCC1=CC(=O)NC(=N1)N2C(=O)C3=C(N2)CCCC3\t135406797\n"
        )
    outfolder = os.path.join(tmp_path, "output_multi")
    # make the config use multiprocessing job manager for testing
    d = make_cfg_dict()
    d["gypsum_dl"]["job_specs"]["job_manager"] = "multiprocessing"
    d["gypsum_dl"]["job_specs"]["num_processors"] = 2
    config = GypsumDLConfig(
        os.path.join(tmp_path, "gypsum_config_multi.yaml"), data_dict=d
    )
    gypsum = Gypsum(
        smiles_file=str(smiles_file),
        outfolder=str(outfolder),
        config_obj=config,
    )
    gypsum.run()
    assert gypsum.smiles_file == str(smiles_file)
    assert gypsum.outfolder == str(outfolder)
    test_output_file = os.path.join(outfolder, "gypsum_dl_success.sdf")
    assert os.path.exists(test_output_file)
    with open(test_output_file, "r", encoding="utf-8") as out_f:
        out_data = out_f.read()
    assert "1043331" in out_data
    assert "135406797" in out_data
    assert np.isclose(len(out_data), 50000, rtol=0.05)
