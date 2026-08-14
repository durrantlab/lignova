# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""testing for yaml class to write configuration files."""

import os
import re
from typing import Any

import pytest
import yaml

from lignova.yaml.config import YamlConfig

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# ---------- Fixtures & helpers ----------


@pytest.fixture
def sample_dict() -> dict[str, Any]:
    """Provide a sample dictionary for testing."""
    return {
        "pdb2pqr": {
            "mandatory_options": {
                "ff": "PARSE",
            },
            "general_options": {
                "keep-chain": True,
                "assign-only": False,
            },
            "pka_options": {
                "titration-state-method": "propka",
                "with-ph": 7.0,
            },
            "propka_options": {
                "reference": "neutral",
                "chain": "A",
                "pH": 7.0,
                "log-level": "INFO",
            },
        }
    }


yaml_file_path = "tests/tmp/test_config.yaml"


@pytest.fixture
def yaml_file(tmp_path: pytest.TempPathFactory, sample_dict: dict[str, Any]) -> str:
    """Create a real YAML file to simulate an existing config on disk."""
    p = tmp_path / "config.yaml"
    with open(p, "w", encoding="utf-8") as f:
        # Write sample dict to yaml file without yaml package
        f.write(str(sample_dict).replace("'", '"'))
    return p


def read_yaml(file_path: str) -> str:
    """Helper to read a YAML file into a dictionary."""
    with open(file_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


# ---------- Tests ----------
def test_init(yaml_file, sample_dict):
    """Test initialization of YamlConfig with existing file."""
    cfg = YamlConfig(str(yaml_file))
    assert cfg.file_path == str(yaml_file)
    assert cfg.data_dict == sample_dict


def test_init_missing_param(tmp_path, sample_dict):
    r"""Test initialization creates missing file from provided dictionary."""
    target = tmp_path / "new_config.yaml"
    assert not target.exists()
    cfg = YamlConfig(file_path=target, data_dict=sample_dict)
    assert target.exists()
    assert read_yaml(target) == sample_dict
    assert cfg.file_path == target


def test_init_no_inputs(tmp_path):
    r"""Test initialization raises error when no inputs provided."""
    target = tmp_path / "missing.yaml"
    with pytest.raises(
        ValueError, match="Either file_path or dictionary must be provided"
    ):
        YamlConfig(str(target))


def test_read_config_returns_dict(yaml_file, sample_dict):
    r"""Test reading configuration returns correct dictionary."""
    cfg = YamlConfig(str(yaml_file))
    data = cfg.read_config()
    assert isinstance(data, dict)
    assert data == sample_dict


def test_write_config_overwrites_file(yaml_file, sample_dict):
    r"""Test writing configuration overwrites existing file."""
    cfg = YamlConfig(str(yaml_file))

    cfg.write_config(sample_dict)
    assert cfg.read_config() == sample_dict

    # Now overwrite with a different values
    new_data = {"a": 1, "b": [1, 2, 3]}
    cfg.write_config(new_data)
    assert cfg.read_config() == new_data


def test_update_config(yaml_file, sample_dict):
    cfg = YamlConfig(str(yaml_file))
    cfg.write_config(sample_dict)
    cfg.update_config(
        updates={
            "pdb2pqr": {"general_options": {"keep-chain": False}},
            "new_key": "added",
        }
    )
    expected = {
        "pdb2pqr": {"general_options": {"keep-chain": False}},
        "new_key": "added",
    }
    assert cfg.read_config() == expected


def test_update_config_nested(yaml_file, sample_dict):
    cfg = YamlConfig(str(yaml_file))
    # add nested update test
    cfg.write_config(sample_dict)  # reset
    cfg.update_config(
        updates={"keep-chain": False},
        nested=True,
        parent_key="pdb2pqr",
    )
    expected_nested = sample_dict.copy()
    expected_nested["pdb2pqr"]["keep-chain"] = False
    assert cfg.read_config() == expected_nested


def test_delete_key(yaml_file, sample_dict):
    cfg = YamlConfig(str(yaml_file))
    cfg.write_config(sample_dict)

    before = cfg.read_config()
    assert "pdb2pqr" in before

    cfg.delete_key("pdb2pqr")
    after = cfg.read_config()

    assert "pdb2pqr" not in after


def test_to_cli(yaml_file, sample_dict):
    cfg = YamlConfig(str(yaml_file))
    cfg.write_config(sample_dict)

    cli_data = cfg.to_cli()
    cli_str = " ".join(cli_data)
    expect_pattern = r"--\w+ [^\s]+"
    matches = re.findall(expect_pattern, cli_str)
    expected_data = " --ff PARSE --keep-chain --titration-state-method propka --with-ph 7.0 --reference neutral --chain A --pH 7.0 --log-level INFO"
    assert matches
    # Check that certain expected substrings are in the CLI string
    assert "--ff PARSE" in cli_str
    assert "--keep-chain " in cli_str
    assert "--titration-state-method propka" in cli_str
    assert "--reference neutral" in cli_str
    assert cli_str.strip() == expected_data.strip()
