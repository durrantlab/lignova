# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Tests for receptor to PDBQT conversion with Meeko and MGLTools."""

import os
import shutil

import pytest

from lignova.preparation.meeko import Meeko
from lignova.preparation.mgltools import MglTools
from lignova.yaml.meeko_config import MeekoConfig
from lignova.yaml.mgltools_config import MglToolsConfig

os.chdir(os.path.dirname(os.path.realpath(__file__)))

_FILES_DIR = os.path.join(os.path.dirname(__file__), "files", "6oav")
PDB_FILE = os.path.join(_FILES_DIR, "6oav.pdb")
PQR_FILE = os.path.join(_FILES_DIR, "6oav_reference.pqr")
_WIDE_PQR_FILE = os.path.join(os.path.dirname(__file__), "files", "receptor.pqr")

LIGAND_RESNAME = "M3A"


def _mgltools_available() -> bool:
    """Whether the isolated mgltools pixi environment is installed.

    Never triggers an install.
    """
    try:
        from lignova.io import mgltools_env_exists
    except ImportError:
        return False
    try:
        return mgltools_env_exists()
    except Exception:
        return False


requires_mgltools = pytest.mark.skipif(
    not _mgltools_available(),
    reason="mgltools pixi environment is not installed "
    "(run `pixi install -e mgltools` to enable)",
)


def copy_wide_pqr(directory: str) -> str:
    """Copy the four-digit-numbered PQR"""
    dest = os.path.join(directory, os.path.basename(_WIDE_PQR_FILE))
    shutil.copy2(_WIDE_PQR_FILE, dest)
    return dest


def copy_pdb(directory: str) -> str:
    """Copy the 6oav PDB into a working directory."""
    dest = os.path.join(directory, os.path.basename(PDB_FILE))
    shutil.copy2(PDB_FILE, dest)
    return dest


def copy_pqr(directory: str) -> str:
    """Copy the 6oav reference PQR into a working directory."""
    dest = os.path.join(directory, os.path.basename(PQR_FILE))
    shutil.copy2(PQR_FILE, dest)
    return dest


def _build_meeko(directory: str) -> tuple[type[Meeko], str, MeekoConfig]:
    """Construct a Meeko object over the shared PDB."""
    config = MeekoConfig(os.path.join(directory, "meeko_config.yaml"))
    perception = config.data_dict["meeko"]["receptor_perception"]
    perception["charge_model"] = "gasteiger"
    perception["allow_bad_res"] = True
    return Meeko, copy_pdb(directory), config


def _build_mgltools(directory: str) -> tuple[type[MglTools], str, MglToolsConfig]:
    """Construct an MGLTools object over the shared PDB."""
    config = MglToolsConfig(os.path.join(directory, "mgltools_config.yaml"))
    perception = config.data_dict["mgltools"]["receptor_perception"]
    perception["preserve_charges"] = False
    return MglTools, copy_pdb(directory), config


toolS = [
    pytest.param("meeko"),
    pytest.param("mgltools", marks=requires_mgltools),
]
BUILDERS = {"meeko": _build_meeko, "mgltools": _build_mgltools}


@pytest.fixture
def converted(request, tmp_path) -> tuple[str, str]:
    """Run the requested tool and return its input and output paths."""
    tool = request.param
    cls, input_path, config = BUILDERS[tool](str(tmp_path))
    prep = cls(input_path, os.path.join(str(tmp_path), "rec"), config)
    out = prep.run()
    assert out == prep.pdbqt_file
    return input_path, out


def _atom_lines(path: str) -> list[str]:
    """Read the ATOM/HETATM records of a structure file."""
    return [
        line
        for line in open(path).read().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]


def _charges(lines: list[str]) -> list[float]:
    """Extract the partial charges from PDBQT atom lines."""
    return [float(line[66:76]) for line in lines]


@pytest.mark.parametrize("converted", toolS, indirect=True)
def test_valid_pdbqt(converted):
    """Every tool must produce a valid PDBQT"""
    _, out = converted
    assert os.path.exists(out)
    assert _atom_lines(out), "PDBQT contains no atom records"
    untyped = [line for line in _atom_lines(out) if not line[77:79].strip()]
    assert not untyped, f"{len(untyped)} atoms have no AutoDock type"
    for line in _atom_lines(out)[:50]:
        for start, stop in ((30, 38), (38, 46), (46, 54)):
            float(line[start:stop])
    charges = _charges(_atom_lines(out))
    assert charges, "no atoms to check"
    assert any(abs(q) > 1e-6 for q in charges), "every charge is zero"
    assert len({round(q, 3) for q in charges}) > 5, "charges look like a constant"
    names_per_key: dict[tuple[str, str], set[str]] = {}
    for line in _atom_lines(out):
        resseq = line[22:26].strip()
        assert resseq.lstrip("-").isdigit(), f"non-numeric residue number {resseq!r}"
        key = (line[21], resseq)
        names_per_key.setdefault(key, set()).add(line[17:20].strip())

    collapsed = {k: v for k, v in names_per_key.items() if len(v) > 1}
    assert not collapsed, f"residue numbers shared by different residues: {collapsed}"


@pytest.mark.parametrize("converted", toolS, indirect=True)
def test_expected_behavior(converted):
    """The default parameters should produce the expected output file."""
    input_path, out = converted
    src = {line[22:26].strip() for line in _atom_lines(input_path)}
    dst = {line[22:26].strip() for line in _atom_lines(out)}
    assert dst <= src, f"output has residue numbers absent from the input: {dst - src}"
    src_chains = {line[21] for line in _atom_lines(input_path)}
    out_chains = {line[21] for line in _atom_lines(out)}
    assert out_chains <= src_chains, (
        f"output chains {sorted(out_chains)} are not all present in the input "
        f"{sorted(src_chains)}"
    )


@requires_mgltools
@pytest.mark.parametrize(
    ("cleanup", "expect_ligand"),
    [
        ("nphs_lps_waters_nonstdres", False),
        ("nphs_lps_waters", True),
    ],
)
def test_mgltools_ligand_follows_cleanup(tmp_path, cleanup, expect_ligand):
    """Ligand retention must be driven by the cleanup mode.

    Parametrized rather than sequential because the HETATM strip rewrites the
    input file in place, so the two cases cannot share one copy.
    """
    cls, input_path, config = _build_mgltools(str(tmp_path))
    config.update_config(
        {"cleanup": cleanup}, parent_key=("mgltools", "receptor_perception")
    )
    assert (
        LIGAND_RESNAME in open(input_path).read()
    ), f"6oav.pdb has no {LIGAND_RESNAME} to test with"

    out = cls(input_path, os.path.join(str(tmp_path), "rec"), config).run()

    assert (LIGAND_RESNAME in open(out).read()) is expect_ligand


def test_meeko_default_behavior(tmp_path):
    """Meeko should build a template for the ligand rather than dropping it."""
    cls, input_path, config = _build_meeko(str(tmp_path))
    out = cls(input_path, os.path.join(str(tmp_path), "rec"), config).run()

    assert (
        LIGAND_RESNAME in open(out).read()
    ), f"{LIGAND_RESNAME} absent from the Meeko output"


@requires_mgltools
def test_mgltools_pqr_behavior(tmp_path):
    """A four-digit-numbered PQR must survive the chain blank-and-restore."""
    config = MglToolsConfig(os.path.join(str(tmp_path), "mgltools_config.yaml"))
    input_path = copy_wide_pqr(str(tmp_path))
    original = open(input_path).read()
    src_chains = {line[21] for line in _atom_lines(input_path)}

    out = MglTools(input_path, os.path.join(str(tmp_path), "rec"), config).run()

    for line in _atom_lines(out):
        assert (
            line[22:26].strip().isdigit()
        ), f"residue number {line[22:26]!r} was truncated"
    assert {
        line[21] for line in _atom_lines(out)
    } <= src_chains, "chain identifier was not restored"
    backup = input_path.replace(".pqr", "_org.pqr")
    assert os.path.exists(backup), "no backup written before the input was rewritten"
    assert open(backup).read() == original
