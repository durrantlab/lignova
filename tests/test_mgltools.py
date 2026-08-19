# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""testing for the yaml class that writes MGLTools configuration files."""

import copy
import os
from typing import Any

import numpy as np
import pytest
from loguru import logger

from lignova.preparation.mgltools import (
    MglTools,
    format_pqr_atom_line,
    parse_pqr_atom_line,
    strip_hetatm_lines,
)
from lignova.yaml.mgltools_config import MglToolsConfig

os.chdir(os.path.dirname(os.path.realpath(__file__)))


def _mgltools_available() -> bool:
    """Whether the isolated mgltools pixi environment is already installed."""
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


@pytest.fixture
def defaults() -> dict[str, Any]:
    """Provide a fresh copy of the default MGLTools configuration."""
    return copy.deepcopy(MglToolsConfig.declare_defaults(MglToolsConfig))


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
def bad_extension_file(tmp_path) -> str:
    """Create a tmp file with an extension MGLTools does not accept."""
    p = tmp_path / "receptor.xyz"
    p.write_text("1\ncomment\nN 0.0 0.0 0.0\n")
    return str(p)


@pytest.fixture
def squashed_pqr(tmp_path) -> str:
    """A PQR whose coordinate columns run together, as PDB2PQR can emit."""
    p = tmp_path / "squashed.pqr"
    p.write_text(
        "ATOM      1  N   ALA A   1     -18.713-100.168  12.345 -0.4000 1.5000\n"
        "ATOM      2  CA  ALA A   1       1.000   2.000   3.000 -0.0000 2.0000\n"
        "HETATM    3  O   HOH A   2       5.000   6.000   7.000 -0.8340 1.4000\n"
        "HETATM    4  C1  LIG A 900       8.000   9.000  10.000  0.1000 1.7000\n"
        "TER\n"
        "END\n"
    )
    return str(p)


@pytest.fixture
def hetatm_pdb(tmp_path) -> str:
    """A PDB with waters and a non-standard residue, in valid fixed columns."""
    p = tmp_path / "receptor_het.pdb"
    p.write_text(
        "ATOM   4353  N   SER A 807       8.361 -19.707   7.711  1.00 51.01           N\n"
        "HETATM 4400  O   HOH A 900       1.000   2.000   3.000  1.00 30.00           O\n"
        "HETATM 4401  C1  M3A A 901       4.000   5.000   6.000  1.00 25.00           C\n"
        "END\n"
    )
    return str(p)


def build(
    tmp_path, defaults: dict[str, Any], **sections: dict[str, Any]
) -> MglToolsConfig:
    """Build a MglToolsConfig from the defaults with the given section overrides."""
    config = copy.deepcopy(defaults)
    for section, params in sections.items():
        config["mgltools"][section].update(params)
    return MglToolsConfig(str(tmp_path / "mgltools_config.yaml"), data_dict=config)


def test_default_behavior(tmp_path, defaults):
    """Test a default config behavior."""
    target = tmp_path / "mgltools_config.yaml"
    assert not target.exists()
    cfg = MglToolsConfig(str(target))
    assert target.exists()
    assert set(cfg.data_dict["mgltools"]) == {
        "input_output",
        "receptor_perception",
    }
    cfg.update_config({"mgltools": {"receptor_perception": {"repair": "bonds"}}})
    perception = cfg.data_dict["mgltools"]["receptor_perception"]
    assert perception["repair"] == "bonds"
    cfg = build(tmp_path, defaults)
    perception = cfg.data_dict["mgltools"]["receptor_perception"]
    assert perception["repair"] == "checkhydrogens"
    assert perception["cleanup"] == "nphs_lps_waters_nonstdres"
    assert perception["preserve_charges"] is True


def test_defaults_rule(tmp_path, defaults):
    """Test all default values are correctly initialized."""
    assert isinstance(MglToolsConfig._ALLOWED_INPUT_EXTENSIONS, frozenset)
    assert (
        MglTools._ALLOWED_INPUT_EXTENSIONS is MglToolsConfig._ALLOWED_INPUT_EXTENSIONS
    )
    with pytest.raises(ValueError, match="Invalid repair mode"):
        build(tmp_path, defaults, receptor_perception={"repair": "magic"})

    cfg = build(tmp_path, defaults, receptor_perception={"repair": "bonds_hydrogens"})
    assert cfg.data_dict["mgltools"]["receptor_perception"]["repair"] == (
        "bonds_hydrogens"
    )
    with pytest.raises(ValueError, match="Invalid cleanup mode"):
        build(tmp_path, defaults, receptor_perception={"cleanup": "everything"})

    cfg = build(tmp_path, defaults, receptor_perception={"cleanup": "nphs_lps_waters"})
    assert cfg.data_dict["mgltools"]["receptor_perception"]["cleanup"] == (
        "nphs_lps_waters"
    )


def test_receptor_rules(tmp_path, defaults, pqr_file, bad_extension_file):
    """Test receptor path existence and extension are validated."""
    with pytest.raises(FileNotFoundError, match="does not exist"):
        build(
            tmp_path,
            defaults,
            input_output={"receptor": str(tmp_path / "nope.pqr")},
        )
    with pytest.raises(ValueError, match="Invalid receptor extension"):
        build(tmp_path, defaults, input_output={"receptor": bad_extension_file})

    cfg = build(tmp_path, defaults, input_output={"receptor": pqr_file})
    assert cfg.data_dict["mgltools"]["input_output"]["receptor"] == pqr_file


def test_preserve_charges_warns_on_pdb(tmp_path, defaults, pdb_file):
    """-C is meaningless for PDB input, which carries no charge column."""
    messages: list[str] = []
    handler_id = logger.add(messages.append, format="{message}", level="WARNING")
    try:
        build(
            tmp_path,
            defaults,
            input_output={"receptor": pdb_file},
            receptor_perception={"preserve_charges": True},
        )
        assert any("preserve_charges" in m for m in messages)

        messages.clear()
        build(
            tmp_path,
            defaults,
            input_output={"receptor": pdb_file},
            receptor_perception={"preserve_charges": False},
        )
        assert not any("preserve_charges" in m for m in messages)
    finally:
        logger.remove(handler_id)


def test_to_cli(tmp_path, defaults, pqr_file):
    """Test the configuration renders the expected flags."""
    outfile = str(tmp_path / "rec.pdbqt")
    cfg = build(
        tmp_path,
        defaults,
        input_output={"receptor": pqr_file, "outfile": outfile},
    )
    args = cfg.to_cli()

    assert args[args.index("-r") + 1] == pqr_file
    assert args[args.index("-o") + 1] == outfile
    assert args[args.index("-A") + 1] == "checkhydrogens"
    assert args[args.index("-U") + 1] == "nphs_lps_waters_nonstdres"
    assert "-C" in args
    assert not any(a.startswith("--") for a in args)
    cfg = build(
        tmp_path,
        defaults,
        input_output={"receptor": pqr_file, "outfile": str(tmp_path / "rec.pdbqt")},
        receptor_perception={"preserve_charges": False},
    )
    args = cfg.to_cli()
    assert "-C" not in args
    assert "False" not in args


@pytest.mark.parametrize(
    ("line", "chain", "resseq", "icode"),
    [
        # fixed-column PQR
        (
            "ATOM      1  N   PHE A 537      -8.834   3.827 -18.279 -0.7800 1.5000",
            "A",
            537,
            " ",
        ),
        # chain fused to residue number
        (
            "ATOM   4353  N   SER A1086       8.361 -19.707   7.711 -0.4000 1.5000",
            "A",
            1086,
            " ",
        ),
        # trailing insertion code
        (
            "ATOM   4353  N   SER A 807A      8.361 -19.707   7.711 -0.4000 1.5000",
            "A",
            807,
            "A",
        ),
    ],
)
def test_malformed_pqr(line, chain, resseq, icode):
    """Test mgltools handle mPQR parsing."""
    atom = parse_pqr_atom_line(line)
    assert atom is not None
    assert atom.chain == chain
    assert atom.resseq == resseq
    assert atom.icode == icode


@pytest.mark.parametrize(
    (
        "line",
        "serial",
        "atom_name",
        "resname",
        "chain",
        "resseq",
        "icode",
        "coords",
        "atol",
    ),
    [
        # insertion code present -> lands in column 27
        (
            "ATOM   4353  N   SER A 807A      8.361 -19.707   7.711 -0.4000 1.5000",
            4353,
            "N",
            "SER",
            "A",
            807,
            "A",
            (8.361, -19.707, 7.711),
            0.0,
        ),
        # no insertion code -> column 27 blank, coordinate columns unshifted
        (
            "ATOM      1  N   PHE A 537      -8.834   3.827 -18.279 -0.7800 1.5000",
            1,
            "N",
            "PHE",
            "A",
            537,
            " ",
            (-8.834, 3.827, -18.279),
            0.0,
        ),
        # |x| >= 100 sheds a decimal to keep a leading space (documents current behaviour;
        # -100.168 -> -100.17, ~0.002 A error, see _fmt_float)
        (
            "ATOM      1  N   ALA A   1     -18.713-100.168  12.345 -0.4000 1.5000",
            1,
            "N",
            "ALA",
            "A",
            1,
            " ",
            (-18.713, -100.168, 12.345),
            0.01,
        ),
    ],
)
def test_format_places_fields_in_pdb_columns(
    line, serial, atom_name, resname, chain, resseq, icode, coords, atol
):
    """Every field must land in its fixed PDB column, stable across icode and width."""
    out = format_pqr_atom_line(parse_pqr_atom_line(line))
    assert out[:6].strip() == "ATOM"
    assert int(out[6:11]) == serial
    assert out[12:16].strip() == atom_name
    assert out[17:20].strip() == resname
    assert out[21] == chain
    assert int(out[22:26]) == resseq
    assert out[26] == icode
    x, y, z = coords
    assert np.isclose(float(out[30:38]), x, atol=atol)
    assert np.isclose(float(out[38:46]), y, atol=atol)
    assert np.isclose(float(out[46:54]), z, atol=atol)


@pytest.mark.parametrize(
    ("cleanup", "expect_water", "expect_ligand"),
    [
        ("nphs_lps", True, True),
        ("nphs_lps_waters", False, True),
        ("nphs_lps_waters_nonstdres", False, False),
    ],
)
def test_strip_hetatm_lines_follows_cleanup(cleanup, expect_water, expect_ligand):
    """Stripping works on raw PDB-format lines, not just PQR."""
    lines = [
        "ATOM   4353  N   SER A 807       8.361 -19.707   7.711  1.00 51.01           N",
        "HETATM 4400  O   HOH A 900       1.000   2.000   3.000  1.00 30.00           O",
        "HETATM 4401  C1  M3A A 901       4.000   5.000   6.000  1.00 25.00           C",
    ]
    kept, removed_waters, removed_hetatm = strip_hetatm_lines(lines, cleanup)
    text = "\n".join(kept)
    assert ("HOH" in text) is expect_water
    assert ("M3A" in text) is expect_ligand
    assert removed_waters == (0 if expect_water else 1)
    assert removed_hetatm == (0 if expect_ligand else 1)


@pytest.mark.parametrize(
    ("cleanup", "expected"),
    [
        ("nphs_lps", False),
        ("nphs_lps_waters", True),
        ("nphs_lps_waters_nonstdres", True),
    ],
)
def test_cleanup_drops_hetatm(tmp_path, defaults, pqr_file, cleanup, expected):
    """Only cleanup modes that remove HETATMs should trigger the PDB repair."""
    cfg = build(tmp_path, defaults, receptor_perception={"cleanup": cleanup})
    prep = MglTools(pqr_file, str(tmp_path / "rec"), cfg)
    assert prep._cleanup_drops_hetatm is expected


def test_full_pqr_fixes(tmp_path, defaults, squashed_pqr):
    """Test the full _fix_pqr_spacing method, which is called by run()."""
    cfg = build(tmp_path, defaults)
    prep = MglTools(squashed_pqr, str(tmp_path / "rec"), cfg)
    prep._fix_pqr_spacing(squashed_pqr, cleanup="nphs_lps")

    lines = [
        line
        for line in open(squashed_pqr).read().splitlines()
        if line.startswith("ATOM")
    ]
    first = lines[0]
    # fixed-width coordinate columns must now parse independently
    assert np.isclose(float(first[30:38]), -18.713)
    # atol: the wide-negative shim, see test_format_places_fields_in_pdb_columns
    assert np.isclose(float(first[38:46]), -100.168, atol=0.01)
    assert np.isclose(float(first[46:54]), 12.345)


@pytest.mark.parametrize("passes", [1, 2])
def test_backup_behavior(
    tmp_path, defaults, squashed_pqr, passes
):
    """Test that the original PQR is backed up before any in-place edits, even if run() is called multiple times."""
    cfg = build(tmp_path, defaults)
    prep = MglTools(squashed_pqr, str(tmp_path / "rec"), cfg)
    original = open(squashed_pqr).read()

    for _ in range(passes):
        prep._fix_pqr_spacing(squashed_pqr, cleanup="nphs_lps")

    backup = squashed_pqr.replace(".pqr", "_org.pqr")
    assert os.path.exists(backup)
    assert open(backup).read() == original


@pytest.mark.parametrize(
    ("cleanup", "expect_water", "expect_ligand"),
    [
        ("nphs_lps", True, True),
        ("nphs_lps_waters", False, True),
        ("nphs_lps_waters_nonstdres", False, False),
    ],
)
def test_hetatm_removal(
    tmp_path, defaults, squashed_pqr, cleanup, expect_water, expect_ligand
):
    """Test that HETATM removal should follow the -U cleanup mode."""
    cfg = build(tmp_path, defaults)
    prep = MglTools(squashed_pqr, str(tmp_path / "rec"), cfg)
    prep._fix_pqr_spacing(squashed_pqr, cleanup=cleanup)

    text = open(squashed_pqr).read()
    assert ("HOH" in text) is expect_water
    assert ("LIG" in text) is expect_ligand


@requires_mgltools
def test_run_generates_pdbqt(tmp_path, defaults, pqr_file):
    """End-to-end: prepare_receptor4.py should produce a parseable PDBQT."""
    import shutil

    local_pqr = tmp_path / "receptor.pqr"
    shutil.copy2(pqr_file, local_pqr)

    cfg = build(tmp_path, defaults)
    prep = MglTools(str(local_pqr), str(tmp_path / "rec"), cfg)
    out = prep.run()

    assert out == prep.pdbqt_file
    assert os.path.exists(out)
    atoms = [
        line
        for line in open(out).read().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]
    assert atoms, "PDBQT contains no atom records"
    assert atoms[0][77:79].strip()
