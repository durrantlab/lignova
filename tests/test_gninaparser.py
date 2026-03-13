r"""Test the GNINA_Results SDF parser, DockedPose container, and DockingDataset."""

import glob
import gzip
import io
import os
import textwrap

import numpy as np
import pyarrow.compute as pc
import pytest

from lignova.analysis.gnina_parser import (
    _SCORE_DIRECTIONS,
    DOCKING_SCHEMA,
    DockedPose,
    DockingDataset,
    GNINA_Results,
    as_poses,
)

os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_filepaths = {
    "write_dir": "./tmp/gnina",
    "sdf_dir": "./tmp/gnina/sdf",
    "dataset_dir": "./tmp/gnina/dataset",
    "docking_tree": "./tmp/gnina/docking_tree",
}


def prep_dirs():
    r"""Prepare directories for writing files."""
    for d in context_filepaths.values():
        os.makedirs(d, exist_ok=True)


if not os.path.exists(context_filepaths["write_dir"]):
    prep_dirs()

_BLOCK_TEMPLATE = textwrap.dedent("""\
    {mol_name}


     36 37  0  0  0  0  0  0  0  0999 V2000
        2.7264    2.1529   27.6800 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.2995    1.5277   27.2217 C   0  0  0  0  0  0  0  0  0  0  0  0
       -0.1200    2.5566   27.8173 O   0  0  0  0  0  0  0  0  0  0  0  0
        1.7020    1.2219   27.1975 N   0  0  0  0  0  0  0  0  0  0  0  0
        2.0136    0.3579   26.6975 H   0  0  0  0  0  0  0  0  0  0  0  0
       -0.6813    0.6377   26.5124 C   0  0  0  0  0  0  0  0  0  0  0  0
       -1.3704    1.3593   25.3429 C   0  0  0  0  0  0  0  0  0  0  0  0
       -1.0097   -0.8235   24.2672 N   0  0  0  0  0  0  0  0  0  0  0  0
       -1.3356    0.5795   24.0109 C   0  0  0  0  0  0  0  0  0  0  0  0
       -1.1085   -1.3540   23.3718 H   0  0  0  0  0  0  0  0  0  0  0  0
       -1.7273   -1.2176   24.9182 H   0  0  0  0  0  0  0  0  0  0  0  0
       -0.3182    1.1577   23.0695 C   0  0  0  0  0  0  0  0  0  0  0  0
        0.6127    1.8875   23.5035 O   0  0  0  0  0  0  0  0  0  0  0  0
       -0.3933    0.8795   21.7112 O   0  0  0  0  0  0  0  0  0  0  0  0
        2.7544    2.1410   29.1823 C   0  0  0  0  0  0  0  0  0  0  0  0
        2.2208    3.0838   29.8300 O   0  0  0  0  0  0  0  0  0  0  0  0
        3.4290    1.0759   29.8677 N   0  0  0  0  0  0  0  0  0  0  0  0
        3.2341    0.7233   31.2762 C   0  0  0  0  0  0  0  0  0  0  0  0
        4.5528    0.3996   31.9758 C   0  0  0  0  0  0  0  0  0  0  0  0
        4.7319   -0.8369   32.6302 C   0  0  0  0  0  0  0  0  0  0  0  0
        5.9234   -1.1226   33.3030 C   0  0  0  0  0  0  0  0  0  0  0  0
        6.9414   -0.1718   33.3594 C   0  0  0  0  0  0  0  0  0  0  0  0
        6.7581    1.0801   32.7712 C   0  0  0  0  0  0  0  0  0  0  0  0
        5.5689    1.3715   32.0954 C   0  0  0  0  0  0  0  0  0  0  0  0
        2.2100   -0.3822   31.3736 C   0  0  0  0  0  0  0  0  0  0  0  0
        1.3261   -0.3394   32.2697 O   0  0  0  0  0  0  0  0  0  0  0  0
        2.2022   -1.4270   30.4545 O   0  0  0  0  0  0  0  0  0  0  0  0
        2.4331    3.5640   27.1620 C   0  0  0  0  0  0  0  0  0  0  0  0
        3.9514    4.5835   27.1884 S   0  0  0  0  0  0  0  0  0  0  0  0
        3.4459    5.7134   25.8466 C   0  0  0  0  0  0  0  0  0  0  0  0
        4.6013    5.9682   24.9157 C   0  0  0  0  0  0  0  0  0  0  0  0
        4.3993    6.6825   23.7221 C   0  0  0  0  0  0  0  0  0  0  0  0
        5.4692    6.9166   22.8494 C   0  0  0  0  0  0  0  0  0  0  0  0
        6.7457    6.4367   23.1589 C   0  0  0  0  0  0  0  0  0  0  0  0
        6.9537    5.7248   24.3443 C   0  0  0  0  0  0  0  0  0  0  0  0
        5.8859    5.4906   25.2210 C   0  0  0  0  0  0  0  0  0  0  0  0
      8  9  1  0
      9  7  1  0
      7  6  1  0
      6  2  1  0
      2  3  2  0
      2  4  1  0
      4  1  1  0
      1 28  1  0
     28 29  1  0
     29 30  1  0
     30 31  1  0
     31 32  2  0
     32 33  1  0
     33 34  2  0
     34 35  1  0
     35 36  2  0
      1 15  1  0
     15 16  2  0
     15 17  1  0
     17 18  1  0
     18 25  1  0
     25 26  2  0
     25 27  1  0
     18 19  1  0
     19 20  2  0
     20 21  1  0
     21 22  2  0
     22 23  1  0
     23 24  2  0
      9 12  1  0
     12 13  2  0
     12 14  1  0
     36 31  1  0
     24 19  1  0
      8 10  1  0
      8 11  1  0
      4  5  1  0
    M  END
    > <UniqueID>
    {uid}

    > <SMILES>
    {smiles}

    > <Energy>
    {energy}

    > <minimizedAffinity>
    {vina}

    > <CNNscore>
    {cnn_score}

    > <CNNaffinity>
    {cnn_aff}

    > <CNN_VS>
    {cnn_vs}

    > <CNNaffinity_variance>
    {cnn_var}
    """)

_N_ATOMS_TEMPLATE = 36


def _offsets(table, row):
    """Extract (block_start, block_end) from a table row."""
    return (
        table.column("block_start")[row].as_py(),
        table.column("block_end")[row].as_py(),
    )


def _energy(i):
    return -100.0 + i * 1.0


def _vina(i):
    return -8.0 + i * 0.1


def _cnn(i):
    return 0.50 + i * 0.01


def _cnn_aff(i):
    return 5.0 + i * 0.1


def _cnn_vs(i):
    return 0.30 + i * 0.02


def _cnn_var(i):
    return 0.010 + i * 0.001


def _make_block(mol_name, uid, smiles, idx):
    return _BLOCK_TEMPLATE.format(
        mol_name=mol_name,
        uid=uid,
        smiles=smiles,
        energy=f"{_energy(idx):.4f}",
        vina=f"{_vina(idx):.4f}",
        cnn_score=f"{_cnn(idx):.4f}",
        cnn_aff=f"{_cnn_aff(idx):.4f}",
        cnn_vs=f"{_cnn_vs(idx):.4f}",
        cnn_var=f"{_cnn_var(idx):.6f}",
    )


def _build_protein_sdf(ligands, uids_per_lig, n_conformers, num_modes):
    """Each (lig, uid) group gets n_conformers * num_modes blocks."""
    blocks, idx = [], 0
    for li, lig in enumerate(ligands):
        for uid in range(1, uids_per_lig + 1):
            smiles = f"C{'C' * (li + uid)}O"
            for _ in range(n_conformers * num_modes):
                blocks.append(_make_block(lig, uid, smiles, idx))
                idx += 1
    return "$$$$\n".join(blocks) + "\n$$$$\n", idx


def _write_sdf(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def _write_sdf_gz(path, text):
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(text)


# Protein configs
#   P1: 3 lig × 3 uid × (3 conf × 3 modes) = 81 blocks
#   P2: 2 lig × 2 uid × (2 conf × 2 modes) = 16 blocks
#   P3: 1 lig × 1 uid × (1 conf × 1 mode)  =  1 block
_P1_LIGS, _P1_UIDS, _P1_CONF, _P1_MODES = ["L1A", "L1B", "L1C"], 3, 3, 3
_P1_N = len(_P1_LIGS) * _P1_UIDS * _P1_CONF * _P1_MODES
_P1_PROTOMERS = len(_P1_LIGS) * _P1_UIDS
_P1_CONFORMERS = _P1_PROTOMERS * _P1_CONF
_P1_BLOCKS_PER_GROUP = _P1_CONF * _P1_MODES

_P2_LIGS, _P2_UIDS, _P2_CONF, _P2_MODES = ["L2A", "L2B"], 2, 2, 2
_P2_N = len(_P2_LIGS) * _P2_UIDS * _P2_CONF * _P2_MODES
_P2_BLOCKS_PER_GROUP = _P2_CONF * _P2_MODES

_P3_LIGS, _P3_UIDS, _P3_CONF, _P3_MODES = ["L3A"], 1, 1, 1
_P3_N = 1

_sdf = context_filepaths["sdf_dir"]

for _n, _l, _u, _c, _m in [
    ("P1", _P1_LIGS, _P1_UIDS, _P1_CONF, _P1_MODES),
    ("P2", _P2_LIGS, _P2_UIDS, _P2_CONF, _P2_MODES),
    ("P3", _P3_LIGS, _P3_UIDS, _P3_CONF, _P3_MODES),
]:
    _txt, _ = _build_protein_sdf(_l, _u, _c, _m)
    _write_sdf(os.path.join(_sdf, f"{_n}_docked.sdf"), _txt)
    _write_sdf_gz(os.path.join(_sdf, f"{_n}_docked.sdf.gz"), _txt)

# Edge case: (E1=5 blocks, E2=3 blocks)
_edge_blocks = []
for i in range(5):
    _edge_blocks.append(_make_block("E1", 1, "CCO", i))
for i in range(5, 8):
    _edge_blocks.append(_make_block("E2", 1, "CCCO", i))
_write_sdf(
    os.path.join(_sdf, "EDGE_docked.sdf"),
    "$$$$\n".join(_edge_blocks) + "\n$$$$\n",
)

_p1 = GNINA_Results(
    os.path.join(_sdf, "P1_docked.sdf"), num_modes=_P1_MODES, protein_id="P1"
)
_p2 = GNINA_Results(
    os.path.join(_sdf, "P2_docked.sdf"), num_modes=_P2_MODES, protein_id="P2"
)
_p3 = GNINA_Results(
    os.path.join(_sdf, "P3_docked.sdf"), num_modes=_P3_MODES, protein_id="P3"
)

# Docking tree + dataset at module level
_tree = context_filepaths["docking_tree"]
_ds = context_filepaths["dataset_dir"]
_DDS_P2_N = 9

for _name, _ligs, _uids in [("P1", _P1_LIGS, _P1_UIDS), ("P2", ["L2A"], 1)]:
    _gdir = os.path.join(_tree, _name, "gypsum_out_0")
    os.makedirs(_gdir, exist_ok=True)
    _txt, _ = _build_protein_sdf(_ligs, _uids, 3, 3)
    _write_sdf(os.path.join(_gdir, "lig_docked.sdf"), _txt)
    with open(os.path.join(_tree, _name, f"{_name}_cleaned.pdb"), "w") as _f:
        _f.write("ATOM mock\nEND\n")

_dds = DockingDataset(_ds)
_dds.build_from_docking_tree(_tree, num_modes=3, copy_blocks=True, copy_proteins=True)
_DDS_TOTAL = _P1_N + _DDS_P2_N


def test_init():
    """Defaults, gz variant, and protein_id inference."""
    dr = GNINA_Results(
        os.path.join(_sdf, "P1_docked.sdf"), num_modes=_P1_MODES, protein_id="P1"
    )
    assert dr.protein_id == "P1"
    assert dr.n_blocks == _P1_N
    assert dr.schema == DOCKING_SCHEMA

    dr_gz = GNINA_Results(
        os.path.join(_sdf, "P1_docked.sdf.gz"), num_modes=_P1_MODES, protein_id="P1"
    )
    assert dr_gz.n_blocks == _P1_N

    dr_infer = GNINA_Results(os.path.join(_sdf, "P1_docked.sdf"), num_modes=_P1_MODES)
    assert dr_infer.protein_id == "P1"


def test_init_invalid():
    with pytest.raises(ValueError, match="does not exist"):
        GNINA_Results("/tmp/_no_such_.sdf", protein_id="X")

    bad = os.path.join(context_filepaths["write_dir"], "bad.txt")
    with open(bad, "w") as f:
        f.write("junk")
    with pytest.raises(ValueError, match="not a .sdf"):
        GNINA_Results(bad, protein_id="X")


def test_counts_across_proteins():
    for dr, ligs, uids, n_conf, n_blk in [
        (_p1, _P1_LIGS, _P1_UIDS, _P1_CONF, _P1_N),
        (_p2, _P2_LIGS, _P2_UIDS, _P2_CONF, _P2_N),
        (_p3, _P3_LIGS, _P3_UIDS, _P3_CONF, _P3_N),
    ]:
        assert dr.n_blocks == n_blk
        assert len(dr) == n_blk
        assert dr.n_ligands == len(ligs)
        assert dr.n_protomers == len(ligs) * uids
        assert dr.n_conformers == len(ligs) * uids * n_conf

    assert _p1.n_protomers == 9
    assert _p1.n_conformers == 27


def test_table_content():
    for field in DOCKING_SCHEMA:
        assert field.name in _p1.table.column_names
    assert all(v == "P1" for v in _p1.table.column("protein_id").to_pylist())
    src = set(_p1.table.column("source_file").to_pylist())
    assert len(src) == 1 and "P1_docked.sdf" in src.pop()


def test_score_values():
    t = _p1.table
    for i in range(t.num_rows):
        assert t.column("CNNscore")[i].as_py() == pytest.approx(_cnn(i), abs=1e-6)
        assert t.column("Vina_affinity")[i].as_py() == pytest.approx(_vina(i), abs=1e-6)
        assert t.column("Energy")[i].as_py() == pytest.approx(_energy(i), abs=1e-6)
    for col in ["CNNscore", "Vina_affinity", "Energy"]:
        assert all(v is not None and not np.isnan(v) for v in t.column(col).to_pylist())


def test_n_atoms_and_conformer_assignment():
    t = _p1.table
    assert all(v == _N_ATOMS_TEMPLATE for v in t.column("n_atoms").to_pylist())

    for i in range(9):
        assert t.column("conformer_idx")[i].as_py() == i // _P1_MODES
        assert t.column("pose_rank")[i].as_py() == i % _P1_MODES
        assert t.column("ligand_id")[i].as_py() == "L1A"
        assert t.column("UniqueID")[i].as_py() == 1

    assert t.column("UniqueID")[9].as_py() == 2
    assert t.column("conformer_idx")[9].as_py() == 0


def test_num_modes_variations():
    assert _p1.num_modes == [_P1_MODES]

    nm_map = _p1.num_modes_per_group
    assert all(v == _P1_MODES for v in nm_map.values())
    assert len(nm_map) == _P1_PROTOMERS

    dr_list = GNINA_Results(
        os.path.join(_sdf, "P1_docked.sdf"), num_modes=[3, 9], protein_id="P1"
    )
    nm = dr_list.num_modes
    assert isinstance(nm, list)
    assert all(v in (3, 9) for v in nm)

    dr_auto = GNINA_Results(
        os.path.join(_sdf, "P1_docked.sdf"), num_modes=None, protein_id="P1"
    )
    nm2 = dr_auto.num_modes
    assert isinstance(nm2, list)
    assert _P1_BLOCKS_PER_GROUP in nm2
    dr = GNINA_Results(
        os.path.join(_sdf, "EDGE_docked.sdf"), num_modes=[3], protein_id="EDGE"
    )
    assert dr.n_blocks == 8
    assert 3 in dr.num_modes

    nm_map = dr.num_modes_per_group
    assert nm_map[("E2", 1)] == 3
    assert _p3.num_modes == [_P3_MODES]
    assert _p3.table.column("conformer_idx")[0].as_py() == 0
    assert _p3.table.column("pose_rank")[0].as_py() == 0


# GNINA_Results — filtering
def test_get_ligand():
    assert _p1.get_ligand("L1A").num_rows == _P1_UIDS * _P1_BLOCKS_PER_GROUP
    assert _p2.get_ligand("L2A").num_rows == 2 * _P2_BLOCKS_PER_GROUP
    assert _p1.get_ligand("NOPE").num_rows == 0


def test_get_protomer():
    sub = _p1.get_protomer("L1A", 1)
    assert sub.num_rows == _P1_BLOCKS_PER_GROUP
    assert set(sub.column("conformer_idx").to_pylist()) == {0, 1, 2}
    assert _p1.get_protomer("L1A", 999).num_rows == 0


def test_get_conformer():
    sub = _p1.get_conformer("L1A", 1, conformer_idx=2)
    assert sub.num_rows == _P1_MODES
    assert set(sub.column("conformer_idx").to_pylist()) == {2}
    assert _p1.get_conformer("L1A", 1, conformer_idx=99).num_rows == 0


def test_filter():
    assert _p1.filter(ligand_id="L1C").num_rows == _P1_UIDS * _P1_BLOCKS_PER_GROUP
    assert _p1.filter(ligand_id="L1C", UniqueID=2).num_rows == _P1_BLOCKS_PER_GROUP

    sub_min = _p1.filter(CNNscore_min=0.6)
    assert sub_min.num_rows == 71

    assert _p1.filter(CNNscore_min=0.55, CNNscore_max=0.65).num_rows == 11

    assert _p1.filter(ligand_id="L1A", CNNscore_min=0.6).num_rows == 17

    assert _p1.filter().num_rows == _P1_N


# GNINA_Results — ranking
def test_get_top_poses():
    top = _p1.get_top_poses()
    assert top.num_rows == _p1.n_conformers

    for i in range(top.num_rows):
        lig = top.column("ligand_id")[i].as_py()
        uid = top.column("UniqueID")[i].as_py()
        cidx = top.column("conformer_idx")[i].as_py()
        group = _p1.get_conformer(lig, uid, cidx)
        assert top.column("CNNscore")[i].as_py() == pytest.approx(
            pc.max(group.column("CNNscore")).as_py()
        )
    assert _p1.get_top_poses(n=2).num_rows == _p1.n_conformers * 2
    assert _p1.get_top_poses(per="protomer").num_rows == _p1.n_protomers

    top_lig = _p1.get_top_poses(per="ligand")
    assert top_lig.num_rows == _p1.n_ligands
    expected = {"L1A": _cnn(26), "L1B": _cnn(53), "L1C": _cnn(80)}
    for i in range(top_lig.num_rows):
        lig = top_lig.column("ligand_id")[i].as_py()
        assert top_lig.column("CNNscore")[i].as_py() == pytest.approx(
            expected[lig], abs=1e-6
        )

    top_w = _p1.get_top_poses(
        by=["CNNscore", "Vina_affinity"], weights=[0.7, 0.3], per="ligand"
    )
    assert top_w.num_rows == _p1.n_ligands


def test_get_top_poses_edge():
    for per in ("conformer", "protomer", "ligand"):
        assert _p3.get_top_poses(per=per).num_rows == 1
    assert _p1.get_top_poses(n=100, per="conformer").num_rows == _P1_N
    with pytest.raises(ValueError, match="Length mismatch"):
        _p1.get_top_poses(by=["CNNscore", "Vina_affinity"], ascending=[True])
    with pytest.raises(ValueError, match="Length mismatch"):
        _p1.get_top_poses(by=["CNNscore", "Vina_affinity"], weights=[1.0])


def test_get_best_per_ligand():
    best = _p1.get_best_per_ligand()
    assert best.num_rows == _p1.n_ligands

    best_v = _p1.get_best_per_ligand(by="Vina_affinity")
    expected = {"L1A": _vina(0), "L1B": _vina(27), "L1C": _vina(54)}
    for i in range(best_v.num_rows):
        lig = best_v.column("ligand_id")[i].as_py()
        assert best_v.column("Vina_affinity")[i].as_py() == pytest.approx(
            expected[lig], abs=1e-6
        )


def test_best_direction():
    for score, direction in _SCORE_DIRECTIONS.items():
        assert GNINA_Results.best_direction(score) == direction
    with pytest.raises(ValueError, match="Unknown score"):
        GNINA_Results.best_direction("not_a_score")


def test_coords_and_blocks():
    t = _p1.table
    bs0, be0 = _offsets(t, 0)
    bs1, be1 = _offsets(t, 1)

    coords, elems = _p1.get_coords(bs0, be0)
    assert coords.shape == (_N_ATOMS_TEMPLATE, 3)
    assert len(elems) == _N_ATOMS_TEMPLATE
    assert elems[0] == "C"
    assert coords[0, 0] == 2.7264

    coords2, _ = _p1.get_coords(bs1, be1)
    np.testing.assert_array_almost_equal(coords, coords2)

    offsets_batch = [_offsets(t, i) for i in [0, 5, 80]]
    batch = _p1.get_coords_batch(offsets_batch)
    assert len(batch) == 3

    text = _p1.get_block_by_offsets(bs0, be0)
    assert "V2000" in text and "L1A" in text

    assert GNINA_Results.read_block_from_file(_p1.filepath, bs0, be0) == text


def test_read_block_from_gz():
    gz = os.path.join(_sdf, "P1_docked.sdf.gz")
    dr = GNINA_Results(gz, num_modes=_P1_MODES, protein_id="P1")
    row = dr.table.slice(0, 1)
    bs, be = row.column("block_start")[0].as_py(), row.column("block_end")[0].as_py()
    text = GNINA_Results.read_block_from_file(gz, bs, be)
    assert "L1A" in text and "V2000" in text


def test_get_mol():
    t = _p1.table
    bs, be = _offsets(t, 0)
    mol = _p1.get_mol(bs, be)
    assert mol is not None and mol.GetNumAtoms() == _N_ATOMS_TEMPLATE


def test_summary():
    buf = io.StringIO()
    _p1.summary(per="global", output=buf)
    text = buf.getvalue()
    assert f"Total blocks:   {_P1_N}" in text
    # Verify offset range format
    assert "offset=" in text
    for line in text.split("\n"):
        if "offset=" in line:
            offset_part = line.split("offset=")[1].split(",")[0]
            assert ":" in offset_part

    out_lig = os.path.join(context_filepaths["write_dir"], "summary_lig.txt")
    _p1.summary(per="ligand", output=out_lig)
    with open(out_lig) as f:
        lig_text = f.read()
    assert all(lig in lig_text for lig in _P1_LIGS)

    out_prot = os.path.join(context_filepaths["write_dir"], "summary_prot.txt")
    _p1.summary(per="protomer", output=out_prot)
    with open(out_prot) as f:
        assert "UID=" in f.read()


def test_summary_invalid():
    with pytest.raises(ValueError, match="Unknown per"):
        _p1.summary(per="galaxy")


def test_to_csv():
    out = os.path.join(context_filepaths["write_dir"], "out.csv")
    _p1.to_csv(out)
    import pyarrow.csv as pcsv

    assert pcsv.read_csv(out).num_rows == _P1_N


def test_to_sdf():
    wdir = context_filepaths["write_dir"]

    out = os.path.join(wdir, "all.sdf")
    _p1.to_sdf(out, overwrite=True)
    with open(out) as f:
        assert f.read().count("$$$$") == _P1_N

    out2 = os.path.join(wdir, "subset.sdf")
    subset = _p1.table.take([0, 1])
    _p1.to_sdf(out2, table=subset, overwrite=True)
    with open(out2) as f:
        assert f.read().count("$$$$") == 2

    out3 = os.path.join(wdir, "all.sdf.gz")
    _p1.to_sdf(out3, overwrite=True)
    with gzip.open(out3, "rt") as f:
        assert f.read().count("$$$$") == _P1_N


# GNINA_Results — aggregation, repr, len
def test_poses_per_conformer():
    vc = _p1.poses_per_conformer()
    counts = [e["counts"] for e in vc.to_pylist()]
    assert all(c == _P1_MODES for c in counts)
    assert len(counts) == _p1.n_conformers
    t = _p1.protomer_counts()
    assert t.num_rows == _p1.n_protomers
    for i in range(t.num_rows):
        assert t.column("n_blocks")[i].as_py() == _P1_BLOCKS_PER_GROUP
        assert t.column("n_conformers")[i].as_py() == _P1_CONF


def test_repr_len():
    r = repr(_p1)
    assert "GNINA_Results" in r and f"blocks={_P1_N}" in r
    assert len(_p1) == _P1_N
    assert len(_p3) == 1


# DockedPose
def test_docked_pose():
    pose = DockedPose(0, _p1.table, _p1)
    assert pose.ligand_id == "L1A"
    assert pose.UniqueID == 1
    assert pose.protein_id == "P1"
    assert pose.conformer_idx == 0
    assert pose.pose_rank == 0
    assert pose.CNNscore == pytest.approx(_cnn(0), abs=1e-6)
    assert pose.Vina_affinity == pytest.approx(_vina(0), abs=1e-6)
    assert "V2000" in pose.block_text

    mol = pose.mol
    assert mol is not None and mol.GetNumAtoms() == _N_ATOMS_TEMPLATE

    last = DockedPose(_P1_N - 1, _p1.table, _p1)
    assert last.ligand_id == "L1C"
    assert last.CNNscore == pytest.approx(_cnn(80), abs=1e-6)
    assert "DockedPose" in repr(last)


# as_poses
def test_as_poses():
    poses = as_poses(_p1.table, _p1)
    assert len(poses) == _P1_N
    assert all(isinstance(p, DockedPose) for p in poses)

    sub = _p1.get_ligand("L1A")
    assert len(as_poses(sub, _p1)) == _P1_UIDS * _P1_BLOCKS_PER_GROUP


# DockingDataset
def test_dataset_reads():
    assert "P1" in _dds.protein_ids() and "P2" in _dds.protein_ids()

    s = _dds.stats()
    assert s["n_proteins"] == 2 and s["total_poses"] == _DDS_TOTAL

    t = _dds.read_protein("P1")
    assert t.num_rows == _P1_N
    assert _dds.read_protein("P1", columns=["ligand_id"]).column_names == ["ligand_id"]

    assert _dds.read_proteins(["P1", "P2"]).num_rows == _DDS_TOTAL
    assert _dds.read_proteins(["P1", "MISSING"]).num_rows == _P1_N
    empty = _dds.read_proteins(["X", "Y"])
    assert empty.num_rows == 0 and empty.schema == DOCKING_SCHEMA

    t_all = _dds.read_all()
    assert t_all.num_rows == _DDS_TOTAL
    assert set(_dds.read_all(ligand_id="L1A").column("ligand_id").to_pylist()) == {
        "L1A"
    }
    with pytest.raises(FileNotFoundError):
        _dds.read_protein("NOPE")
    assert "proteins=2" in repr(_dds)


def test_dataset_iter_batches():
    assert (
        sum(b.num_rows for b in _dds.iter_batches(protein_id="P1", batch_size=10))
        == _P1_N
    )
    assert sum(b.num_rows for b in _dds.iter_batches(batch_size=10)) == _DDS_TOTAL


# DockingDataset — write tests
def test_dataset_build_from_tree(tmp_path):
    tree = str(tmp_path / "tree")
    for name, ligs, uids in [("P1", _P1_LIGS, _P1_UIDS), ("P2", ["L2A"], 1)]:
        gdir = os.path.join(tree, name, "gypsum_out_0")
        os.makedirs(gdir)
        txt, _ = _build_protein_sdf(ligs, uids, 3, 3)
        _write_sdf(os.path.join(gdir, "lig_docked.sdf"), txt)
        with open(os.path.join(tree, name, f"{name}_cleaned.pdb"), "w") as f:
            f.write("ATOM\nEND\n")

    ds1 = str(tmp_path / "ds1")
    dds1 = DockingDataset(ds1)
    results = dds1.build_from_docking_tree(
        tree, num_modes=3, copy_blocks=True, copy_proteins=True
    )
    assert results["P1"] == _P1_N and results["P2"] == _DDS_P2_N
    assert os.path.isdir(os.path.join(ds1, "blocks", "P1"))
    assert os.path.isfile(os.path.join(ds1, "proteins", "P1.pdb"))

    ds2 = str(tmp_path / "ds2")
    dds2 = DockingDataset(ds2)
    dds2.build_from_docking_tree(
        tree, num_modes=3, copy_blocks=False, copy_proteins=False
    )
    assert not os.path.isdir(os.path.join(ds2, "blocks"))
    ds = str(tmp_path / "ds_add")
    dds = DockingDataset(ds)
    assert (
        dds.add_protein(
            "P3", [os.path.join(_sdf, "P3_docked.sdf")], num_modes=_P3_MODES
        )
        == _P3_N
    )
    assert (
        dds.add_protein(
            "P2", [os.path.join(_sdf, "P2_docked.sdf.gz")], num_modes=_P2_MODES
        )
        == _P2_N
    )


def test_dataset_batched_parquets(tmp_path):
    out1 = str(tmp_path / "b1")
    assert _dds.to_batched_parquets(out1, proteins_per_batch=1) == 2
    assert len(glob.glob(os.path.join(out1, "batch_*.parquet"))) == 2

    out2 = str(tmp_path / "b2")
    assert _dds.to_batched_parquets(out2, proteins_per_batch=100) == 1
    assert DockingDataset.read_batched_parquets(out2).num_rows == _DDS_TOTAL
    assert set(
        DockingDataset.read_batched_parquets(out2, ligand_id="L1A")
        .column("ligand_id")
        .to_pylist()
    ) == {"L1A"}
    assert (
        sum(b.num_rows for b in DockingDataset.iter_batches_from(out2, batch_size=10))
        == _DDS_TOTAL
    )
    ds = str(tmp_path / "empty")
    os.makedirs(os.path.join(ds, "parquet"))
    with pytest.raises(FileNotFoundError):
        DockingDataset(ds).to_batched_parquets(str(tmp_path / "out"))


# Test against real GNINA output sample
_SAMPLE_SDF = "./files/sample_docked.sdf.gz"


def test_real_sdf_sample():
    """Parse a real GNINA output excerpt and verify exact values."""
    dr = GNINA_Results(_SAMPLE_SDF, num_modes=None, protein_id="10GS")
    assert dr.n_blocks == 52
    assert dr.n_ligands == 3
    assert dr.n_protomers == 6
    assert dr.n_conformers == 6
    assert len(dr) == 52
    assert dr.protein_id == "10GS"
    assert dr.num_modes == [7, 9]
    nm_map = dr.num_modes_per_group
    assert nm_map[("VWW_10GS", 377)] == 9
    assert nm_map[("SAS_13GS", 382)] == 7
    t = dr.table
    assert t.column("ligand_id")[0].as_py() == "VWW_10GS"
    assert t.column("UniqueID")[0].as_py() == 377
    assert t.column("conformer_idx")[0].as_py() == 0
    assert t.column("pose_rank")[0].as_py() == 0
    assert t.column("n_atoms")[0].as_py() == 36
    assert t.column("SMILES")[0].as_py() == (
        "N[C@@H](CCC(=O)N[C@@H](CSCc1ccccc1)"
        "C(=O)[N-][C@@H](C(=O)[O-])c1ccccc1)C(=O)[O-]"
    )
    assert t.column("Energy")[0].as_py() == 42.79254465218859
    assert t.column("Vina_affinity")[0].as_py() == -11.35966
    assert t.column("CNNscore")[0].as_py() == 0.5107182264
    assert t.column("CNNaffinity")[0].as_py() == 4.9588913918
    assert t.column("CNN_VS")[0].as_py() == 2.5325961113
    assert t.column("CNNaffinity_variance")[0].as_py() == 0.1850149035

    assert t.column("ligand_id")[8].as_py() == "VWW_10GS"
    assert t.column("conformer_idx")[8].as_py() == 0
    assert t.column("pose_rank")[8].as_py() == 8
    assert t.column("CNNscore")[8].as_py() == 0.3145772517
    assert t.column("Vina_affinity")[8].as_py() == -8.63333

    assert t.column("ligand_id")[9].as_py() == "0HH_12GS"
    assert t.column("UniqueID")[9].as_py() == 378
    assert t.column("n_atoms")[9].as_py() == 32
    assert t.column("CNNscore")[9].as_py() == 0.3482419252
    assert t.column("Vina_affinity")[9].as_py() == -7.26635

    assert t.column("ligand_id")[18].as_py() == "0HH_12GS"
    assert t.column("UniqueID")[18].as_py() == 379
    assert t.column("conformer_idx")[18].as_py() == 0
    assert t.column("pose_rank")[18].as_py() == 0
    assert t.column("CNNscore")[18].as_py() == 0.3611989319

    assert t.column("ligand_id")[27].as_py() == "SAS_13GS"
    assert t.column("UniqueID")[27].as_py() == 380
    assert t.column("n_atoms")[27].as_py() == 29
    assert t.column("CNNscore")[27].as_py() == 0.5947299004
    assert t.column("Vina_affinity")[27].as_py() == -10.93801

    assert t.column("UniqueID")[36].as_py() == 381
    assert t.column("n_atoms")[36].as_py() == 30
    assert t.column("CNNscore")[36].as_py() == 0.5309475660

    assert t.column("UniqueID")[45].as_py() == 382
    assert t.column("n_atoms")[45].as_py() == 29
    assert t.column("CNNscore")[45].as_py() == 0.5417990685

    assert t.column("ligand_id")[51].as_py() == "SAS_13GS"
    assert t.column("pose_rank")[51].as_py() == 6
    assert t.column("CNNscore")[51].as_py() == 0.3796257079

    assert dr.get_ligand("VWW_10GS").num_rows == 9  # 1 UID × 9 modes
    assert dr.get_ligand("0HH_12GS").num_rows == 18  # 2 UIDs × 9 modes
    assert dr.get_ligand("SAS_13GS").num_rows == 25  # 9 + 9 + 7 (UID 382 truncated)

    p = dr.get_protomer("0HH_12GS", 379)
    assert p.num_rows == 9
    assert set(p.column("UniqueID").to_pylist()) == {379}

    c = dr.get_conformer("SAS_13GS", 380, conformer_idx=0)
    assert c.num_rows == 9

    top = dr.get_top_poses(per="ligand")
    assert top.num_rows == 3
    top_map = {
        top.column("ligand_id")[i].as_py(): top.column("CNNscore")[i].as_py()
        for i in range(3)
    }
    assert top_map["VWW_10GS"] == 0.5107182264
    assert top_map["SAS_13GS"] == 0.5947299004
    assert top_map["0HH_12GS"] == 0.3611989319

    # get_best_per_ligand by Vina (lowest = best)
    best_v = dr.get_best_per_ligand(by="Vina_affinity")
    vina_map = {
        best_v.column("ligand_id")[i].as_py(): best_v.column("Vina_affinity")[i].as_py()
        for i in range(3)
    }
    assert vina_map["VWW_10GS"] == -11.35966
    assert vina_map["SAS_13GS"] == -12.53481
    assert vina_map["0HH_12GS"] == -9.84825

    bs0, be0 = _offsets(t, 0)
    coords, elems = dr.get_coords(bs0, be0)
    assert coords.shape == (36, 3)
    assert len(elems) == 36
    assert elems[0] == "C"
    assert coords[0, 0] == 2.7264
    assert coords[0, 1] == 2.1529
    assert coords[0, 2] == 27.6800

    mol = dr.get_mol(bs0, be0)
    assert mol is not None
    assert mol.GetNumAtoms() == 36

    pose = DockedPose(0, t, dr)
    assert pose.ligand_id == "VWW_10GS"
    assert pose.UniqueID == 377
    assert pose.CNNscore == 0.5107182264
    assert "V2000" in pose.block_text
    assert pose.mol is not None

    # CNNscore >= 0.50
    hi = dr.filter(CNNscore_min=0.50)
    assert hi.num_rows == 5
    scores = hi.column("CNNscore").to_pylist()
    assert all(s >= 0.50 for s in scores)
    assert 0.5107182264 in scores
    assert 0.5947299004 in scores
    assert 0.5309475660 in scores

    vc = dr.poses_per_conformer()
    counts = sorted([e["counts"] for e in vc.to_pylist()])
    assert len(counts) == 6
    assert counts == [7, 9, 9, 9, 9, 9]

    pc_t = dr.protomer_counts()
    assert pc_t.num_rows == 6
    # Build a lookup by UniqueID for flexible checking
    pc_map = {
        pc_t.column("UniqueID")[i].as_py(): (
            pc_t.column("n_blocks")[i].as_py(),
            pc_t.column("n_conformers")[i].as_py(),
        )
        for i in range(pc_t.num_rows)
    }
    # All protomers have 1 conformer; UID 382 has 7 blocks, others have 9
    for uid in (377, 378, 379, 380, 381):
        assert pc_map[uid] == (9, 1)
    assert pc_map[382] == (7, 1)

    buf = io.StringIO()
    dr.summary(per="global", output=buf)
    txt = buf.getvalue()
    assert "Total blocks:   52" in txt
    assert "10GS" in txt
    assert "offset=" in txt
    for line in txt.split("\n"):
        if "offset=" in line:
            offset_part = line.split("offset=")[1].split(",")[0]
            assert ":" in offset_part

    out_csv = os.path.join(context_filepaths["write_dir"], "sample.csv")
    dr.to_csv(out_csv)
    import pyarrow.csv as pcsv

    assert pcsv.read_csv(out_csv).num_rows == 52

    out_sdf = os.path.join(context_filepaths["write_dir"], "sample_out.sdf")
    subset = t.take([0, 27])
    dr.to_sdf(out_sdf, table=subset, overwrite=True)
    with open(out_sdf) as f:
        assert f.read().count("$$$$") == 2

    poses = as_poses(t, dr)
    assert len(poses) == 52
    assert poses[0].ligand_id == "VWW_10GS"
    assert poses[51].ligand_id == "SAS_13GS"
