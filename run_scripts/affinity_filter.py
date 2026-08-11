#!/usr/bin/env python
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""
Filter Gypsum-DL output SDFs to keep only compounds with usable affinity
data AND a rotatable-bond count at or below the configured threshold.
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass

from loguru import logger

from lignova.hdf5 import ParquetParser

DEFAULT_CONVERTIBLE_TYPES = {
    "Kd",
    "Ki",
    "Kb",
    "Kieq",
    "Kic",
    "IC50",
    "EC50",
    "fIC50",
    "fEC50",
}

DEFAULT_MAX_ROTATABLE_BONDS = 10

PDB_LIGAND_ID_RE = re.compile(r"^([A-Za-z0-9]{3})_([A-Za-z0-9]{4})(?:_copy_\d+)?$")


@dataclass(frozen=True)
class FilterSets:
    """Identifier sets used to keep/drop compounds in an SDF.

    Attributes:
        usable:      Identifiers that passed affinity + RB filters.
        known:       Identifiers that appeared in the affinity parquet
                     for at least one of this PDB's genes, regardless of
                     whether they passed. Used so PDB-derived ligands
                     are filtered by affinity ONLY when data exists.
        rb_excluded: Identifiers with RotatableBonds > threshold and no
                     acceptable row anywhere.
    """

    usable: set[str]
    known: set[str]
    rb_excluded: set[str]


def pdb_ligand_id_candidates(ligand_id: str) -> set[str]:
    r"""
    Given a ligand ID from a gypsum SDF title, return the set of candidate the
    PDB-derived ligand codes that could match it.
    Args:
        ligand_id: The ligand ID from the SDF title.
    Returns:
        A set of candidate PDB-derived ligand codes.
    """
    m = PDB_LIGAND_ID_RE.match(ligand_id)
    if not m:
        return set()
    code, pdb = m.group(1), m.group(2)
    return {code, f"{code}_{pdb}"}


def _load_parquet_with_columns(path: str, required: list[str]) -> ParquetParser:
    """Open a parquet, validate required columns, return the read table.
    Args:
        path : Path to the parquet file.
        required : List of required column names.
    Returns:
        The ParquetParser object containing the read table.
    """
    if not os.path.exists(path):
        logger.error(f"The file {path} does not exist")
        raise FileNotFoundError(f"The file {path} does not exist")
    parser = ParquetParser(path)
    missing = [c for c in required if c not in parser.schema.names]
    if missing:
        logger.error(f"The file {path} is missing required columns: {missing}")
        raise ValueError(f"The file {path} is missing required columns: {missing}")
    return parser.read(columns=required)


def build_pdb_to_gene(file_path: str) -> dict[str, list[str]]:
    """Map each PDB ID to its list of associated gene IDs.
    Args:
        file_path : Path to the parquet file containing "Represenatives" and "members" columns.
    Returns:
        A dictionary mapping PDB IDs to lists of gene IDs.
    """
    table = _load_parquet_with_columns(file_path, ["Represenatives", "members"])
    pdb_to_genes: dict[str, list[str]] = defaultdict(list)
    for pdb, gene in zip(
        table.column("Represenatives").to_pylist(),
        table.column("members").to_pylist(),
    ):
        if gene not in pdb_to_genes[pdb]:
            pdb_to_genes[pdb].append(gene)
    return dict(pdb_to_genes)


def build_cid_filter_sets(
    affinity_path: str,
    gene_ids: set[str],
    convertible_types: set[str],
    max_rotatable_bonds: int | None = DEFAULT_MAX_ROTATABLE_BONDS,
) -> FilterSets:
    """Build identifier sets used to filter SDFs.
    Args:
        affinity_path: Path to the experimental affinity parquet.
        gene_ids: Set of gene IDs associated with the current PDB.
        convertible_types: Set of affinity types considered convertible.
        max_rotatable_bonds: Maximum number of rotatable bonds allowed.Defaults is DEFAULT_MAX_ROTATABLE_BONDS =10
    Returns:
        A FilterSets object containing the usable, known, and rb_excluded identifier sets.
    """
    required_cols = ["actives", "gene_id", "affinity", "type"]
    if max_rotatable_bonds is not None:
        required_cols.append("RotatableBonds")
    raw = _load_parquet_with_columns(affinity_path, required_cols)

    cids = raw.column("actives").to_pylist()
    genes = raw.column("gene_id").to_pylist()
    affs = raw.column("affinity").to_pylist()
    types = raw.column("type").to_pylist()
    rbs = (
        raw.column("RotatableBonds").to_pylist()
        if max_rotatable_bonds is not None
        else [None] * len(cids)
    )

    usable: set[str] = set()
    known: set[str] = set()
    rb_high: set[str] = set()
    rb_ok: set[str] = set()
    n_match_rows = 0

    for cid, gene, aff, tp, rb in zip(cids, genes, affs, types, rbs):
        if cid is None:
            continue
        cid_s = str(cid)
        if max_rotatable_bonds is not None and rb is not None:
            (rb_high if rb > max_rotatable_bonds else rb_ok).add(cid_s)

        if gene is None or str(gene) not in gene_ids:
            continue

        known.add(cid_s)

        if aff is None or aff <= 0 or tp not in convertible_types:
            continue
        if cid_s in rb_high and cid_s not in rb_ok:
            continue
        usable.add(cid_s)
        n_match_rows += 1

    rb_excluded = rb_high - rb_ok

    if max_rotatable_bonds is not None:
        logger.info(
            f"  {len(usable):,} usable out of {len(known):,} known CIDs from "
            f"{n_match_rows:,} affinity rows; "
            f"{len(rb_excluded):,} IDs flagged by RotatableBonds > "
            f"{max_rotatable_bonds} (applies to PubChem and PDB ligands) "
            f"(genes: {sorted(gene_ids)})"
        )
    else:
        logger.info(
            f"  {len(usable):,} usable out of {len(known):,} known CIDs from "
            f"{n_match_rows:,} affinity rows (genes: {sorted(gene_ids)})"
        )
    return FilterSets(usable=usable, known=known, rb_excluded=rb_excluded)


def filter_sdf_by_cid(
    sdf_path: str,
    sets: FilterSets,
    output_path: str,
    pdb_affinity_check: bool = True,
) -> dict:
    """Write an SDF containing only molecules that pass the configured filters.
    Args:
        sdf_path : Path to the input gypsum SDF file.
        sets : Identifier sets (usable, known, rb_excluded) built by build_cid_filter_sets.
        output_path : Path where the filtered SDF will be written.
        pdb_affinity_check: If True, PDB-derived ligands whose identifiers appear in sets.known are kept only when also in sets.usable; if False, PDB-derived ligands are filtered by RotatableBonds only. Default is True.
    Returns:
        A Counter dict with keys total, kept, kept_affinity, kept_pdb, drop_pdb_rb, drop_pdb_aff, dropped, wrote.
    """
    with open(sdf_path, "r", encoding="utf-8") as f:
        content = f.read()

    molecules = content.split("$$$$\n")
    kept_blocks: list[str] = []
    n_total = 0
    n_kept_affinity = 0
    n_kept_pdb = 0
    n_drop_pdb_rb = 0
    n_drop_pdb_aff = 0

    for mol in molecules:
        if not mol.strip():
            continue
        n_total += 1
        first_line = next((ln for ln in mol.splitlines() if ln.strip()), "")
        title = first_line.strip()
        candidates = pdb_ligand_id_candidates(title)

        if candidates:
            if candidates & sets.rb_excluded:
                n_drop_pdb_rb += 1
                continue
            if pdb_affinity_check and (candidates & sets.known):
                if not (candidates & sets.usable):
                    n_drop_pdb_aff += 1
                    continue
            kept_blocks.append(mol)
            n_kept_pdb += 1
        else:
            cid = title.split("_", 1)[0]
            if cid in sets.usable:
                kept_blocks.append(mol)
                n_kept_affinity += 1

    n_kept = n_kept_affinity + n_kept_pdb
    n_dropped = n_total - n_kept

    counters = {
        "total": n_total,
        "kept": n_kept,
        "kept_affinity": n_kept_affinity,
        "kept_pdb": n_kept_pdb,
        "drop_pdb_rb": n_drop_pdb_rb,
        "drop_pdb_aff": n_drop_pdb_aff,
        "dropped": n_dropped,
        "wrote": 0,
    }

    if n_kept == 0:
        logger.warning(
            f"  {os.path.basename(sdf_path)}: 0/{n_total} compounds passed. Skipping output"
        )
        return counters

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("$$$$\n".join(kept_blocks))
        f.write("$$$$\n")
    counters["wrote"] = 1
    logger.info(
        f"  {os.path.basename(sdf_path)}: kept {n_kept:,}/{n_total:,} "
        f"(affinity={n_kept_affinity:,}, pdb={n_kept_pdb:,}, "
        f"pdb dropped: RB={n_drop_pdb_rb:,}, aff={n_drop_pdb_aff:,}, "
        f"total dropped={n_dropped:,})  {os.path.basename(output_path)}"
    )
    return counters


def filter_protein_directory(
    pdb_dir: str,
    affinity_path: str,
    file_path: str,
    convertible_types: set[str],
    overwrite: bool = False,
    max_rotatable_bonds: int | None = DEFAULT_MAX_ROTATABLE_BONDS,
    pdb_affinity_check: bool = True,
) -> dict:
    """Build filter sets for one PDB and process all of its gypsum SDFs.
    Args:
        pdb_dir: The PDB protein directory containing gypsum_out_*/ subdirs.
        affinity_path : Path to the experimental affinity parquet.
        file_path : Path to protein_clustered_data.parquet for PDB to gene resolution.
        convertible_types : Set of affinity types considered convertible.
        overwrite : If True, re-filter even when a *_filtered.sdf already exists. Default is False.
        max_rotatable_bonds: Drop identifiers with RotatableBonds greater than this value. Set to None to disable. Default is DEFAULT_MAX_ROTATABLE_BONDS = 10.
        pdb_affinity_check: If True, PDB-derived ligands whose identifiers appear in the parquet are also filtered by affinity. Default is True.
    Returns:
        Aggregate counters for this PDB, including the protein_id field.
    """
    protein_id = os.path.basename(os.path.normpath(pdb_dir))
    logger.info(f"Filtering {protein_id}...")

    pdb_to_gene = build_pdb_to_gene(file_path)
    gene_ids = set(pdb_to_gene.get(protein_id, []))

    if not gene_ids:
        logger.warning(
            f"{protein_id}: no gene mapping. "
            f"Only PDB-derived compounds (passing RB filter) will be kept"
        )

    sets = build_cid_filter_sets(
        affinity_path,
        gene_ids,
        convertible_types,
        max_rotatable_bonds=max_rotatable_bonds,
    )

    if not sets.usable:
        logger.warning(
            f"{protein_id}: 0 CIDs with usable affinity. "
            f"Only PDB-derived compounds (passing RB + affinity-when-known) will be kept"
        )

    sdf_files = sorted(glob.glob(os.path.join(pdb_dir, "gypsum_out_*", "*.sdf")))
    sdf_files = [
        f
        for f in sdf_files
        if not f.endswith("_filtered.sdf") and not f.endswith("_docked.sdf")
    ]

    if not sdf_files:
        logger.warning(f"{protein_id}: no gypsum SDFs found")
        return {"protein_id": protein_id, "sdfs": 0}

    totals: dict[str, int] = defaultdict(int)
    for sdf in sdf_files:
        out_path = sdf.replace(".sdf", "_filtered.sdf")
        if not overwrite and os.path.isfile(out_path):
            logger.debug(f"  {os.path.basename(sdf)}: filtered exists, skipping")
            continue
        counters = filter_sdf_by_cid(
            sdf,
            sets,
            out_path,
            pdb_affinity_check=pdb_affinity_check,
        )
        for k, v in counters.items():
            totals[k] += v
        totals["sdfs"] += 1

    kept_pct = 100.0 * totals["kept"] / max(totals["total"], 1)
    dropped_pct = 100.0 * totals["dropped"] / max(totals["total"], 1)
    logger.info(f"{protein_id} summary: {totals['sdfs']} SDFs processed")
    logger.info(
        f"  Compounds: {totals['total']:,} total  "
        f"{totals['kept']:,} kept ({kept_pct:.1f}%) "
        f"[affinity={totals['kept_affinity']:,}, pdb={totals['kept_pdb']:,} "
        f"(pdb dropped: RB={totals['drop_pdb_rb']:,}, aff={totals['drop_pdb_aff']:,})], "
        f"{totals['dropped']:,} dropped ({dropped_pct:.1f}%)"
    )
    logger.info(f"  Output files written: {totals['wrote']}")
    return {"protein_id": protein_id, **dict(totals)}


def run_cli():
    """Command line interface for filtering Gypsum-DL SDFs by affinity and rotatable-bond count."""
    parser = argparse.ArgumentParser(
        description=(
            "Filter Gypsum-DL SDFs by experimental affinity AND rotatable-bond "
            "count. PDB-derived compounds pass unless they fail the rotatable-bond "
            "filter; when affinity data exists for them in the parquet, they are "
            "also filtered by affinity (unless --no-pdb-affinity-filter is set)."
        )
    )
    parser.add_argument(
        "-d",
        "--directory",
        type=str,
        required=True,
        help="PDB protein directory containing gypsum_out_*/ subdirs.",
    )
    parser.add_argument(
        "-a",
        "--affinity",
        type=str,
        required=True,
        help="Path to experimental affinity parquet "
        "(columns: actives, gene_id, affinity, type, RotatableBonds).",
    )
    parser.add_argument(
        "-c",
        "--clustered",
        type=str,
        required=True,
        help="Path to protein_clustered_data.parquet (PDBgene).",
    )
    parser.add_argument(
        "-t",
        "--convertible-types",
        type=str,
        required=False,
        default=None,
        help="Comma-separated override for convertible types. Default: "
        + ",".join(sorted(DEFAULT_CONVERTIBLE_TYPES)),
    )
    parser.add_argument(
        "-m",
        "--max-rotatable-bonds",
        type=int,
        required=False,
        default=DEFAULT_MAX_ROTATABLE_BONDS,
        help="Exclude identifiers with RotatableBonds > this value. Applies to "
        "PubChem CIDs and PDB-derived ligand codes. Pass a negative number "
        f"to disable. Default: {DEFAULT_MAX_ROTATABLE_BONDS}.",
    )
    parser.add_argument(
        "--no-pdb-affinity-filter",
        action="store_true",
        help="Force legacy behavior: PDB-derived ligands are never filtered by "
        "affinity, even when affinity data exists for them.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-filter even when a *_filtered.sdf already exists.",
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    if args.convertible_types:
        convertible = {
            t.strip() for t in args.convertible_types.split(",") if t.strip()
        }
        logger.info(f"Using user-provided convertible types: {sorted(convertible)}")
    else:
        convertible = DEFAULT_CONVERTIBLE_TYPES

    max_rb = args.max_rotatable_bonds if args.max_rotatable_bonds >= 0 else None
    if max_rb is None:
        logger.info("Rotatable-bond filter disabled (negative threshold)")
    else:
        logger.info(f"Rotatable-bond filter: RotatableBonds <= {max_rb}")

    pdb_affinity = not args.no_pdb_affinity_filter
    logger.info(
        "PDB-ligand affinity filter: "
        + ("ON when data is present" if pdb_affinity else "OFF (legacy)")
    )

    filter_protein_directory(
        pdb_dir=args.directory,
        affinity_path=args.affinity,
        file_path=args.clustered,
        convertible_types=convertible,
        overwrite=args.overwrite,
        max_rotatable_bonds=max_rb,
        pdb_affinity_check=pdb_affinity,
    )


if __name__ == "__main__":
    run_cli()
