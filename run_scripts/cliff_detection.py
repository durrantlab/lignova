# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Butina Clustering and activity-cliff detection over a chosen set of targets from the `data_extraction.py` parquets."""

import argparse
import math
import os
import random
import time
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator

import pyarrow as pa
from loguru import logger

from lignova.activity import (
    CONVERTIBLE_TYPES,
    DEFAULT_SPREAD_GATE,
    aggregate_activity,
    cid_crosswalk,
    group_measurements,
)
from lignova.clustering import (
    ButinaClustering,
    ButinaParams,
    CliffParams,
    MorganFeaturizer,
    SaliUndefined,
    SeverityMetric,
    compute_pairwise,
    find_activity_cliffs,
    high_confidence_cliffs,
    label_cliff_severity,
    resolve_smiles,
    same_assay_cliffs,
)
from lignova.hdf5 import ParquetParser

# Column names in the enriched parquet.
GENE_COL = "target_geneid"
"""The column name for the target gene ID."""

KEY_COL = "InChIKey"
"""The column name for the compound InChIKey (used as the canonical compound ID)."""

SMILES_COL = "SMILES"
"""The column name for the compound SMILES (used for featurization)."""

CID_COL = "cid"
"""The column name for the PubChem CID (used for crosswalks and reporting)."""

TYPE_COL = "activity_name"
"""The column name for the measurement type (e.g. IC50, Ki, etc.)."""

VALUE_COL = "activity_value"
"""The column name for the measurement value (numeric)."""

QUALIFIER_COL = "activity_qualifier"
"""The column name for the measurement qualifier (e.g. '=', '<', etc.) (used for filtering to only exact measurements '=' )."""

OUTCOME_COL = "outcome"
"""The column name for the measurement outcome (Active / Inactive)."""

AID_COL = "aid"
"""The column name for the assay ID (used for assay-gate filtering of inactives ensuring inactives are only included if the same experiment reported at least one active)."""


UNIT_COL = "activity_unit"
"""The column name for the measurement unit (e.g. nM, uM, etc.)."""


READ_COLS = [
    GENE_COL,
    KEY_COL,
    SMILES_COL,
    CID_COL,
    TYPE_COL,
    VALUE_COL,
    QUALIFIER_COL,
    AID_COL,
    OUTCOME_COL,
    UNIT_COL,
]


TABLES = ("cliffs", "clusters", "metrics")
"""The three output tables: per-cliff rows, per-compound cluster assignments, and per-(target, fingerprint) metrics summary."""


_ACTIVE_AID_CACHE: dict[str, frozenset[str]] = {}
"""Dictionary mapping the parquet path to the frozenset of assay IDs that report at least one Active row (used for gating inactives.)"""

SPARSE_ABOVE = 30000
""" The threshold above which the Butina clustering will use a sparse representation. """


def _fp_key(radius: int, fp_size: int) -> str:
    """Generate a fingerprint label to identify the featurization config in every output row.

    Args:
        radius: Morgan radius.
        fp_size: Fingerprint bit size.

    Returns:
        A string like "morgan_r2_n2048" that can be used as a column value in
    """
    return f"morgan_r{radius}_n{fp_size}"


def _active_aid_set(path: str) -> frozenset[str]:
    """Return the set of aids that have at least one Active outcome in the parquet.
    Args:
        path: Path to the enriched parquet file.

    Returns:
        The frozenset of aids (as strings) carrying one or more Active outcomes.
    """
    if path in _ACTIVE_AID_CACHE:
        return _ACTIVE_AID_CACHE[path]
    active_aids: set[str] = set()
    scanner = ParquetParser(path).read(columns=[AID_COL, OUTCOME_COL], lazy=True)
    for batch in scanner.to_batches():
        for row in batch.to_pylist():
            if str(row.get(OUTCOME_COL) or "").strip().lower() == "active":
                aid = str(row.get(AID_COL) or "").strip()
                if aid:
                    active_aids.add(aid)
    frozen = frozenset(active_aids)
    _ACTIVE_AID_CACHE[path] = frozen
    logger.info("There are {n} aids report >=1 active.", n=len(frozen))
    return frozen


def _iter_enriched_rows(path: str, columns: list[str]) -> Iterator[dict]:
    """Stream the requested columns of the enriched parquet as row dictionaries where inactives
    are only yielded if their aid appears in the set of aids that report at least one Active outcome.

    Args:
        path: Path to the enriched parquet file.
        columns: The columns to read.

    Yields:
        One row dictionary at a time (inactive rows failing the assay gate are skipped).
    """
    active_aids = _active_aid_set(path)
    read_columns = list(columns)
    gate_only = [c for c in (OUTCOME_COL, AID_COL) if c not in read_columns]
    read_columns += gate_only
    scanner = ParquetParser(path).read(columns=read_columns, lazy=True)
    for batch in scanner.to_batches():
        for row in batch.to_pylist():
            if str(row.get(OUTCOME_COL) or "").strip().lower() == "inactive":
                aid = str(row.get(AID_COL) or "").strip()
                if aid not in active_aids:
                    continue
            if gate_only:
                for c in gate_only:
                    row.pop(c, None)
            yield row


def _write_table(rows: list[dict] | None, path: str) -> None:
    """Write a list of row dictionaries to a parquet file.

    Args:
        rows: The rows to write. If empty, nothing is written and a warning is logged.
        path: The path to write the output parquet file.
    """
    if not rows:
        logger.warning("No rows to write for {p}; skipping.", p=path)
        return
    table = pa.Table.from_pylist(rows)
    ParquetParser(path, schema=table.schema).write(table)
    logger.info("Wrote {n} rows to {p}", n=len(rows), p=path)


def _per_target_dir(out_dir: str, table: str) -> str:
    """Return the directory that holds the per-gene parquets for one table.

    Args:
        out_dir: The run output directory.
        table: One of TABLES.

    Returns:
        The parquet subdirectory path.
    """
    return os.path.join(out_dir, "parquets", table)


def _done_dir(out_dir: str) -> str:
    """Return the directory that holds the per-gene completion markers.

    Args:
        out_dir: The run output directory.

    Returns:
        The marker subdirectory path.
    """
    return os.path.join(out_dir, "parquets", "_done")


def _done_marker(out_dir: str, gene: str) -> str:
    """Return the `.done` marker path for one gene.

    Args:
        out_dir: The output directory of the run.
        gene: The gene id used for the marker.

    Returns:
        The marker file path for that gene.
    """
    return os.path.join(_done_dir(out_dir), f"{gene}.done")


def _is_done(out_dir: str, gene: str) -> bool:
    """Return True if this gene has a completion marker so it was done before.

    Args:
        out_dir: The output directory of the run.
        gene: The gene id.

    Returns:
        True if the gene's `.done` marker exists.
    """
    return os.path.exists(_done_marker(out_dir, gene))


def _write_parquet(rows: list[dict], out_dir: str, table: str, gene: str) -> None:
    """Write one gene's rows for one table to its per-gene parquet parquet.

    Args:
        rows: The rows for this (gene, table).
        out_dir: The output directory of the run.
        table: One of TABLES.
        gene: The gene id to use for the parquet filename.
    """
    parquet_dir = _per_target_dir(out_dir, table)
    os.makedirs(parquet_dir, exist_ok=True)
    _write_table(rows, os.path.join(parquet_dir, f"{gene}.parquet"))


def _mark_done(out_dir: str, gene: str) -> None:
    """Write the `.done` marker for a gene, indicating that all its parquets are done before.

    Args:
        out_dir: The output directory of the run.
        gene: The gene id to use for the marker filename.
    """
    done_dir = _done_dir(out_dir)
    os.makedirs(done_dir, exist_ok=True)
    with open(_done_marker(out_dir, gene), "w", encoding="utf-8") as fh:
        fh.write("")


def _merge_parquets(out_dir: str, table: str) -> None:
    """Merge every per-gene parquet for one table into the combined out_dir/<table>.parquet.

    Args:
        out_dir: The output directory of the run.
        table: One of TABLES.
    """
    parquet_dir = _per_target_dir(out_dir, table)
    if not os.path.isdir(parquet_dir) or not any(
        f.endswith(".parquet") for f in os.listdir(parquet_dir)
    ):
        logger.warning("No {t} parquets to merge; skipping {t}.parquet.", t=table)
        return
    dataset = ParquetParser.open_dataset(parquet_dir)
    merged = dataset.to_table()
    out_path = os.path.join(out_dir, f"{table}.parquet")
    ParquetParser(out_path, schema=merged.schema).write(merged)
    logger.info("Merged {n} {t} rows into {p}", n=merged.num_rows, t=table, p=out_path)


def extract_genes(rows: Iterable[dict]) -> list[str]:
    """Return the distinct gene ids in first-seen order.

    Args:
        rows: The rows to scan for gene ids (must contain the GENE_COL).

    Returns:
        The unique gene ids, in the order first encountered.
    """
    seen: set[str] = set()
    order: list[str] = []
    for row in rows:
        gene = row.get(GENE_COL)
        if gene and str(gene) not in seen:
            seen.add(str(gene))
            order.append(str(gene))
    return order


def count_usable_by_gene(rows: Iterable[dict]) -> dict[str, int]:
    """Count, per gene, the distinct compounds that carry at least one usable measurement
        based on the rules in `[CONVERTIBLE_TYPES][lignova.activity.CONVERTIBLE_TYPES]`.

    Args:
        rows: Enriched rows with at least the gene, key, type, value and qualifier columns.

    Returns:
        A dictionary with the gene ids as the keys and the count of distinct usable compounds
        as the values.
    """
    by_gene: dict[str, set[str]] = {}
    for row in rows:
        gene = row.get(GENE_COL)
        key = row.get(KEY_COL)
        if (
            (not gene or not key)
            or str(row.get(TYPE_COL) or "") not in CONVERTIBLE_TYPES
            or str(row.get(QUALIFIER_COL) or "").strip() != "="
        ):
            continue
        try:
            if float(row.get(VALUE_COL)) <= 0:
                continue
        except (TypeError, ValueError):
            continue
        by_gene.setdefault(str(gene), set()).add(str(key))
    return {gene: len(keys) for gene, keys in by_gene.items()}


def stratified_sample(counts: dict[str, int], n_targets: int, seed: int) -> list[str]:
    """Sample targets in thirds by usable-compound count from top, middle, and bottom used for testing.

    Args:
        counts: A dictionary with Gene id as keys and usable-compound count (from `count_usable_by_gene`) as values.
        n_targets: Total number of targets to sample.
        seed: Random seed for reproducible target selection.

    Returns:
        The sampled gene ids. Fewer than `n_targets` are returned if a stratum is too small.
    """
    ranked = [g for g, _ in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]
    if n_targets >= len(ranked):
        return ranked

    third = len(ranked) // 3
    stratafied_count = [ranked[:third], ranked[third : 2 * third], ranked[2 * third :]]
    base, rem = divmod(n_targets, 3)
    sizes = [base + (1 if i < rem else 0) for i in range(3)]

    rng = random.Random(seed)
    chosen: list[str] = []
    for stratum, size in zip(stratafied_count, sizes):
        chosen.extend(rng.sample(stratum, min(size, len(stratum))))
    logger.info(
        "Stratified sample: {p} targets ({s}) from {t} ranked genes.",
        p=len(chosen),
        s=sizes,
        t=len(ranked),
    )
    return chosen


def _apply_bounds(
    genes: list[str], counts: dict[str, int], lo: int, hi: float
) -> list[str]:
    """Keep only genes whose usable-compound count is within [lo, hi], preserving order.
        genes with no usable compounds (count = 0) are kept only if lo == 0.

    Args:
        genes: The gene ids to filter.
        counts: A dictionary with gene ids as keys and usable-compound counts as values.
        lo: Minimum usable-compound count to keep.
        hi: Maximum usable-compound count to keep.

    Returns:
        The filtered gene ids, in the same order as the input.

    """
    kept = [g for g in genes if lo <= counts.get(g, 0) <= hi]
    dropped = len(genes) - len(kept)
    if dropped:
        logger.info(
            "Usable-size bounds [{lo}, {hi}]: kept {k}, dropped {d} of {t} genes.",
            lo=lo,
            hi=hi,
            k=len(kept),
            d=dropped,
            t=len(genes),
        )
    return kept


def select_genes(path: str, args: argparse.Namespace) -> list[str]:
    """Select which gene ids to process, per user selection mode and any usable-compound count bounds.

    Args:
        path: Path to the enriched parquet.
        args: Parsed CLI args

    Returns:
        The selected gene ids.
    """
    lo = args.min if args.min is not None else 0
    hi = args.max if args.max is not None else math.inf
    bounded = not (lo == 0 and hi == math.inf)

    if args.select == "all":
        genes = extract_genes(_iter_enriched_rows(path, [GENE_COL]))
        if bounded:
            counts = count_usable_by_gene(_iter_enriched_rows(path, READ_COLS))
            genes = _apply_bounds(genes, counts, lo, hi)
        return genes

    if args.select == "range":
        start, end = args.range
        genes = extract_genes(_iter_enriched_rows(path, [GENE_COL]))
        if start > len(genes):
            logger.info(
                "Range selected start {s} exceeds {t} genes. Nothing to do.",
                s=start,
                t=len(genes),
            )
            return 3
        if bounded:
            counts = count_usable_by_gene(_iter_enriched_rows(path, READ_COLS))
            genes = _apply_bounds(genes, counts, lo, hi)
        selected = genes[start - 1 : end]
        logger.info(
            "Range select: genes {a}-{b} of {t} leading to {n} targets.",
            a=start,
            b=end,
            t=len(genes),
            n=len(selected),
        )
        return selected

    counts = count_usable_by_gene(_iter_enriched_rows(path, READ_COLS))
    if bounded:
        counts = {g: c for g, c in counts.items() if lo <= c <= hi}
        logger.info(
            "Usable-size bounds [{lo}, {hi}]: {n} genes eligible for stratification.",
            lo=lo,
            hi=hi,
            n=len(counts),
        )
    return stratified_sample(counts, args.n_targets, args.seed)


def run_target(
    gene: str,
    rows: list[dict],
    params: CliffParams,
    fp_sizes: list[int] = [2048],
    cutoff: float = 0.55,
    radius: int = 2,
    spread_gate: float = DEFAULT_SPREAD_GATE,
    standardize: bool = True,
) -> tuple[list[dict], list[dict], list[dict]]:
    """Run clustering and cliff detection for one target across every user specified fingerprint size.

    Args:
        gene: The gene id being processed.
        rows: The enriched rows for this gene only.
        params: The cliff detection and severity configuration (metric, min_delta, thresholds, min_similarity).
        fp_sizes: The fingerprint bit sizes to evaluate.Default is [2048].
        cutoff: Butina similarity cutoff (used in clustering). Default is 0.55.
        radius: Morgan radius (shared across the fingerprint sizes). Default is 2.
        spread_gate: The quality gate fails when the winning-type spread exceeds this (in log units). Default is DEFAULT_SPREAD_GATE =1.
        standardize: Whether to standardize the SMILES to neutral parent forms when handling multiple SMILES variants for the same InChIKey. Default is True.

    Returns:
        A tuple of (cliff_rows, cluster_rows, metric_rows), each tagged with the gene id and the
        fingerprint key.
    """
    grouped = group_measurements(rows)
    activity = aggregate_activity(grouped, spread_gate=spread_gate)
    xref = cid_crosswalk(grouped)
    variants: dict[str, set[str]] = defaultdict(set)
    for r in rows:
        k, s = r.get(KEY_COL), r.get(SMILES_COL)
        if k and s:
            variants[str(k)].add(str(s))
    resolved = {
        k: resolve_smiles(v, standardize=standardize) for k, v in variants.items()
    }
    smiles = {k: chosen for k, (chosen, _) in resolved.items()}

    n_conflicts = sum(1 for v in variants.values() if len(v) > 1)
    if n_conflicts:
        methods = Counter(m for k, (_, m) in resolved.items() if len(variants[k]) > 1)
        logger.info(
            "Resolved {n} multi-SMILES keys for {g}: {m}",
            n=n_conflicts,
            g=gene,
            m=dict(methods),
        )

    cliff_rows: list[dict] = []
    cluster_rows: list[dict] = []
    metric_rows: list[dict] = []

    for fp_size in fp_sizes:
        fp_key = _fp_key(radius, fp_size)
        t0 = time.perf_counter()
        feats = MorganFeaturizer(radius=radius, fp_size=fp_size).featurize(smiles)
        dense = len(feats.items) <= SPARSE_ABOVE
        sims = compute_pairwise(
            feats.items, min_sim=min(params.min_similarity, cutoff), dense=dense
        )
        clusters = ButinaClustering(
            ButinaParams(similarity_cutoff=cutoff, sparse_above=SPARSE_ABOVE)
        ).cluster(sims)
        result = label_cliff_severity(find_activity_cliffs(sims, activity, params))

        cliff_rows.extend(_cliff_rows(gene, fp_key, result, xref, smiles))
        cluster_rows.extend(_cluster_rows(gene, fp_key, clusters))
        metric_rows.append(
            _metric_row(
                gene,
                fp_key,
                feats,
                sims,
                clusters,
                result,
                params,
                time.perf_counter() - t0,
                spread_gate,
                standardize_conflicts=standardize,
                n_smiles_conflicts=n_conflicts,
            )
        )

    return cliff_rows, cluster_rows, metric_rows


def _cliff_rows(
    gene: str, fp_key: str, result, xref: dict[str, set[str]], smiles: dict[str, str]
) -> list[dict]:
    """Flatten a labelled CliffResult into per-cliff rows for cliffs.parquet keyed by InChIKey.

    Args:
        gene: The gene id being processed.
        fp_key: The fingerprint key (e.g. "morgan_r2_n2048").
        result: The labelled CliffResult.
        xref: The compound key to source CIDs map (from `cid_crosswalk`).
        smiles: The compound key to chosen SMILES map (from `resolve_smiles`).

    Returns:
        A list of dictionaries, one per cliff, with the gene and fp_key added.
    """
    out: list[dict] = []
    for c in result.cliffs:
        out.append(
            {
                "gene": gene,
                "fp": fp_key,
                "id_a": c.id_a,
                "id_b": c.id_b,
                "smiles_a": smiles.get(c.id_a),
                "smiles_b": smiles.get(c.id_b),
                "cids_a": ";".join(sorted(xref.get(c.id_a, set()))),
                "cids_b": ";".join(sorted(xref.get(c.id_b, set()))),
                "similarity": c.similarity,
                "pAct_diff": c.pAct_diff,
                "landscape_index": c.landscape_index,
                "severity": c.severity_label.value if c.severity_label else None,
                "is_cross_type": c.is_cross_type,
                "is_low_quality": c.is_low_quality,
                "is_same_assay": c.is_same_assay,
                "involves_inactive": c.involves_inactive,
            }
        )
    return out


def _cluster_rows(gene: str, fp_key: str, clusters) -> list[dict]:
    """Flatten a ClusterResult into per-compound rows for clusters.parquet.

    Args:
        gene: The gene id being processed.
        fp_key: The fingerprint key (e.g. "morgan_r2_n2048").
        clusters: The ClusterResult.

    Returns:
        A list of dictionaries, one per compound, with the gene and fp_key added.In addition
        a `cluster_uid` where gene, fp_key, and cluster_id are concatenated to make a globally unique cluster identifier.
    """
    reps = set(clusters.representatives.values())
    return [
        {
            "gene": gene,
            "fp": fp_key,
            "inchikey": key,
            "cluster_id": cid,
            "cluster_uid": f"{gene}:{fp_key}:{cid}",
            "is_representative": key in reps,
        }
        for key, cid in clusters.labels.items()
    ]


def _metric_row(
    gene,
    fp_key,
    feats,
    sims,
    clusters,
    result,
    params,
    seconds,
    spread_gate,
    standardize_conflicts,
    n_smiles_conflicts,
) -> dict:
    """Assemble one per-(target, fingerprint) summary row for metrics.parquet.

    Args:
        gene: The gene id being processed.
        fp_key: The fingerprint key (e.g. "morgan_r2_n2048").
        feats: The MorganFeaturizer result.
        sims: The TanimotoSimilarities result.
        clusters: The ButinaClustering result.
        result: The labelled CliffResult.
        params: The CliffParams used for the cliff detection and severity labelling.
        seconds: Wall-clock seconds for this one (target, fp) pass; for spotting slow configs at scale.
        spread_gate: The quality gate threshold for the winning-type spread.
        standardize_conflicts: Whether SMILES were standardized to neutral parent forms when resolving multiple SMILES variants for the same InChIKey.
        n_smiles_conflicts: The number of InChIKeys that had multiple SMILES variants that were resolved to a single representative SMILES.

    Returns:
        A dictionary with the metrics for this (gene, fp_key) pass.
    """
    sizes = [len(v) for v in clusters.clusters().values()]
    return {
        "gene": gene,
        "fp": fp_key,
        "metric": params.metric.value,
        "standardize_conflicts": standardize_conflicts,
        "n_smiles_conflicts": n_smiles_conflicts,
        "spread_gate": spread_gate,
        "min_similarity": params.min_similarity,
        "butina_cutoff": clusters.params.similarity_cutoff,
        "min_delta": params.min_delta,
        "n_featurized": len(feats.items),
        "n_skipped": len(feats.skipped),
        "n_edges": len(sims.edges),
        # These are the pairs that make SALI undefined
        "n_identical_pairs": sum(1 for *_, s in sims.edges if s >= 1.0),
        "sali_undefined": (
            params.sali_undefined.value
            if params.metric is SeverityMetric.SALI
            else None
        ),
        "n_clusters": clusters.n_clusters,
        "n_singletons": sum(1 for s in sizes if s == 1),
        "n_cliffs": result.n_cliffs,
        "n_undefined": result.n_undefined,
        "n_high_confidence": len(high_confidence_cliffs(result.cliffs)),
        "n_same_assay": len(same_assay_cliffs(result.cliffs)),
        "n_involves_inactive": sum(1 for c in result.cliffs if c.involves_inactive),
        "n_extreme": _band(result, "extreme"),
        "n_strong": _band(result, "strong"),
        "n_moderate": _band(result, "moderate"),
        # wall-clock seconds for this one (target, fp) pass; for spotting slow configs at scale.
        "runtime_seconds": round(seconds, 3),
    }


def _band(result, name: str) -> int:
    """Count cliffs in a named severity band."""
    return sum(
        1 for c in result.cliffs if c.severity_label and c.severity_label.value == name
    )


def main() -> int:
    """Parse arguments, select targets, run the evaluation, and write the three output tables."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "-p", "--parquet", required=True, help="Path to the enriched parquet."
    )
    ap.add_argument(
        "-o", "--out", required=True, help="Output directory for the three tables."
    )
    ap.add_argument(
        "-s",
        "--select",
        choices=["all", "range", "stratified"],
        default="all",
        help="Which targets to process. Default: all. use range for parallelization, stratified for a small test sample.",
    )
    ap.add_argument(
        "-r",
        "--range",
        type=int,
        nargs=2,
        metavar=("START", "END"),
        help="For --select range: 1-based inclusive slice of genes in file order.",
    )
    ap.add_argument(
        "-n",
        "--n-targets",
        type=int,
        default=9,
        help="For --select stratified: total targets to sample (in thirds).",
    )
    ap.add_argument(
        "-se",
        "--seed",
        type=int,
        default=0,
        help="For --select stratified: seed for reproducible sampling.",
    )
    ap.add_argument(
        "-m",
        "--min",
        type=int,
        default=None,
        help="Exclude genes with fewer than this many usable compounds "
        "(applies to all selection modes). Default: no lower bound.",
    )
    ap.add_argument(
        "-ma",
        "--max",
        type=int,
        default=None,
        help="Exclude genes with more than this many usable compounds. Useful for a "
        "sweep to skip a giant target whose O(n^2) pass would dominate runtime. "
        "Default: no upper bound.",
    )
    ap.add_argument(
        "-f",
        "--fp-sizes",
        type=int,
        nargs="+",
        default=[2048],
        help="One or more fingerprint bit sizes to evaluate (the sweep).",
    )
    ap.add_argument(
        "--radius", type=int, default=2, help="Morgan radius for all fp sizes."
    )
    ap.add_argument(
        "-fl", "--floor", type=float, default=0.55, help="Similarity floor for edges."
    )
    ap.add_argument(
        "-c", "--cutoff", type=float, default=0.55, help="Butina similarity cutoff."
    )
    ap.add_argument(
        "-d",
        "--min-delta",
        type=float,
        default=2.0,
        help="Minimum pActivity gap for a cliff.",
    )
    ap.add_argument(
        "-up",
        "--undefined-sali",
        choices=[p.value for p in SaliUndefined],
        default=SaliUndefined.NEXT_LARGEST.value,
        help="SALI policy for Tanimoto==1 (undefined) pairs: 'next_largest' (default), "
        "'max_severity' (rank above all finite cliffs), or 'exclude' (drop them).",
    )

    ap.add_argument(
        "-sp",
        "--spread-gate",
        type=float,
        default=DEFAULT_SPREAD_GATE,
        help="Quality gate: winning-type pActivity spread (log units) above which a compound "
        "fails the gate. Default: aggregate_activity's built-in default.",
    )
    ap.add_argument(
        "-i",
        "--metric",
        choices=[m.value for m in SeverityMetric],
        default=SeverityMetric.TS_SALI.value,
    )
    ap.add_argument(
        "--no-resume",
        action="store_true",
        help="Reprocess every selected gene even if a .done marker exists so existing parquets are overwritten. "
        "Default: skip genes already marked done.",
    )
    ap.add_argument(
        "--merge-only",
        action="store_true",
        help="Skip all per-target work and just merge whatever parquets already exist into the "
        "three combined tables.",
    )

    ap.add_argument(
        "--no-merge",
        action="store_true",
        help="Do not merge parquets after processing is complete. Useful for parallel runs where "
        "each task writes its own per-gene parquets and a single final merge is done later ",
    )
    ap.add_argument(
        "--standardize",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Neutralize and salt-strip InChIKeys that carry multiple SMILES. "
        "Default: on. Use --no-standardize to pick deterministically instead.",
    )
    args = ap.parse_args()

    if args.select == "range" and not args.range:
        ap.error("--select range requires --range START END")
    if args.select == "range" and (args.range[0] < 1 or args.range[1] < args.range[0]):
        ap.error("--range needs 1 <= START <= END")
    if args.merge_only and args.no_merge:
        ap.error("--merge-only and --no-merge are mutually exclusive")
    params = CliffParams(
        min_delta=args.min_delta,
        metric=SeverityMetric(args.metric),
        sali_undefined=SaliUndefined(args.undefined_sali),
        min_similarity=args.floor,
    )
    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    if args.merge_only:
        for table in TABLES:
            _merge_parquets(out_dir, table)
        logger.info("merge-only complete.")
        return 0

    selected = select_genes(args.parquet, args)
    if not selected:
        logger.error("No targets selected; nothing to do.")
        return 1

    if args.no_resume:
        todo = list(selected)
    else:
        todo = [g for g in selected if not _is_done(out_dir, g)]
        skipped = len(selected) - len(todo)
        if skipped:
            logger.info(
                "Resume: {s} of {t} selected genes already done; {n} to process.",
                s=skipped,
                t=len(selected),
                n=len(todo),
            )

    if not todo:
        if args.no_merge:
            logger.info(
                "All selected genes already processed; skipping merge (--no-merge)."
            )
        else:
            logger.info("All selected genes already processed; merging parquets.")
            for table in TABLES:
                _merge_parquets(out_dir, table)
        return 0

    if args.select == "all":
        logger.warning(
            "--select all buffers rows for every gene ({n}) in memory. For a very large parquet, "
            "prefer a gene-sorted stream or process in ranges.",
            n=len(todo),
        )

    todo_set = set(todo)
    per_gene: dict[str, list[dict]] = {g: [] for g in todo}
    for row in _iter_enriched_rows(args.parquet, READ_COLS):
        g = str(row.get(GENE_COL) or "")
        if g in todo_set:
            per_gene[g].append(row)

    processed = 0
    n_cliffs = 0
    for gene, rows in per_gene.items():
        logger.info("Processing gene {g} and its {n} rows.", g=gene, n=len(rows))
        cliffs, clusters, metrics = run_target(
            gene=gene,
            rows=rows,
            fp_sizes=args.fp_sizes,
            cutoff=args.cutoff,
            params=params,
            radius=args.radius,
            spread_gate=args.spread_gate,
            standardize=args.standardize,
        )
        _write_parquet(cliffs, out_dir, "cliffs", gene)
        _write_parquet(clusters, out_dir, "clusters", gene)
        _write_parquet(metrics, out_dir, "metrics", gene)
        _mark_done(out_dir, gene)
        processed += 1
        n_cliffs += len(cliffs)

    if args.no_merge:
        logger.info(
            "Processed {n} gene(s) this run; skipping merge (--no-merge).", n=processed
        )
    else:
        logger.info("Processed {n} gene(s) this run so merging parquets.", n=processed)
        for table in TABLES:
            _merge_parquets(out_dir, table)
    logger.info(
        "Done: {t} target(s) this run x {f} fp sizes produced {c} total cliffs.",
        t=processed,
        f=len(args.fp_sizes),
        c=n_cliffs,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
