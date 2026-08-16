# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Bulk-download PubChem bioactivities, then enrich compound properties for
every CID of a chosen outcome across the WHOLE parquet (all actives or all
inactives), writing one combined Parquet and caching per-CID properties.

It offers two modes:
  * Assay mode  (--aids ...) : enrich the given assays, filtered by outcome.
      Output: one row per assay record, WITH activity data (activity_value,
      activity_name, activity_unit, aid, target, pubmed_id) plus properties.
  * Global mode (--all)      : scan the whole parquet for distinct CIDs of the
      chosen outcome and enrich them all.
      Output: one row per CID, compound descriptors and activity data.
"""

import argparse
import asyncio
import os
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from loguru import logger

from lignova.APIs.pubchem import DEFAULT_PROPERTIES, AssayInfo, PubChemAPI, PubChemBulk
from lignova.hdf5 import ParquetParser

BIOACTIVITIES_GZ = "bioactivities.tsv.gz"
BIOACTIVITIES_PARQUET = "bioactivities.parquet"
_BATCH = 150  # CIDs per PubChem property request


DEFAULT_CACHE_NAME = "properties.parquet"
DEFAULT_OUT_TEMPLATE = "enriched_{filter}.parquet"

# Values in possible to be extracted as string properties from pubchem API
_STRING_PROPERTIES = frozenset(
    {
        "SMILES",
        "ConnectivitySMILES",
        "InChI",
        "InChIKey",
        "IUPACName",
        "MolecularFormula",
        "Title",
    }
)


# map from the eaw bioactivities tsv to the output parquet columns
_DUMP_TO_OUT = {
    "CID": "cid",
    "Activity Outcome": "outcome",
    "Activity Value": "activity_value",
    "Activity Name": "activity_name",
    "Activity Unit": "activity_unit",
    "Activity Qualifier": "activity_qualifier",
    "AID": "aid",
    "Protein Accession": "target_accession",
    "Gene ID": "target_geneid",
    "PMID": "pubmed_id",
}

# Output activity columns
_ACTIVITY_OUT = list(_DUMP_TO_OUT.values())

_NUMERIC_ACTIVITY = {"cid", "activity_value"}
_STRING_ACTIVITY = [c for c in _ACTIVITY_OUT if c not in _NUMERIC_ACTIVITY]


def _to_float(v: Any) -> float | None:
    """Coerce a cell to float, mapping blanks/None/unparseable to None.

    Args:
        v: The value to coerce.

    Returns:
        The value as a float, or None if it is blank, None, or unparseable.
    """
    if v is None:
        return None
    if isinstance(v, str):
        v = v.strip()
        if not v:
            return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


class PropertyCache:
    """Class to hold a persistent CID to properties mapping that is backed up to a directory of parquet
    parts on disk holding by default 50 batches of 150 CIDs each.
    """

    def __init__(self, path: str, required: list[str] = DEFAULT_PROPERTIES) -> None:
        """Initialize the cache, loading any existing parts from disk.

        Args:
            path: Directory holding the property-cache parts.
            required: Property names that define the on-disk column schema. Defaults to [DEFAULT_PROPERTIES][lignova.APIs.pubchem.DEFAULT_PROPERTIES].
        """
        self.path = str(path)
        self.required = required
        self._rows: dict[int, dict[str, Any]] = {}
        self._pending: list[dict[str, Any]] = []
        self._part_seq = 0
        self._load()

    def _normalize(self, props: dict[str, Any]) -> dict[str, Any]:
        """Normalize the property dict to consistent per-column types.

        Args:
            props: Raw property mapping, possibly including a CID key.

        Returns:
            A new dict with string properties as str/"" and the rest as float/None.
        """
        clean: dict[str, Any] = {}
        for k, v in props.items():
            if k == "CID":
                continue
            if k in _STRING_PROPERTIES:
                clean[k] = "" if v is None else str(v)
            else:
                clean[k] = _to_float(v)
        return clean

    def _load_df(self, df: pd.DataFrame) -> None:
        """Load a DataFrame's rows into the in-memory dictionary.

        Args:
            df: A cache part with a CID column.
        """
        if "CID" not in df.columns:
            return
        for row in df.to_dict("records"):
            cid = row.get("CID")
            if cid is None or (isinstance(cid, float) and pd.isna(cid)):
                continue
            cid = int(cid)
            props = self._normalize(row)
            self._rows[cid] = props

    def _load(self) -> None:
        """Populate the in-memory index from the cache directory."""
        if os.path.isdir(self.path):
            parts = sorted(f for f in os.listdir(self.path) if f.endswith(".parquet"))
            self._part_seq = len(parts)
            for name in parts:
                self._load_df(
                    ParquetParser(os.path.join(self.path, name)).convert_to_pandas()
                )
            logger.info("Loaded property cache: {n} CIDs.", n=len(self._rows))
            return
        logger.info("No property cache at {p}; starting empty.", p=self.path)

    def cached_cids(self) -> set[int]:
        """Return all CIDs currently in the cache.

        Returns:
            A set of cached CIDs
        """
        return set(self._rows)

    def get(self, cid: int) -> dict[str, Any] | None:
        """Return a CID's cached properties, or None if not present.

        Args:
            cid: The compound ID to look up.

        Returns:
            The property dict for the CID, or None if it isn't cached.
        """
        return self._rows.get(int(cid))

    def add_many(self, fetched: dict[int, dict[str, Any]]) -> None:
        """Add fetched CID with its properties to the cache and queue them for saving.

        Args:
            fetched: Mapping of CID to its raw property dict from the API.
        """
        for cid, props in fetched.items():
            cid = int(cid)
            clean = self._normalize(props)
            self._rows[cid] = clean
            self._pending.append({"CID": cid, **clean})

    def save(self) -> None:
        """Split the pending rows into a new parquet part and write it to disk."""
        if not self._pending:
            return
        os.makedirs(self.path, exist_ok=True)
        part = pd.DataFrame(self._pending).reindex(columns=["CID", *self.required])
        part["CID"] = part["CID"].astype("int64")
        for col in self.required:
            if col in _STRING_PROPERTIES:
                part[col] = part[col].astype(object).where(part[col].notna(), "")
            else:
                part[col] = pd.to_numeric(part[col], errors="coerce")
        out = os.path.join(self.path, f"part-{self._part_seq:06d}.parquet")
        part.to_parquet(out, index=False)
        self._part_seq += 1
        self._pending = []
        logger.info("Saved property cache: {n} CIDs.", n=len(self._rows))


def distinct_cids_for_outcome(parquet_path: str, outcome: str) -> list[int]:
    """Scan the whole parquet for distinct CIDs of one outcome.

    Streams only the CID + Activity Outcome columns and deduplicates in memory

    Args:
        parquet_path: The bioactivities parquet.
        outcome: e.g. "active" or "inactive" (case-insensitive).

    Returns:
        Sorted list of distinct CIDs with that outcome.
    """
    target = outcome.strip().lower()
    parser = ParquetParser(parquet_path)
    scanner = parser.read(columns=["CID", "Activity Outcome"], lazy=True)
    seen: set[int] = set()
    for batch in scanner.to_batches():
        cids = batch.column("CID").to_pylist()
        outs = batch.column("Activity Outcome").to_pylist()
        for cid, oc in zip(cids, outs):
            if cid in (None, "") or oc is None:
                continue
            if oc.strip().lower() == target:
                try:
                    seen.add(int(cid))
                except (TypeError, ValueError):
                    continue
    result = sorted(seen)
    logger.info(
        "Found {n} distinct '{outcome}' CIDs in the parquet.",
        n=len(result),
        outcome=outcome,
    )
    return result


async def ensure_parquet(
    bulk: PubChemBulk, workdir: str, force: bool, verify_md5: bool = False
) -> str:
    """Ensure the bioactivities parquet exists, downloading if needed.

    Args:
        bulk: The bulk client used to download and convert the bioactivities.
        workdir: Directory holding the .tsv.gz and its parquet.
        force: Re-download and re-convert to parquet even if the files already exist.
        verify_md5: Verify the downloaded file's MD5 after download. Defaults to False.

    Returns:
        Path to the bioactivities parquet.
    """
    gz = os.path.join(workdir, BIOACTIVITIES_GZ)
    pq = os.path.join(workdir, BIOACTIVITIES_PARQUET)
    await bulk.download(gz, force=force, verify_md5=verify_md5)
    if os.path.exists(pq) and not force:
        logger.info("Parquet present; skipping ingest.")
    else:
        bulk.to_parquet(gz, pq)
    return pq


async def enrich_cids(
    api: PubChemAPI,
    cids: list[int],
    properties: list[str],
    cache: PropertyCache,
    checkpoint_every: int = 50,
) -> None:
    """Fetch and cache properties for uncached CIDs, concurrently.

    Args:
        api: The PubChem API client.
        cids: CIDs to enrich where those already cached are skipped.
        properties: Property tags to request per CID.
        cache: The property cache to read from and write into.
        checkpoint_every: Save the cache every N completed batches. Defaults to 50.
    """
    have = cache.cached_cids()
    todo = [c for c in cids if c not in have]
    logger.info(
        "{tot} CIDs requested, {have} cached, {new} to fetch.",
        tot=len(cids),
        have=len(cids) - len(todo),
        new=len(todo),
    )
    batches = [todo[i : i + _BATCH] for i in range(0, len(todo), _BATCH)]
    done = 0

    tasks = [asyncio.create_task(api.get_properties(b, properties)) for b in batches]
    for coro in asyncio.as_completed(tasks):
        fetched = await coro
        cache.add_many(fetched)
        done += 1
        if done % max(checkpoint_every, 1) == 0:
            cache.save()
            logger.info(
                "Checkpoint: {done}/{total} batches.", done=done, total=len(batches)
            )
    cache.save()


def stream_measurements_for_outcome(
    parquet_path: str,
    outcome: str,
    cache: PropertyCache,
    out_path: str,
    allowed_cids: set[int],
    chunk_rows: int = 500000,
) -> int:
    """Generate one output row per measurement of `outcome`, joined with properties.

    Args:
        parquet_path: The created bioactivities parquet.
        outcome: Outcome to keep (case-insensitive), e.g. "active".
        cache: Property cache queried for each CID's descriptors.
        out_path: Destination parquet for the enriched per-measurement table.
        allowed_cids: Only rows whose CID is in this set are emitted.
        chunk_rows: Rows buffered before each incremental write.

    Returns:
        Total number of rows written.
    """
    target = outcome.strip().lower()
    prop_cols = list(cache.required)
    out_columns = [*_ACTIVITY_OUT, *prop_cols]

    parser = ParquetParser(parquet_path)
    scanner = parser.read(columns=list(_DUMP_TO_OUT), lazy=True)

    writer = None
    buffer: list[dict[str, Any]] = []
    written = 0

    def flush() -> None:
        nonlocal writer, buffer, written
        if not buffer:
            return
        frame = pd.DataFrame(buffer).reindex(columns=out_columns)
        frame["cid"] = pd.to_numeric(frame["cid"], errors="coerce").astype("int64")
        frame["activity_value"] = pd.to_numeric(
            frame["activity_value"], errors="coerce"
        )
        for col in _STRING_ACTIVITY:
            frame[col] = frame[col].astype(object).where(frame[col].notna(), "")
        for col in prop_cols:
            if col in _STRING_PROPERTIES:
                frame[col] = frame[col].astype(object).where(frame[col].notna(), "")
            else:
                frame[col] = pd.to_numeric(frame[col], errors="coerce")
        table = pa.Table.from_pandas(frame, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(out_path, table.schema)
        writer.write_table(table)
        written += len(frame)
        buffer = []

    for batch in scanner.to_batches():
        d = batch.to_pydict()
        cids = d["CID"]
        outs = d["Activity Outcome"]
        for i, raw_cid in enumerate(cids):
            oc = outs[i]
            if raw_cid in (None, "") or oc is None:
                continue
            if oc.strip().lower() != target:
                continue
            try:
                cid = int(raw_cid)
            except (TypeError, ValueError):
                continue
            if cid not in allowed_cids:
                continue
            av = d["Activity Value"][i]
            # keep only rows with a non-blank activity value
            if av in (None, ""):
                continue
            row = {
                "cid": cid,
                "outcome": oc,
                "activity_value": av,
                "activity_name": d["Activity Name"][i],
                "activity_unit": d["Activity Unit"][i],
                "activity_qualifier": d["Activity Qualifier"][i],
                "aid": d["AID"][i],
                "target_accession": d["Protein Accession"][i],
                "target_geneid": d["Gene ID"][i],
                "pubmed_id": d["PMID"][i],
            }
            props = cache.get(cid)
            if props:
                row.update(props)
            buffer.append(row)
            if len(buffer) >= chunk_rows:
                flush()
    flush()
    if writer is not None:
        writer.close()
    return written


def _filter_by_outcome(assay: AssayInfo, outcome: str | None) -> AssayInfo:
    """Return a copy of the assay keeping only records of one outcome.

    Args:
        assay: The assay to filter.
        outcome: Outcome to keep (case-insensitive); None keeps everything.

    Returns:
        An AssayInfo with the filtered records.
    """
    if outcome is None:
        return assay
    t = outcome.strip().lower()
    kept = [r for r in assay.records if r.outcome.strip().lower() == t]
    return AssayInfo(aid=assay.aid, records=kept)


def build_cid_table(
    cids: list[int], outcome: str, cache: PropertyCache
) -> pd.DataFrame:
    """Build a one-row-per-CID table of the outcome label plus cached properties.

    Args:
        cids: CIDs to emit, in order.
        outcome: Outcome label written into every row (e.g. "inactive").
        cache: Property cache queried for each CID's descriptors.

    Returns:
        A DataFrame with columns cid, outcome, and one column per property.
    """
    rows: list[dict[str, Any]] = []
    for cid in cids:
        row: dict[str, Any] = {"cid": cid, "outcome": outcome}
        props = cache.get(cid)
        if props:
            row.update({k.lower(): v for k, v in props.items()})
        rows.append(row)
    return pd.DataFrame(rows)


async def main() -> int:
    """Parse args, build the cache, and run the enrichment pipeline.

    Returns:
        Process exit code: 0 on success, 1 if no rows were produced.
    """
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--aids", type=int, nargs="+", help="Enrich these assays.")
    mode.add_argument(
        "--all",
        action="store_true",
        help="Enrich all CIDs of the specified outcome across the whole parquet.",
    )
    p.add_argument(
        "--filter",
        "-f",
        type=str,
        choices=["active", "inactive", "all"],
        default="active",
        help="Which outcome to filter on (default: active).",
    )
    p.add_argument(
        "--workdir",
        "-w",
        type=str,
        required=True,
        help="Directory to hold the bioactivities parquet and property cache.",
    )
    p.add_argument(
        "--out",
        "-o",
        type=str,
        default=None,
        help=(
            "Output parquet for the enriched data "
            f"(default: <workdir>/{DEFAULT_OUT_TEMPLATE.format(filter='<filter>')})."
        ),
    )
    p.add_argument(
        "--cache",
        "-c",
        type=str,
        default=None,
        help=(
            "Directory for the property cache "
            f"(default: <workdir>/{DEFAULT_CACHE_NAME})."
        ),
    )
    p.add_argument(
        "--extra-properties",
        "-e",
        nargs="*",
        default=None,
        help="Additional PubChem properties to fetch beyond the defaults.",
    )
    p.add_argument(
        "--limit",
        "-l",
        type=int,
        default=0,
        help="Limit for the CIDs to process where 0 means no limit.",
    )
    p.add_argument(
        "--checkpoint",
        type=int,
        default=50,
        help="Save the cache every N batches.",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Force re-download of bioactivities and re-convert to parquet.",
    )
    p.add_argument(
        "--verify-md5",
        action="store_true",
        default=False,
        help="Verify the MD5 checksum of the downloaded bioactivities file.",
    )
    args = p.parse_args()

    if args.all and args.filter == "all":
        p.error("--all needs a concrete --filter (active or inactive).")

    outcome = None if args.filter == "all" else args.filter
    out_path = args.out
    if out_path is None:
        out_path = os.path.join(
            args.workdir, DEFAULT_OUT_TEMPLATE.format(filter=args.filter)
        )
    cache_path = args.cache
    if cache_path is None:
        cache_path = os.path.join(args.workdir, DEFAULT_CACHE_NAME)
    properties = list(dict.fromkeys(DEFAULT_PROPERTIES + (args.extra_properties or [])))
    cache = PropertyCache(cache_path, required=properties)

    async with PubChemBulk() as bulk:
        parquet_path = await ensure_parquet(
            bulk, args.workdir, args.force, args.verify_md5
        )

        async with PubChemAPI() as api:
            if args.all:
                cids = distinct_cids_for_outcome(parquet_path, args.filter)
                if args.limit:
                    cids = cids[: args.limit]
                    logger.info("Limited to first {n} CIDs.", n=len(cids))
                await enrich_cids(api, cids, properties, cache, args.checkpoint)
                parent = os.path.dirname(out_path)
                if parent:
                    os.makedirs(parent, exist_ok=True)
                n = stream_measurements_for_outcome(
                    parquet_path,
                    args.filter,
                    cache,
                    out_path,
                    allowed_cids=set(cids),
                )
                if n == 0:
                    logger.warning("No rows produced; nothing written.")
                    return 1
                logger.info("Wrote {n} rows to {p}", n=n, p=out_path)
                return 0
            else:
                assays = []
                wanted = set()
                for aid in args.aids:
                    a = _filter_by_outcome(bulk.load_assay(aid, parquet_path), outcome)
                    if a.records:
                        assays.append(a)
                        wanted.update(a.unique_cids)
                cids = sorted(wanted)
                if args.limit:
                    cids = cids[: args.limit]
                await enrich_cids(
                    api,
                    cids,
                    properties,
                    cache,
                    args.checkpoint,
                )
                rows = []
                for a in assays:
                    for r in a.records:
                        row = {
                            "aid": a.aid,
                            "cid": r.cid,
                            "outcome": r.outcome,
                            "activity_value": r.activity_value,
                            "activity_name": r.activity_name,
                            "activity_unit": r.activity_unit,
                            "activity_qualifier": r.activity_qualifier,
                            "target_accession": r.target_accession,
                            "target_geneid": r.target_geneid,
                            "pubmed_id": r.pubmed_id,
                        }
                        if r.cid is not None and (pr := cache.get(r.cid)):
                            row.update(pr)
                        rows.append(row)
                df = pd.DataFrame(rows)

    if df.empty:
        logger.warning(
            "No rows produced for the given assays and outcome thus nothing written."
        )
        return 1
    parent = os.path.dirname(out_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    df.to_parquet(out_path)
    logger.info("Wrote {n} rows to {p}", n=len(df), p=out_path)
    return 0


if __name__ == "__main__":
    asyncio.run(main())
