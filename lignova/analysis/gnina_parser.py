# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Parser for GNINA docked SDF files and dataset builder."""

# NOTE: this is for SDF files output ONLY
import glob
import gzip
import os
import re
import shutil
import tempfile
from collections.abc import Iterator

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pcsv
import pyarrow.dataset as pds
from loguru import logger

from ..hdf5.parquet import ParquetParser

DOCKING_SCHEMA = pa.schema(
    [
        pa.field("block_start", pa.int64()),
        pa.field("block_end", pa.int64()),
        pa.field("protein_id", pa.string()),
        pa.field("ligand_id", pa.string()),
        pa.field("UniqueID", pa.int32()),
        pa.field("conformer_idx", pa.int32()),
        pa.field("pose_rank", pa.int32()),
        pa.field("SMILES", pa.string()),
        pa.field("n_atoms", pa.int32()),
        pa.field("Vina_affinity", pa.float64()),
        pa.field("CNNscore", pa.float64()),
        pa.field("CNNaffinity", pa.float64()),
        pa.field("CNN_VS", pa.float64()),
        pa.field("CNNaffinity_variance", pa.float64()),
        pa.field("Energy", pa.float64()),
        pa.field("source_file", pa.string()),
    ]
)

_GNINA_FLOAT_PROPS = [
    "CNNscore",
    "CNNaffinity",
    "CNN_VS",
    "CNNaffinity_variance",
]

# Defines the best direction for each score:
SCORE_DIRECTIONS: dict[str, str] = {
    "CNNscore": "descending",
    "CNNaffinity": "descending",
    "CNN_VS": "descending",
    "Vina_affinity": "ascending",
    "Energy": "ascending",
    "CNNaffinity_variance": "ascending",
}

_PROP_RE = re.compile(r">\s+<([^>]+)>\n([^\n]*)")
_VALID_EXTENSIONS = (".sdf", ".sdf.gz")


class TruncatedSDFError(Exception):
    """Raised when an SDF file appears truncated and allow_truncated=False."""

    def __init__(self, filepath: str, reason: str):
        self.filepath = filepath
        self.reason = reason
        super().__init__(f"Truncated SDF: {filepath} — {reason}")


class GNINA_Results:
    r"""Lazy-indexed parser for a single GNINA docked SDF file."""

    def __init__(
        self,
        filepath: str,
        num_modes: int | list[int] | None = 9,
        protein_id: str | None = None,
        keep_raw: bool = False,
        allow_truncated: bool = True,
    ):
        r"""Initialize the GNINA_Results parser.
        Args:
            filepath: Path to the .sdf or .sdf.gz file containing the docking results.
            num_modes: The number of modes per ligand/protomer. Can be specified as:
                - int: A single value applied to all groups.
                - list[int]: A list of known mode counts to match against group block counts.
                - None : Auto-detect from block count patterns.
                Default is 9, which is the standard for GNINA docking outputs.
            protein_id: Optional identifier for the target protein. If not provided,
                it will be inferred from the file name.
            keep_raw: Whether to keep the full decompressed text in memory after
                parsing. Default is False, the raw text is released after the table is built
                and a seekable file on disk is used for subsequent block retrieval.
            allow_truncated: Whether to allow files that appear truncated (e.g. missing some expected blocks).
        """
        if not os.path.isfile(filepath):
            raise ValueError(f"File {filepath} does not exist.")
        if not any(filepath.endswith(ext) for ext in _VALID_EXTENSIONS):
            raise ValueError(f"File {filepath} is not a .sdf or .sdf.gz file.")
        self._filepath = filepath
        if protein_id is None:
            logger.warning(
                "No protein_id provided, attempting to infer from file path."
            )
        self._protein_id = protein_id or self._infer_protein_id()
        self._offsets: list[tuple[int, int]] = []
        self._raw: str | None = ""
        self._seekable_path: str | None = None
        self._owns_seekable: bool = False
        self._allow_truncated = allow_truncated
        self._load_offsets()

        self._group_num_modes: int | dict[tuple[str, int], int]
        self._table, self._group_num_modes = self._build_table(num_modes)
        if not keep_raw:
            self._ensure_seekable()
            self._raw = None
            logger.debug(
                f"Released raw text from {os.path.basename(filepath)} into {self._seekable_path}"
            )

    def _ensure_seekable(self) -> None:
        r"""Ensure there is a non-gzipped .sdf file on disk so we can seek into for on-demand block retrieval."""
        if self._seekable_path is not None:
            return

        if not self._filepath.endswith(".gz"):
            self._seekable_path = self._filepath
            self._owns_seekable = False
        else:
            fd, tmp_path = tempfile.mkstemp(suffix=".sdf", prefix="gnina_")
            os.close(fd)
            with open(tmp_path, "w", encoding="utf-8") as out:
                out.write(self._raw)
            self._seekable_path = tmp_path
            self._owns_seekable = True
            logger.debug(
                f"Decompressed {os.path.basename(self._filepath)} to {tmp_path}"
            )

    def _build_table(
        self,
        num_modes: int | list[int] | None,
    ) -> tuple[pa.Table, int | dict[tuple[str, int], int]]:
        r"""Parse all blocks, resolve num_modes, and build the PyArrow table in
        a single pass over the block data.
        Args:
            num_modes: number of poses per protormer (int, list[int], or None).
        Returns:
            A tuple of:
            - The built PyArrow Table with all docking data.
            - The resolved num_modes (int if uniform, per-group dict otherwise).
        """
        n = len(self._offsets)

        metas: list[dict[str, str | float | int]] = []
        group_counts: dict[tuple[str, int], int] = {}

        for idx in range(n):
            block = self._get_block(idx)
            raw = self._parse_block_meta(block)
            metas.append(raw)

            if not isinstance(num_modes, int):
                lig_id = raw.get("mol_name", "")
                try:
                    uid = int(float(raw.get("UniqueID", "0")))
                except (ValueError, TypeError):
                    uid = 0
                key = (lig_id, uid)
                group_counts[key] = group_counts.get(key, 0) + 1

        if isinstance(num_modes, int):
            resolved: int | dict[tuple[str, int], int] = num_modes
        else:
            resolved = self._resolve_modes(
                group_counts,
                modes=num_modes,
            )
        uniform = resolved if isinstance(resolved, int) else None

        col: dict[str, np.ndarray | list[str]] = {
            "block_start": np.array([s for s, _ in self._offsets], dtype=np.int64),
            "block_end": np.array([e for _, e in self._offsets], dtype=np.int64),
            "protein_id": [self._protein_id] * n,
            "ligand_id": [""] * n,
            "UniqueID": np.zeros(n, dtype=np.int32),
            "conformer_idx": np.zeros(n, dtype=np.int32),
            "pose_rank": np.zeros(n, dtype=np.int32),
            "SMILES": [""] * n,
            "n_atoms": np.zeros(n, dtype=np.int32),
            "Vina_affinity": np.full(n, np.nan, dtype=np.float64),
            "CNNscore": np.full(n, np.nan, dtype=np.float64),
            "CNNaffinity": np.full(n, np.nan, dtype=np.float64),
            "CNN_VS": np.full(n, np.nan, dtype=np.float64),
            "CNNaffinity_variance": np.full(n, np.nan, dtype=np.float64),
            "Energy": np.full(n, np.nan, dtype=np.float64),
            "source_file": [self._filepath] * n,
        }

        group_seq: dict[tuple[str, int], int] = {}

        for idx, raw in enumerate(metas):
            lig_id = raw.get("mol_name", "")
            try:
                uid = int(float(raw.get("UniqueID", "0")))
            except (ValueError, TypeError):
                uid = 0

            key = (lig_id, uid)
            seq = group_seq.get(key, 0)
            group_seq[key] = seq + 1

            nm = uniform if uniform is not None else resolved[key]

            col["ligand_id"][idx] = lig_id
            col["UniqueID"][idx] = uid
            col["conformer_idx"][idx] = seq // nm
            col["pose_rank"][idx] = seq % nm
            col["SMILES"][idx] = raw.get("SMILES", "")
            col["n_atoms"][idx] = raw.get("n_atoms", 0)

            # minimizedAffinity → Vina_affinity
            try:
                col["Vina_affinity"][idx] = float(raw["minimizedAffinity"])
            except (KeyError, ValueError, TypeError):
                pass

            for prop in _GNINA_FLOAT_PROPS:
                try:
                    col[prop][idx] = float(raw[prop])
                except (KeyError, ValueError, TypeError):
                    pass
            try:
                col["Energy"][idx] = float(raw.get("Energy", "nan"))
            except (ValueError, TypeError):
                pass

        if uniform is not None:
            for key, cnt in group_seq.items():
                if cnt < uniform:
                    logger.warning(
                        f"{key[0]}:UID{key[1]} has {cnt}/{uniform} blocks (truncated)"
                    )
                elif cnt % uniform != 0:
                    logger.warning(
                        f"{key[0]}:UID{key[1]} has {cnt} blocks, not divisible by num_modes={uniform}"
                    )

        arrays = [pa.array(col[f.name], type=f.type) for f in DOCKING_SCHEMA]
        table = pa.Table.from_arrays(arrays, schema=DOCKING_SCHEMA)
        return table, resolved

    def _resolve_modes(
        self,
        gc: dict[tuple[str, int], int],
        modes: list[int] | None = None,
    ) -> dict[tuple[str, int], int]:
        r"""Resolve per-group num_modes by finding the largest candidate that
        evenly divides each group's block count, or falling back to the closest
        candidate by smallest remainder.
        Args:
            gc: A dictionary mapping (ligand_id, UniqueID) to block counts.
            modes: Optional list of known mode counts. If None, the distinct
                values in gc are used as candidates (auto-detect).
        Returns:
            A dictionary mapping:
            - Keys: Tuples of (ligand_id, UniqueID) representing each group.
            - Values: The resolved num_modes for each group.
        """
        if not gc:
            raise ValueError("No groups found for num_modes resolution.")

        counts = list(gc.values())

        if modes is None and len(set(counts)) == 1:
            nm = counts[0]
            logger.info(f"Auto-detected num_modes={nm}")
            return {k: nm for k in gc}

        candidates = sorted(
            set(modes) if modes is not None else set(counts),
            reverse=True,
        )
        result: dict[tuple[str, int], int] = {}

        for key, cnt in gc.items():
            assigned = False
            for m in candidates:
                if cnt % m == 0:
                    result[key] = m
                    assigned = True
                    break

            if not assigned:
                best = min(candidates, key=lambda m: cnt % m)
                result[key] = best
                if cnt < best:
                    logger.warning(
                        f"{key[0]}:UID{key[1]} has {cnt}/{best} blocks (truncated)"
                    )
                else:
                    logger.warning(
                        f"{key[0]}:UID{key[1]} has {cnt} blocks, not evenly divisible by num_modes={best} "
                    )

        self._log_mode_summary(result)

        if modes is None:
            detected = sorted(set(result.values()))
            if len(detected) > 1:
                logger.warning(f"Auto-detected multiple num_modes: {detected}. ")

        return result

    def _log_mode_summary(
        self,
        mapping: dict[tuple[str, int], int],
    ) -> None:
        r"""Log a summary of how many groups were assigned to each num_modes value.
        Args:
            mapping: A dictionary mapping (ligand_id, UniqueID) to assigned num_modes.
        """
        by_mode: dict[int, list[tuple[str, int]]] = {}
        for key, m in mapping.items():
            by_mode.setdefault(m, []).append(key)
        for m in sorted(by_mode):
            keys = by_mode[m]
            logger.info(
                f"num_modes={m}: {len(keys)} group(s) — "
                f"{', '.join(f'{k[0]}:UID{k[1]}' for k in keys)}"
            )

    def _read_raw(self) -> str:
        r"""Read the entire decompressed text of the SDF file into memory."""
        if self._filepath.endswith(".sdf.gz"):
            with gzip.open(self._filepath, "rb") as f:
                return f.read().decode("utf-8")
        with open(self._filepath, "r", encoding="utf-8") as f:
            return f.read()

    def _get_block(self, block_idx: int) -> str:
        r"""Retrieve a block by slicing into the raw string (if available)
        or by seeking into the seekable file on disk."""
        start, end = self._offsets[block_idx]
        if self._raw is not None:
            return self._raw[start:end]
        # Fall back to disk seek
        with open(self._seekable_path, "r", encoding="utf-8") as f:
            f.seek(start)
            return f.read(end - start)

    def cleanup(self) -> None:
        r"""Remove the temporary seekable file if we created one."""
        if getattr(self, "_owns_seekable", False) and self._seekable_path is not None:
            try:
                os.remove(self._seekable_path)
                logger.debug(f"Cleaned up temp file: {self._seekable_path}")
            except OSError:
                pass
            self._seekable_path = None
            self._owns_seekable = False

    def __del__(self) -> None:
        r"""clean up if in case cleanup() wasn't called explicitly."""
        self.cleanup()

    def __enter__(self):
        r"""Support `with GNINA_Results(...) as dr:` usage."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        r"""Clean up temp file when exiting the context."""
        self.cleanup()
        return False

    def release_raw(self) -> None:
        r"""Manually release the raw text to free memory.
        Subsequent block access will read from disk.
        """
        if self._raw is None:
            return
        self._ensure_seekable()
        self._raw = None

    @property
    def is_raw_loaded(self) -> bool:
        r"""Whether the full decompressed text is currently held in memory."""
        return self._raw is not None

    def _infer_protein_id(self) -> str:
        r"""Infer the protein_id from the file name by taking the substring before the first underscore."""
        name = os.path.basename(self._filepath)
        if "_" in name:
            return name.split("_")[0]
        return ""

    def _load_offsets(self) -> None:
        r"""Read raw file and index block boundaries."""
        self._raw = self._read_raw()
        delimiter = "$$$$"
        start = 0
        while True:
            end = self._raw.find(delimiter, start)
            if end == -1:
                trailing = self._raw[start:].strip()
                if trailing:
                    self._offsets.append((start, start + len(trailing)))
                break
            block_text = self._raw[start:end].strip()
            if block_text:
                actual_start = start
                while actual_start < end and self._raw[actual_start] in (
                    "\n",
                    "\r",
                    " ",
                ):
                    actual_start += 1
                actual_end = end
                while actual_end > actual_start and self._raw[actual_end - 1] in (
                    "\n",
                    "\r",
                    " ",
                ):
                    actual_end -= 1
                if actual_start < actual_end:
                    self._offsets.append((actual_start, actual_end))
            start = end + len(delimiter)
        if not self._allow_truncated:
            stripped = self._raw.rstrip()
            if not self._offsets:
                raise TruncatedSDFError(self._filepath, "No SDF blocks found")
            if not stripped.endswith("$$$$"):
                raise TruncatedSDFError(
                    self._filepath,
                    "File does not end with $$$$ delimiter (last block incomplete)",
                )

    @staticmethod
    def _parse_block_meta(block_text: str) -> dict[str, str | float | int]:
        r"""Parse the metadata properties from an SDF block string.
        Args:
            block_text: The raw text of a single SDF block.
        Returns:
            A dictionary containing the molecule name, atom count, and any properties defined in the block.
        """
        lines = block_text.split("\n")
        mol_name = lines[0].strip() if lines else ""
        n_atoms = 0
        if len(lines) > 3:
            try:
                n_atoms = int(lines[3][:3].strip())
            except (ValueError, IndexError):
                pass
        props = {}
        for m in _PROP_RE.finditer(block_text):
            props[m.group(1)] = m.group(2).strip()
        return {"mol_name": mol_name, "n_atoms": n_atoms, **props}

    @staticmethod
    def _parse_coords_from_block(
        block_text: str,
    ) -> tuple[np.ndarray, list[str]]:
        r"""Parse atomic coordinates and element symbols from an SDF block string.
        Args:
            block_text: The raw text of a single SDF block.
        Returns:
            A tuple containing:
            - A NumPy array of shape (n_atoms, 3) with the atomic coordinates.
            - A list of element symbols corresponding to each atom.
        """
        lines = block_text.split("\n")
        try:
            n_atoms = int(lines[3][:3].strip())
        except (ValueError, IndexError):
            raise ValueError("Cannot parse atom count from block")
        coords = np.empty((n_atoms, 3), dtype=np.float64)
        elems = []
        for i in range(n_atoms):
            parts = lines[4 + i].split()
            coords[i, 0] = float(parts[0])
            coords[i, 1] = float(parts[1])
            coords[i, 2] = float(parts[2])
            elems.append(parts[3])
        return coords, elems

    @staticmethod
    def _parse_mol_from_block(block_text: str):
        r"""Parse an SDF block string into an RDKit molecule object.
        Args:
            block_text: The raw text of a single SDF block.
        Returns:
            An RDKit molecule object representing the structure in the block.
        """
        from rdkit import Chem

        text = block_text + "\n$$$$\n"
        suppl = Chem.SDMolSupplier()
        suppl.SetData(text, removeHs=False)
        mols = [m for m in suppl if m is not None]
        if not mols:
            raise ValueError("RDKit failed to parse block")
        return mols[0]

    @staticmethod
    def best_direction(score: str) -> str:
        r"""Return the sort direction ('ascending' or 'descending') that
        corresponds to 'best' for a given score column.
        Args:
            score: The column name (e.g. 'CNNscore', 'Vina_affinity').
        Returns:
            'ascending' or 'descending'.
        """
        if score not in SCORE_DIRECTIONS:
            raise ValueError(
                f"Unknown score {score!r}. " f"Known scores: {sorted(SCORE_DIRECTIONS)}"
            )
        return SCORE_DIRECTIONS[score]

    def get_block_by_offsets(self, block_start: int, block_end: int) -> str:
        r"""Retrieve a block by byte offsets into the decompressed raw string or from disk if raw has been released.
        Args:
            block_start: Start byte offset of the block.
            block_end: End byte offset of the block.
        Returns:
            The raw text of the SDF block.
        """
        if self._raw is not None:
            return self._raw[block_start:block_end]
        with open(self._seekable_path, "r", encoding="utf-8") as f:
            f.seek(block_start)
            return f.read(block_end - block_start)

    @staticmethod
    def read_block_from_file(
        source_file: str,
        block_start: int,
        block_end: int,
    ) -> str:
        r"""Read a single SDF block directly from a file using stored byte offsets.
        For .sdf files this is an O(1) seek. For .sdf.gz files the entire file
        must be decompressed first (offsets are into decompressed text).
        Args:
            source_file: Path to the .sdf or .sdf.gz file.
            block_start: Start byte offset of the block in the decompressed text.
            block_end: End byte offset of the block in the decompressed text.
        Returns:
            The raw text of the SDF block.
        """
        if source_file.endswith(".sdf.gz"):
            with gzip.open(source_file, "rt", encoding="utf-8") as f:
                raw = f.read()
            return raw[block_start:block_end]
        with open(source_file, "r", encoding="utf-8") as f:
            f.seek(block_start)
            return f.read(block_end - block_start)

    @property
    def num_modes(self) -> list[int]:
        r"""The unique num_modes values found across all groups, sorted.

        Returns a sorted list of distinct mode counts.
        E.g. [9] if uniform, [7, 9] if one group was truncated.
        """
        if isinstance(self._group_num_modes, int):
            return [self._group_num_modes]
        return sorted(set(self._group_num_modes.values()))

    @property
    def num_modes_per_group(self) -> dict[tuple[str, int], int]:
        r"""Per-group num_modes mapping.

        Returns a dictionary mapping (ligand_id, UniqueID): num_modes
        for every group in the file.
        """
        if isinstance(self._group_num_modes, int):
            # Build full dict from the table when a uniform int was used
            t = self._table
            keys = set()
            for i in range(t.num_rows):
                keys.add(
                    (t.column("ligand_id")[i].as_py(), t.column("UniqueID")[i].as_py())
                )
            return {k: self._group_num_modes for k in keys}
        return dict(self._group_num_modes)

    @property
    def table(self) -> pa.Table:
        r"""The full PyArrow Table containing all parsed docking data."""
        return self._table

    @property
    def schema(self) -> pa.Schema:
        r"""The PyArrow Schema for the docking results table."""
        return DOCKING_SCHEMA

    @property
    def filepath(self) -> str:
        r"""The original file path of the SDF or SDF.GZ file that was parsed."""
        return self._filepath

    @property
    def protein_id(self) -> str:
        r"""The identifier for the target protein, either provided or inferred from the file name."""
        return self._protein_id

    @property
    def blocks(self) -> list[str]:
        r"""All raw SDF blocks as a list of strings. This will load the entire raw text into memory if not already loaded."""
        return [self._get_block(i) for i in range(len(self._offsets))]

    @property
    def n_blocks(self) -> int:
        r"""The total number of SDF blocks (poses) parsed from the file."""
        return len(self._offsets)

    @property
    def n_ligands(self) -> int:
        r"""Count of distinct ligand_id values in the table."""
        return pc.count_distinct(self._table.column("ligand_id")).as_py()

    @property
    def n_protomers(self) -> int:
        r"""Combine ligand_id and UniqueID to count distinct protomers"""
        lig = self._table.column("ligand_id")
        uid = self._table.column("UniqueID")
        combined = pc.binary_join_element_wise(
            lig, pc.cast(uid, pa.string()), pa.scalar("_")
        )
        return pc.count_distinct(combined).as_py()

    @property
    def n_conformers(self) -> int:
        r"""Combine ligand_id, UniqueID, and conformer_idx to count distinct conformers per protomer
        (e.g., if num_modes=9, expect 9 conformers per protomer)."""
        t = self._table
        combined = pc.binary_join_element_wise(
            t.column("ligand_id"),
            pc.cast(t.column("UniqueID"), pa.string()),
            pc.cast(t.column("conformer_idx"), pa.string()),
            pa.scalar("_"),
        )
        return pc.count_distinct(combined).as_py()

    def _filter(self, expr: pc.Expression) -> pa.Table:
        r"""Apply a PyArrow compute expression filter to the table and return the filtered result.
        Args:
            expr: A PyArrow compute expression that evaluates to a boolean mask.
        Returns:
            A PyArrow Table containing only the rows where the expression is True.
        """
        return self._table.filter(expr)

    def get_ligand(self, ligand_id: str) -> pa.Table:
        r"""Get all poses for a specific ligand_id (across all protomers and conformers).
        Args:
            ligand_id: The ligand identifier (e.g., "M3A") to filter by.
        Returns:
            A PyArrow Table containing all poses for the specified ligand_id.
        """
        return self._filter(pc.field("ligand_id") == ligand_id)

    def get_protomer(self, ligand_id: str, unique_id: int) -> pa.Table:
        r"""Get all poses for a specific protomer (ligand_id + UniqueID).
        Args:
            ligand_id: The ligand identifier (e.g., "M3A") to filter by.
            unique_id: The UniqueID to filter by, which distinguishes different protomers of the same ligand.
        Returns:
            A PyArrow Table containing all poses for the specified protomer.
        """
        return self._filter(
            (pc.field("ligand_id") == ligand_id) & (pc.field("UniqueID") == unique_id)
        )

    def get_conformer(
        self,
        ligand_id: str,
        unique_id: int,
        conformer_idx: int = 0,
    ) -> pa.Table:
        r"""Get all poses for a specific conformer/pose.
        Args:
            ligand_id: The ligand identifier (e.g., "M3A") to filter by.
            unique_id: The UniqueID to filter by, which distinguishes different protomers of the same ligand.
            conformer_idx: The conformer index to filter by (default 0).
        Returns:
            A PyArrow Table containing the pose info.
        """
        return self._filter(
            (pc.field("ligand_id") == ligand_id)
            & (pc.field("UniqueID") == unique_id)
            & (pc.field("conformer_idx") == conformer_idx)
        )

    def get_top_poses(
        self,
        by: str | list[str] = "CNNscore",
        ascending: bool | list[bool] | None = None,
        weights: list[float] | None = None,
        n: int = 1,
        per: str = "conformer",
    ) -> pa.Table:
        r"""Get the top N poses sorted by one or more score columns
        Args:
            by: Column name or list of column names to sort/rank by.
            ascending: Sort direction(s). Can be:
                - None (default): auto-detect from SCORE_DIRECTIONS based on the 'by' column.
                - bool: applied to all columns in by.
                - list[bool]: one per column in by.
            weights: Optional list of floats (one per column in `by`). When
                provided, a weighted composite score is computed. Weights are
                normalised to sum to 1. Only valid when `by` is a list.
            n: The number of top poses to return per group. Default is 1.
            per: The grouping level for selecting top poses. Options are
                'conformer', 'protomer', or 'ligand'.
        Returns:
            A PyArrow Table containing the top N poses per specified group.
        """
        group_cols = {
            "conformer": ["ligand_id", "UniqueID", "conformer_idx"],
            "protomer": ["ligand_id", "UniqueID"],
            "ligand": ["ligand_id"],
        }
        keys = group_cols.get(per, group_cols["conformer"])

        if isinstance(by, str):
            by = [by]
        if ascending is None:
            ascending = [self.best_direction(c) == "ascending" for c in by]
        elif isinstance(ascending, bool):
            ascending = [ascending] * len(by)
        if len(ascending) != len(by):
            raise ValueError(
                f"Length mismatch: {len(by)} columns but {len(ascending)} ascending flags"
            )

        if len(by) == 1 and weights is None:
            sort_order = "ascending" if ascending[0] else "descending"
            indices = pc.sort_indices(self._table, sort_keys=[(by[0], sort_order)])
            sorted_t = self._table.take(indices)
            return self._top_n_per_group(sorted_t, keys, n)

        if weights is None:
            weights = [1.0] * len(by)
        if len(weights) != len(by):
            raise ValueError(
                f"Length mismatch: {len(by)} columns but {len(weights)} weights"
            )

        w_sum = sum(weights)
        w_norm = [w / w_sum for w in weights]

        composite = np.zeros(self._table.num_rows, dtype=np.float64)
        for col_name, asc, w in zip(by, ascending, w_norm):
            arr = self._table.column(col_name).to_numpy().astype(np.float64)
            lo = np.nanmin(arr)
            hi = np.nanmax(arr)
            if hi - lo > 0:
                normed = (arr - lo) / (hi - lo)
            else:
                normed = np.zeros_like(arr)
            # Flip so higher = better for all columns
            if asc:
                normed = 1.0 - normed
            composite += w * normed

        idx = np.argsort(-composite)
        sorted_t = self._table.take(idx)
        return self._top_n_per_group(sorted_t, keys, n)

    @staticmethod
    def _top_n_per_group(
        sorted_table: pa.Table,
        group_keys: list[str],
        n: int,
    ) -> pa.Table:
        r"""Select the first N rows per group from a pre-sorted table.
        Args:
            sorted_table: A PyArrow Table already sorted by the desired score.
            group_keys: Column names that define the grouping.
            n: Number of rows to keep per group.
        Returns:
            A PyArrow Table with at most n rows per group.
        """
        seen: dict[tuple, int] = {}
        keep = []
        for i in range(sorted_table.num_rows):
            row_key = tuple(sorted_table.column(k)[i].as_py() for k in group_keys)
            count = seen.get(row_key, 0)
            if count < n:
                keep.append(i)
                seen[row_key] = count + 1
        return sorted_table.take(keep)

    def get_best_per_ligand(
        self,
        by: str | list[str] = "CNNscore",
        ascending: bool | list[bool] | None = None,
        weights: list[float] | None = None,
    ) -> pa.Table:
        r"""Convenience method to get the single best pose per ligand_id.
        Args:
            by: Column name or list of column names to rank by.
            ascending: Sort direction(s). None = auto-detect from score type.
            weights: Optional weights for composite scoring (one per column in `by`).
        Returns:
            A PyArrow Table containing the best pose per ligand_id."""
        return self.get_top_poses(
            by=by, ascending=ascending, weights=weights, n=1, per="ligand"
        )

    def filter(self, **kwargs) -> pa.Table:
        expr = None
        for k, v in kwargs.items():
            if k.endswith("_min"):
                e = pc.field(k[:-4]) >= v
            elif k.endswith("_max"):
                e = pc.field(k[:-4]) <= v
            else:
                e = pc.field(k) == v
            expr = e if expr is None else (expr & e)
        if expr is None:
            return self._table
        return self._filter(expr)

    def get_coords(
        self, block_start: int, block_end: int
    ) -> tuple[np.ndarray, list[str]]:
        r"""Parse atomic coordinates and element symbols from a block.
        Args:
            block_start: Start byte offset of the block.
            block_end: End byte offset of the block.
        Returns:
            A tuple containing:
            - A NumPy array of shape (n_atoms, 3) with the atomic coordinates.
            - A list of element symbols corresponding to each atom.
        """
        return self._parse_coords_from_block(
            self.get_block_by_offsets(block_start, block_end)
        )

    def get_coords_batch(
        self,
        offsets: list[tuple[int, int]],
    ) -> list[tuple[np.ndarray, list[str]]]:
        r"""Get coordinates for a batch of blocks.
        Args:
            offsets: List of (block_start, block_end) tuples.
        Returns:
            A list of (coords_array, elements_list) tuples.
        """
        return [self.get_coords(bs, be) for bs, be in offsets]

    def get_mol(self, block_start: int, block_end: int):
        r"""Parse a block into an RDKit molecule object.
        Args:
            block_start: Start byte offset of the block.
            block_end: End byte offset of the block.
        Returns:
            An RDKit molecule object.
        """
        return self._parse_mol_from_block(
            self.get_block_by_offsets(block_start, block_end)
        )

    def summary(
        self,
        per: str = "global",
        output: str | None = None,
    ) -> None:
        r"""write a human-readable summary of the docking results
        Args:
            per: The level of detail for the summary. Options are:
                - 'global': Overall summary for the entire dataset.
                - 'ligand': Summary broken down by ligand_id, showing score ranges and best pose per ligand.
                - 'protomer': Summary broken down by protomer (ligand_id + UniqueID), showing score ranges and best pose per protomer.
            output: The output destination for the summary. Options are:
                - None (default): Print to console using loguru.
                - str: A file path to write the summary to.
        """
        t = self._table
        score_cols = list(SCORE_DIRECTIONS.keys())
        lines = [
            f"File:           {os.path.basename(self._filepath)}",
            f"Protein:        {self._protein_id}",
            f"Total blocks:   {self.n_blocks}",
            f"Unique ligands: {self.n_ligands}",
            f"Protomers:      {self.n_protomers}",
            f"Conformers:     {self.n_conformers}",
        ]

        def _ranges(tgt, sub, indent="  "):
            for c in score_cols:
                lo = pc.min(sub.column(c)).as_py()
                hi = pc.max(sub.column(c)).as_py()
                if lo is not None:
                    tgt.append(f"{indent}{c:25s} {lo:.4f}  to  {hi:.4f}")

        def _best(tgt, sub, label="Best pose (by CNNscore)"):
            mx = pc.max(sub.column("CNNscore"))
            row = pc.index(sub.column("CNNscore"), mx).as_py()
            if row >= 0:
                r = sub.slice(row, 1)
                bs = r.column("block_start")[0].as_py()
                be = r.column("block_end")[0].as_py()
                tgt.append(f"  {label}:")
                tgt.append(
                    f"    offset={bs}:{be}, "
                    f"UID={r.column('UniqueID')[0].as_py()}, "
                    f"CNNscore={r.column('CNNscore')[0].as_py():.4f}, "
                    f"CNNaffinity={r.column('CNNaffinity')[0].as_py():.3f}, "
                    f"CNN_VS={r.column('CNN_VS')[0].as_py():.3f}, "
                    f"Energy={r.column('Energy')[0].as_py():.3f}, "
                    f"Vina={r.column('Vina_affinity')[0].as_py():.2f}"
                )

        if per == "global":
            lines += ["", "Score ranges (all poses):"]
            _ranges(lines, t)
            lines.append("")
            _best(lines, t)
        elif per == "ligand":
            ligs = pc.unique(t.column("ligand_id")).to_pylist()
            lines += ["", f"Score ranges per ligand ({len(ligs)} ligand(s)):"]
            for lig in ligs:
                sub = self._filter(pc.field("ligand_id") == lig)
                n_uid = pc.count_distinct(sub.column("UniqueID")).as_py()
                lines += ["", f"  [{lig}]  poses={sub.num_rows}  protomers={n_uid}"]
                _ranges(lines, sub, "    ")
                _best(lines, sub)
        elif per == "protomer":
            ligs = pc.unique(t.column("ligand_id")).to_pylist()
            lines += ["", "Score ranges per protomer:"]
            for lig in ligs:
                lig_sub = self._filter(pc.field("ligand_id") == lig)
                uids = pc.unique(lig_sub.column("UniqueID")).to_pylist()
                lines += ["", f"  Ligand: {lig}"]
                for uid in sorted(uids):
                    sub = self._filter(
                        (pc.field("ligand_id") == lig) & (pc.field("UniqueID") == uid)
                    )
                    smi = sub.column("SMILES")[0].as_py()
                    lines += ["", f"    UID={uid}  poses={sub.num_rows}  SMILES={smi}"]
                    _ranges(lines, sub, "      ")
                    _best(lines, sub)
        else:
            raise ValueError(
                f"Unknown per={per!r}. Choose 'global', 'ligand', or 'protomer'."
            )

        text = "\n".join(lines)
        if output is None:
            logger.info(text)
        elif isinstance(output, str):
            with open(output, "w", encoding="utf-8") as fh:
                fh.write(text + "\n")
            logger.info(f"Summary written to {output}")
        else:
            output.write(text + "\n")

    def to_csv(self, output_path: str, overwrite: bool = False) -> None:
        r"""Export the docking table to a CSV file.
        Args:
            output_path: Path to the output CSV file.
            overwrite: If False (default), skip if the file already exists.
        """
        if not overwrite and os.path.isfile(output_path):
            logger.info(f"Skipping CSV export, file exists: {output_path}")
            return
        pcsv.write_csv(self._table, output_path)
        logger.info(f"Exported {self.n_blocks} rows to {output_path}")

    def to_sdf(
        self,
        output_path: str,
        table: pa.Table | None = None,
        overwrite: bool = False,
    ) -> None:
        r"""Export SDF blocks to a file.
        Args:
            output_path: Path to the output .sdf or .sdf.gz file.
            table: A PyArrow Table with block_start/block_end columns.
                None = export all blocks.
            overwrite: If False (default), skip if the file already exists.
        """
        if not overwrite and os.path.isfile(output_path):
            logger.info(f"Skipping SDF export, file exists: {output_path}")
            return
        if table is None:
            table = self._table
        blocks = []
        for i in range(table.num_rows):
            bs = table.column("block_start")[i].as_py()
            be = table.column("block_end")[i].as_py()
            blocks.append(self.get_block_by_offsets(bs, be))
        text = "\n$$$$\n".join(blocks) + "\n$$$$\n"
        if output_path.endswith(".sdf.gz"):
            with gzip.open(output_path, "wt", encoding="utf-8") as f:
                f.write(text)
        else:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(text)
        logger.info(f"Wrote {len(blocks)} blocks to {output_path}")

    def poses_per_conformer(self) -> pa.Array:
        r"""Count the number of poses per unique conformer (ligand_id + UniqueID + conformer_idx).
        Returns:
            A PyArrow Array where each element corresponds to the count of poses for a unique conformer.
            should be 9 for each conformer if num_modes=9."""
        t = self._table
        combined = pc.binary_join_element_wise(
            t.column("ligand_id"),
            pc.cast(t.column("UniqueID"), pa.string()),
            pc.cast(t.column("conformer_idx"), pa.string()),
            pa.scalar("_"),
        )
        return pc.value_counts(combined)

    def protomer_counts(self) -> pa.Table:
        r"""Count the number of poses and conformers per protomer (ligand_id + UniqueID)."""
        t = self._table.group_by(["ligand_id", "UniqueID", "SMILES"]).aggregate(
            [("block_start", "count"), ("conformer_idx", "count_distinct")]
        )
        return t.rename_columns(
            ["ligand_id", "UniqueID", "SMILES", "n_blocks", "n_conformers"]
        )

    def __repr__(self) -> str:
        r"""Return a string representation of the GNINA_Results object, including file name and summary stats."""
        return (
            f"GNINA_Results({os.path.basename(self._filepath)!r}, "
            f"blocks={self.n_blocks}, ligands={self.n_ligands}, "
            f"protomers={self.n_protomers}, "
            f"conformers={self.n_conformers})"
        )

    def __len__(self) -> int:
        r"""override len() to return the number of blocks (poses) in the dataset."""
        return self.n_blocks


class DockedPose:
    r"""Lightweight container for a single docked pose."""

    __slots__ = [
        "block_start",
        "block_end",
        "ligand_id",
        "protein_id",
        "UniqueID",
        "conformer_idx",
        "pose_rank",
        "SMILES",
        "Energy",
        "Vina_affinity",
        "CNNscore",
        "CNNaffinity",
        "CNN_VS",
        "CNNaffinity_variance",
        "_dr",
    ]

    def __init__(self, row_idx: int, table: pa.Table, dr: GNINA_Results):
        r"""Initialize a DockedPose by reading the specified row from the table and storing
        a reference to the GNINA_Results for block access.
        Args:
            row_idx: The index of the row in the table to read.
            table: The PyArrow Table containing the docking data.
            dr: The GNINA_Results instance to reference for block access.
        """
        for name in self.__slots__:
            if name == "_dr":
                continue
            try:
                setattr(self, name, table.column(name)[row_idx].as_py())
            except (KeyError, IndexError):
                setattr(self, name, None)
        self._dr = dr

    @property
    def block_text(self) -> str:
        r"""Get the raw SDF block text for this pose using stored byte offsets.
        Returns:
            The raw text of the SDF block for this pose.
        """
        return self._dr.get_block_by_offsets(self.block_start, self.block_end)

    @property
    def coords(self) -> tuple[np.ndarray, list[str]]:
        r"""Get the atomic coordinates and element symbols for this pose.
        Returns:
            A tuple containing:
            - A NumPy array of shape (n_atoms, 3) with the atomic coordinates.
            - A list of element symbols corresponding to each atom."""
        return GNINA_Results._parse_coords_from_block(self.block_text)

    @property
    def mol(self):
        r"""Get the RDKit molecule object for this pose.
        Returns:
            An RDKit molecule object representing the structure of this pose."""
        return GNINA_Results._parse_mol_from_block(self.block_text)

    def __repr__(self) -> str:
        cnn = f"{self.CNNscore:.4f}" if self.CNNscore is not None else "None"
        return (
            f"DockedPose(lig={self.ligand_id!r}, UID={self.UniqueID}, "
            f"conf={self.conformer_idx}, rank={self.pose_rank}, "
            f"CNN={cnn})"
        )


def as_poses(table: pa.Table, dr: GNINA_Results) -> list[DockedPose]:
    r"""Convert a PyArrow Table of docking results into a list of DockedPose objects.
    Args:
        table: A PyArrow Table containing the docking data, with required columns.
        dr: The GNINA_Results instance to reference for block access.
    Returns:
            A list of DockedPose objects, one for each row in the table.
    """
    return [DockedPose(i, table, dr) for i in range(table.num_rows)]


class DockingDataset:
    r"""Manages per protein docking dataset and converts to batched Parquet

    one parquet per protein layout::

        dataset_dir/
        ├── parquet/{PROTEIN_ID}.parquet
        ├── blocks/{PROTEIN_ID}/*.sdf   (optional, only if copy_blocks=True)
        └── proteins/{PROTEIN_ID}.pdb   (optional copy of cleaned receptor PDBs, only if copy_proteins=True)

    Parameters
    ----------
    dataset_dir : str
        Root of the dataset directory.
    """

    def __init__(self, dataset_dir: str):
        self._root = dataset_dir
        self._parquet_dir = os.path.join(self._root, "parquet")
        self._blocks_dir = os.path.join(self._root, "blocks")
        self._proteins_dir = os.path.join(self._root, "proteins")

    def _parser_for(self, protein_id: str) -> ParquetParser:
        """Return a ParquetParser pointed at a protein's parquet file."""
        pq_path = os.path.join(self._parquet_dir, f"{protein_id}.parquet")
        return ParquetParser(pq_path)

    @staticmethod
    def _copy_sdf_decompressed(src_path: str, dest_dir: str) -> str:
        r"""Copy an SDF file into dest_dir, decompressing .gz files so that
        byte offsets stored in parquet seeks on the copy.
        Args:
            src_path: Path to the source .sdf or .sdf.gz file.
            dest_dir: Destination directory for the copy.
        Returns:
            The basename of the written file (always .sdf, never .gz).
        """
        os.makedirs(dest_dir, exist_ok=True)
        sdf_name = os.path.basename(src_path)
        if sdf_name.endswith(".gz"):
            sdf_name = sdf_name[:-3]
        dest = os.path.join(dest_dir, sdf_name)
        if not os.path.exists(dest):
            if src_path.endswith(".sdf.gz"):
                with gzip.open(src_path, "rb") as fin, open(dest, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
            else:
                shutil.copy2(src_path, dest)
        return sdf_name

    def build_from_docking_tree(
        self,
        docking_root: str,
        num_modes: int = 9,
        copy_blocks: bool = False,
        copy_proteins: bool = False,
        protein_glob: str = "*_cleaned.pdb",
        overwrite: bool = False,
    ) -> dict[str, int]:
        r"""Scan a raw docking tree and build the per protein dataset.

        Args:
            docking_root :Root containing the sdf files
            num_modes : Poses per conformer (GNINA default 9).
            copy_blocks :If True, copy  docked SDF files into
                blocks/<PDB_ID>/.  Gzipped files are decompressed
            copy_proteins : If True, copy cleaned PDB files into proteins/.
            protein_glob : Glob pattern for the cleaned receptor PDB inside each
                protein directory.
            overwrite : If False (default), skip proteins that already have a parquet file.
        Returns:
            A dictionary mapping {protein_id: n_poses} for every protein processed.
        """
        os.makedirs(self._parquet_dir, exist_ok=True)
        if copy_proteins:
            os.makedirs(self._proteins_dir, exist_ok=True)

        protein_dirs = sorted(
            os.path.join(docking_root, d)
            for d in os.listdir(docking_root)
            if os.path.isdir(os.path.join(docking_root, d)) and not d.startswith(".")
        )

        results: dict[str, int] = {}

        for prot_dir in protein_dirs:
            protein_id = os.path.basename(prot_dir)

            sdf_files = sorted(
                glob.glob(os.path.join(prot_dir, "gypsum_out_*", "*_docked.sdf*"))
            )
            if not sdf_files:
                logger.debug(f"Skipping {protein_id}: no docked SDFs found")
                continue
            pq_path = os.path.join(self._parquet_dir, f"{protein_id}.parquet")
            if not overwrite and os.path.isfile(pq_path):
                parser = self._parser_for(protein_id)
                meta = parser.read_metadata()
                results[protein_id] = meta.num_rows
                logger.debug(f"{protein_id}: already exists ({meta.num_rows} poses)")
                continue

            tables: list[pa.Table] = []
            for sdf_path in sdf_files:
                try:
                    dr = GNINA_Results(
                        sdf_path,
                        num_modes=num_modes,
                        protein_id=protein_id,
                    )
                except Exception as e:
                    logger.error(f"Failed to parse {sdf_path}: {e}")
                    continue

                if copy_blocks:
                    dest_dir = os.path.join(self._blocks_dir, protein_id)
                    sdf_name = self._copy_sdf_decompressed(sdf_path, dest_dir)
                    # Relative path within the dataset
                    rel = f"blocks/{protein_id}/{sdf_name}"
                else:
                    # Store the full absolute path to the original file
                    rel = os.path.abspath(sdf_path)

                n = dr.n_blocks
                src = pa.array([rel] * n, type=pa.string())
                col_idx = dr.table.schema.get_field_index("source_file")
                t = dr.table.set_column(col_idx, "source_file", src)
                tables.append(t)

            if not tables:
                continue

            combined = pa.concat_tables(tables)

            parser = self._parser_for(protein_id)
            parser.write(combined)
            results[protein_id] = combined.num_rows
            logger.info(
                f"{protein_id}: {combined.num_rows} poses → {protein_id}.parquet"
            )

            if copy_proteins:
                for pdb in glob.glob(os.path.join(prot_dir, protein_glob)):
                    dest = os.path.join(self._proteins_dir, f"{protein_id}.pdb")
                    if not os.path.exists(dest):
                        shutil.copy2(pdb, dest)
                    break

        logger.info(
            f"Built per protein dataset: {len(results)} proteins, "
            f"{sum(results.values()):,} total poses in {self._root}"
        )
        return results

    def _resolve_source(self, source_file: str) -> str:
        r"""Resolve a source_file value from parquet to an actual file path.

        If source_file is a relative path (starts with 'blocks/'), it's
        resolved relative to the dataset root.  Otherwise it's treated
        as an absolute path to the original file.
        """
        if source_file.startswith("blocks/"):
            return os.path.join(self._root, source_file)
        return source_file

    def add_protein(
        self,
        protein_id: str,
        sdf_files: list[str],
        num_modes: int = 9,
        copy_blocks: bool = False,
        overwrite: bool = False,
    ) -> int:
        r"""Add a single protein's docking results from given SDF files.
        Args:
            protein_id: The identifier for the protein (e.g., "1abc").
            sdf_files: List of paths to the docked SDF files for this protein.
            num_modes: Number of poses per conformer (GNINA default 9).
            copy_blocks: If True, copy SDF files into blocks/<protein_id>/, decompressing if needed.
            overwrite: If False (default), skip if a parquet for this protein already exists.
        Returns:
            The number of poses added for this protein."""
        os.makedirs(self._parquet_dir, exist_ok=True)
        tables: list[pa.Table] = []
        pq_path = os.path.join(self._parquet_dir, f"{protein_id}.parquet")
        if not overwrite and os.path.isfile(pq_path):
            parser = self._parser_for(protein_id)
            meta = parser.read_metadata()
            logger.info(f"{protein_id}: already exists ({meta.num_rows} poses)")
            return meta.num_rows

        for sdf_path in sdf_files:
            dr = GNINA_Results(
                sdf_path,
                num_modes=num_modes,
                protein_id=protein_id,
            )

            if copy_blocks:
                dest_dir = os.path.join(self._blocks_dir, protein_id)
                sdf_name = self._copy_sdf_decompressed(sdf_path, dest_dir)
                rel = f"blocks/{protein_id}/{sdf_name}"
            else:
                rel = os.path.abspath(sdf_path)

            # deleted the duplicate line that was here

            n = dr.n_blocks
            src = pa.array([rel] * n, type=pa.string())
            col_idx = dr.table.schema.get_field_index("source_file")
            t = dr.table.set_column(col_idx, "source_file", src)
            tables.append(t)

        combined = pa.concat_tables(tables)
        parser = self._parser_for(protein_id)
        parser.write(combined)
        logger.info(f"Added {protein_id}: {combined.num_rows} poses")
        return combined.num_rows

    def _as_dataset(self) -> pds.Dataset:
        r"""Open all per-protein parquets as a unified lazy dataset."""
        return ParquetParser.open_dataset(self._parquet_dir)

    def read_protein(
        self,
        protein_id: str,
        columns: list[str] | None = None,
    ) -> pa.Table:
        r"""Read one protein's parquet."""
        parser = self._parser_for(protein_id)
        if not os.path.isfile(parser.file_path):
            raise FileNotFoundError(
                f"No parquet for protein {protein_id} at {parser.file_path}"
            )
        return parser.read(columns=columns)

    def read_block(
        self,
        protein_id: str,
        row_idx: int,
    ) -> str:
        r"""Read a single SDF block from the dataset by protein and row index.

        Works whether blocks were copied into the dataset or not.
        For copied blocks (decompressed .sdf), this is a fast seek.
        For original .sdf.gz files, this requires full decompression.

        Args:
            protein_id: The protein ID.
            row_idx: Row index within that protein's parquet.

        Returns:
            The raw SDF block text.
        """
        t = self.read_protein(
            protein_id, columns=["source_file", "block_start", "block_end"]
        )
        source = t.column("source_file")[row_idx].as_py()
        bs = t.column("block_start")[row_idx].as_py()
        be = t.column("block_end")[row_idx].as_py()
        path = self._resolve_source(source)
        return GNINA_Results.read_block_from_file(path, bs, be)

    def read_proteins(
        self,
        protein_ids: list[str],
        columns: list[str] | None = None,
    ) -> pa.Table:
        r"""Read a subset of proteins into one table.
        Args:
            protein_ids: List of protein IDs to read.
            columns: Optional list of columns to read. None = all columns.
        Returns:
            A PyArrow Table concatenating the specified proteins. If a protein ID is not found, it is skipped with a warning.
        """
        tables = []
        for pid in protein_ids:
            try:
                tables.append(self.read_protein(pid, columns=columns))
            except FileNotFoundError:
                logger.warning(f"Protein {pid} not found, skipping")
        if not tables:
            return pa.table(
                {f.name: [] for f in DOCKING_SCHEMA},
                schema=DOCKING_SCHEMA,
            )
        return pa.concat_tables(tables)

    def read_all(
        self,
        columns: list[str] | None = None,
        **filters,
    ) -> pa.Table:
        r"""Cross-protein query via pyarrow.dataset (lazy scan).

        Args:
        columns : list[str], optional
            Column subset.
        **filters
            Equality filters pushed down to row groups,
            e.g. ligand_id="M3A".
        """
        ds = self._as_dataset()
        expr = None
        for k, v in filters.items():
            e = pds.field(k) == v
            expr = e if expr is None else (expr & e)
        return ds.to_table(filter=expr, columns=columns)

    def iter_batches(
        self,
        protein_id: str | None = None,
        batch_size: int = 200,
        columns: list[str] | None = None,
    ) -> Iterator[pa.RecordBatch]:
        r"""Stream RecordBatches for loops.

        If *protein_id* is given, streams from that single parquet.
        Otherwise streams across the whole dataset.
        Args:
            protein_id: Optional single protein ID to stream from.
            batch_size: Number of rows per batch when streaming across the whole dataset.
            columns: Optional list of columns to read. None = all columns.
        Returns:
            An iterator of PyArrow RecordBatch objects.
        """
        if protein_id is not None:
            parser = self._parser_for(protein_id)
            scanner = parser.read(columns=columns, lazy=True)
        else:
            ds = self._as_dataset()
            scanner = ds.scanner(columns=columns, batch_size=batch_size)
        yield from scanner.to_batches()

    def protein_ids(self) -> list[str]:
        r"""List all protein IDs that have parquet files."""
        return sorted(
            os.path.splitext(f)[0]
            for f in os.listdir(self._parquet_dir)
            if f.endswith(".parquet")
        )

    def stats(self) -> dict[str, int]:
        r"""Quick counts without reading data (parquet metadata only)."""
        total_rows = 0
        n_proteins = 0
        for f in os.listdir(self._parquet_dir):
            if not f.endswith(".parquet"):
                continue
            protein_id = os.path.splitext(f)[0]
            parser = self._parser_for(protein_id)
            meta = parser.read_metadata()
            total_rows += meta.num_rows
            n_proteins += 1
        return {"n_proteins": n_proteins, "total_poses": total_rows}

    def to_batched_parquets(
        self,
        output_dir: str,
        proteins_per_batch: int = 200,
        row_group_size: int = 50_000,
        overwrite: bool = False,
    ) -> int:
        r"""
        Groups proteins into batch files, each containing *proteins_per_batch* proteins sorted by protein_id.
        Batch files are named by the first and last PDB ID within e.g. bt_0_10GS_4HVP.parquet.
        Args:
            output_dir :Destination directory for batch_000.parquet, etc.
            proteins_per_batch : How many proteins per batch file.
            row_group_size : Target rows per row group within each batch file.
            overwrite : If False (default), skip writing a batch file if it already exists.
        Returns:
            Number of batch files written.
        """
        os.makedirs(output_dir, exist_ok=True)

        pq_files = sorted(glob.glob(os.path.join(self._parquet_dir, "*.parquet")))
        if not pq_files:
            raise FileNotFoundError(f"No parquet files found in {self._parquet_dir}")

        n_batches = 0
        total_rows = 0

        for batch_start in range(0, len(pq_files), proteins_per_batch):
            chunk = pq_files[batch_start : batch_start + proteins_per_batch]

            first_pid = os.path.splitext(os.path.basename(chunk[0]))[0]
            last_pid = os.path.splitext(os.path.basename(chunk[-1]))[0]
            batch_name = f"bt_{n_batches}_{first_pid}_{last_pid}.parquet"
            batch_path = os.path.join(output_dir, batch_name)

            if not overwrite and os.path.isfile(batch_path):
                logger.debug(f"Skipping batch, file exists: {batch_name}")
                n_batches += 1
                continue

            tables = []
            for f in chunk:
                pid = os.path.splitext(os.path.basename(f))[0]
                tables.append(self._parser_for(pid).read())
            combined = pa.concat_tables(tables)

            batch_parser = ParquetParser(batch_path)
            batch_parser.write(combined, group_size=row_group_size)

            total_rows += combined.num_rows
            n_batches += 1
            logger.debug(
                f"  Batch {n_batches}: {len(chunk)} proteins, has {combined.num_rows:,} rows {batch_name}"
            )

        logger.info(
            f"Wrote {total_rows:,} rows across {n_batches} batch files in {output_dir}"
        )
        return n_batches

    @staticmethod
    def read_batched_parquets(
        directory: str,
        columns: list[str] | None = None,
        **filters,
    ) -> pa.Table:
        r"""Read from an batch directory with row-group pruning.
        Args:
            directory: Path to the batch directory containing batch_*.parquet files.
            columns: Optional list of columns to read. None = all columns.
            **filters: Equality filters for pruning, e.g. ligand_id="M3A".
        Returns:
            A PyArrow Table containing the filtered results from all batch files.
        """
        ds = ParquetParser.open_dataset(directory)
        expr = None
        for k, v in filters.items():
            e = pds.field(k) == v
            expr = e if expr is None else (expr & e)
        return ds.to_table(filter=expr, columns=columns)

    @staticmethod
    def iter_batches_from(
        path: str,
        batch_size: int = 200,
        columns: list[str] | None = None,
        **filters,
    ) -> Iterator[pa.RecordBatch]:
        r"""Stream batches from any parquet file or directory.
        Args:
            path: Path to a single parquet file or a directory containing batch_*.parquet files.
            batch_size: Number of rows per batch when streaming.
            columns: Optional list of columns to read. None = all columns.
            **filters: Equality filters for pruning, e.g. ligand_id="M3A".
        Returns:
            An iterator of PyArrow RecordBatch objects matching the filters.
        """
        ds = ParquetParser.open_dataset(path)
        expr = None
        for k, v in filters.items():
            e = pds.field(k) == v
            expr = e if expr is None else (expr & e)
        scanner = ds.scanner(filter=expr, columns=columns, batch_size=batch_size)
        yield from scanner.to_batches()

    def __repr__(self) -> str:
        s = self.stats()
        return (
            f"DockingDataset({self._root!r}, "
            f"proteins={s['n_proteins']}, poses={s['total_poses']:,})"
        )

    @property
    def schema(self) -> pa.Schema:
        r"""Aggregated schema across all per-protein parquets in this dataset."""
        return self._as_dataset().schema

    @staticmethod
    def topn_per_pair(
        table: pa.Table,
        rank_by: str = "CNNscore",
        n: int = 5,
        add_rank_col: bool = True,
        group_keys: tuple[str, ...] = ("protein_id", "ligand_id"),
        rank_col: str = "pair_rank",
    ) -> pa.Table:
        r"""Return the top-N rows per protein-lifand with a rank column added if requested.
        Args:
            table: PyArrow table.
            rank_by: Score column used for ranking within each group.
            n: Keep only top n rowa based on the score.
            group_keys: Columns defining the group.
            rank_col: Name of the 0-indexed rank column added to the output.
        Returns:
            A PyArrow Table, pre-sorted by (group_keys, rank_by), with rank column added if requested.
        """
        if add_rank_col and rank_col in table.column_names:
            logger.error(f"Rank column {rank_col!r} already exists in the input table")
            raise ValueError(
                f"Column {rank_col!r} already exists in table pass a different rank_col or set add_rank_col=False."
            )
        if table.num_rows == 0:
            if add_rank_col:
                return table.append_column(rank_col, pa.array([], type=pa.int32()))
            return table

        direction = GNINA_Results.best_direction(rank_by)
        sort_keys = [(k, "ascending") for k in group_keys]
        sort_keys.append((rank_by, direction))
        sorted_t = table.sort_by(sort_keys)

        m = sorted_t.num_rows
        group_start = np.zeros(m, dtype=bool)
        group_start[0] = True
        if m > 1:
            changes = np.zeros(m - 1, dtype=bool)
            for k in group_keys:
                col = sorted_t.column(k).to_numpy(zero_copy_only=False)
                changes |= col[1:] != col[:-1]
            group_start[1:] = changes

        start_positions = np.where(group_start)[0]
        last_start = np.searchsorted(start_positions, np.arange(m), side="right") - 1
        rank_in_group = (np.arange(m) - start_positions[last_start]).astype(np.int32)

        keep = rank_in_group < n
        out = sorted_t.filter(pa.array(keep, type=pa.bool_()))
        if add_rank_col:
            out = out.append_column(
                rank_col,
                pa.array(rank_in_group[keep], type=pa.int32()),
            )
        return out
