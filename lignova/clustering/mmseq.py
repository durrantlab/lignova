# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Sequence-identity clustering with MMseqs2 easy-cluster."""

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from typing import ClassVar

from loguru import logger

from .base import (
    Clusterer,
    ClusterMethod,
    ClusterParams,
    ClusterResult,
)


@dataclass(frozen=True, slots=True)
class MMseqsParams(ClusterParams):
    """Parameters for MMseqs2 sequence-identity clustering."""

    method: ClassVar[ClusterMethod] = ClusterMethod.MMSEQS2
    """The clustering method used, which is MMseqs2."""

    min_seq_id: float = 0.9
    """Minimum sequence identity for two sequences to cluster together (`--min-seq-id`).
    Must be between 0 and 1; higher values are stricter."""

    sensitivity: float = 7.0
    """Search sensitivity (`-s`). Higher values find more remote homologs at
    higher cost. MMseqs2's scale runs 1.0 (faster), 4.0 (fast), to 7.5
    (sensitive); accepted range is 1.0–7.5."""

    coverage_mode: int = 0
    """Coverage mode (`--cov-mode`), one of {0, 1, 2, 3}."""

    cluster_mode: int = 0
    """Cluster mode (`--cluster-mode`), one of {0, 1, 2, 3}."""

    self_match: bool = True
    """Whether to add self-matches (`--add-self-matches`)."""

    cluster_reassign: bool = True
    """Whether to reassign members to their best cluster after cascaded clustering
    (`--cluster-reassign`). Kept as a param because it changes cluster assignment;
    if you never toggle it, it could instead be an internal constant."""

    sort_results: bool = True
    """Whether to sort alignment results by sequence identity (`--sort-results`).
    It sets the TSV row order that determines integer cluster-id assignment. Defaults is True."""

    output_dir: str | None = None
    """Directory to write MMseqs2 outputs to.
    Default is None, which uses a temporary directory that is deleted after clustering.
    If a path is given, outputs are written there and kept after the run."""

    def __post_init__(self) -> None:
        if not (0.0 <= self.min_seq_id <= 1.0):
            raise ValueError(
                f"min_seq_id must be a float between 0 and 1, got {self.min_seq_id}"
            )
        if not (1.0 <= self.sensitivity <= 7.5):
            raise ValueError(
                f"sensitivity must be between 1.0 and 7.5, got {self.sensitivity}"
            )
        if self.coverage_mode not in (0, 1, 2, 3):
            raise ValueError(
                f"coverage_mode must be one of {{0, 1, 2, 3}}, got {self.coverage_mode}"
            )
        if self.cluster_mode not in (0, 1, 2, 3):
            raise ValueError(
                f"cluster_mode must be one of {{0, 1, 2, 3}}, got {self.cluster_mode}"
            )
        if self.output_dir is not None and not isinstance(self.output_dir, str):
            raise TypeError(
                f"output_dir must be a str path or None, got {type(self.output_dir).__name__}"
            )


class MMseqsClustering(Clusterer[MMseqsParams]):
    """MMseqs2 `easy-cluster` sequence-identity clustering."""

    # `easy-cluster` forwards createdb options; dbtype 1 = amino-acid sequences.
    _DBTYPE_AMINO_ACID: ClassVar[str] = "1"

    def __init__(self, params: MMseqsParams) -> None:
        """Initialize the MMseqs2 clusterer.

        Args:
            params: the `MMseqsParams` object containing the configuration for the clustering algorithm and output location.
        """
        super().__init__(params)

    def cluster(self, sequences: dict[str, str]) -> ClusterResult:
        """Cluster sequences by sequence identity with MMseqs2.

        Args:
            sequences: a dictionary with the keys being stable id of the sequence (used for fasta headers) and the values being the sequence string.

        Returns:
            `ClusterResult` object whose `labels` and `representatives` are
            keyed by the keys in `sequences`.
        """
        ids = list(sequences)
        n = len(ids)

        for sid in ids:
            if not sid or any(ch.isspace() for ch in sid):
                raise ValueError(
                    f"Sequence id {sid!r} is empty or contains whitespace. MMseqs2 truncates FASTA headers at whitespace."
                )

        if n == 0:
            return ClusterResult(labels={}, representatives={}, params=self.params)
        if n == 1:
            return ClusterResult(
                labels={ids[0]: 0}, representatives={0: ids[0]}, params=self.params
            )

        root, cleanup = self._resolve_root()
        try:
            query_fasta = os.path.join(root, "input.fasta")
            self._write_fasta(sequences, query_fasta)
            labels, representatives = self._run_easy_cluster([query_fasta], root)
        finally:
            if cleanup:
                shutil.rmtree(root, ignore_errors=True)

        logger.info(
            "MMseqs2 produced {rep_len} clusters from {n} sequences "
            "at min_seq_id {min_seq_id}",
            rep_len=len(representatives),
            n=n,
            min_seq_id=self.params.min_seq_id,
        )
        return ClusterResult(
            labels=labels, representatives=representatives, params=self.params
        )

    def cluster_fasta(self, *fasta_paths: str) -> ClusterResult:
        """Cluster sequences already on disk as FASTA file(s) instead of in-memory dictionary.

        Args:
            *fasta_paths: One or more paths to FASTA files.

        Returns:
            A `ClusterResult` object keyed by the FASTA headers (as MMseqs2 parses
            them, i.e. truncated at the first whitespace).
        """
        if not fasta_paths:
            raise ValueError("cluster_fasta requires at least one FASTA path.")
        for path in fasta_paths:
            if not os.path.exists(path):
                raise FileNotFoundError(f"FASTA file {path} not found.")

        self._warn_truncation(list(fasta_paths))
        root, cleanup = self._resolve_root()
        try:
            labels, representatives = self._run_easy_cluster(list(fasta_paths), root)
        finally:
            if cleanup:
                shutil.rmtree(root, ignore_errors=True)

        logger.info(
            "MMseqs2 produced {rep_len} clusters from "
            "{fasta_len} FASTA input(s) at min_seq_id {min_seq_id}",
            rep_len=len(representatives),
            fasta_len=len(fasta_paths),
            min_seq_id=self.params.min_seq_id,
        )
        return ClusterResult(
            labels=labels, representatives=representatives, params=self.params
        )

    @staticmethod
    def _warn_truncation(fasta_paths: list[str]) -> None:
        """Warn if any FASTA header will be altered or collapse under MMseqs2's
        whitespace headers trunctions.

        Args:
            fasta_paths: List of FASTA file paths to check.

        """
        seen: dict[str, str] = {}
        for path in fasta_paths:
            with open(path, encoding="utf-8") as fh:
                for line in fh:
                    if not line.startswith(">"):
                        continue
                    full = line[1:].rstrip("\n")
                    truncated = full.split(None, 1)[0] if full else ""

                    if truncated != full:
                        logger.warning(
                            "FASTA header {full!r} in {path} will be truncated to "
                            "id {truncated!r} by MMseqs2.",
                            full=full,
                            path=path,
                            truncated=truncated,
                        )

                    if truncated in seen:
                        if seen[truncated] == full:
                            logger.warning(
                                "Duplicate FASTA header {full!r} produces id "
                                "{truncated!r} more than once; entries will "
                                "collide in the clustering result.",
                                full=full,
                                truncated=truncated,
                            )
                        else:
                            logger.warning(
                                "FASTA headers {a!r} and {b!r} both reduce to "
                                "{truncated!r} after truncation; they will collide "
                                "in the clustering result.",
                                a=seen[truncated],
                                b=full,
                                truncated=truncated,
                            )
                    else:
                        seen[truncated] = full

    def _resolve_root(self) -> tuple[str, bool]:
        """Ensure the output directory exists, or create a temporary one.
        Returns:
            A tuple of root_path and cleanup where
            - root_path is the directory to use
            - cleanup is True if the directory would be deleted after use.
        """
        if self.params.output_dir is not None:
            os.makedirs(self.params.output_dir, exist_ok=True)
            return self.params.output_dir, False
        return tempfile.mkdtemp(prefix="mmseqs_"), True

    @classmethod
    def from_cluster_tsv(cls, tsv_path: str, params: MMseqsParams) -> ClusterResult:
        """Build a `ClusterResult` object from an existing MMseqs2 `*_cluster.tsv` result.

        Args:
            tsv_path: Path to an MMseqs2 cluster TSV.
            params: The params used to produce it (recorded on the result).

        Returns:
            A `ClusterResult` object keyed by the ids in the TSV.
        """
        if not os.path.exists(tsv_path):
            raise FileNotFoundError(f"TSV file {tsv_path} not found.")
        labels, representatives = cls._parse_cluster_tsv(tsv_path)
        logger.info(
            "Parsed {rep_len} clusters from {tsv_path}",
            rep_len=len(representatives),
            tsv_path=tsv_path,
        )
        return ClusterResult(
            labels=labels, representatives=representatives, params=params
        )

    def _run_easy_cluster(
        self, input_fastas: list[str], tmp_root: str
    ) -> tuple[dict[str, int], dict[int, str]]:
        """Run `easy-cluster` inside `tmp_root` and parse its cluster table.
        Args:
            input_fastas: List of FASTA file paths to cluster.
            tmp_root: Path to a temporary directory for MMseqs2 to use.
        Returns:
            A tuple of labels and representatives)as produced by `_parse_cluster_tsv`.
        """
        result_prefix = os.path.join(tmp_root, "result")
        mmseqs_tmp = os.path.join(tmp_root, "tmp")
        command = self._build_command(input_fastas, result_prefix, mmseqs_tmp)

        logger.debug("Running MMseqs2: {}", " ".join(command))
        proc = subprocess.run(command, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            logger.error(
                "MMseqs2 sequence-identity clustering failed (exit {returncode})",
                returncode=proc.returncode,
            )
            raise RuntimeError(
                f"MMseqs2 failed (exit {proc.returncode}):\n{proc.stderr.strip()}"
            )
        logger.info("MMseqs2 sequence-identity clustering completed")

        tsv_path = f"{result_prefix}_cluster.tsv"
        if not os.path.exists(tsv_path):
            raise RuntimeError(
                f"MMseqs2 exited 0 but produced no cluster table at {tsv_path}."
            )
        return self._parse_cluster_tsv(tsv_path)

    def _build_command(
        self, input_fastas: list[str], result_prefix: str, mmseqs_tmp: str
    ) -> list[str]:
        """Build the command to run `mmseqs easy-cluster`.
        Args:
            input_fastas: List of FASTA file paths to cluster.
            result_prefix: Prefix for the output files.
            mmseqs_tmp: Path to a temporary directory for MMseqs2 to use.
        Returns:
            A list of command-line arguments to pass to `subprocess.run`.
        """
        p = self.params
        command = ["mmseqs", "easy-cluster", *input_fastas]
        command.extend(
            [
                result_prefix,
                mmseqs_tmp,
                "-s",
                str(p.sensitivity),
                "--add-self-matches",
                "1" if p.self_match else "0",
                "--cov-mode",
                str(p.coverage_mode),
                "--min-seq-id",
                str(p.min_seq_id),
                "--cluster-mode",
                str(p.cluster_mode),
                "--cluster-reassign",
                "1" if p.cluster_reassign else "0",
                "--sort-results",
                "1" if p.sort_results else "0",
                "--dbtype",
                self._DBTYPE_AMINO_ACID,
            ]
        )
        return command

    @staticmethod
    def _write_fasta(sequences: dict[str, str], path: str) -> None:
        """Write `sequences` to a FASTA file at `path`."""
        with open(path, "w", encoding="utf-8") as fh:
            for sid, seq in sequences.items():
                fh.write(f">{sid}\n{seq}\n")

    @staticmethod
    def _parse_cluster_tsv(
        tsv_path: str,
    ) -> tuple[dict[str, int], dict[int, str]]:
        """Parse an MMseqs2 two-column cluster TSV into labels and representatives.
        with representatives keyed by cluster id and labels keyed by sequence id.

        Args:
            tsv_path: Path to an MMseqs2 cluster TSV file.

        Returns:
            A tuple of (labels, representatives) where:
            - `labels` is a dictionary mapping sequence ids to their integer cluster ids.
            - `representatives` is a dictionary mapping integer cluster ids to their representative sequence ids.
        """
        labels: dict[str, int] = {}
        representatives: dict[int, str] = {}
        rep_to_cid: dict[str, int] = {}

        with open(tsv_path, encoding="utf-8") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if not line:
                    continue
                rep, member = line.split("\t")
                cid = rep_to_cid.get(rep)
                if cid is None:
                    cid = len(rep_to_cid)
                    rep_to_cid[rep] = cid
                    representatives[cid] = rep
                labels[member] = cid

        # Ensure that every representative is present in labels with its own cluster id
        for cid, rep in representatives.items():
            labels.setdefault(rep, cid)

        return labels, representatives
