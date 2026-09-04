# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Test protein clustering module."""

import os

import pytest

from lignova.clustering import MMseqsClustering, MMseqsParams

os.chdir(os.path.dirname(os.path.realpath(__file__)))

QUERY_FASTA = "./files/query.fasta"
REFERENCE_FASTA = "./files/reference.fasta"
WRITE_DIR = "./tmp/clustering"


EXPECTED_REPRESENTATIVES = {"2PAM_1|Chains", "4Q0Q_1|Chain"}


if not os.path.exists(WRITE_DIR):
    os.makedirs(WRITE_DIR)


def _read_fasta(path: str) -> dict[str, str]:
    """Read a FASTA file into a dictionary of sequence IDs and sequences."""
    seqs: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks)
                current = line[1:].split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        seqs[current] = "".join(chunks)
    return seqs


@pytest.fixture(scope="module")
def params() -> MMseqsParams:
    """Return a set of parameters for MMseqs2 clustering that will produce two clusters with the expected representatives, writing outputs to WRITE_DIR so the TSV persists."""
    return MMseqsParams(
        min_seq_id=0.9,
        sensitivity=7.0,
        coverage_mode=0,
        cluster_mode=0,
        self_match=True,
        sort_results=True,
        output_dir=WRITE_DIR,
    )


def test_cluster_fasta(params: MMseqsParams) -> None:
    """Test that the FASTA file path path produces the expected partition."""
    result = MMseqsClustering(params).cluster_fasta(QUERY_FASTA, REFERENCE_FASTA)

    assert result.n_clusters == 2
    reps = {rep for rep in result.representatives.values()}
    assert reps == EXPECTED_REPRESENTATIVES
    assert set(result.labels.values()) == set(result.representatives)


def test_cluster_dict(params: MMseqsParams) -> None:
    """Test that the in-memory dict path produces the same partition as the file path."""
    sequences = {**_read_fasta(QUERY_FASTA), **_read_fasta(REFERENCE_FASTA)}
    result = MMseqsClustering(params).cluster(sequences)

    assert result.n_clusters == 2
    reps = {rep for rep in result.representatives.values()}
    assert reps == EXPECTED_REPRESENTATIVES


def test_output_dir_persists_tsv(params: MMseqsParams) -> None:
    """Test that setting output_dir leaves MMseqs2's outputs on disk under WRITE_DIR."""
    MMseqsClustering(params).cluster_fasta(QUERY_FASTA, REFERENCE_FASTA)

    assert os.path.exists(os.path.join(WRITE_DIR, "result_cluster.tsv"))
    assert os.path.exists(os.path.join(WRITE_DIR, "result_rep_seq.fasta"))
    assert os.path.exists(os.path.join(WRITE_DIR, "result_all_seqs.fasta"))


def test_from_cluster_tsv(params: MMseqsParams) -> None:
    """Test that populating the ClusterResult from the tsv file remain consistent"""
    run_result = MMseqsClustering(params).cluster_fasta(QUERY_FASTA, REFERENCE_FASTA)

    tsv_path = os.path.join(WRITE_DIR, "result_cluster.tsv")
    assert os.path.exists(tsv_path)

    parsed_result = MMseqsClustering.from_cluster_tsv(tsv_path, params)

    assert parsed_result.n_clusters == run_result.n_clusters == 2
    assert parsed_result.labels == run_result.labels
    assert parsed_result.representatives == run_result.representatives
    reps = {rep for rep in parsed_result.representatives.values()}
    assert reps == EXPECTED_REPRESENTATIVES
