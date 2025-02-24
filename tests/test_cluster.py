r"""Test protein clustering module."""

import os

from lignova.clustering import mmseqs_cluster, mmseqs_parser

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

context_filepaths = {
    "reference_filepath": "./files/reference.fasta",
    "query_filepath": "./files/query.fasta",
    "write_dir": "./tmp/clustering",
}


def prep_dirs():
    r"""Prepare directories for writing files."""
    os.makedirs(context_filepaths["write_dir"])


if not os.path.exists(context_filepaths["write_dir"]):
    prep_dirs()


def test_mmseqs_clustering():
    """Test mmseqs clustering."""
    mmseqs_cluster(
        context_filepaths["query_filepath"],
        context_filepaths["reference_filepath"],
        outfile_name_suffix=os.path.join(context_filepaths["write_dir"], "clusters"),
        tmp_dir=os.path.join(context_filepaths["write_dir"], "tmp"),
        sort=True,
        coverage_mode=0,
        sensitivity=7.0,
        cluster_mode=0,
        self_match=True,
    )
    assert os.path.exists(
        os.path.join(context_filepaths["write_dir"], "clusters_rep_seq.fasta")
    )
    with open(
        os.path.join(context_filepaths["write_dir"], "clusters_rep_seq.fasta"),
        encoding="utf-8",
    ) as f:
        lines = f.readlines()
    assert os.path.exists(
        os.path.join(context_filepaths["write_dir"], "clusters_all_seqs.fasta")
    )
    assert os.path.exists(
        os.path.join(context_filepaths["write_dir"], "clusters_cluster.tsv")
    )
    assert len([line for line in lines if line.startswith(">")]) == 2
    assert [
        line.strip(">").split("|")[0] for line in lines if line.startswith(">")
    ] == ["2PAM_1", "4Q0Q_1"]


def test_mmseqs_parser():
    r"""Test mmseqs parser."""
    mmseqs_cluster(
        context_filepaths["query_filepath"],
        context_filepaths["reference_filepath"],
        outfile_name_suffix=os.path.join(context_filepaths["write_dir"], "clusters"),
        tmp_dir=os.path.join(context_filepaths["write_dir"], "tmp"),
        sort=True,
        coverage_mode=0,
        sensitivity=7.0,
        cluster_mode=0,
        self_match=True,
    )
    clusters = mmseqs_parser(
        os.path.join(context_filepaths["write_dir"], "clusters_cluster.tsv"), save=True
    )

    assert len(clusters) == 2
    assert os.path.exists(
        os.path.join(context_filepaths["write_dir"], "clusters_cluster_parsed.csv")
    )
