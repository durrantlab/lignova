# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"Implementation of the MMSeq2 clustering algorithm. https://github.com/soedinglab/MMseqs2"

import os
import shutil
import subprocess
import time
from typing import TextIO

import pandas as pd
from loguru import logger


def mmseqs_cluster(
    query_fasta: str | TextIO,
    reference_fasta: str | None = None,
    cluster_threshold: float | None = 0.9,
    sort: bool = True,
    coverage_mode: None | int = 0,
    sensitivity: float = 7.0,
    outfile_name_suffix: str | None = "clusters",
    tmp_dir: None | str = "/tmp",
    cluster_mode: None | int = 0,
    self_match: bool = True,
) -> None:
    """
    Cluster sequences using MMSeqs2 https://mmseqs.com/latest/userguide.pdf

    Args:
        query_fasta : Path to the query FASTA file.
        reference_fasta : Path to the reference FASTA file. Default is None.
        cluster_threshold : Cluster threshold. Default is 0.9.
        sort : Sort the output. Default is True.
        coverage_mode : Coverage mode. Default is 0. See MMSeqs2 documentation
            for more details.
        sensitivity : Sensitivity. Default is 7.0.
        outfile_name_suffix : Suffix for the output file and path.
            Default is "clusters" saved in the current working directory.
        tmp_dir : Temporary directory. Default is `"/tmp"` saved in the current
            working directory.
        cluster_mode : Cluster mode. Default is 0. See MMSeqs2 documentation for
            more details.
        self_match : Include self-matches. Default is True.
    """
    # check if the query and reference fasta files exist
    if not os.path.exists(query_fasta):
        raise FileNotFoundError(f"Query fasta file {query_fasta} not found.")
    if reference_fasta:
        if not os.path.exists(reference_fasta):
            raise FileNotFoundError(
                f"Reference fasta file {reference_fasta} not found."
            )
    if coverage_mode not in [0, 1, 2, 3]:
        raise ValueError(f"Coverage mode {coverage_mode} is not valid.")
    if cluster_mode not in [0, 1, 2, 3]:
        raise ValueError(f"Cluster mode {cluster_mode} is not valid.")
    command = [
        "mmseqs",
        "easy-cluster",
        query_fasta,
    ]
    if reference_fasta:
        command.extend([reference_fasta])
    command.extend(
        [
            outfile_name_suffix,
            tmp_dir,
            "-s",
            str(sensitivity),
            "--add-self-matches",
            "1" if self_match else "0",
            "--cov-mode",
            str(coverage_mode),
            "--min-seq-id",
            str(cluster_threshold),
            "--cluster-mode",
            str(cluster_mode),
            "--cluster-reassign",
            "1",
            "--dbtype",
            "1",
            "--sort-results",
            "1" if sort else "0",
        ]
    )
    with subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ) as process:
        stdout, stderr = process.communicate()
    save_directory = os.path.dirname(outfile_name_suffix)
    if process.returncode == 0:
        logger.info("Sequence Identity based clustering is completed")
        # from the path in outfile_name_suffix variable extract the directory path
        # save the standard output to a log file
        with open(
            os.path.join(save_directory, "mmseqs2.log"), "w", encoding="utf-8"
        ) as log_file:
            log_file.write(stdout.decode("utf-8"))
        cluster_output_file = outfile_name_suffix + "_rep_seq.fasta"
        while not os.path.exists(cluster_output_file):
            logger.info(f"Waiting for the {cluster_output_file} to be written")
            time.sleep(2)  # Wait for 1 second before checking again
        logger.info(f"Output files are written to {cluster_output_file}")
        # delete the temporary directory and its contents
        if os.path.exists(tmp_dir):
            logger.info(f"Deleting temporary directory {tmp_dir}")
            shutil.rmtree(tmp_dir)
    else:
        logger.error("Sequence Identity based clustering failed")
        # save the standard error to a log file
        with open(
            os.path.join(save_directory, "mmseqs2.log"), "w", encoding="utf-8"
        ) as log_file:
            log_file.write(stderr.decode("utf-8"))


def mmseqs_parser(
    tsv_filename: str, save: bool = False, save_filename: None | str = None
) -> pd.DataFrame:
    """
    Parse MMSeq2 output TSV file.

    Args:
        tsv_filename : Path to the MMSeq2 output TSV file.
        save : Save the parsed clusters to a file. Default is False.
        save_filename : File Path to save the parsed clusters. Default is None.

    Returns:
        DataFrame containing the parsed MMSeq2 output.

    """
    if not os.path.exists(tsv_filename):
        raise FileNotFoundError(f"TSV file {tsv_filename} not found.")
    # Read the TSV file
    clusters = pd.read_csv(
        tsv_filename, sep="\t", header=None, names=["cluster", "members"]
    )
    logger.info(f"Read {len(clusters)} clusters from {tsv_filename}")
    # Group by cluster and aggregate members into lists
    grouped_clusters = clusters.groupby("cluster")["members"].apply(list).reset_index()
    grouped_clusters.rename(columns={"cluster": "representatives"}, inplace=True)
    logger.info(f"Found {len(grouped_clusters)} unique clusters")
    if save:
        logger.info("Saving parsed clusters")
        if save_filename is None:
            save_filename = tsv_filename.replace(".tsv", "_parsed.csv")
            logger.info(f"Saving parsed clusters to {save_filename}")
        grouped_clusters.to_csv(save_filename, index=False)
    return grouped_clusters
