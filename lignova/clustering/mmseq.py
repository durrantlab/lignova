r" Implementation of the MMSeq2 clustering algorithm. https://github.com/soedinglab/MMseqs2"

from typing import TextIO, Union

import os
import shutil
import subprocess
import time

import pandas as pd
from loguru import logger


def mmseqs_cluster(
    query_fasta: str,
    reference_fasta: str,
    cluster_threshold: Union[float, None] = 0.9,
    sort: bool = True,
    coverage_mode: Union[None, int] = 0,
    sensitivity: float = 7.0,
    outfile_name_suffix: Union[str, None] = "clusters",
    tmp_dir: Union[None, str] = "/tmp",
    cluster_mode: Union[None, int] = 0,
    self_match: bool = True,
) -> str:
    """
    Cluster sequences using MMSeqs2 https://mmseqs.com/latest/userguide.pdf
    Parameters
    ----------
    query_fasta : Union[str, TextIO]
        Path to the query FASTA file.
    reference_fasta : str
        Path to the reference FASTA file.
    cluster_threshold : float
        Cluster threshold. Default is 0.9.
    sort : bool
        Sort the output. Default is True.
    coverage_mode : int
        Coverage mode. Default is 0. See MMSeqs2 documentation for more details.
    sensitivity : float
        Sensitivity. Default is 7.0.
    outfile_name_suffix : str
        Suffix for the output file and path. Default is "clusters" saved in the current working directory.
    tmp_dir : str
        Temporary directory. Default is "/tmp" saved in the current working directory.
    cluster_mode : int
        Cluster mode. Default is 0. See MMSeqs2 documentation for more details.
    self_match : bool
        Include self-matches. Default is True.
    """
    # check if the query and reference fasta files exist
    if not os.path.exists(query_fasta):
        raise FileNotFoundError(f"Query fasta file {query_fasta} not found.")
    if not os.path.exists(reference_fasta):
        raise FileNotFoundError(f"Reference fasta file {reference_fasta} not found.")
    if coverage_mode not in [0, 1, 2, 3]:
        raise ValueError(f"Coverage mode {coverage_mode} is not valid.")
    if cluster_mode not in [0, 1, 2, 3]:
        raise ValueError(f"Cluster mode {cluster_mode} is not valid.")
    command = [
        "mmseqs",
        "easy-cluster",
        query_fasta,
        reference_fasta,
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
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode == 0:
        logger.info("Sequence Identity based clustering is completed")
        # from the path in outfile_name_suffix variable extract the directory path
        save_directory = os.path.dirname(outfile_name_suffix)
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
        with open("mmseqs2.log", "w", encoding="utf-8") as log_file:
            log_file.write(stderr)


def mmseqs_parser(
    tsv_filename: str, save: bool = False, save_filename: Union[None, str] = None
) -> pd.DataFrame:
    """
    Parse MMSeq2 output TSV file.
    Parameters
    ----------
    tsv_filename : str
        Path to the MMSeq2 output TSV file.
    save : bool
        Save the parsed clusters to a file. Default is False.
    save_filename : str
       File Path to save the parsed clusters. Default is None.
    Returns
    -------
    pd.DataFrame
        DataFrame containing the parsed MMSeq2 output.

    """
    if not os.path.exists(tsv_filename):
        raise FileNotFoundError(f"TSV file {tsv_filename} not found.")
    # Read the TSV file
    clusters = pd.read_csv(
        tsv_filename, sep="\t", header=None, names=["cluster", "members"]
    )
    logger.info(f"Read {len(clusters)} clusters from {tsv_filename}")
    # Group by cluster
    logger.info(f"found {len(clusters.groupby('cluster'))} unique clusters")
    if save:
        if save_filename is None:
            save_filename = tsv_filename.replace(".tsv", "_parsed.csv")
        with open(save_filename, "w", encoding="utf-8") as f:
            for name, group in clusters.groupby("cluster"):
                f.write(f"Cluster {name}:\n")
                for member in group["members"]:
                    f.write(f"{member}\n")
    return clusters.groupby("cluster")
