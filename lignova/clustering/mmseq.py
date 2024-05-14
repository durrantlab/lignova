r" Implementation of the MMSeq2 clustering algorithm. https://github.com/soedinglab/MMseqs2"

from typing import TextIO, Union

import os
import subprocess


# NOTE: WORKING PROGRESS
def mmseqs_cluster(
    query_fasta: str,
    reference_fasta: str,
    cluster_threshold: Union[float, None] = 0.9,
    sort: bool = True,
    coverage_mode: Union[None, int] = 0,
    sensitivity: float = 7.0,
    outfile_name: Union[str, None] = "clusters",
    tmp_dir: Union[None, str] = "/tmp",
    self_match: bool = True,
) -> str:
    """
    Cluster sequences using MMSeqs2.
    Parameters
    ----------
    query_fasta : Union[str, TextIO]
        Path to the query FASTA file.
    reference_fasta : str
        Path to the reference FASTA file.
    cluster_threshold : float
        Cluster threshold.
    Returns
    -------
    str
        Path to the output file.
    """
    # check if the query and reference fasta files exist
    if not os.path.exists(query_fasta):
        raise FileNotFoundError(f"Query fasta file {query_fasta} not found.")
    if not os.path.exists(reference_fasta):
        raise FileNotFoundError(f"Reference fasta file {reference_fasta} not found.")

    # Run MMSeqs2
    subprocess.run(
        [
            "mmseqs",
            "easy-cluster",
            query_fasta,
            reference_fasta,
            outfile_name,
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
            "0",
            "--cluster-reassign",
            "1",
            "--dbtype",
            "1",
            "--sort-results",
            "1" if sort else "0",
        ]
    )
