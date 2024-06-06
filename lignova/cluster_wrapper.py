r" Implemtnation for a wrapper for MMseqs2 clustering"
from typing import TextIO, Union

import ast
import os
import time

import pandas as pd
from loguru import logger

from lignova.clustering.mmseq import mmseqs_cluster, mmseqs_parser
from lignova.structure.utils import get_rcsb_data, validate_ligands, validate_pdb


def fasta_parser(fasta: str | TextIO, delimiter: str | None = None) -> list:
    r"""Parse FASTA files to get the protein ids
    Parameters
    ----------
    fasta : Union[str, TextIO]
        Path to the FASTA file.
    delimiter : Union[str, None]
        Delimiter to split the protein id from the FASTA header. Default is None.
    Returns
    -------
    list
        List of protein ids.
    """
    # check if the fasta file exists and has .fasta extension
    if not os.path.exists(fasta):
        raise FileNotFoundError(f"FASTA file {fasta} not found.")
    if not fasta.endswith(".fasta"):
        raise ValueError(f"FASTA file {fasta} must have .fasta extension.")
    # read the fasta file
    with open(fasta, "r", encoding="utf-8") as file:
        lines = file.readlines()
    # get the protein ids by filtering the lines starting with ">" a
    # and splitting them by the delimiter
    if delimiter is not None:
        protein_ids = [
            line.split(delimiter)[0].strip(">")
            for line in lines
            if line.startswith(">")
        ]
    else:
        protein_ids = [line.strip(">") for line in lines if line.startswith(">")]
    return protein_ids


def pdb_validations(csvfilenames: str) -> dict:
    r"""Validate the PDB files for the proteins in the CSV file.
    Parameters
    ----------
    csvfilenames : str
        Path to the CSV file containing the protein ids.
    """
    start_time = time.time()
    # check if the CSV file exists
    if not os.path.exists(csvfilenames):
        raise FileNotFoundError(f"CSV file {csvfilenames} not found.")
    # read the CSV file using pandas
    protein_ids = pd.read_csv(csvfilenames)
    # get the PDB column
    pdb_ids = protein_ids["PDB"]
    # check if new_pdb_files exists and if so read it as a dictionary
    if os.path.exists("new_pdb_files"):
        new_file = pd.read_csv("new_pdb_files")
        new_file = new_file.to_dict()
    else:
        new_file = {}
    # iterate over the PDB ids and validate them while keeping track of the values of the From column of the CSV file
    for gene_id, pdb_id in zip(protein_ids["From"], pdb_ids):
        # split each pdb id by ; and skip the last one
        pdb_id = pdb_id.split(";")[:-1]
        logger.debug(f"Gene id: {gene_id} PDB id: {pdb_id}")
        # find the value of the From column
        if gene_id not in new_file:
            new_file[gene_id] = []
        else:
            continue
        for pdb in pdb_id:
            # validate the pdb id
            if len(get_rcsb_data(pdb)) == 0:
                logger.error(f"Invalid PDB id {pdb} for protein {gene_id}")
                continue
            if not validate_pdb(pdb) or not validate_ligands(pdb):
                logger.error(f"Invalid PDB id {pdb} for protein {gene_id}")
                continue
            logger.info(f"Valid PDB id {pdb} for gene_id {gene_id}")
            new_file[gene_id].append(pdb)
            # save this new file to a csv file every 1hr
            if time.time() - start_time > 60 * 60:
                df = pd.DataFrame(new_file.items(), columns=["Gene_id", "PDB"])
                df.to_csv("new_pdb_files", index=False)
                start_time = time.time()

    return new_file


def fasta_filter(
    fasta: str | TextIO,
    outfile_name: str,
    csvfilenames: str,
    delimiter: str | None = "|",
) -> TextIO:
    r"""filter the proteins in the fasta file that are not in the csv file
    Parameters
    ----------
    fasta : Union[str, TextIO]
        Path to the FASTA file.
    outfile_name : str
        Path to the new FASTA file.
    csvfilenames : str
        Path to the CSV file containing the protein ids.
    delimiter : Union[str, None]
        Delimiter to split the protein id from the FASTA header. Default is |.
    Returns
    -------
    TextIO
        New FASTA file with the proteins in the csv file.
    """
    # Check if the fasta file exists and has .fasta extension
    if not os.path.exists(fasta) or not fasta.endswith(".fasta"):
        raise FileNotFoundError(f"FASTA file {fasta} not found or not a valid file.")

    # Read the CSV file using pandas
    protein_ids_df = pd.read_csv(csvfilenames)
    # Get the PDB column and make it a list
    pdb_ids = protein_ids_df["pdb_id"].tolist()

    # Read the fasta file and filter the sequences
    new_fasta = []
    keep_sequence = False

    with open(fasta, "r", encoding="utf-8") as file:
        for line in file:
            if line.startswith(">"):
                protein_id = line.split(delimiter)[0].strip(">").split("_")[0]
                keep_sequence = protein_id in pdb_ids
                if keep_sequence:
                    new_fasta.append(line)
            elif keep_sequence:
                new_fasta.append(line)

    # Write the new fasta file
    with open(outfile_name, "w", encoding="utf-8") as file:
        file.writelines(new_fasta)


def clean_cluster_files(file_path: str, delim: list = ["[", "]"]) -> list:
    """
    This function takes the mmseqs2 cluster file and cleans it by removing clusters with no representatives
    Parameters
    ----------
    file_path : str
        The path to the file containing the PDB IDs.
    delim : str, optional
        The delimiter used in the file. The default is ":".
    Returns
    -------
    pdb_ids : list
        A list of PDB IDs.
    """
    with open(file_path, "r", encoding="utf-8") as file:
        pdb_ids = []
        for line in file:
            if line.startswith("Cluster: []"):
                continue
            if line.startswith("Cluster "):
                cluster_id = (
                    line.split(" ")[1]
                    .strip()
                    .replace(delim[0], "")
                    .replace(delim[1], "")
                    .replace("\n", "")
                )
                pdb_ids.extend(cluster_id.split(", "))
        pdb_ids = [x for x in pdb_ids if x not in (" ", "''")]
        return pdb_ids


if __name__ == "__main__":
    PUBCHEM_FASTA = "../PUBCHEM_HDF5.fasta"
    """
    #Parse the fasta files to get the protein ids
    PubChem_protein_ids = fasta_parser(PUBCHEM_FASTA)
    logger.info("Number of proteins in PubChem: {}", len(PubChem_protein_ids))
    valid_pubchem = pdb_validations( "../gene-id_with_pdb.csv")
    df = pd.DataFrame(valid_pubchem.items(), columns=['Gene_id','PDB'])
    df.to_csv("../valid_pubchem.csv", index=False)
    """

    # cluster the PubChem fasta file
    mmseqs_cluster(
        PUBCHEM_FASTA,
        outfile_name_suffix="../clusters",
        tmp_dir="../tmp",
    )
    # parse the cluster file
    mmseqs_parser("../clusters_cluster.tsv", save=True)

    valid_pubchem = pd.read_csv("../valid_pubchem.csv")
    with open("../clusters_cluster_parsed.csv", "r", encoding="utf-8") as clust_file:
        lines = clust_file.readlines()
    i = 0
    new_lines = []
    count = len(clean_cluster_files("../clusters_cluster_parsed.csv"))
    while i < len(lines):
        line = lines[i]
        if line.startswith("Cluster"):
            cluster_id = line.split("Cluster")[1].split(":")[0].strip()
            # check if this cluster_id has | in it and if so skip it
            if "|" in cluster_id:
                logger.info(f"Cluster {cluster_id} has | in it")
                line = f"Cluster : {cluster_id}\n"
                new_lines.append(line)
                i += 1
                count -= 1
                continue
            elif (
                int(cluster_id) in valid_pubchem["Gene_id"].tolist()
                and valid_pubchem[valid_pubchem["Gene_id"] == int(cluster_id)][
                    "PDB"
                ].tolist()[0]
                != "[]"
            ):
                # new line is Cluster with the value of the PDB column in the valid_pubchem file
                new_line = f"Cluster : {ast.literal_eval(valid_pubchem[valid_pubchem['Gene_id']==int(cluster_id)]['PDB'].to_list()[0])}\n"
                new_lines.append(new_line)
                i += 1
                continue
            else:
                tmp_rep = cluster_id
                # check the members of the clusters and if the members are in the valid_pubchem file or had | in the name then they are the representatives
                i += 1
                cluster_members = []
                while i < len(lines) and not lines[i].startswith("Cluster"):
                    cluster_members.append(lines[i].strip())
                    i += 1
                if len(cluster_members) == 1:
                    continue
                else:
                    # check if any cluster members is in the valid_pubchem file or has | in the name
                    for member in cluster_members:
                        if "|" in member:
                            tmp_rep = member
                            break
                        elif (
                            int(member) in valid_pubchem["Gene_id"].tolist()
                            and valid_pubchem[valid_pubchem["Gene_id"] == int(member)][
                                "PDB"
                            ].tolist()[0]
                            != "[]"
                        ):
                            tmp_rep = member
                            break
                        else:
                            continue
                    if tmp_rep != cluster_id:
                        logger.info(f"new representative: {tmp_rep}")
                        if "|" not in tmp_rep:
                            tmp_rep = valid_pubchem[
                                valid_pubchem["Gene_id"] == int(tmp_rep)
                            ]["PDB"].tolist()[0]
                        new_line = f"Cluster : {tmp_rep}\n"
                        new_lines.append(new_line)
                        # add to new_lines the cluster members
                        for member in cluster_members:
                            tmp = member + "\n"
                            new_lines.append(tmp)

                    else:
                        continue
        else:
            i += 1
            new_lines.append(line)

    # write the new_lines to a new file
    with open("../new_clusters_cluster_parsed.csv", "w", encoding="utf-8") as f:
        f.writelines(new_lines)
    # NOTE:EXTRA SANITY CHECK
    with open("../new_clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    """
    # exclude the clusters with only | in their members
    new_lines = []
    i = 0
    member = []
    while i < len(lines):
        # add the member list to new_lines
        member = []
        if lines[i].startswith("Cluster"):
            cluster_id = lines[i].split(":")[1].strip()
            if "|" in cluster_id:
                new_lines.append(lines[i])
                i += 1
                # save line to new file
                while i < len(lines) and not lines[i].startswith("Cluster"):
                    # save the line in members list
                    member.append(lines[i])
                    i += 1
                # check if all members have | in them and if so add members to new_lines
                # and if not delete te last new_lines added
                if all("|" in mem for mem in member):
                    new_lines.extend(member)
                else:
                    new_lines.pop()
                    continue
            else:
                i += 1
                continue
        else:
            i += 1
            continue

    with open("../no_pubchem.csv", "w", encoding="utf-8") as f:
        f.writelines(new_lines)
    # count the number of clusters with no PubChem members
    with open("../no_pubchem.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    logger.info(
        "Number of clusters with no PubChem members: {}",
        len([line for line in lines if line.startswith("Cluster")]),
    )
    # compare lines list and new lines list and delete common lines between them
    with open("../new_clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    with open("../no_pubchem.csv", "r", encoding="utf-8") as f:
        new_lines = f.readlines()
    # delete common lines between lines and new_lines
    pubchem_only = [line for line in lines if line not in new_lines]
    with open("../pubchem_only_clusters.csv", "w", encoding="utf-8") as f:
        f.writelines(pubchem_only)
    logger.info(
        "Number of clusters with only PubChem members: {}",
        len([line for line in pubchem_only if line.startswith("Cluster")]),
    )
    """
    logger.info("Number of PubChem clusters: {}", count)
    # find the number of representatives in the new_clusters_cluster_parsed.csv file
    representatives = clean_cluster_files("../new_clusters_cluster_parsed.csv")
    logger.info("Number of representatives after parsing: {}", len(representatives))
    # find the number of representatives using lines starting with Cluster
    with open("../new_clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    clusters = [line for line in lines if line.startswith("Cluster")]
    logger.info("Number of clusters: {}", len(clusters))
    # find the number of representatives in the original clusters_cluster_parsed.csv file
    representatives = clean_cluster_files("../clusters_cluster_parsed.csv")
    logger.info(
        "Number of representatives in the original clusters_cluster_parsed.csv file: {}",
        len(representatives),
    )
    # find the number of clusters in the original clusters_cluster_parsed.csv file
    with open("../clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    clusters = [line for line in lines if line.startswith("Cluster")]
    logger.info(
        "Number of clusters in the original clusters_cluster_parsed.csv file: {}",
        len(clusters),
    )
    # find the number of proteins in ../PUBCHEM_HDF5.fasta file
    with open("../PUBCHEM_HDF5.fasta", "r", encoding="utf-8") as f:
        lines = f.readlines()
    proteins = [line for line in lines if line.startswith(">")]
    logger.info("Number of proteins in PUBCHEM_FASTA: {}", len(proteins))

    # Find out the number of gene ids i have left
    with open("../new_clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        if line.startswith("Cluster") or "|" in line:
            continue
        else:
            new_lines.append(line.strip())
    # find length of new_lines
    logger.info("Length of new_lines: {}", len(new_lines))
    logger.info("new_lines: {}", new_lines)
