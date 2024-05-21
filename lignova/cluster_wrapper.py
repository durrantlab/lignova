r" Implemtnation for a wrapper for MMseqs2 clustering"
from typing import TextIO, Union

import ast
import os
import time

import pandas as pd
from loguru import logger

from lignova.clustering.mmseq import mmseqs_cluster, mmseqs_parser
from lignova.structure.utils import get_rcsb_data, validate_ligands, validate_pdb


def fasta_parser(fasta: Union[str, TextIO], delimiter: Union[str, None] = None) -> list:
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


def binding_moad_validation(csvfilenames: str) -> list:
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
    pdb_ids = protein_ids["pdb_id"]
    # make a list to store the valid pdb ids
    valid_pdb_ids = []
    # iterate over the PDB ids and validate them
    for pdb_id in pdb_ids:
        # validate the pdb id
        if len(get_rcsb_data(pdb_id)) == 0:
            logger.error(f"Invalid PDB id {pdb_id}")
            continue
        if not validate_pdb(pdb_id) or not validate_ligands(pdb_id):
            logger.error(f"Invalid PDB id {pdb_id}")
            continue
        logger.info(f"Valid PDB id {pdb_id}")
        valid_pdb_ids.append(pdb_id)

        if time.time() - start_time > 60 * 60:
            df = pd.DataFrame(valid_pdb_ids, columns=["pdb_id"])
            df.to_csv("valid_binding_moad.csv", index=False)
            start_time = time.time()
    return valid_pdb_ids


def fasta_filter(
    fasta: Union[str, TextIO],
    outfile_name: str,
    csvfilenames: str,
    delimiter: Union[str, None] = "|",
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
    with open(outfile_name, "w", encoding="utf-8") as f:
        f.writelines(new_fasta)


def clean_cluster_files(
    file_path: str, delim: list = ["[", "]"]
):  # NOTE:change this delim when dealing with filtered to ( and )
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
    BIND_MOAD_FASTA = "/home/mma121/PubChem_small/try_schrodinger/MOAD_sequences.fasta"
    PUBCHEM_FASTA = "/home/mma121/PubChem_small/try_schrodinger/PUBCHEM_HDF5.fasta"
    """

    #Parse the fasta files to get the protein ids
    bind_moad_protein_ids = fasta_parser(BIND_MOAD_FASTA, delimiter="|")
    PubChem_protein_ids = fasta_parser(PUBCHEM_FASTA)
    logger.info("Number of proteins in bind MOAD: {}", len(bind_moad_protein_ids))
    logger.info("Number of proteins in PubChem: {}", len(PubChem_protein_ids))
    # save the bind MOAD protein ids to a file as csv file
    #split the protein ids by _ and take the first element
    bind_moad_protein_ids = [protein_id.split("_")[0] for protein_id in bind_moad_protein_ids]
    df = pd.DataFrame(bind_moad_protein_ids, columns=["pdb_id"])
    df.to_csv("bind_moad_protein_ids.csv", index=False)
    

    # validate the PDB files for both the bind MOAD and PubChem proteins
    new_file = binding_moad_validation(
        "/home/mma121/PubChem_small/try_schrodinger/bind_moad_protein_ids.csv"
    )
    # save the new file which is a list to a csv file
    df = pd.DataFrame(new_file, columns=["pdb_id"])
    df.to_csv("../valid_binding_moad.csv", index=False)
    valid_pubchem = pdb_validations( "/home/mma121/PubChem_small/try_schrodinger/gene-id_with_pdb.csv")
    df = pd.DataFrame(valid_pubchem.items(), columns=['Gene_id','PDB'])
    df.to_csv("../valid_pubchem.csv", index=False)


    # Delete the invalid Binding MOAD proteins from the fasta file and save the new fasta file
    fasta_filter(
        BIND_MOAD_FASTA, "../valid_binding_moad.csv", outfile_name="../valid_MOAD.fasta"
    )
    NEW_BINDING_MOAD_FASTA = "../valid_MOAD.fasta"
   

    # cluster the valid Binding MOAD fasta file and the PubChem fasta file
    mmseqs_cluster(
        PUBCHEM_FASTA,
        NEW_BINDING_MOAD_FASTA,
        outfile_name_suffix="../clusters",
        tmp_dir="../tmp",
    )
    # parse the cluster file
    mmseqs_parser("../clusters_rep_seq.fasta.tsv", save=True)

    """
    valid_pubchem = pd.read_csv("../valid_pubchem.csv")
    # NOTE: this is a working code to assess the loss of data in the clustering process NOT DONE YET
    # read the ../clusters_cluster_parsed.csv file
    with open("../clusters_cluster_parsed.csv", "r", encoding="utf-8") as clust_file:
        lines = clust_file.readlines()
    i = 0
    new_lines = []
    count = 7014
    pubcount = 2247
    while i < len(lines):
        line = lines[i]
        if line.startswith("Cluster"):
            cluster_id = line.split("Cluster")[1].split(":")[0].strip()
            logger.debug(cluster_id)
            logger.info(f"Cluster {cluster_id}")
            # check if this cluster_id has | in it and if so skip it
            if "|" in cluster_id:
                logger.info(f"Cluster {cluster_id} has | in it")
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
                logger.info(f"Cluster {cluster_id} has a valid PDB ID")
                logger.debug(
                    valid_pubchem[valid_pubchem["Gene_id"] == int(cluster_id)][
                        "PDB"
                    ].tolist()
                )
                # new line is Cluster with the value of the PDB column in the valid_pubchem file
                new_line = f"Cluster {ast.literal_eval(valid_pubchem[valid_pubchem['Gene_id']==int(cluster_id)]['PDB'].to_list()[0])}\n"
                new_lines.append(new_line)
                i += 1
                pubcount -= 1
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
                        if member in valid_pubchem["Gene_id"].tolist() or "|" in member:
                            tmp_rep = member
                            break
                    if tmp_rep != cluster_id:
                        logger.info(f"new representative: {tmp_rep}")
                        if "|" not in tmp_rep:
                            tmp_rep = valid_pubchem[
                                valid_pubchem["Gene_id"] == tmp_rep
                            ]["Gene_id"].tolist()[0]
                        new_line = f"Cluster {tmp_rep}:\n"
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
    with open("new_clusters_cluster_parsed.csv", "w", encoding="utf-8") as f:
        f.writelines(new_lines)

    logger.info("Number of PubChem clusters: {}", count)
    logger.info("Number of invalid PubChem clusters: {}", pubcount)
    # find the number of representatives in the new_clusters_cluster_parsed.csv file
    representatives = clean_cluster_files("new_clusters_cluster_parsed.csv")
    logger.info("Number of representatives: {}", len(representatives))
    # find the number of representatives using lines starting with Cluster
    with open("new_clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
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
    """
    with open("../clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    clusters = [line for line in lines if line.startswith("Cluster")]
    logger.info("Number of clusters: {}", len(clusters))
    #clean the cluster files using the clean_cluster_files function
    representatives=clean_cluster_files("../clusters_cluster_parsed.csv")
    logger.info("Number of representatives: {}", len(representatives))
    #save the representatives to a file
    with open("../representative_prot.txt", "w", encoding="utf-8") as f:
        for rep in representatives:
            f.write(rep + "\n")
    """
    """
    #read this file '../representative_prot.txt' using pandas and find the number of representatives without | in the name
    representatives=pd.read_csv("../representative_prot.txt", header=None)
    representatives=representatives[0].tolist()
    #find the representatives without | in the name
    representatives_pubchem=[int(rep.rstrip(":")) for rep in representatives if "|" not in rep]
    logger.info("Number of representatives without | in the name: {}", len(representatives_pubchem))
    #find theses representatives in the ../valid_pubchem.csv file
    valid_pubchem=pd.read_csv("../valid_pubchem.csv")
    #find the representatives in the valid_pubchem file
    valid_pubchem=valid_pubchem[valid_pubchem["Gene_id"].isin(representatives_pubchem)]
    logger.info("Number of representatives in the valid_pubchem file: {}", len(valid_pubchem))
    #what are the not found representatives
    not_found=[rep for rep in representatives_pubchem if rep not in valid_pubchem["Gene_id"].tolist()]
    logger.info("Number of representatives not found in the valid_pubchem file: {}", len(not_found))
    #loop through the found representatives and check if the PDB column is []
    count=0
    need_to_check_gp=[]
    for rep in valid_pubchem["Gene_id"].tolist():
        if valid_pubchem[valid_pubchem["Gene_id"]==rep]["PDB"].tolist()[0]=="[]":
            count=count+1
            need_to_check_gp.append(rep)
    logger.info("Number of representatives with no PDB IDs: {}", count)
    #save the need_to_check_gp to a file
    with open("need_to_check_gp.txt", "w", encoding="utf-8") as f:
        for rep in need_to_check_gp:
            f.write(str(rep) + "\n")
    # read the need_to_check_gp file and find the clusters with members more than 1
    with open("need_to_check_gp.txt", "r", encoding="utf-8") as f:
        need_to_check_gp = f.readlines()
    need_to_check_gp = [rep.strip() for rep in need_to_check_gp]
    # find the need_to_check_gp in the clusters_cluster_parsed.csv file and
    # check if the cluster has members more than 1
    # by reading line by line and checking if the line starts with Cluster and the
    # number of lines till the next Cluster line is more than 1
    with open("../clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    clusters = {}
    cluster_members = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("Cluster"):
            tmp = line.split("Cluster")[1].split(":")[0].strip()
            logger.info(tmp)
            if tmp in need_to_check_gp:
                i += 1
                while i < len(lines) and not lines[i].startswith("Cluster"):
                    cluster_members.append(lines[i].strip())
                    i += 1
                if len(cluster_members) > 1:
                    clusters[tmp] = cluster_members
                    logger.info(f"Cluster members: {len(cluster_members)}")
                else:
                    cluster_members = []
                clusters[tmp] = cluster_members
                logger.info(f"Cluster members: {len(cluster_members)}")
                cluster_members = []
            else:
                i += 1
        else:
            i += 1
    #save the clusters dictionary to a txt file using pandas
    df = pd.DataFrame(clusters.items(), columns=["Gene_id", "PDB"])
    df.to_csv("clusters_members.csv", index=False)
    # read the ../valid_pubchem.csv file
    valid_pubchem = pd.read_csv("../valid_pubchem.csv")
    # get a list of the Gene_id column
    gene_ids = valid_pubchem["Gene_id"].tolist()
    # find the number of clusters where the value is []
    count = 0
    for key, value in clusters.items():
        if value == []:
            count += 1
    logger.info("Number of clusters with no members: {}", count)
    # loop over the clusters_cluster_parsed.csv lines and if the clusters value is [] delete the lines till the next Cluster line
    with open("../clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    new_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        new_lines.append(line)
        if line.startswith("Cluster"):
            tmp = line.split("Cluster")[1].split(":")[0].strip()
            if tmp in clusters and clusters[tmp] == []:
                i += 1
                while i < len(lines) and not lines[i].startswith("Cluster"):
                    i += 1
                # delete the last item in the new_lines list
                new_lines.pop()
            elif tmp not in gene_ids:
                new_lines.pop()
                i += 1
            else:
                new_lines.append(line)
                i += 1
        else:
            new_lines.append(line)
            i += 1
    # write the new lines to a new file
    with open("new_clusters_cluster_parsed.csv", "w", encoding="utf-8") as f:
        f.writelines(new_lines)
"""
