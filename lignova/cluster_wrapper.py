r" Implemtnation for a wrapper for MMseqs2 clustering"
from typing import TextIO, Union

import ast
import os
import pickle
import time

import matplotlib
import matplotlib.pyplot as plt
import matplotlib_inline
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from loguru import logger

from lignova.analysis.utils import obabel_convert
from lignova.clustering.mmseq import mmseqs_cluster, mmseqs_parser
from lignova.clustering.tanimoto import TanimotoClustering
from lignova.hdf5.parquet import ParquetParser
from lignova.hdf5.parser import HDF5Parser
from lignova.structure import Protein
from lignova.structure.utils import (
    convert_cif2pdb,
    get_ligand_names,
    get_rcsb_data,
    get_smiles,
    validate_ligands,
    validate_pdb,
)


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
                pdb_ids.append(cluster_id.split(", "))
        pdb_ids = [x for x in pdb_ids if x not in (" ", "''")]
        return pdb_ids


def get_protein_clusters(cluster_file: str) -> dict[tuple | list]:
    r"""Get the clusters from the cluster file.
    Parameters
    ----------
    cluster_file : str
        Path to the cluster file.
    Returns
    -------
    dict
        Dictionary of clusters. The keys are the cluster representative and the values are the cluster members.
    """

    # loop through the cluster file and write the data to the new file
    with open(cluster_file, "r", encoding="utf-8") as file:
        lines = file.readlines()
    # find lines starting with Cluster and get the number of lines
    clusters = [line for line in lines if line.startswith("Cluster")]
    cluster_number = 0
    members = []
    representatives = []
    tmp = []
    cluster_dict = {}
    for line in lines:
        if line.startswith("Cluster"):
            if tmp != []:
                members.append(tmp)
                cluster_dict[representatives[cluster_number - 1]] = members
                tmp = []
                members = []
            # split the line by : and append the value to the representatives
            cluster_number += 1
            representatives.append(line.split(":")[1].strip())
            continue
        else:
            tmp.append(line.strip())
        # check if this is the last line and if so append the members to the members list in the cluster_dict
        if line == lines[-1]:
            members.append(tmp)
            cluster_dict[representatives[cluster_number - 1]] = members
    return cluster_dict


# NOTE:THIS FUNCTION IS ONLY IMPLEMENTED FOR PARQUET FILES
def make_protein_cluster_file(new_file: str, cluster_csv: str) -> None:
    r"""Make a new cluster file with the new lines.
    Parameters
    ----------
    new_file : str
        Path to the new cluster file. if the file not found, it will be created.
    cluster_csv : str
        Path to the cluster file.
    """
    # check if the cluster file exists
    if not os.path.exists(cluster_csv) or not cluster_csv.endswith(".csv"):
        raise FileNotFoundError(
            f"Cluster file {cluster_csv} not found or not a valid file."
        )
    # check the extension of the new file if it is h5 or parquet
    if new_file.endswith(".h5") or new_file.endswith(".hdf5"):
        if not os.path.exists(new_file):
            hdf5 = HDF5Parser(new_file)
            hdf5.create()
    elif new_file.endswith(".parquet"):
        schema = pa.schema(
            [
                ("Cluster number", pa.int64()),
                ("Represenatives", pa.string()),
                ("members", pa.string()),
                ("member_compound", pa.list_(pa.string())),
            ]
        )
        clusters = get_protein_clusters(cluster_csv)
        logger.info(f"Number of clusters: {len(clusters)}")
        data = []
        for cluster_number, (representatives, members_list) in enumerate(
            clusters.items(), start=1
        ):
            for rep in ast.literal_eval(representatives):
                for member in members_list[0]:
                    data.append((cluster_number, rep, member, []))
        logger.debug(data[:5])
        logger.debug(data[25])
        parquet = ParquetParser(new_file, schema)
        if not os.path.exists(new_file):
            parquet.create()
        parquet.write(data, schema)
    else:
        raise ValueError(f"Invalid file extension {new_file}")


def hdf5_raw_file_parser(hdf5_file: str) -> tuple[dict, dict]:
    r"""Parse the HDF5 file to get the protein ids and the ligand ids.
    Parameters
    ----------
    hdf5_file : str
        Path to the HDF5 file.
    Returns
    -------
    tuple
        Tuple of dictionaries containing the protein ids and the ligand ids.
    """
    # check if the HDF5 file exists
    if not os.path.exists(hdf5_file) or not hdf5_file.endswith(".hdf5"):
        raise FileNotFoundError(f"HDF5 file {hdf5_file} not found or not a valid file.")
    # read the HDF5 file
    hdf5 = HDF5Parser(hdf5_file)
    # read the hdf5 file
    aids = hdf5.read("aids")
    aid_2_target = {}
    aid_2_cids = {}
    logger.info(f"Number of aids: {len(aids)}")
    for aid in aids:
        gene_id = hdf5.read(f"/aids/{aid}/targets_gene_id")
        cids = hdf5.read(f"/aids/{aid}/cids")
        logger.debug(f"cids: {cids}")
        aid_2_target[aid] = gene_id[0]
        aid_2_cids[aid] = cids
    logger.debug(f"Number of aids: {len(aid_2_target)}")
    logger.debug(f"Number of aids: {len(aid_2_cids)}")
    return aid_2_target, aid_2_cids


def add_compounds(
    original_file: ParquetParser,
    data_source: tuple(dict[str | str]),
    overwrite: bool = True,
) -> None:
    r"""Add the compounds to the parquet file.
    Parameters
    ----------
    original_file : ParquetParser
        ParquetParser object of the original file.
    data_source : tuple
        Tuple containing the dictionaries of the aids and cids.
    overwrite : bool
        Whether to overwrite the original file. Default is True.
        if false, the new file will be created and names as original_file + "_new".
    """
    # read the data_source as a dataframe
    aid_2_target, aid_2_cids = data_source
    # convert the dictionanries to a pandas dataframe
    aid_2_target = pd.DataFrame(aid_2_target.items(), columns=["AID", "Gene_id"])
    aid_2_cids = pd.DataFrame(aid_2_cids.items(), columns=["AID", "CIDs"])
    old_data = original_file.convert_to_pandas()
    # loop through the old_data and update the member_compound column
    for index, row in old_data.iterrows():
        gene_id = row["members"]
        member = row["member_compound"]
        # get the aids for each gene_id
        aids = list(
            aid_2_target[aid_2_target["Gene_id"].astype(str) == str(gene_id)]["AID"]
        )
        # loop through the aids and get the cids for each aid
        cids = []
        for aid in aids:
            # add the cids to the cids list knowing that the cids are ['5291'] sting of list
            cids.extend(aid_2_cids[aid_2_cids["AID"] == aid]["CIDs"].to_list()[0])
            logger.debug(f"cids: {cids}")
        # update the member_compound column
        old_data.at[index, "member_compound"] = cids
        logger.info(old_data.at[index, "member_compound"])
    if overwrite:
        original_file.write(old_data)
        logger.info(f"Data written to Parquet file at {original_file.file_path}")
    else:
        new_file = "new_" + original_file.file_path
        new_parser = ParquetParser(new_file, original_file.schema)
        new_parser.write(old_data, original_file.schema)
        logger.info(f"New file created at {new_file}")


def make_ligand_cluster_file(
    old_parquet: ParquetParser,
    new_file_name: ParquetParser,
    aid_2_cids: pd.DataFrame,
    aid_2_target: pd.DataFrame,
    hdf5_file: str,
) -> None:
    r"""Make a new Parquet file for the ligand clustering.
    Parameters
    ----------
    old_parquet : ParquetParser
        ParquetParser object of the old file.
    new_file_name : ParquetParser
        ParquetParser object of the new file.
    aid_2_cids : pd.DataFrame
        DataFrame containing the aids and cids.
    aid_2_target : pd.DataFrame
        DataFrame containing the aids and gene ids.
    hdf5_file : str
        Path to the HDF5 file.
    """
    initial_data = []
    hdf5 = HDF5Parser("../PubChem_data_edited_400k.hdf5")
    old_data = old_parquet.convert_to_pandas().groupby("Cluster number")
    # check if the progress_cache.pkl file exists and if so read it
    if os.path.exists("progress_cache.pkl"):
        with open("progress_cache.pkl", "rb") as f:
            initial_data = pickle.load(f)

    # loop through the data and get the gene ids and the cids
    for cluster_number, data in old_data:
        logger.debug(f'Cluster represenatives: {data["Represenatives"]}')
        for pdb_id in data["Represenatives"]:
            # check if the data["Represenatives"] is in the initial_data pdb_id
            if any(pdb_id == x[1] for x in initial_data):
                logger.debug(f"pdb_id: {pdb_id} already in initial_data")
                continue
            protein = Protein()
            # check if the pdb file exists
            if not os.path.exists(f"../representatives/{pdb_id.lower()}.pdb"):
                # then use the protein.load function to get the pdb file
                protein.load(
                    pdb_id=pdb_id,
                    write=True,
                    write_path=f"../representatives/{pdb_id.lower()}.pdb",
                )
            else:
                protein.load(file_path=f"../representatives/{pdb_id.lower()}.pdb")
            logger.debug(protein._pdb_file_path)
            ligands = get_ligand_names(protein._pdb_file_path)
            if len(ligands) > 1:
                for ligand in ligands:
                    smiles = get_smiles(ligand)
                    initial_data.append(
                        (cluster_number, pdb_id, ligand, smiles["stereo_smiles"], None)
                    )
            elif len(ligands) == 1:
                smiles = get_smiles(ligands[0])
                initial_data.append(
                    (cluster_number, pdb_id, ligands[0], smiles["stereo_smiles"], None)
                )
            else:
                logger.error(f"No ligands found for {pdb_id}")
                continue
        for member in data["members"]:
            # check if the data["members"] is in the initial_data member
            if any(member == x[2] for x in initial_data):
                logger.debug(f"member: {member} already in initial_data")
                continue
            # find the aids for each gene id
            aids = list(
                aid_2_target[aid_2_target["Gene_id"].astype(str) == str(member)]["AID"]
            )
            # loop through the aids and get their cids from the aid_2_cids dataframe
            for aid in aids:
                logger.debug(f"aid: {aid}")
                cids = ast.literal_eval(
                    aid_2_cids[aid_2_cids["AID"] == aid]["CIDs"].to_list()[0]
                )
                for cid in cids:
                    # get the smiles for each cid using hdf5.read(f'aids/{aid}/cids/{cid}/smiles')
                    try:
                        smiles = hdf5.read(f"aids/{aid}/cids/{cid}/smiles").astype(str)
                        logger.debug(f"smiles: {smiles}")
                        initial_data.append(
                            (cluster_number, member, cid, smiles[0], None)
                        )
                    except Exception as e:
                        logger.error(f"Error: {e}")
                        initial_data.append((cluster_number, member, cid, None, None))
                        continue
            # save the progress to cache file
            cache_file = "progress_cache.pkl"
            logger.warning(f"Saving progress to {cache_file}")
            with open(cache_file, "wb") as f:
                pickle.dump(initial_data, f)
    new_parquet.write(initial_data, new_schema)


def add_smiles_cluster(
    old_parquet: ParquetParser, similarity_cutoff: float, new_parquet: ParquetParser
):
    r"""Clustering smiles and Adding the smiles cluster to the parquet file.
    Parameters
    ----------
    old_parquet : ParquetParser
        ParquetParser object of the old file.
    similarity_cutoff : float
        The cutoff value for the tanimoto clustering similarity.
    new_parquet : ParquetParser
        ParquetParser object of the new file.
    """
    # READ THE DATA FROM THE compound_clustered_pubchem.parquet file
    data = parquet.convert_to_pandas().groupby("Protein Cluster number")
    TanimotoClustering = TanimotoClustering()
    for cluster, rest in data:
        compounds_count = []
        logger.debug(f"Length of compounds with duplication: {len(rest)}")
        cluster_number = cluster
        # get the rows with unique values in the compound id column
        unique_rows = rest.drop_duplicates(
            subset="Compound ID", ignore_index=True
        ).dropna()
        logger.debug(f"unique rows: {unique_rows}")
        logger.debug(f"Length of unique rows: {len(unique_rows)}")
        smiles = unique_rows["Smiles"]
        compound_ids = unique_rows["Compound ID"]
        logger.debug(f"Length of smiles: {len(smiles)}")
        original_length = len(smiles)
        # Get fingerprints for each molecule and tanimoto similarity as one liner
        fingerprints = [
            TanimotoClustering.get_morgan_fingerprint(smile) for smile in smiles
        ]
        logger.debug(f"unique fingerprints: {len(fingerprints)}")
        similarity = []
        # loop over the fingerprints to get the similarity for each molecule with the previous molecules
        for i in range(1, len(fingerprints)):
            similarity.append(
                TanimotoClustering.tanimoto_similarity(
                    fingerprints[:i], fingerprints[i]
                )
            )
        logger.debug(f"Length of similarity: {len(similarity)}")
        # run TanimotoClustering.cluster_tanimoto(similarity, compound_ids, cutoff)
        clusters = TanimotoClustering.cluster_tanimoto(
            similarity, compound_ids, similarity_cutoff
        )
        # loop through the clusters and add the cluster number to the dataframe Ligand Cluster number
        for cluster_number, cluster in enumerate(clusters, start=1):
            for compound in cluster:
                rest.loc[rest["Compound ID"] == compound][
                    "Ligand Cluster number"
                ] = cluster_number

        # save the data to the new parquet file
        new_parquet.write(rest, new_parquet.schema)


if __name__ == "__main__":
    """
        PUBCHEM_FASTA = "../PUBCHEM_HDF5.fasta"
        #Parse the fasta files to get the protein ids
        PubChem_protein_ids = fasta_parser(PUBCHEM_FASTA)
        logger.info("Number of proteins in PubChem: {}", len(PubChem_protein_ids))
        valid_pubchem = pdb_validations( "../gene-id_with_pdb.csv")
        df = pd.DataFrame(valid_pubchem.items(), columns=['Gene_id','PDB'])
        df.to_csv("../valid_pubchem.csv", index=False)

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
                cluster_id = line.split("Cluster")[1].split(":")[1].strip()
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

                        continue
        else:
            i += 1
            new_lines.append(line)

    # write the new_lines to a new file
    with open("../new_clusters_cluster_parsed.csv", "w", encoding="utf-8") as f:
        f.writelines(new_lines)
    # find the number of representatives in the new_clusters_cluster_parsed.csv file
    with open("../clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    clusters = [line for line in lines if line.startswith("Cluster")]
    logger.info("Number of clusters before parsing: {}", len(clusters))
    # find the number of representatives using lines starting with Cluster
    with open("../new_clusters_cluster_parsed.csv", "r", encoding="utf-8") as f:
        lines = f.readlines()
    clusters = [line for line in lines if line.startswith("Cluster")]
    logger.info("Number of clusters after parsing: {}", len(clusters))
    # find the number of proteins in ../PUBCHEM_HDF5.fasta file
    with open("../PUBCHEM_HDF5.fasta", "r", encoding="utf-8") as f:
        lines = f.readlines()
    proteins = [line for line in lines if line.startswith(">")]
    logger.info("Number of gene ids in PUBCHEM_FASTA: {}", len(proteins))

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
    logger.info("Length of gene ids after parsing: {}", len(new_lines))

    #make_protein_cluster_file("../clustered_pubchem.parquet", "../new_clusters_cluster_parsed.csv")
    # read the parquet file

    logger.debug(data.scanner(filter=ds.field("Cluster number") == 1).to_table())
    # read the csv files
    aid_2_target = pd.read_csv("../aid_2_target.csv")
    aid_2_cids = pd.read_csv("../aid_2_cids.csv")
    # loop through the old_data and update the member_compound column
    for index, row in old_data.iterrows():
        gene_id = row["members"]
        member = row["member_compound"]
        # get the aids for each gene_id
        aids = list(
            aid_2_target[aid_2_target["Gene_id"].astype(str) == str(gene_id)]["AID"]
        )
        # loop through the aids and get the cids for each aid
        cids = []
        for aid in aids:
            # add the cids to the cids list knowing that the cids are ['5291'] sting of list
            cids.extend(
                ast.literal_eval(
                    aid_2_cids[aid_2_cids["AID"] == aid]["CIDs"].to_list()[0]
                )
            )
        # update the member_compound column
        old_data.at[index, "member_compound"] = cids
        logger.info(old_data.at[index, "member_compound"])
    logger.info(old_data)
    # write the updated data to the parquet file
    ParquetParser("../full_clustered_pubchem.parquet").write(old_data, schema)

    # read the parquet file
    parquet = ParquetParser("../full_clustered_pubchem.parquet")
    schema = pa.schema(
        [
            ("Cluster number", pa.int64()),
            ("Represenatives", pa.string()),
            ("members", pa.string()),
            ("member_compound", pa.list_(pa.string())),
        ]
    )
    data = parquet.convert_to_pandas(schema)
    logger.info(data.tail())
    origi_parquet = ParquetParser("../clustered_pubchem.parquet")
    origi_data = origi_parquet.convert_to_pandas(schema)
    logger.info(origi_data.tail())

    schema = pa.schema(
        [
            ("Cluster number", pa.int64()),
            ("Represenatives", pa.string()),
            ("members", pa.string()),
            ("member_compound", pa.list_(pa.string())),
        ]
    )

    make_protein_cluster_file("clustered_pubchem.parquet", "../new_clusters_cluster_parsed.csv")
    hdf5_result=hdf5_raw_file_parser('../PubChem_data_edited_400k.hdf5')
    add_compounds(ParquetParser("clustered_pubchem.parquet",schema=schema), hdf5_result,overwrite=False)
    """
    # read all .cif files in the ../representatives folder and convert them to .pdb files
    cif_files = [
        file for file in os.listdir("../representatives") if file.endswith(".cif")
    ]
    logger.info(f"Number of cif files: {len(cif_files)}")
    for file in cif_files:
        if os.path.exists(f"../representatives/{file.replace('.cif','.pdb')}"):
            continue
        logger.info(f"Converting {file} to .pdb file")
        convert_cif2pdb(
            f"../representatives/{file}",
            f"../representatives/{file.replace('.cif','.pdb')}",
        )
        logger.info(f"Converted {file} to .pdb file")
    # read aid_2_target and aid_2_cids from the csv files
    aid_2_target = pd.read_csv("../aid_2_target.csv")
    aid_2_cids = pd.read_csv("../aid_2_cids.csv")
    hdf5 = HDF5Parser("../PubChem_data_edited_420k.hdf5")
    # read the parquet file
    schema = pa.schema(
        [
            ("Cluster number", pa.int64()),
            ("Represenatives", pa.string()),
            ("members", pa.string()),
            ("member_compound", pa.list_(pa.string())),
        ]
    )
    old_parquet = ParquetParser("../protein_clustered_data.parquet", schema)
    # get cluster number 1
    # data=parquet.convert_to_table().group_by("Cluster number")
    new_schema = pa.schema(
        [
            ("Protein Cluster number", pa.int64()),
            ("PDB/Gene ID", pa.string()),
            ("Compound ID", pa.string()),
            ("Smiles", pa.string()),
            ("Ligand Cluster number", pa.int64()),
        ]
    )
    # make a new parquet file with the new schema and the data from the trial_data
    new_parquet = ParquetParser("compounds_clustered_pubchem.parquet", new_schema)

    make_ligand_cluster_file(
        old_parquet,
        new_parquet,
        aid_2_cids,
        aid_2_target,
        hdf5_file="../PubChem_data_edited_400k.hdf5",
    )
    """
    # group the data by the cluster number
    initial_data = []
    trial_data = trial_data.groupby("Cluster number")
    # loop through the data and get the gene ids and the cids
    for cluster_number, data in trial_data:
        for pdb_id in data["Represenatives"]:
            protein = Protein()
            protein.load(file_path=f"../representatives/{pdb_id.lower()}.pdb")
            logger.debug(protein._pdb_file_path)
            ligands = get_ligand_names(protein._pdb_file_path)
            if len(ligands) > 1:
                for ligand in ligands:
                    smiles = get_smiles(ligand)
                    initial_data.append(
                        (cluster_number, pdb_id, ligand, smiles["stereo_smiles"], None)
                    )
            else:
                smiles = get_smiles(ligands[0])
                initial_data.append(
                    (cluster_number, pdb_id, ligands[0], smiles["stereo_smiles"], None)
                )
        for member in data["members"]:
            # find the aids for each gene id
            aids = list(
                aid_2_target[aid_2_target["Gene_id"].astype(str) == str(member)]["AID"]
            )
            # loop through the aids and get their cids from the aid_2_cids dataframe
            for aid in aids:
                logger.debug(f"aid: {aid}")
                cids = ast.literal_eval(
                    aid_2_cids[aid_2_cids["AID"] == aid]["CIDs"].to_list()[0]
                )
                for cid in cids:
                    # get the smiles for each cid using hdf5.read(f'aids/{aid}/cids/{cid}/smiles')
                    try:
                        smiles = hdf5.read(f"aids/{aid}/cids/{cid}/smiles").astype(str)
                        logger.debug(f"smiles: {smiles}")
                        initial_data.append(
                            (cluster_number, member, cid, smiles[0], None)
                        )
                        logger.info(f"Data: {initial_data}")
                    except Exception as e:
                        logger.error(f"Error: {e}")
                        continue
    new_parquet.write(initial_data, new_schema)
    #read the protein_clustered_data.parquet file and find the length of the data
    schema = pa.schema(
        [
            ("Cluster number", pa.int64()),
            ("Represenatives", pa.string()),
            ("members", pa.string()),
            ("member_compound", pa.list_(pa.string())),
        ]
    )
    parquet = ParquetParser("../protein_clustered_data.parquet", schema)
    logger.info(parquet.convert_to_pandas().head())
    logger.info(parquet.convert_to_pandas().tail())
    #find if there any rows with missing values in member_compound column
    logger.info(f'The length of empty columns\n:{parquet.convert_to_pandas().isnull().sum()}')
    #length of the data in the parquet file
    logger.info(f"Length of entire data in parquaet file, {len(parquet.convert_to_pandas())}")
    #extract the member_compound column from the parquet file and extract all the unique values
    member_compound = parquet.convert_to_pandas()["member_compound"]
    #append all the values in the member_compound column to a list
    all_compounds = []
    for compounds in member_compound:
        all_compounds.extend(compounds)
    #logger.info(all_compounds)
    logger.info(f"Length of all compounds left after protein clustering,{len(list(set(all_compounds)))}")
    #get the length of member_compound column in the parquet file where the members == 25
    cids=parquet.convert_to_pandas()[parquet.convert_to_pandas()['members'] == '25']
    logger.debug(cids)
    logger.info(f"{set(cids['member_compound'][0])}")
    logger.info(f"Length of member_compound column where members == 25: {len(list(set((cids['member_compound'][0]))))}")

    matplotlib.use('Agg')  # or another appropriate backend like 'Qt5Agg', 'Agg', etc.
    #read the compound_clustered_pubchem.parquet file
    schema = pa.schema(
        [
            ("Protein Cluster number", pa.int64()),
            ("PDB/Gene ID", pa.string()),
            ("Compound ID", pa.string()),
            ("Smiles", pa.string()),
            ("Ligand Cluster number", pa.int64()),
        ]
    )
    parquet = ParquetParser("../compounds_clustered_pubchem.parquet", schema)
    logger.info(parquet.convert_to_pandas().head())
    logger.info(parquet.convert_to_pandas().tail())
    #find if there any rows with missing values in Smiles column
    logger.info(f'The length of empty columns\n:{parquet.convert_to_pandas().isnull().sum()}')
    #length of the data in the parquet file
    logger.info(f"Length of entire data in parquaet file, {len(parquet.convert_to_pandas())}")
    #READ THE DATA FROM THE compound_clustered_pubchem.parquet file
    data = parquet.convert_to_pandas().groupby("Protein Cluster number")
    TanimotoClustering = TanimotoClustering()
    # loop through each group and get the smiles for all the compounds in the group
    cutoff_values = range(10, 100, 5)

    for cluster, rest in data:
        compounds_count = []
        logger.debug(f"Length of compounds with duplication: {len(rest)}")
        cluster_number=cluster
        # get the rows with unique values in the compound id column
        unique_rows = rest.drop_duplicates(subset='Compound ID',ignore_index=True)
        logger.debug(f'unique rows: {unique_rows}')
        logger.debug(f"Length of unique rows: {len(unique_rows)}")
        smiles = unique_rows['Smiles']
        compound_ids = unique_rows['Compound ID']
        logger.debug(f"Length of smiles: {len(smiles)}")
        original_length=len(smiles)
        # Get fingerprints for each molecule and tanimoto similarity as one liner
        fingerprints = [TanimotoClustering.get_morgan_fingerprint(smile) for smile in smiles]
        logger.debug(f'unique fingerprints: {len(fingerprints)}')
        similarity = []
        # loop over the fingerprints to get the similarity for each molecule with the previous molecules
        for i in range(1, len(fingerprints)):
            similarity.append(
                TanimotoClustering.tanimoto_similarity(fingerprints[:i], fingerprints[i])
            )
        logger.debug(f'Length of similarity: {len(similarity)}')
        for cutoff in cutoff_values:
            compounds_number = []
            cutoff = cutoff / 100
            clusters = TanimotoClustering.cluster_tanimoto(similarity, compound_ids, cutoff)
            count = 0
            # read each cluster and find the clusters where it has a member with whose id has letters in it
            for cluster in clusters:
                for member in cluster:
                    if any(char.isalpha() for char in member) and len(cluster) > 1:
                        compounds_number.extend(cluster)
                        count += 1
                        break
            compounds_count.append(len(compounds_number))
            logger.info(f"Number of clusters with members with letters for cutoff {cutoff}: {count}")
            logger.info(f"Number of compounds with letters for cutoff {cutoff}: {len(compounds_number)}")
        #get the cutoff values as a list
        logger.debug(list(cutoff_values))
        logger.debug(compounds_count)
        logger.debug(f'Plotting the graph for cluster number: {cluster_number}')
        # Plotting the graph
        plt.plot(cutoff_values, compounds_count)
        plt.xticks(range(10, 105, 5))
        plt.xlabel('Similarity Cutoff')
        plt.ylabel('Number of Compounds')
        plt.title('Effect of Similarity Cutoff on Number of Compounds')
        #add the number of original compounds on the plot as a horizontal line with red color and dashed line style
        # and label it as Original Compounds with the number of original compounds
        plt.axhline(y=original_length, color='red', linestyle='--', label=f'Original Compounds:{original_length}')
        #add a grid to the plot
        #save the plot to a file with the cluster number as the name with similarity_cutoff.png
        plt.savefig(f'../{cluster_number}_similarity_cutoff.png')
        plt.clf()

"""
