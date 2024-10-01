r" Implemtnation for a wrapper for MMseqs2 clustering"
from typing import TextIO

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
    get_ligand_names,
    get_rcsb_data,
    get_smiles,
    map_genid_to_pdb,
    validate_ligands,
    validate_pdb,
)

# OUTLINE
# 1. Read the PubChem HDF5 file and extract the fasta sequences
# 2. map the gene id cluster representatives to the PDB ids using ID mapping webservices
# 3. Validate the PDB ids and write the valid PDB ids to a new file
# 5. Cluster the fasta sequences using MMseqs2
# 3. Parse the MMseqs2 output to find representatives and members
# 4. Validate the PDB ids and write the valid PDB ids to a new file
# and write the clusters to a new file
# 6. Write the clusters data into protein_cluster.parquet file
# 7. Parse the hdf5 file to get the aids and cids
# 8. Add the compounds to the parquet file with the protein clusters from the protein_cluster.parquet file
# 9. Run tanimoto clustering on the compounds and add the ligand clusters to the parquet file
# 10. Write the ligand clusters data into ligand_cluster.parquet file


def create_fasta_file(hdf5_file: str, fasta_file: str) -> None:
    r"""Create a FASTA file from the HDF5 file.
    Parameters
    ----------
    hdf5_file : str
        Path to the HDF5 file.
    fasta_file : str
        Path to the new FASTA file.
    """
    # check if the HDF5 file exists
    if not os.path.exists(hdf5_file) or not hdf5_file.endswith(".hdf5"):
        raise FileNotFoundError(f"HDF5 file {hdf5_file} not found or not a valid file.")
    # read the HDF5 file
    hdf5 = HDF5Parser(hdf5_file)
    # read the hdf5 file
    aids = hdf5.read("aids")
    logger.info(f"Number of aids: {len(aids)}")
    # loop through the aids and get the gene ids
    with open(fasta_file, "w", encoding="utf-8") as file:
        for aid in aids:
            if "protein_sequence" in hdf5.read(f"aids/{aid}"):
                gene_id = hdf5.read(f"aids/{aid}/targets_gene_id")[0]
                sequence = hdf5.read(f"aids/{aid}/protein_sequence")
                sequence = "".join([byte.decode("utf-8") for byte in sequence])
                file.write(f">{gene_id}\n{sequence}\n")
            else:
                logger.error(f"No protein sequence found for aid {aid}")
                continue


def fasta_parser(fasta: str | TextIO, delimiter: str | None = None) -> list:
    r"""Parse FASTA files to get the protein ids
    Parameters
    ----------
    fasta : str |TextIO
        Path to the FASTA file.
    delimiter : str| None
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


def pdb_validations(csvfilenames: str) -> None:
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
    fasta : str| TextIO
        Path to the FASTA file.
    outfile_name : str
        Path to the new FASTA file.
    csvfilenames : str
        Path to the CSV file containing the protein ids.
    delimiter : str| None (default =|)
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
    data_source: tuple[dict[str, str]],
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
    initial_data: list | None = None,
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
    start_time = time.time()
    if os.path.exists("backup.parquet"):
        done_data = ParquetParser(
            "backup.parquet", new_file_name.schema
        ).convert_to_pandas()
    else:
        done_data = pd.DataFrame()
    initial_data = initial_data if initial_data is not None else []
    hdf5 = HDF5Parser("../PubChem_data_edited.hdf5")
    old_data = old_parquet.convert_to_pandas().groupby("Cluster number")
    # check if the progress_cache.pkl file exists and if so read it
    if os.path.exists("cache.pkl"):
        try:
            with open("cache.pkl", "rb") as f:
                initial_data = pickle.load(f)
        except Exception as e:
            logger.error(f"Error: {e}")
    # loop through the data and get the gene ids and the cids
    for cluster_number, data in old_data:
        logger.debug(f'Cluster represenatives: {data["Represenatives"]}')
        for pdb_id in data["Represenatives"]:
            # check if the pdb_id is in the done_data
            if any(pdb_id == x[1] for x in done_data):
                logger.debug(f"pdb_id: {pdb_id} already in done_data. Skipping")
                continue
            # check if the data["Represenatives"] is in the initial_data pdb_id
            if any(pdb_id == x[1] for x in initial_data):
                logger.debug(f"pdb_id: {pdb_id} already in initial_data. Skipping")
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
                try:
                    smiles = get_smiles(ligands[0])
                    initial_data.append(
                        (
                            cluster_number,
                            pdb_id,
                            ligands[0],
                            smiles["stereo_smiles"],
                            None,
                        )
                    )
                except Exception as e:
                    logger.error(f"Error: {e}")
                    initial_data.append(
                        (cluster_number, pdb_id, ligands[0], None, None)
                    )
                    continue
            else:
                logger.error(f"No ligands found for {pdb_id}")
                continue
        for member in data["members"]:
            # check if the data["members"] is in the initial_data member
            if any(member == x[2] for x in initial_data):
                logger.debug(f"member: {member} already in initial_data. Skipping")
                continue
            if any(member == x[2] for x in done_data):
                logger.debug(f"member: {member} already in done_data. Skipping")
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
            cache_file = "cache.pkl"
            logger.warning(f"Saving progress to {cache_file}")
            with open(cache_file, "wb") as f:
                pickle.dump(initial_data, f)
            # write the parquet file every 3hr
            if time.time() - start_time > 60 * 60 * 3:
                # make a parquet file with the initial_data named backup
                backup = ParquetParser("backup.parquet", new_file_name.schema)
                backup.write(initial_data, new_file_name.schema)
                logger.info(
                    f"Data written to backup Parquet file at {backup.file_path}"
                )
                start_time = time.time()
                # read the backup file as a pandas dataframe
                done_data = backup.convert_to_pandas()

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
    data = old_parquet.convert_to_pandas().groupby("Protein Cluster number")
    tanimoto = TanimotoClustering()
    all_results = []
    for cluster, rest in data:
        logger.debug(f"Length of compounds: {len(rest)}")
        # get the rows with unique values in the compound id column
        unique_rows = rest.drop_duplicates(subset=["Compound ID"])
        logger.debug(f"Length of unique rows: {len(unique_rows)}")
        smiles = unique_rows["Smiles"]
        compound_ids = unique_rows["Compound ID"].tolist()
        logger.debug(f"Length of smiles: {len(smiles)}")
        original_length = len(smiles)
        # if length of the smiles is less than 2, skip the cluster
        if original_length < 2:
            logger.error(f"Length of smiles: {original_length} is less than 2")
            continue

        # Get fingerprints for each molecule and tanimoto similarity as one liner
        fingerprints = []
        valid_compound_ids = []
        for smile, compound_id in zip(smiles, compound_ids):
            if smile is not None:
                fingerprint = tanimoto.get_morgan_fingerprint(smile)
                if fingerprint is not None:
                    fingerprints.append(fingerprint)
                    valid_compound_ids.append(compound_id)

        # Update compound_ids to only include valid ones
        compound_ids = valid_compound_ids
        """
        # Get fingerprints for each molecule and tanimoto similarity as one liner
        fingerprints = [tanimoto.get_morgan_fingerprint(smile) for smile in smiles]

        # remove the None values from the fingerprints and their corresponding compound_ids
        fingerprints = [
            fingerprint for fingerprint in fingerprints if fingerprint is not None
        ]
        compound_ids = [
            compound_id
            for compound_id, fingerprint in zip(compound_ids, fingerprints)
            if fingerprint is not None
        ]
        """
        # remove the corresponding compound_ids from the rest dataframe and reset the index
        rest = rest[rest["Compound ID"].isin(compound_ids)]
        logger.debug(f"unique fingerprints: {len(fingerprints)}")
        similarity = []
        # Loop over the fingerprints to get the similarity for each molecule with the previous molecules
        for i in range(1, len(fingerprints)):
            similarity.append(
                tanimoto.tanimoto_similarity(fingerprints[:i], fingerprints[i])
            )
        logger.debug(f"Length of similarity: {len(similarity)}")
        if len(similarity) == 0:
            logger.error("Length of similarity is 0")
            continue
        # Run TanimotoClustering.cluster_tanimoto(similarity, compound_ids, cutoff)
        clusters = tanimoto.cluster_tanimoto(
            similarity, compound_ids, similarity_cutoff
        )
        # Add a column for Ligand Cluster number if it doesn't exist
        if "Ligand Cluster number" not in rest.columns:
            rest["Ligand Cluster number"] = None

        # Loop through the clusters and add the cluster number to the dataframe Ligand Cluster number
        for cluster_number, cluster in enumerate(clusters, start=1):
            for compound in cluster:
                rest.loc[rest["Compound ID"] == compound, "Ligand Cluster number"] = (
                    cluster_number
                )

        all_results.append(rest)

    # Concatenate all results and save the data to the new parquet file
    final_result = pd.concat(all_results)
    new_parquet.write(final_result, new_parquet.schema)


if __name__ == "__main__":
    HDF5_FILE = "../PubChem_data_edited.hdf5"
    FASTA_FILE = "../protein_sequences.fasta"
    GENE_ID2PDB_ID_FILE = "gene_id2pdb_id.csv"
    # create_fasta_file(HDF5_FILE, FASTA_FILE)
    mmseqs_cluster(FASTA_FILE, outfile_name_suffix="../clusters", tmp_dir="../tmp")
    PubChem_protein_ids = fasta_parser(FASTA_FILE)
    logger.info(f"Number of protein ids: {len(set(PubChem_protein_ids))}")
    gene2pdb = {}
    batch_size = 400
    unique_protein_ids = list(set(PubChem_protein_ids))
    for i in range(0, len(unique_protein_ids), batch_size):
        batch = [
            protein_id.strip() for protein_id in unique_protein_ids[i : i + batch_size]
        ]
        logger.info(f"Processing batch {i // batch_size + 1}: {batch}")
        batch_results = map_genid_to_pdb(batch)
        for result in batch_results:
            gene_id = result["Gene ID"]
            pdb_id = result["PDB IDs"]
            gene2pdb[gene_id] = pdb_id
    # write the gene2pdb dictionary to a csv file
    gene2pdb_df = pd.DataFrame(gene2pdb.items(), columns=["From", "PDB"])
    gene2pdb_df.to_csv(GENE_ID2PDB_ID_FILE, index=False)
    # validate the PDB files
    # pdb_validations("gene_id2pdb_id.csv")
    # NOTE:THESE FILES ARE NOT REAL
    protein_cluster_file = "../../protein_cluster.parquet"
    ligand_cluster_file = "../../ligand_cluster.parquet"
    # TODO:Continue rewiting the code from here forward
