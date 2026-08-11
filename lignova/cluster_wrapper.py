# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"Implementation for a wrapper for MMseqs2 clustering"

import ast
import os
import time
from typing import TextIO

import pandas as pd
import pyarrow as pa
from loguru import logger

from lignova.clustering.mmseq import mmseqs_cluster, mmseqs_parser
from lignova.clustering.tanimoto import TanimotoClustering
from lignova.hdf5.parquet import ParquetParser
from lignova.hdf5.parser import HDF5Parser
from lignova.structure import Protein
from lignova.structure.editing import convert_cif2pdb
from lignova.structure.utils import (
    get_ligand_names,
    get_smiles,
    map_genid_to_pdb,
    validate_ligands,
    validate_pdb,
)

# OUTLINE
# 1. Read the PubChem HDF5 file and extract the fasta sequences
# 2. Cluster the fasta sequences using MMseqs2
# 3. map the gene id to the PDB ids using ID mapping webservices
# 4. Validate the PDB ids and write the valid PDB ids to a new file
# 5. Parse the MMseqs2 output to find representatives and members
# of the clusters, clear the clusters with no representatives
# 6. Write the clusters data into protein_cluster.parquet file
# 7. Parse the hdf5 file to get the aids and cids
# 8. Add the compounds to the parquet file with the protein clusters
# from the protein_cluster.parquet file
# 9. Run tanimoto clustering on the compounds and add the ligand clusters to the parquet file
# 10. Write the ligand clusters data into parquet file


def create_fasta_file(hdf5_file: str, fasta_file: str) -> None:
    r"""Create a FASTA file from the HDF5 file.

    Args:
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
    done_gene_ids = []
    # loop through the aids and get the gene ids
    with open(fasta_file, "w", encoding="utf-8") as file:
        for aid in aids:
            if "protein_sequence" not in hdf5.read(f"aids/{aid}"):
                logger.error(f"No protein sequence found for aid {aid}")
                continue
            gene_id = hdf5.read(f"aids/{aid}/targets_gene_id")[0]
            if gene_id in done_gene_ids:
                continue
            sequence = hdf5.read(f"aids/{aid}/protein_sequence")
            sequence = "".join([byte.decode("utf-8") for byte in sequence])
            file.write(f">{gene_id}\n{sequence}\n")
            done_gene_ids.append(gene_id)


def fasta_parser(fasta: str | TextIO, delimiter: str | None = None) -> list:
    r"""Parse FASTA files to get the protein ids

    Args:
    fasta : str |TextIO
        Path to the FASTA file.
    delimiter : str| None
        Delimiter to split the protein id from the FASTA header. Default is None.
    Returns:
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

    Args:
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
    new_file = {}
    # check if f"valid_{csvfilenames}" exists and if so read it as a dictionary
    if os.path.exists(f"valid_{csvfilenames}"):
        file_data = pd.read_csv(f"valid_{csvfilenames}")
        # make the new_file as a dictionary
        # Create the dictionary with Gene_id as keys and PDB as values
        for gene_id, pdb_list in zip(file_data["Gene_id"], file_data["PDB"]):
            new_file[gene_id] = ast.literal_eval(pdb_list)
    for gene_id, pdb_id in zip(protein_ids["From"], pdb_ids):
        # read the pdb_id as a list
        pdb_id = ast.literal_eval(pdb_id)
        logger.debug(f"Gene id: {gene_id} PDB id: {pdb_id}")
        # find the value of the From column
        if gene_id not in new_file:
            new_file[gene_id] = []
        for pdb in pdb_id:
            # validate the pdb id
            # check if the pdb id is valid is in the new_file[gene_id] values
            if new_file[gene_id] is not None:
                if pdb in new_file[gene_id]:
                    logger.debug(
                        f"PDB id {pdb} already in new_file for gene_id {gene_id}"
                    )
                    continue
            if not validate_pdb(pdb) or not validate_ligands(pdb):
                logger.error(f"Invalid PDB id {pdb} for protein {gene_id}")
                continue
            logger.info(f"the current new_file dict: {new_file[gene_id]}")
            logger.info(f"the type of new_file dict: {type(new_file[gene_id])}")
            logger.info(f"Valid PDB id {pdb} for gene_id {gene_id}")
            tmp = new_file[gene_id]
            tmp.append(pdb)
            new_file[gene_id] = tmp
            # save this new file to a csv file every 10 minutes
            if time.time() - start_time > 60 * 10:
                new_csv_df = pd.DataFrame(new_file.items(), columns=["Gene_id", "PDB"])
                new_csv_df.to_csv(f"valid_{csvfilenames}", index=False)
                start_time = time.time()
    # save the new_file to a csv file
    new_csv_df = pd.DataFrame(new_file.items(), columns=["Gene_id", "PDB"])
    new_csv_df.to_csv(f"valid_{csvfilenames}", index=False)


def clean_cluster_files(cluster_filepath: str, valid_csvpath: str) -> None:
    """
    Clean the cluster file by removing the clusters with no raw
    no members and the clusters with invalid PDB ids.

    Args:
    cluster_filepath : str
        Path to the parsed mmseq cluster file.
    valid_csvpath : str
        Path to the csv mapping the gene ids to the valid PDB ids.
    Returns:
    None
    """
    # check if the cluster file exists
    if not os.path.exists(cluster_filepath):
        raise FileNotFoundError(f"Cluster file {cluster_filepath} not found.")
    if not os.path.exists(valid_csvpath):
        raise FileNotFoundError(f"Valid CSV file {valid_csvpath} not found.")

    # Read the valid PDB ids from the CSV file
    valid_pdbs = pd.read_csv(valid_csvpath)
    valid_pdbs_dict = dict(zip(valid_pdbs["Gene_id"], valid_pdbs["PDB"]))

    # Read the cluster file
    cluster_data = pd.read_csv(cluster_filepath)
    logger.info(f"Number of clusters: {len(cluster_data)}")
    # Ensure the cluster file has the required columns
    if (
        "representatives" not in cluster_data.columns
        or "members" not in cluster_data.columns
    ):
        raise ValueError(
            "Cluster file must contain 'representatives' and 'members' columns."
        )

    # Filter clusters based on valid representatives or members
    valid_clusters = []
    for _, row in cluster_data.iterrows():
        representatives = int(row["representatives"])
        members = ast.literal_eval(row["members"])
        # Check if any representative is valid
        new_reps = []
        if representatives in valid_pdbs_dict.keys():
            new_reps.extend(ast.literal_eval(valid_pdbs_dict[representatives]))
        # Check if any member maps to a valid PDB id
        for member in members:
            if int(member) in valid_pdbs_dict.keys():
                new_reps.extend(ast.literal_eval(valid_pdbs_dict[member]))
        row["representatives"] = list(set(new_reps))
        row["members"] = members
        valid_clusters.append(row)

    # Create a DataFrame from the valid clusters and save it
    valid_clusters_df = pd.DataFrame(valid_clusters)
    # delete the rows with empty lists in the representatives column
    valid_clusters_df = valid_clusters_df[
        valid_clusters_df["representatives"].apply(lambda x: len(x) > 0)
    ]
    valid_clusters_df.to_csv(f"valid_{os.path.basename(cluster_filepath)}", index=False)
    logger.info(
        f"Cleaned cluster file saved to valid_{os.path.basename(cluster_filepath)}"
    )


def make_protein_cluster_file(new_file: ParquetParser, cluster_csv: str) -> None:
    r"""Make a  parquet cluster file with the new lines.

    Args:
    new_file : ParquetParser
        Parquet object. if the file not found, it will be created.
    cluster_csv : str
        Path to the cluster file.
    """
    # check if the cluster file exists
    if not os.path.exists(cluster_csv) or not cluster_csv.endswith(".csv"):
        raise FileNotFoundError(
            f"Cluster file {cluster_csv} not found or not a valid file."
        )
    # check the new file exists
    if not os.path.exists(new_file.file_path):
        logger.info(f"New file {new_file.file_path} not found. Creating it.")
        new_file.create()
    # Read the cluster CSV file using pandas
    cluster_data = pd.read_csv(cluster_csv)
    logger.info(f"Number of clusters: {len(cluster_data)}")

    # Ensure the cluster file has the required columns
    if (
        "representatives" not in cluster_data.columns
        or "members" not in cluster_data.columns
    ):
        raise ValueError(
            "Cluster file must contain 'representatives' and 'members' columns."
        )

    data = []
    for cluster_number, row in cluster_data.iterrows():
        representatives = ast.literal_eval(row["representatives"])
        members_list = ast.literal_eval(row["members"])
        for rep in representatives:
            for member in members_list:
                data.append((cluster_number + 1, str(rep), str(member), []))
    new_file.write(data)
    logger.info(f"Data written to Parquet file at {new_file.file_path}")


def hdf5_raw_file_parser(hdf5_file: str) -> tuple[dict, dict]:
    r"""Parse the HDF5 file to get the protein ids and the ligand ids.

    Args:
    hdf5_file : str
        Path to the HDF5 file.

    Returns:
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
    aid_2_target_dict = {}
    aid_2_cids_dict = {}
    logger.info(f"Number of aids: {len(aids)}")
    for aid in aids:
        gene_id = hdf5.read(f"/aids/{aid}/targets_gene_id")
        cids = hdf5.read(f"/aids/{aid}/cids")
        logger.debug(f"cids: {cids}")
        aid_2_target_dict[aid] = gene_id[0]
        aid_2_cids_dict[aid] = cids
    logger.debug(f"Number of aids: {len(aid_2_target_dict)}")
    logger.debug(f"Number of aids: {len(aid_2_cids_dict)}")
    return aid_2_target_dict, aid_2_cids_dict


def add_compounds(
    original_file: ParquetParser,
    data_source: tuple[dict[str, str]],
    overwrite: bool = True,
) -> None:
    r"""Add the compounds to the parquet file.

    Args:
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
        # get the aids for each gene_id
        aids = list(
            aid_2_target[aid_2_target["Gene_id"].astype(str) == str(gene_id)]["AID"]
        )
        # loop through the aids and get the cids for each aid
        cids = []
        for aid in aids:
            # add the cids to the cids list knowing that the cids are ['5291'] sting of list
            cids.extend(aid_2_cids[aid_2_cids["AID"] == aid]["CIDs"].to_list()[0])
        # update the member_compound column
        old_data.at[index, "member_compound"] = cids
    if overwrite:
        original_file.write(old_data)
        logger.info(f"Data written to Parquet file at {original_file.file_path}")
    else:
        new_file = "new_" + original_file.file_path
        new_parser = ParquetParser(new_file, original_file.schema)
        new_parser.write(old_data)
        logger.info(f"New file created at {new_file}")


def make_ligand_cluster_file(
    old_parquet: ParquetParser,
    new_file_name: ParquetParser,
    aid_2_cids_df: pd.DataFrame,
    aid_2_target_df: pd.DataFrame,
    hdf5_file: str,
    initial_data: list | None = None,
) -> None:
    r"""Make a new Parquet file for the ligand clustering.

    Args:
    old_parquet : ParquetParser
        ParquetParser object of the old file.
    new_file_name : ParquetParser
        ParquetParser object of the new file.
    aid_2_cids_df : pd.DataFrame
        DataFrame containing the aids and cids.
    aid_2_target_df : pd.DataFrame
        DataFrame containing the aids and gene ids.
    hdf5_file : str
        Path to the HDF5 file.
    intial_data : list
        List of initial data to be added to the new file. Default is None.
    """
    start_time = time.time()
    if os.path.exists("backup.parquet"):
        done_data = ParquetParser(
            "backup.parquet", new_file_name.schema
        ).convert_to_pandas()
    else:
        done_data = pd.DataFrame()
    initial_data = initial_data if initial_data is not None else []
    hdf5 = HDF5Parser(hdf5_file)
    old_data = old_parquet.convert_to_pandas().groupby("Cluster number")
    # loop through the data and get the gene ids and the cids
    for cluster_number, data in old_data:
        logger.debug(f'Cluster represenatives: {data["Represenatives"]}')
        for pdb_id in data["Represenatives"]:
            logger.info(f"pdb_id: {pdb_id}")
            # check if the pdb_id is in the done_data
            if any(pdb_id == x[1] for x in done_data) or any(
                pdb_id == x[1] for x in initial_data
            ):
                logger.debug(f"pdb_id: {pdb_id} already in done_data. Skipping")
                continue
            protein = Protein()
            # check if the pdb file exists
            if not os.path.exists(f"raw/{pdb_id.lower()}.pdb"):
                # then use the protein.load function to get the pdb file
                file_ext = (
                    "pdb"
                    if protein.get_pdb_from_rcsb(pdb_id).startswith("HEADER")
                    else "cif"
                )
                logger.info(f"Downloading {pdb_id} from RCSB with extension {file_ext}")
                protein.load(
                    pdb_id=pdb_id,
                    write=True,
                    write_path=f"raw/{pdb_id.lower()}.{file_ext}",
                )
                if file_ext == "cif":
                    logger.info(f"Converting {pdb_id} to pdb format")
                    convert_cif2pdb(
                        f"raw/{pdb_id.lower()}.cif", f"raw/{pdb_id.lower()}.pdb"
                    )
            protein = Protein(file_path=f"raw/{pdb_id.lower()}.pdb")
            ligands = get_ligand_names(protein.file_path)
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
            if any(member == x[1] for x in done_data) or any(
                member == x[1] for x in initial_data
            ):
                logger.debug(f"member: {member} already in done_data. Skipping")
                continue
            # find the aids for each gene id
            aids = list(
                aid_2_target_df[aid_2_target_df["Gene_id"].astype(str) == str(member)][
                    "AID"
                ]
            )
            # loop through the aids and get their cids from the aid_2_cids_df dataframe
            for aid in aids:
                logger.debug(f"aid: {aid}")
                cids = ast.literal_eval(
                    aid_2_cids_df[aid_2_cids_df["AID"] == aid]["CIDs"].to_list()[0]
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
            # save the progress to backup parquet file every 20 minutes
            if time.time() - start_time > 60 * 20:
                # make a parquet file with the initial_data named backup
                backup = ParquetParser("backup.parquet", new_file_name.schema)
                backup.write(initial_data)
                logger.info(
                    f"Data written to backup Parquet file at {backup.file_path}"
                )
                start_time = time.time()
                # read the backup file as a pandas dataframe
                done_data = backup.convert_to_pandas()

    new_file_name.write(initial_data)


def add_smiles_cluster(
    old_parquet: ParquetParser, similarity_cutoff: float, new_parquet: ParquetParser
):
    r"""Clustering smiles and Adding the smiles cluster to the parquet file.

    Args:
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
        smiles = [
            smile.decode("utf-8") if isinstance(smile, bytes) else smile
            for smile in unique_rows["Smiles"]
        ]
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
        # remove the corresponding compound_ids from the rest dataframe and reset the index
        rest = rest[rest["Compound ID"].isin(compound_ids)]
        logger.debug(f"unique fingerprints: {len(fingerprints)}")
        similarity = []
        # Loop over the fingerprints to get the similarity for each molecule
        # with the previous molecules
        for index in range(1, len(fingerprints)):
            similarity.append(
                tanimoto.tanimoto_similarity(fingerprints[:index], fingerprints[index])
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

        # Loop through the clusters and add the cluster number to the Ligand Cluster number
        for cluster_number, cluster in enumerate(clusters, start=1):
            for compound in cluster:
                rest.loc[rest["Compound ID"] == compound, "Ligand Cluster number"] = (
                    cluster_number
                )

        all_results.append(rest)
    # Concatenate all results and save the data to the new parquet file
    final_result = pd.concat(all_results)
    new_parquet.write(final_result)


# TODO: Turn into a run function for CLI
if __name__ == "__main__":
    HDF5_FILE = "../PubChem_data_edited.hdf5"
    FASTA_FILE = "../protein_sequences.fasta"
    GENE_ID2PDB_ID_FILE = "gene_id2pdb_id.csv"
    RAW_ID_MAPPING_FILE = "id_mapping.csv"
    VALID_PDBS_FILE = "valid_gene_id2pdb_id.csv"
    PROTEIN_CLUSTER_FILE = "protein_cluster_4_10.parquet"
    AID_2_CIDS_FILE = "../aid_2_cids.csv"
    AID_2_TARGET_FILE = "../aid_2_target.csv"
    FINAL_FULL_PARQUET = "final_ligand_cluster.parquet"

    # READ THE protein_cluster_4_10.parquet file as a ParquetParser object
    schema = pa.schema(
        [
            ("Cluster number", pa.int64()),
            ("Represenatives", pa.string()),
            ("members", pa.string()),
            ("member_compound", pa.list_(pa.string())),
        ]
    )
    prot_parq = ParquetParser("../protein_clustered_data.parquet", schema)
    # head the protein_cluster_4_10.parquet file as a pandas dataframe
    prot_data = prot_parq.convert_to_pandas()
    logger.info(
        f"Number of rows in the protein_cluster_4_10.parquet file: {len(prot_data)}"
    )
    print(prot_data.head())
    # read the filtered_version_2_9_1_2025.csv file as a pandas dataframe and extract the from column
    filtered_data = pd.read_csv("./filtered_version_6_9_1_2025.csv")
    represenatives_id = filtered_data["From"]
    # find the represenatives_id that are not in the protein_cluster_4_10.parquet file
    missing_represenatives = [
        rep_id
        for rep_id in represenatives_id
        if rep_id not in prot_data["Represenatives"].to_list()
    ]
    logger.info(f"Number of missing represenatives: {len(missing_represenatives)}")
    # ingnore the missing represenatives and find the length of the member_compound column for each represenatives
    sum = 0
    data_no = []
    for rep_id in represenatives_id:
        if rep_id in missing_represenatives:
            continue
        member_compound = prot_data[prot_data["Represenatives"] == rep_id][
            "member_compound"
        ].to_list()[0]
        logger.info(
            f"Length of member_compound for represenatives {rep_id}: {len(member_compound)}"
        )
        data_no.append(tuple(member_compound))
        sum += len(member_compound)
    # print length of the data_no without duplication
    print(len(data_no))
    # remove the duplication in the data_no
    data_no = list(set(data_no))
    print(len(data_no))
    print(sum)
    # 1. Read the PubChem HDF5 file and extract the fasta sequences
    create_fasta_file(HDF5_FILE, FASTA_FILE)
    # 2. Cluster the fasta sequences using MMseqs2
    mmseqs_cluster(FASTA_FILE, outfile_name_suffix="../clusters", tmp_dir="../tmp")

    # 3. map the gene id to the PDB ids using ID mapping webservices
    PubChem_protein_ids = fasta_parser(FASTA_FILE)
    logger.info(f"Number of protein ids: {len(PubChem_protein_ids)}")
    BATCH_SIZE = 400
    mapped_geneids = []
    unique_protein_ids = list(set(PubChem_protein_ids))
    logger.info(f"Number of unique protein ids: {len(unique_protein_ids)}")
    for i in range(0, len(unique_protein_ids), BATCH_SIZE):
        batch = [
            protein_id.strip() for protein_id in unique_protein_ids[i : i + BATCH_SIZE]
        ]
        logger.info(f"Processing batch {i // BATCH_SIZE + 1}: {batch}")
        batch_results = map_genid_to_pdb(batch)
        for result in batch_results:
            flattened_result = {
                "Gene ID": result["Gene ID"],
                "UniprotID": result["UniprotID"],
                "Organism": result["Organism"],
                "Protein Name": result["Protein Name"],
                "Gene Name": result["Gene Name"],
                "PDB IDs": ",".join(result["PDB IDs"]),
                "AlphaFold IDs": ",".join(result["AlphaFold IDs"]),
            }
            mapped_geneids.append(flattened_result)

    mapped_geneids_df = pd.DataFrame(mapped_geneids)

    gene2pdb = mapped_geneids_df[["Gene ID", "PDB IDs"]]
    # remove rows with nan values in the PDB column
    gene2pdb = gene2pdb.dropna(subset=["PDB IDs"])
    gene2pdb.columns = ["From", "PDB"]
    # make the pdb column a list of strings
    gene2pdb["PDB"] = gene2pdb["PDB"].apply(lambda x: x.split(","))
    mapped_geneids_df.to_csv(RAW_ID_MAPPING_FILE, index=False)
    gene2pdb.to_csv(GENE_ID2PDB_ID_FILE, index=False)

    # 4. Validate the PDB ids and write the valid PDB ids to a new file
    # read the gene2pdb csv file and delete any rows and empty list in the PDB column
    gene2pdb = pd.read_csv(GENE_ID2PDB_ID_FILE)
    logger.info(f"Number of gene2pdb: {len(gene2pdb)}")
    pdb_validations(GENE_ID2PDB_ID_FILE)
    # read the valid_gene_id2pdb_id csv file and remove gene ids with empty PDB column
    valid_gene2pdb = pd.read_csv(VALID_PDBS_FILE)
    valid_gene2pdb = valid_gene2pdb[valid_gene2pdb["PDB"] != "[]"]
    valid_gene2pdb.to_csv(VALID_PDBS_FILE, index=False)
    # 5. Parse the MMseqs2 output to find representatives and members
    mmseqs_parser("../clusters_cluster.tsv", save=True)
    # 6. Write the clusters data into protein_cluster.parquet file
    schema = pa.schema(
        [
            ("Cluster number", pa.int64()),
            ("Represenatives", pa.string()),
            ("members", pa.string()),
            ("member_compound", pa.list_(pa.string())),
        ]
    )
    prot_parq = ParquetParser(PROTEIN_CLUSTER_FILE, schema)

    clean_cluster_files("../clusters_cluster_parsed.csv", VALID_PDBS_FILE)
    make_protein_cluster_file(prot_parq, "valid_clusters_cluster_parsed.csv")
    # 7. Parse the hdf5 file to get the aids and cids
    aids_2_target, aid_2_cid = hdf5_raw_file_parser(HDF5_FILE)
    # save the aids_2_target and aid_2_cid as csv files
    aids_2_target_df = pd.DataFrame(aids_2_target.items(), columns=["AID", "Gene_id"])
    aid_2_cid_df = pd.DataFrame(aid_2_cid.items(), columns=["AID", "CIDs"])
    aids_2_target_df.to_csv("../aid_2_target.csv", index=False)
    aid_2_cid_df.to_csv("../aid_2_cids.csv", index=False)

    # 8. Add the compounds to the parquet file with the protein clusters
    add_compounds(prot_parq, (aids_2_target, aid_2_cid), overwrite=False)
    prot_with_compounds = ParquetParser("new_protein_cluster_4_10.parquet", schema)
    # 9. Run tanimoto clustering on the compounds and add the ligand clusters to the parquet file
    # 10. Write the ligand clusters data into ligand_cluster.parquet file

    LIGAND_CLUST_FILE = "ligand_cluster.parquet"
    ligand_cluster_schema = pa.schema(
        [
            ("Protein Cluster number", pa.int64()),
            ("PDB/Gene ID", pa.string()),
            ("Compound ID", pa.string()),
            ("Smiles", pa.string()),
            ("Ligand Cluster number", pa.int64()),
        ]
    )
    lig_parq = ParquetParser(LIGAND_CLUST_FILE, ligand_cluster_schema)

    # read the lig_parq data as a pandas dataframe then convet it to a list
    lig_data = lig_parq.convert_to_pandas()
    # read the ligand_cluster_file as a ParquetParser object
    aid_2cids_df = pd.read_csv(AID_2_CIDS_FILE)
    aid_2target_df = pd.read_csv(AID_2_TARGET_FILE)

    make_ligand_cluster_file(
        prot_with_compounds,
        lig_parq,
        aid_2cids_df,
        aid_2target_df,
        HDF5_FILE,
        lig_data.values.tolist(),
    )
    final_parq = ParquetParser(FINAL_FULL_PARQUET, ligand_cluster_schema)
    add_smiles_cluster(lig_parq, 0.7, final_parq)

    final_data = final_parq.convert_to_pandas()
    # find if there is duplicate data in the final_parq file with full exaxt rows
    logger.info(f"Number of rows in the final_parq file: {len(final_data)}")
    # check the duplicates in the final_parq file interms of Compound ID and Smiles
    duplicates = final_data[final_data.duplicated(subset=["Compound ID"])]
    # save the unique data to a new parquet file named unique_final_ligand_cluster.parquet
    unique_data = final_data.drop_duplicates(subset=["Compound ID"])
    unique_parq = ParquetParser("unique_final_ligand_cluster.parquet", schema)
    unique_parq.write(unique_data)
    logger.info(
        f"Number of unique rows in the unique_final_parq file: {len(unique_data)}"
    )
    # how many unique ligand clusters are there
    logger.info(
        f"Number of unique ligand clusters: {len(final_data['Ligand Cluster number'].unique())}"
    )
    # how many unique protein clusters are there
    logger.info(
        f"Number of unique protein clusters: {len(final_data['Protein Cluster number'].unique())}"
    )
    # how many unique ligands are there
    logger.info(f"Number of unique ligands: {len(final_data['Compound ID'].unique())}")
    logger.info(f"Number of unique proteins: {len(final_data['PDB/Gene ID'].unique())}")
    # how many unique smiles are there
    logger.info(f"Number of unique smiles: {len(final_data['Smiles'].unique())}")
    # how many pdb/gene ids are there that are only numbers and not a mix of numbers and letters
    pdb_ids = final_data["PDB/Gene ID"]
    pdb_ids = pdb_ids[pdb_ids.str.isdigit()]
    logger.info(
        f"Number of pdb/gene ids that are only numbers: {len(pdb_ids.unique())}"
    )
    # get the number of the compound ids for these pdb_ids
    logger.info(
        f"Number of compound ids for these pdb/gene ids: {len(final_data[final_data['PDB/Gene ID'].isin(pdb_ids)]['Compound ID'].unique())}"
    )
