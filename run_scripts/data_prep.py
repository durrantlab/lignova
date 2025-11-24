#!/usr/bin/env python
"""
Prepare proteins for PDB2PQR from a parquet file containing protein-ligand data.
"""

import argparse
import os

import pandas as pd
import pyarrow as pa
from loguru import logger

from lignova.hdf5.parquet import ParquetParser
from lignova.structure.editing import convert_cif2pdb, write_mda_universe

from lignova.structure.protein import Protein
from lignova.structure.utils import (
    separate_protein_ligand,
    validate_ligands,
    validate_pdb,
)


def get_pdb_ids_from_parquet(
    file_path: str, schema: pa.Schema | None = None
) -> list:
    r"""
    Get the pdb ids from the parquet file

    Args:
        file_path : The path to the parquet file

    Returns:
        The list of pdb ids.
    """
    if not os.path.exists(file_path):
        logger.error(f"The file {file_path} does not exist")
        raise FileNotFoundError(f"The file {file_path} does not exist")
    if schema is None:
        schema = pa.schema(
            [
                ("Protein Cluster number", pa.int64()),
                ("PDB/Gene ID", pa.string()),
                ("Compound ID", pa.string()),
                ("Smiles", pa.string()),
                ("Ligand Cluster number", pa.int64()),
            ]
        )
    new_parquet = ParquetParser(file_path, schema)
    raw_prot_ids = new_parquet.convert_to_pandas()["PDB/Gene ID"].unique()
    pdb_ids = [
        raw_prot_ids[i]
        for i in range(len(raw_prot_ids))
        if any(char.isalpha() for char in raw_prot_ids[i])
    ]
    return pdb_ids


def get_pdb_coordinates(pdb_id: str, work_dir: str) -> None:
    """
    Download/prepare PDB coordinates if they pass validation.

    Args:
        pdb_id : The PDB id to be downloaded
        work_dir : The working directory where the PDB file will be downloaded.
    """
    current_dir = os.getcwd()
    protein = Protein()
    if not os.path.exists(work_dir):
        logger.info("Output_Dir not found,Creating it in working directory")
        os.mkdir(os.path.join(current_dir, work_dir))
    if (
        not os.path.exists(os.path.join(work_dir, pdb_id.lower() + ".pdb"))
        and validate_pdb(pdb_id)
        and validate_ligands(pdb_id)
    ):
        logger.info(f"Downloading PDB file for {pdb_id}")
        file_ext = (
            "pdb" if protein.get_pdb_from_rcsb(pdb_id).startswith("HEADER") else "cif"
        )
        protein.load(
            pdb_id=pdb_id,
            write=True,
            write_path=os.path.join(work_dir, pdb_id.lower() + "." + file_ext),
        )
        if file_ext == "cif":
            logger.info(f"Converting {pdb_id} to pdb format")
            convert_cif2pdb(
                os.path.join(work_dir, pdb_id.lower() + ".cif"),
                os.path.join(work_dir, pdb_id.lower() + ".pdb"),
            )
    elif os.path.exists(os.path.join(work_dir, pdb_id.lower() + ".pdb")):
        logger.info(f"{pdb_id} already exists in the directory")
    else:
        logger.warning(f"{pdb_id} failed validation test")

def prep_proteins(
    pdb_file: str,
    output_dir: str,
) -> str:
    """
    Prepare the protein for PDB2PQR
    Args:
        pdb_file : The path to the pdb file
        output_dir : The output directory
    Returns:
        The path to the prepared protein file
    """
    temp_prot = Protein(file_path=pdb_file)
    if not os.path.exists(output_dir):
        logger.warning(f"The directory {output_dir} does not exist.Creating it")
        os.makedirs(output_dir)
    
    if not os.path.exists(os.path.join(output_dir, temp_prot.file_id + "_cleaned.pdb")):
        protein, ligand = separate_protein_ligand(
            pdb_file,
            keep_het_chain="A",
        )
        write_mda_universe(
            protein, os.path.join(output_dir, f"{temp_prot.file_id}_cleaned.pdb")
        )
        write_mda_universe(
            ligand, os.path.join(output_dir, f"{temp_prot.file_id}_ligand.pdb")
        )
    else:
        logger.info(f"Cleaned protein file for {temp_prot.file_id}_cleaned.pdb already exists.")
    return os.path.join(output_dir, f"{temp_prot.file_id}_cleaned.pdb")
        
def run_cli(
    input_dir:str,
    parquet_path: str,
    output_dir: str,
    num_proteins: list[int] | None = None,
) -> list[str]:
    """
    CLI to prepare proteins pre-PDB2PQR from a parquet file.
    Args:
        input_dir : The input directory containing raw pdb files
        parquet_path : The path to the clustered parquet file containing PDB IDs OR name of a single pdb file
        output_dir : The base output directory
        num_proteins : range of proteins to prepare
    Returns:
        The list of paths to the prepared protein files.
    """
    logger.info(f"Reading PDB IDs from parquet: {parquet_path}")
    #check that parquet file has a valid extension i.e parquet or .pq or pdb
    if not parquet_path.endswith((".parquet", ".pq", ".pdb")):
        raise ValueError("parquet_path must be a parquet file with .parquet or .pq extension or a single pdb id with .pdb extension.")
    if parquet_path.endswith(".parquet") or parquet_path.endswith(".pq"):
        pdb_ids = get_pdb_ids_from_parquet(parquet_path)
        if num_proteins is not None:
            #validate num_proteins is a list of two integers
            if not (isinstance(num_proteins, list) and len(num_proteins) == 2 and all(isinstance(i, int) for i in num_proteins)):
                raise ValueError("num_proteins must be a list of two integers representing the range of proteins to prepare.")
            #check if the end index is greater than the length of pdb_ids
            if num_proteins[1] > len(pdb_ids):
                logger.warning(f"End index {num_proteins[1]} is greater than the number of available PDB IDs {len(pdb_ids)}. Adjusting to {len(pdb_ids)}.")
                num_proteins[1] = len(pdb_ids)
            
            #knowing that num_proteins is a list define the range of proteins to prepare
            pdb_ids = pdb_ids[num_proteins[0] : num_proteins[1]]
            # Example: if num_proteins = [0, 50], prepare first 50 proteins
            logger.info(f"Preparing proteins in the range: {num_proteins[0]} to {num_proteins[1]}")
            
        if not pdb_ids:
            raise RuntimeError("No PDB IDs found in parquet file.")
        #Ensure output directory exists
    else:
        #single pdb file provided
        if not os.path.exists(os.path.join(input_dir, parquet_path)):
            raise FileNotFoundError(f"The file {parquet_path} does not exist.")
        #run the preparation on the single pdb file directly
        pdb_output_dir = os.path.join(output_dir, os.path.basename(parquet_path).split(".pdb")[0])
        prepped = prep_proteins(
            pdb_file=os.path.join(input_dir, parquet_path),
            output_dir=pdb_output_dir)
        logger.info(
            f"Cleaned protein for {os.path.basename(parquet_path)} saved to {prepped}"
        )
        return [prepped]
    os.makedirs(output_dir, exist_ok=True)
    
    prepared_prot_path=[]
    for pdb_id in pdb_ids:
        
        logger.info(f"Processing PDB: {pdb_id}")
        get_pdb_coordinates(pdb_id, input_dir)

        pdb_filepath = os.path.join(input_dir, f"{pdb_id.lower()}.pdb")
        if not os.path.exists(pdb_filepath):
            logger.warning(f"Skipping {pdb_id}: PDB file not present after download.")
            continue
        pdb_output_dir = os.path.join(output_dir, pdb_id)
        os.makedirs(pdb_output_dir, exist_ok=True)
        # 3. prepare protein
        try:
            prepped = prep_proteins(
                pdb_file=pdb_filepath,
                output_dir=pdb_output_dir)
            prepared_prot_path.append(prepped)
            logger.info(
                f"Cleaned protein for {pdb_id} saved to {prepped}"
            )
        except Exception as exc:
            logger.error(f"Failed to prepare protein {pdb_id}: {exc}")
    return prepared_prot_path


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare proteins for running PBQ2PQR from a parquet file."
    )
    parser.add_argument("-i", "--input-dir", required=True, help="Input directory for raw PDB files.")
    parser.add_argument(
        "-p",
        "--parquet",
        required=True,
        help="Path to the parquet file containing PDB IDs. Can also be a single PDB file name.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Base output directory.",
    )
    parser.add_argument(
        "-s",
        "--start-index",
        type=int,
        default=0,
        help="Start index (inclusive) of proteins to prepare (default: 0).",
    )
    parser.add_argument(
        "-e",
        "--end-index",
        type=int,
        default=None,
        help="End index (exclusive) of proteins to prepare. "
             "If not provided, computed as start-index + num-proteins, "
             "or all proteins if both are omitted.",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    
    cleaned_file=run_cli(
        input_dir=args.input_dir,
        parquet_path=args.parquet,
        output_dir=args.output_dir,
        num_proteins=[args.start_index, args.end_index])
    
    output_file=os.path.join(args.output_dir, f"cleaned_pdblist_{args.start_index}_{args.end_index}.txt")
    with open(output_file, "w") as f:
        for item in cleaned_file:
            f.write(f"{item}\n")
    logger.info(f"List of cleaned PDB files saved to {output_file}")

# Print the help message when running the script without arguments
if not any(arg in os.sys.argv for arg in ("-p", "--parquet")):
    build_arg_parser().print_help()
    os.sys.exit(1)


if __name__ == "__main__":
    main()
