#!/usr/bin/env python
"""
Prepare proteins for PDB2PQR from a parquet file containing protein-ligand data.
Prepare ligand using gypsum-dl.
"""

import argparse
import os

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

#declare default sources 
clustered_parquet_source = "../../lignova_parquets/protein_clustered_data.parquet"
pdb_ligand_parquet_source = "../../lignova_parquets/final_ligand_cluster_0.7_Tc.parquet"

def get_pdb_ids_from_parquet(
    file_path: str, schema: pa.Schema | None = None
) -> list[str]:
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
    r"""
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

def extract_ligand_id(parquet_file: str, pdb_id:str) -> list[int]:
    r"""
    Extract ligand information from the parquet file.
    Args:
        parquet_file : The path to the parquet file
        pdb_id : The PDB id to extract ligands for
    Returns:
        The list of ligands.
    """
    if not os.path.exists(parquet_file):
        logger.error(f"The file {parquet_file} does not exist")
        raise FileNotFoundError(f"The file {parquet_file} does not exist")
    schema = pa.schema(
        [
            ("Cluster number", pa.int64()),
            ("Representatives", pa.string()),
            ("memberd", pa.string()),
            ("member_compound", list(pa.string())),
        ]
    )
    new_parquet = ParquetParser(parquet_file, schema)
    member_compound_data= new_parquet.filter_data(condition=lambda x:x == pdb_id, column="Representatives")["member_compound"][0]
    if len(member_compound_data) == 0:
        logger.warning(f"No ligands found for PDB ID {pdb_id} in the parquet file.")
        return []
    ligand_ids = [int(ligand) for ligand in member_compound_data]
    return ligand_ids

def make_smi_file_from_parquet(parquet:str,ligand_source:list[int], output_smi: str) -> None:
    r"""
    Create a smiles file from the ligand data in the parquet file.
    Args:
        parquet : The path to the parquet file containing the smiles for ligands
        ligand_source : The parquet file containing ligand information
        output_smi : The path to the output smiles file
    """
    if len(ligand_source) == 0:
        logger.warning("No ligand IDs provided to create smiles file.")
        return 
    if not os.path.exists(parquet):
        logger.error(f"The file {parquet} does not exist")
        raise FileNotFoundError(f"The file {parquet} does not exist")
    schema = pa.schema(
        [
            ("aid", pa.int32()),
            ("gene_id", pa.int32()),
            ("pubmed", pa.int32()),
            ("actives", pa.int32()),
            ("type", pa.string()),
            ("affinity", pa.float64()),
            ("smiles", pa.string()),
        ]
    )
    new_parquet = ParquetParser(parquet, schema)
    with open(output_smi, "w",encoding="utf-8") as smi_file:
        for ligand_id in ligand_source:
            smiles_data = new_parquet.filter_data(condition=lambda x:x == ligand_id, column="actives")
            if smiles_data["smiles"][0] is None or smiles_data["smiles"][0] == "" or smiles_data["pubmed"][0] is None:
                logger.warning(f"Missing smiles or pubmed data for ligand ID {ligand_id}. Skipping.")
                continue
            smi_file.write(f"{smiles_data['smiles'][0]} {smiles_data['actives'][0]}_{smiles_data['pubmed'][0]}\n")
    

def add_pdb_ligand_to_smi(parquet_file: str, pdb_id: str, output_smi: str) -> None:
    r"""
    Add ligands for a specific PDB ID from the parquet file to a smiles file.
    Args:
        parquet_file : The path to the parquet file containing ligand information
        pdb_id : The PDB id to extract ligands for
        output_smi : The path to the output smiles file
    """
    if not os.path.exists(parquet_file):
        logger.error(f"The file {parquet_file} does not exist")
        raise FileNotFoundError(f"The file {parquet_file} does not exist")
    schema= pa.schema(
        [
            ("Protein Cluster number", pa.int64()),
            ("PDB/Gene ID", pa.string()),
            ("Compound ID", pa.string()),
            ("Smiles", pa.string()),
            ("Ligand Cluster number", pa.int64()),
        ]
    )
    new_parquet = ParquetParser(parquet_file, schema)
    #get the protein_cluster number for the given pdb_id
    filtered_data = new_parquet.filter_data(condition=lambda x: x == pdb_id, column="PDB/Gene ID")
    if filtered_data.num_rows == 0:
        logger.warning(f"No entries found for PDB ID {pdb_id} in the parquet file.")
        raise ValueError(f"No entries found for PDB ID {pdb_id} in the parquet file.")
    protein_cluster_number = filtered_data["Protein Cluster number"][0]
    #get all ligands for the protein_cluster_number that only have a length of 3
    ligands_data = new_parquet.filter_data(
        condition=lambda x: x == protein_cluster_number,
        column="Protein Cluster number"
    )
    with open(output_smi, "a",encoding="utf-8") as smi_file:
        for i in range(ligands_data.num_rows):
            compound_id = ligands_data["Compound ID"][i]
            smiles = ligands_data["Smiles"][i]
            pdb_entry = ligands_data["PDB/Gene ID"][i]
            if len(compound_id) == 3:
                smi_file.write(f"{smiles} {compound_id}_{pdb_entry}\n")
    
def prep_proteins(
    pdb_file: str,
    output_dir: str,
) -> str:
    r"""
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

def run_cli_prot(
    input_dir:str,
    pdb_parq_source: str,
    output_dir: str,
    num_proteins: list[int] | None = None,
) -> list[str]:
    r"""
    CLI to prepare proteins pre-PDB2PQR from a parquet file 
    Args:
        input_dir : The input directory containing raw pdb files
        pdb_parq_source : The path to the clustered parquet file containing PDB IDs OR name of a single pdb file
        output_dir : The base output directory
        num_proteins : range of proteins to prepare
    Returns:
        The list of paths to the prepared protein files.
    """
    logger.info(f"Reading PDB IDs from parquet: {pdb_parq_source}")
    #check that parquet file has a valid extension i.e parquet or .pq or pdb
    if not pdb_parq_source.endswith((".parquet", ".pq", ".pdb")):
        raise ValueError("pdb_parq_source must be a parquet file with .parquet or .pq extension or a single pdb id with .pdb extension.")
    if pdb_parq_source.endswith(".parquet") or pdb_parq_source.endswith(".pq"):
        pdb_ids = get_pdb_ids_from_parquet(pdb_parq_source)
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
        if not os.path.exists(os.path.join(input_dir, pdb_parq_source)):
            raise FileNotFoundError(f"The file {pdb_parq_source} does not exist.")
        #run the preparation on the single pdb file directly
        pdb_output_dir = os.path.join(output_dir, os.path.basename(pdb_parq_source).split(".pdb")[0])
        prepped = prep_proteins(
            pdb_file=os.path.join(input_dir, pdb_parq_source),
            output_dir=pdb_output_dir)
        logger.info(
            f"Cleaned protein for {os.path.basename(pdb_parq_source)} saved to {prepped}"
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


def make_mangable_files(file_list: list[str],split:int=250) -> list[str]:
    r""" Split a large .smi files into smaller manageable files.
    Args:
        file_list : The list of .smi files to be split
        split : The number of lines per split file
    """
    new_paths = []
    for file in file_list:
        if not os.path.exists(file):
            logger.error(f"The file {file} does not exist")
            raise FileNotFoundError(f"The file {file} does not exist")
        with open(file, "r",encoding="utf-8") as f:
            lines = f.readlines()
        #ensure that the file has more lines than the split size
        if len(lines) <= split:
            logger.info(f"The file {file} has {len(lines)} lines which is less than or equal to the split size of {split}. No splitting needed.")
            new_paths.append(file)
            continue
        base_name = os.path.splitext(file)[0]
        for i in range(0, len(lines), split):
            chunk = lines[i:i + split]
            chunk_file = f"{base_name}_part{i//split + 1}.smi"
            with open(chunk_file, "w",encoding="utf-8") as chunk_f:
                chunk_f.writelines(chunk)
            logger.info(f"Created chunk file: {chunk_file}")
            new_paths.append(chunk_file)
    return new_paths

def run_cli_ligand(
    pdb_id_source: str,
    output_dir: str,
    ligand_source: str,
    output_smi: str,
    num_proteins: list[int] | None = None,
    add_pdb_ligand: bool = True,
    chunk_size: int | None = 250,
) -> list[str]:
    r"""
    CLI to prepare ligands pre-gypsum from a parquet file 
    Args:
        pdb_id_source : The path to the clustered parquet file containing PDB IDs OR name of a single pdb file
        output_dir : The base output directory
        ligand_source : The parquet file containing ligand information
        output_smi : The suffix for the output smiles file 
        num_proteins : range of proteins to prepare. If None, prepare all proteins.
        add_pdb_ligand : Whether to add ligands from the pdb_ligand_parquet_source parquet files. Defaults to True.
        chunk_size : The number of lines per split file when making manageable smi files. Default is 250. if chunk_size is None no splitting is done.
    Returns:
        The list of paths to the prepared ligand smiles files.
    """
    outputpaths = []
    if output_smi.endswith(".smi"):
        output_smi = output_smi.replace(".smi","")
    if not os.path.exists(ligand_source):
        logger.error(f"The file {ligand_source} does not exist")
        raise FileNotFoundError(f"The file {ligand_source} does not exist")
    if add_pdb_ligand:
        #validate that the default files exists
        if not os.path.exists(pdb_ligand_parquet_source) or not os.path.isfile(clustered_parquet_source):
            logger.error(f"The file {pdb_ligand_parquet_source} or {clustered_parquet_source} does not exist")
            raise FileNotFoundError(f"The file {pdb_ligand_parquet_source} or {clustered_parquet_source} does not exist")
    #check if pdb_id_source is a parquet file or a single pdb file
    if pdb_id_source.endswith(".parquet") or pdb_id_source.endswith(".pq"):
        pdb_ids = get_pdb_ids_from_parquet(pdb_id_source)
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
            logger.info(f"Preparing ligands for proteins in the range: {num_proteins[0]} to {num_proteins[1]}")
        if not pdb_ids:
            raise RuntimeError("No PDB IDs found in parquet file.")
    else:
        #single pdb file provided
        pdb_id = os.path.basename(pdb_id_source).split(".pdb")[0]
        pdb_ids = [pdb_id]
    os.makedirs(output_dir, exist_ok=True)
    for pdb_id in pdb_ids:
        output_dir=os.path.join(output_dir, pdb_id)
        os.makedirs(output_dir, exist_ok=True)
        pathfile=os.path.join(output_dir, f"{pdb_id}_{output_smi}.smi")
        if os.path.exists(pathfile):
            logger.info(f"Smiles file for {pdb_id} already exists at {pathfile}. Skipping generation.")
            outputpaths.append(pathfile)
            continue
        if add_pdb_ligand:
            ligand_ids = extract_ligand_id(clustered_parquet_source, pdb_id)
            make_smi_file_from_parquet(ligand_source, ligand_ids,pathfile)
            if os.path.exists(pathfile):
                outputpaths.append(pathfile)
            else:
                logger.warning(f"No ligands found for PDB ID {pdb_id} in the ligand source parquet file.")
            add_pdb_ligand_to_smi(pdb_ligand_parquet_source, pdb_id, pathfile)
        else:
            make_smi_file_from_parquet(ligand_source, extract_ligand_id(ligand_source, pdb_id),pathfile)
            if os.path.exists(pathfile):
                outputpaths.append(pathfile)
            else:
                logger.warning(f"No ligands found for PDB ID {pdb_id} in the ligand source parquet file.")
    if chunk_size is not None:
        outputpaths = make_mangable_files(outputpaths,split=chunk_size)
    return outputpaths

    
def build_arg_parser() -> argparse.ArgumentParser:
    r""" Build the argument parser for the CLI."""
    parser = argparse.ArgumentParser(
        description="Prepare proteins for running PBQ2PQR from a parquet file and prepare ligands using gypsum-dl."
    )
    parser.add_argument("-m","--mode", choices=["protein","ligand"], required=True, help="Mode to run the script in: 'protein' to prepare proteins, 'ligand' to prepare ligands.")
    
    parser.add_argument("-i", "--input-dir", required=False, help="Input directory for raw PDB files.")
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
        required=True,
        help="Start index (inclusive) of proteins to prepare (default: 0).",
    )
    parser.add_argument(
        "-e",
        "--end-index",
        type=int,
        default=None,
        required=True,
        help="End index (exclusive) of proteins to prepare. "
             "If not provided, computed as start-index + num-proteins, "
             "or all proteins if both are omitted.",)
    parser.add_argument(
        "-ch","--chunk-size",
        type=int,
        required=False,
        help="Number of lines per split smiles file when making manageable files (default: 250). If not provided, no splitting is done.",
    )
    parser.add_argument(
            "-l",
            "--ligand-source",
            required=False,
            help="Path to the parquet file containing ligand information.",
        )
    parser.add_argument(
            "-osmi",
            "--output-smi",
            required=False,
            help="Suffix for the output smiles file.",
        )
    parser.add_argument("--add-pdb-ligand", 
                            action="store_true", 
                            default=True,
                            help="Whether to add ligands from the pdb_ligand_parquet_source parquet files.(Default:True)")
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    #check if mode is protein and if -p is provided as .parquet or .pq file that the input-dir is also provided
    if args.mode == "protein":
        if args.parquet.endswith((".parquet", ".pq")) and args.input_dir is None:
            parser.error("The --input-dir argument is required when --parquet is a parquet file in protein mode.")
        cleaned_file=run_cli_prot(
            input_dir=args.input_dir,
            pdb_parq_source=args.parquet,
            output_dir=args.output_dir,
            num_proteins=[args.start_index, args.end_index])
        
        output_file=os.path.join(args.output_dir, f"cleaned_pdblist_{args.start_index}_{args.end_index}.txt")
        with open(output_file, "w",encoding="utf-8") as f:
            for item in cleaned_file:
                f.write(f"{item}\n")
        logger.info(f"List of cleaned PDB files saved to {output_file}")
    else:
        if args.ligand_source is None or args.output_smi is None:
            parser.error("The --ligand-source and --output-smi arguments are required in ligand mode.")
        if args.chunk_size is None:
            smiles_files=run_cli_ligand(
                pdb_id_source=args.parquet,
                output_dir=args.output_dir,
                ligand_source=args.ligand_source,
                output_smi=args.output_smi,
                num_proteins=[args.start_index, args.end_index],
                add_pdb_ligand=args.add_pdb_ligand)
        else:
            smiles_files=run_cli_ligand(
                pdb_id_source=args.parquet,
                output_dir=args.output_dir,
                ligand_source=args.ligand_source,
                output_smi=args.output_smi,
                num_proteins=[args.start_index, args.end_index],
                add_pdb_ligand=args.add_pdb_ligand,
                chunk_size=args.chunk_size)
        output_file=os.path.join(args.output_dir, f"smiles_filelist_{args.start_index}_{args.end_index}.txt")
        with open(output_file, "w",encoding="utf-8") as f:
            for item in smiles_files:
                f.write(f"{item}\n")
        

# Print the help message when running the script without arguments
if not any(arg in os.sys.argv for arg in ("-p", "--parquet")):
    build_arg_parser().print_help()
    os.sys.exit(1)


if __name__ == "__main__":
    main()
