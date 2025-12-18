#!/usr/bin/env python
"""
run gypsum-dl to protonate proteins
"""

import argparse
import sys
import os
from lignova.preparation.gypsumdl import Gypsum
from lignova.yaml.ligprep_config import GypsumDLConfig
from loguru import logger

def conformer_gen(
    smi_file: str,
    output_folderpath: str,
    config_path: str,
) -> None:
    """
    Protonate the protein using Gypsum-dl.

    Args:
        smi_file : The path to the smiles file
        output_folderpath : The path to the output folder
        config_path : The path to the gypsum configuration file if not provided default config will be used

    Returns:
        The path to the protonated protein file
    """
    #check if the output directory exists if not create it
    output_dir = os.path.dirname(output_folderpath)
    if not os.path.exists(output_dir):
        logger.warning(f"The directory {output_dir} does not exist.Creating it")
        os.makedirs(output_dir)
    protonation_config = GypsumDLConfig(config_path,data_dict=None)
    dict_config = protonation_config.data_dict
    current_num_procs = dict_config["gypsum_dl"]["job_specs"].get("num_processors")
    #get the number of available processors on the machine and set it to the config
    available_procs = int(os.environ.get("GYPSUM_NUM_PROCS", "1"))
    if available_procs is not None and current_num_procs != available_procs:
        logger.warning(f"Setting number of processors to {available_procs}")
        protonation_config.data_dict["gypsum_dl"]["job_specs"][
            "num_processors"
        ] = available_procs
    gpys_obj = Gypsum(smiles_file=smi_file, outfolder=output_folderpath, config_obj=protonation_config)
    gpys_obj.run()
    
def run_cli():
    """Command line interface for protonation using PDB2PQR."""
    parser = argparse.ArgumentParser(
        description="Protonate proteins using PDB2PQR."
    )
    parser.add_argument("-s",
        "--smi_file",
        type=str,
        required=True,
        help="Path to the input smi file.",
    )
    parser.add_argument("-o",
        "--output_folderpath",
        type=str,
        required=True,
        help="Path to the output sdf file.",
    )
    parser.add_argument("-c",
        "--config_path",
        type=str,
        required=False,
        help="Path to the ligand conformer generation configuration file.",
    )
    
    #if no arguments are provided print help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    conformer_gen(
        smi_file=args.smi_file,
        output_folderpath=args.output_folderpath,
        config_path=args.config_path,
    )
    

if __name__ == "__main__":
    run_cli()
    