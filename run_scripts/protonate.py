#!/usr/bin/env python
"""
run PDB2PQR to protonate proteins
"""

import argparse
import os
import sys

from loguru import logger

from lignova.preparation.pdb2pqr import PDB2PQR
from lignova.yaml.protonation_contig import ProtonationConfig


def protonate(
    pdb_file: str,
    output_filepath: str,
    config_path: str,
) -> None:
    """
    Protonate the protein using PDB2PQR.

    Args:
        pdb_file : The path to the pdb file
        output_filepath : The path to the output file
        config_path : The path to the protonation configuration file if not provided default config will be used
    """
    # check if the output directory exists if not create it
    output_dir = os.path.dirname(output_filepath)
    if not os.path.exists(output_dir):
        logger.warning(f"The directory {output_dir} does not exist.Creating it")
        os.makedirs(output_dir)
    protonation_config = ProtonationConfig(config_path, data_dict=None)
    pdb2pqr = PDB2PQR(
        pdb_file=pdb_file, outfile=output_filepath, config_obj=protonation_config
    )
    pdb2pqr.run()


def run_cli():
    """Command line interface for protonation using PDB2PQR."""
    parser = argparse.ArgumentParser(description="Protonate proteins using PDB2PQR.")
    parser.add_argument(
        "-p",
        "--pdb_file",
        type=str,
        required=True,
        help="Path to the input pdb file.",
    )
    parser.add_argument(
        "-o",
        "--output_filepath",
        type=str,
        required=True,
        help="Path to the output protonated pdb file.",
    )
    parser.add_argument(
        "-c",
        "--config_path",
        type=str,
        required=False,
        help="Path to the protonation configuration file.",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    protonate(
        pdb_file=args.pdb_file,
        output_filepath=args.output_filepath,
        config_path=args.config_path,
    )


if __name__ == "__main__":
    run_cli()
