# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Implementation for protein preparation using PDB2PQR pre-docking."""

import os
import subprocess

from loguru import logger

from lignova.yaml.protonation_config import ProtonationConfig


class PDB2PQR:
    """Class to handle protein preparation using PDB2PQR."""

    def __init__(
        self, pdb_file: str, outfile: str, config_obj: ProtonationConfig
    ) -> None:
        """Initialize PDB2PQR with a given configuration file.

        Args:
            pdb_file (str): Path to the input PDB file.
            outfile (str): Path to the output PDB file.
            config_obj (ProtonationConfig): Configuration object for PDB2PQR.
        """
        self.pdb_file = pdb_file
        # check the input file exists
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"Input PDB file {pdb_file} does not exist.")
        self.config = config_obj
        self.outfile = outfile

    def run(self) -> None:
        """Run the PDB2PQR preparation process."""
        pdb2pqr_config = self.config.to_cli()
        cmd = ["pdb2pqr"]
        outdir = os.path.dirname(self.outfile)
        if outdir and not os.path.exists(outdir):
            logger.debug(f"Output directory {outdir} does not exist. Creating it.")
            os.makedirs(outdir)
        if pdb2pqr_config:
            cmd.extend(pdb2pqr_config)
        cmd.append(self.pdb_file)
        cmd.append(self.outfile)
        try:
            logger.debug(f"Running PDB2PQR with command: {' '.join(cmd)}")
            cmd_str = " ".join(cmd)
            process = subprocess.run(
                cmd_str, capture_output=True, text=True, shell=True
            )
            if process.returncode == 0:
                logger.info(
                    f"PDB2PQR completed successfully. Output written to {self.outfile}"
                )
            else:
                logger.error(f"PDB2PQR failed with error: {process.stderr}")
                raise RuntimeError(f"PDB2PQR failed with error: {process.stderr}")

        except Exception as e:
            logger.error(f"An error occurred while running PDB2PQR: {str(e)}")
            raise e
