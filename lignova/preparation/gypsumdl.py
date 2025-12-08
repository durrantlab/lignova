r"""Implementation for ligand preparation using Gypsum-dl pre-docking."""

import glob
import os
import subprocess

from loguru import logger

from lignova.yaml.ligprep_config import GypsumDLConfig


class Gypsum:
    """Class to handle protein preparation using Gypsum-dl."""

    def __init__(
        self,
        smiles_file: str,
        outfolder: str,
        config_obj: GypsumDLConfig,
    ) -> None:
        """Initialize gypsum-dl with a given configuration file.

        Args:
            smiles_file (str): Path to the input PDB file.
            outfolder (str): Path to the output folder.
            config_obj (GypsumDLConfig): Configuration object for gypsum-dl.
            outfile (str): name of the output file
        """
        if not smiles_file.endswith(".smi"):
            raise ValueError("Input file must have a .smi extension.")
        self.smiles_file = smiles_file
        # check the input file exists
        if not os.path.exists(smiles_file):
            raise FileNotFoundError(f"Input smiles file {smiles_file} does not exist.")
        self.config = config_obj
        self.outfolder = outfolder

    def run(self) -> None:
        """Run the gypsum-dl preparation process."""
        gypsum_config = self.config.to_cli()
        job_manager = self.config.data_dict["gypsum_dl"]["job_specs"].get("job_manager")
        ntasks = self.config.data_dict["gypsum_dl"]["job_specs"].get(
            "tasks_per_processor"
        )
        gypsum_config = [
            arg for arg in gypsum_config if not arg.startswith("--tasks_per_processor")
        ]
        cmd = ["gypsum-dl"]
        if not os.path.exists(self.outfolder):
            logger.debug(
                f"Output directory {self.outfolder} does not exist. Creating it."
            )
            os.makedirs(self.outfolder)
        if job_manager == "mpi":
            cmd = [
                "mpirun",
                "-n",
                str(ntasks),
                "python",
                "-m",
                "mpi4py",
                "run_gypsum_dl.py",
                "--source",
                self.smiles_file,
                "--output_folder",
                self.outfolder,
            ]
            cmd.extend(gypsum_config)
        else:
            cmd.extend(
                ["--source", self.smiles_file, "--output_folder", self.outfolder]
            )
            cmd.extend(gypsum_config)
        logger.debug(f"Running gypsum-dl with command: {' '.join(cmd)}")
        cmd_str = " ".join(cmd)
        try:
            process = subprocess.run(
                cmd_str,
                capture_output=True,
                text=True,
                shell=True,
            )
            if process.returncode == 0:
                logger.info(
                    f"gypsum-dl completed successfully. Output written to {self.outfolder}"
                )
            else:
                logger.error(
                    "gypsum-dl failed with return code "
                    f"{process.returncode}. Stderr:\n{process.stderr}"
                )
                raise RuntimeError(f"gypsum-dl failed with error: {process.stderr}")

        except Exception as e:
            logger.error(f"An error occurred while running gypsum-dl: {str(e)}")
            raise e
