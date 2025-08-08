r"""Implementation of RMSD using OpenBabel."""

from typing import TextIO, override

import subprocess

from loguru import logger

from .base import RMSDBase


class obabelRMSD(RMSDBase):
    r"""Class to calculate RMSD using OpenBabel."""

    @override
    def calculate(
        self,
        firstonly: bool = True,
        save: bool = False,
        minimize: bool = False,
        output_filename: str | TextIO | None = None,
    ) -> float:
        """Calculate RMSD between reference and target file using OpenBabel,
        Args:
            firstonly :if True, calculate the RMSD for the first structure
                only in the reference file. Default is True.
            save :if True Write the RMSD to a txt file. Default is False.
            minimize :if True Compute minimum RMSD. Default is False.
            output_filename : Output file name if save is true. Default is None.
        """
        command = ["obrms"]
        if firstonly:
            command.append("-f")
        if minimize:
            command.append("-m")
        if not self.reference.file_path or not self.target.file_path:
            raise ValueError("File paths for reference or target cannot be None")
        command.extend([self.reference.file_path, self.target.file_path])
        try:
            with subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            ) as process:
                stdout, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(f"RMSD calculation completed for {self.target.file_id}")
                logger.info(f"Output:\n{stdout}")
            else:
                logger.error(f"RMSD calculation failed for {self.target.file_id}")
                logger.error(f"Error Output:\n{stderr}")
                raise subprocess.CalledProcessError(
                    process.returncode, " ".join(command)
                )
        except subprocess.CalledProcessError as e:
            logger.error(f"An error occurred during RMSD calculation: {str(e)}")
            raise e

        # Parse the RMSD from the output
        rmsd_result = stdout.strip()

        # If an output filename is provided, write the RMSD to the file
        if save:
            if output_filename is not None:
                if isinstance(output_filename, str):
                    with open(output_filename + ".txt", "w", encoding="utf-8") as file:
                        _ = file.write(f"RMSD: {rmsd_result}\n")
                else:
                    raise ValueError("Output filename must be a string if save is True")
            else:
                raise ValueError("Output filename is required if csv is True")
        # get the last value of the line
        return float(rmsd_result.split()[-1])
