"""Class for calculating RMSD using Spyrmsd tool."""

from typing import TextIO, override

import os
import subprocess

from loguru import logger

from lignova.docking.contexts import GlideContext
from lignova.structure.ligand import DockedLigand
from lignova.structure.protein import Protein

from ..utils import obabel_convert
from .base import RMSDBase


class spyrmsdRMSD(RMSDBase):
    r"""
    This class uses the Spyrmsd tool to calculate RMSD between a reference ligand
    and a target ligand. It can handle symmetry, hydrogens, and superimposition."""

    def __init__(self, target: DockedLigand, reference: Protein, context: GlideContext):
        r"""Initialize the spyrmsdRMSD class."""
        super().__init__(target, reference, context)
        state: int = self._validate_file_format()
        self._state: int = state
        if state != 0:
            self.fix_file_format()
        if not os.path.exists(self.reference.file_path):
            raise FileNotFoundError(
                f"Reference file {self.reference.file_path} not found."
            )
        if self.target.file_path is None or not os.path.exists(self.target.file_path):
            raise FileNotFoundError(f"Target file {self.target.file_path} not found.")

    def _validate_file_format(self) -> int:
        r"""Validate the file format for Spyrmsd."""
        if self.reference.file_ext in [
            ".sdf",
            ".pdb",
        ] and self.target.file_ext in [".sdf", ".pdb"]:
            return 0
        else:
            return 1

    @override
    def calculate(
        self,
        symmetry: bool = True,
        hydrogens: bool = False,
        superimpose: bool = False,
        save: bool = False,
        output_filename: str | TextIO | None = None,
    ) -> float:
        r"""Calculate RMSD between reference and target ligand using Spyrmsd,
        taking into account the symmetry of the molecules.
        **Note:** that this can only be used for small molecules and not for proteins

        Args:
            symmetry :
                Use symmetry information. Default is True.
            hydrogens :
                Include hydrogens in the calculation. Default is False.
            superimpose :
                Superimpose the molecules. Default is False. (i.e perform in-place RMSD)
        """
        if self._state != 0:
            self.fix_file_format(target_format="sdf")
            # make new ligand and reference objects
            self.reference: Protein = self.reference.__class__(
                str(self.reference.file_path).replace(self.reference.file_ext, ".sdf")
            )
            self.target: DockedLigand = self.target.__class__(
                str(self.target.file_path).replace(self.target.file_ext, ".sdf")
            )
        command = ["python", "-m", "spyrmsd"]
        if not symmetry:
            command.append("-n")
        if hydrogens:
            command.append("--h")
        if superimpose:
            command.append("-m")
        if self.reference.file_path is None or self.target.file_path is None:
            raise ValueError("File paths for reference or target cannot be None.")
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
        rmsd_value = float(stdout.strip())
        # If an output filename is provided, write the RMSD to the file
        if save:
            if output_filename is not None:
                if isinstance(output_filename, str):
                    with open(output_filename + ".txt", "w", encoding="utf-8") as file:
                        _ = file.write(f"RMSD: {rmsd_value}\n")
                else:
                    raise ValueError("Output filename must be a string if save is True")
            else:
                raise ValueError("Output filename is required if save is True")
        return rmsd_value

    def fix_file_format(
        self,
        target_format: str | None = "sdf",
    ) -> None:
        r"""Handles file format conversion for the reference and target files to work with Spyrmsd.
        This method ensures that the reference and target files are in the correct format
        for RMSD calculations using Spyrmsd. It can convert files to 'sdf' or 'pdb' format.
        Default is 'sdf'.
        Args:
            target_format : The desired file format for the reference and target files.
                Supported formats are 'sdf' and 'pdb'. Default is 'sdf'.
        Returns:
            None
        """
        # check if the target format is valid
        if target_format not in ["sdf", "pdb"]:
            raise ValueError("Target format must be either 'sdf' or 'pdb'.")
        # Convert the reference file if necessary
        if self.reference.file_ext not in [".sdf", ".pdb"]:
            output_file = self.reference.file_path.replace(
                self.reference.file_ext, f".{target_format}"
            )
            obabel_convert(self.reference.file_path, output_file)
        # Convert the target file if necessary
        if self.target.file_ext not in [".sdf", ".pdb"]:
            output_file = str(self.target.file_path).replace(
                self.target.file_ext, f".{target_format}"
            )
            if self.target.file_path is None:
                raise ValueError("Target file path cannot be None.")
            obabel_convert(self.target.file_path, output_file)
