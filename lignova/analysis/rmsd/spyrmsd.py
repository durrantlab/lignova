"""Class for calculating RMSD using Spyrmsd tool."""

import os
from typing import override

from lignova.structure.ligand import DockedLigand, Ligand

from ..utils import obabel_convert
from .base import RMSDBase


class spyrmsdRMSD(RMSDBase):
    r"""Calculate RMSD using spyrmsd, accounting for molecular symmetry.

    Note: This is only suitable for small molecules, not proteins.
    """

    def __init__(self, target: DockedLigand, reference: Ligand):
        r"""Initialize the spyrmsdRMSD class.

        Validates that input files exist and are in a supported format.
        Converts to SDF automatically if needed.
        """
        super().__init__(target, reference)

        if not os.path.exists(self.reference.file_path):
            raise FileNotFoundError(
                f"Reference file {self.reference.file_path} not found."
            )
        if not self.target.file_path or not os.path.exists(self.target.file_path):
            raise FileNotFoundError(f"Target file {self.target.file_path} not found.")

        self._ensure_compatible_format()

    def _ensure_compatible_format(self, target_format: str = "sdf") -> None:
        r"""Convert reference and target to a spyrmsd-compatible format if needed.

        Args:
            target_format: Desired format. Must be 'sdf' or 'pdb'.
        """
        supported = ["sdf", "pdb"]

        if self.reference.file_ext not in supported:
            ref_out = self.reference.file_path.replace(
                self.reference.file_ext, f".{target_format}"
            )
            obabel_convert(self.reference.file_path, ref_out)
            self.reference = self.reference.__class__(ref_out)

        if self.target.file_ext not in supported:
            tgt_out = str(self.target.file_path).replace(
                self.target.file_ext, f".{target_format}"
            )
            obabel_convert(self.target.file_path, tgt_out)
            self.target = self.target.__class__(tgt_out)

    @override
    def calculate(
        self,
        symmetry: bool = True,
        hydrogens: bool = False,
        superimpose: bool = False,
        save: bool = False,
        output_filename: str | None = None,
    ) -> list[float]:
        r"""Calculate RMSD between reference and target using spyrmsd.

        Args:
            symmetry: Use symmetry-corrected RMSD. Default is True.
            hydrogens: Include hydrogens. Default is False.
            superimpose: Superimpose before calculation. Default is False.
                Set to False for in-place RMSD of docked poses.
            save: If True, write the result to a text file. Default is False.
            output_filename: Output file path (without extension) if save
                is True.

        Returns:
            List of RMSD values.
        """
        command = ["python3", "-m", "spyrmsd"]
        if not symmetry:
            command.append("-n")
        if hydrogens:
            command.append("--hydrogens")
        if superimpose:
            command.append("-m")
        command.extend([self.reference.file_path, self.target.file_path])

        result = self._run_command(command)

        values = [
            float(line.strip()) for line in result.strip().splitlines() if line.strip()
        ]

        if save:
            self._save_result(values, output_filename)

        return values
