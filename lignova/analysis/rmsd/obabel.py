# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Implementation of RMSD using OpenBabel."""

import os
from typing import override

from .base import RMSDBase


class obabelRMSD(RMSDBase):
    r"""Class to calculate RMSD using OpenBabel."""

    @override
    def calculate(
        self,
        firstonly: bool = True,
        save: bool = False,
        minimize: bool = False,
        output_filename: str | None = None,
    ) -> list[float]:
        r"""Calculate RMSD between reference and target file using OpenBabel,
        Args:
            firstonly :if True, calculate the RMSD for the first structure
                only in the reference file. Default is True.
            save :if True Write the RMSD to a txt file. Default is False.
            minimize :if True Compute minimum RMSD. Default is False.
            output_filename : Output file name if save is true. Default is None.
        Returns:
            list of RMSD values
        """
        if not self.reference.file_path or not self.target.file_path:
            raise ValueError("File paths for reference or target cannot be None.")
        if not os.path.exists(self.reference.file_path):
            raise FileNotFoundError(
                f"Reference file {self.reference.file_path} not found."
            )
        if not os.path.exists(self.target.file_path):
            raise FileNotFoundError(f"Target file {self.target.file_path} not found.")
        command = ["obrms"]
        if firstonly:
            command.append("-f")
        if minimize:
            command.append("-m")
        command.extend([self.reference.file_path, self.target.file_path])

        result = self._run_command(command)

        values = []
        for line in result.strip().splitlines():
            tokens = line.split()
            if tokens:
                values.append(float(tokens[-1]))

        if save:
            self._save_result(values, output_filename)

        return values
