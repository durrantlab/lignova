# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Base class for RMSD calculators."""

import subprocess
from abc import ABC, abstractmethod

from loguru import logger

from ...structure.ligand import DockedLigand, Ligand


class RMSDBase(ABC):
    """
    Base class for RMSD calculators.
    """

    def __init__(
        self,
        target: DockedLigand,
        reference: Ligand,
    ):
        r"""Initialize RMSD class.

        Args:
            target:
                Docked ligand(s) Object that will be analyzed.
            reference:
                Reference ligand(s) in a Ligand object that will be used for comparison.
        """
        assert isinstance(target, DockedLigand), "Ligand must be a DockedLigand object."
        assert isinstance(reference, Ligand), "Reference must be a Ligand object."
        self.target: DockedLigand = target
        self.reference: Ligand = reference

    @abstractmethod
    def calculate(self) -> list[float]:
        r"""Calculate RMSD."""
        raise NotImplementedError()

    def _run_command(self, command: list[str]) -> str:
        r"""Run a subprocess command and return stdout.

        Args:
            command: Command and arguments to run.

        Returns:
            The stdout output as a string.
        """
        logger.info(f"Running: {' '.join(command)}")

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            logger.error(f"Command failed for {self.target.file_id}")
            logger.error(f"stderr:\n{result.stderr}")
            raise subprocess.CalledProcessError(
                result.returncode, " ".join(command), result.stdout, result.stderr
            )

        logger.info(f"Completed for {self.target.file_id}")
        return result.stdout

    def _save_result(
        self,
        values: list[float],
        output_filename: str | None,
    ) -> str:
        r"""Write RMSD values to a text file.

        Args:
            values: RMSD values to write.
            output_filename: Path without extension.

        Returns:
            Path to the written file.
        """
        if output_filename is None:
            raise ValueError("output_filename is required when save=True.")

        path = f"{output_filename}.txt"
        with open(path, "w", encoding="utf-8") as f:
            for v in values:
                f.write(f"{v}\n")

        logger.info(f"Saved RMSD values to {path}")
        return path
