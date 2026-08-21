# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Implementation for receptor preparation using Meeko pre-docking."""

import os
import subprocess

from loguru import logger

from lignova.yaml.meeko_config import MeekoConfig


class Meeko:
    """Class to handle receptor preparation using Meeko."""

    _SCRIPT = "mk_prepare_receptor.py"
    """The script name for Meeko's receptor preparation tool."""

    _ALLOWED_INPUT_EXTENSIONS = {".pqr", ".pdb", ".cif", ".mmcif"}
    """Allowed input file extensions for Meeko receptor preparation."""

    def __init__(
        self, input_file: str, output_basename: str | None, config_obj: MeekoConfig
    ) -> None:
        """Initialize Meeko with a given configuration object.

        Args:
            input_file : Path to the input receptor file (PQR from PDB2PQR).
            output_basename : Basename used for the generated receptor files. if None, defaults to same as input file without extension.
            config_obj : Configuration object for mk_prepare_receptor.

        """
        self.input_file = os.path.abspath(input_file)
        # check the input file exists
        if not os.path.exists(self.input_file):
            raise FileNotFoundError(f"Input file {input_file} does not exist.")
        extension = os.path.splitext(self.input_file)[1].lower()
        if extension not in self._ALLOWED_INPUT_EXTENSIONS:
            raise ValueError(
                f"Input file {input_file} must be one of {sorted(self._ALLOWED_INPUT_EXTENSIONS)}."
            )
        self.config = config_obj
        if output_basename is None:
            self.output_basename = os.path.splitext(self.input_file)[0]
        else:
            self.output_basename = (
                os.path.splitext(output_basename)[0]
                if output_basename.endswith(".pdbqt")
                else output_basename
            )
        self._write_paths()

    def _write_paths(self) -> None:
        """Write the input and output paths into the configuration and revalidate.
        Warns if the configuration already held a different input or basename.

        """
        input_cfg = self.config.data_dict.get("meeko", {}).get("input_output", {})
        reader = self.config._EXTENSION_READERS[
            os.path.splitext(self.input_file)[1].lower()
        ]

        existing_input = input_cfg.get(reader)
        if existing_input is not None and existing_input != self.input_file:
            logger.warning(
                "Config '{reader}' is '{existing_input}' but '{input_file}' was provided. Overwriting with provided value.",
                reader=reader,
                existing_input=existing_input,
                input_file=self.input_file,
            )
        existing_basename = input_cfg.get("output_basename")
        if existing_basename is not None and existing_basename != self.output_basename:
            logger.warning(
                "Config 'output_basename' is '{existing_basename}' but '{output_basename}' was provided. Overwriting with provided value.",
                existing_basename=existing_basename,
                output_basename=self.output_basename,
            )

        for key in self.config._INPUT_PARAMS:
            input_cfg[key] = self.input_file if key == reader else None
        input_cfg["output_basename"] = self.output_basename

        self.config.validate()

    @property
    def pdbqt_file(self) -> str:
        """Path to the rigid PDBQT that mk_prepare_receptor writes."""
        return f"{self.output_basename}.pdbqt"

    def run(self) -> str:
        """Run the Meeko receptor preparation process.

        Returns:
            Path to the generated PDBQT file.
        """
        meeko_config = self.config.to_cli()
        cmd = [self._SCRIPT]
        outdir = os.path.dirname(self.output_basename)
        if outdir and not os.path.exists(outdir):
            logger.debug(
                "Output directory {outdir} does not exist. Creating it.",
                outdir=outdir,
            )
            os.makedirs(outdir)
        if meeko_config:
            cmd.extend(meeko_config)
        try:
            logger.debug("Running Meeko with command: {cmd}", cmd=" ".join(cmd))
            cmd_str = " ".join(cmd)
            process = subprocess.run(
                cmd_str, capture_output=True, text=True, shell=True
            )
            if process.returncode == 0:
                logger.info(
                    "Meeko completed successfully. Output written to {output_basename}",
                    output_basename=self.output_basename,
                )
            else:
                logger.error("Meeko failed with error: {stderr}", stderr=process.stderr)
                raise RuntimeError(f"Meeko failed with error: {process.stderr}")

        except Exception as e:
            logger.error("An error occurred while running Meeko: {error}", error=str(e))
            raise e

        # NOTE: mk_prepare_receptor appends _rigid when flexible residues are requested
        rigid_file = f"{self.output_basename}_rigid.pdbqt"
        if os.path.exists(rigid_file):
            logger.info(
                "Rigid PDBQT receptor written to {rigid_file}", rigid_file=rigid_file
            )
            return rigid_file
        if not os.path.exists(self.pdbqt_file):
            raise RuntimeError(
                f"Meeko completed but output not found at {self.pdbqt_file}"
            )
        logger.info(
            "PDBQT receptor written to {pdbqt_file}", pdbqt_file=self.pdbqt_file
        )
        return self.pdbqt_file
