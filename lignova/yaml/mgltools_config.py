# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Implementation for the yaml class to write configuration for MGLTools."""

import os
from typing import Any, override

from loguru import logger

from .config import YamlConfig


class MglToolsConfig(YamlConfig):
    """Class to handle YAML configuration files for MGLTools receptor preparation."""

    # Schema for prepare_receptor4.py configuration
    _REPAIR_ALLOWED = {
        "bonds_hydrogens",
        "bonds",
        "hydrogens",
        "checkhydrogens",
        "None",
    }
    _CLEANUP_ALLOWED = {
        "nphs",
        "lps",
        "waters",
        "nonstdres",
        "deleteAltB",
        "nphs_lps",
        "nphs_lps_waters",
        "nphs_lps_waters_nonstdres",
    }
    _ALLOWED_INPUT_EXTENSIONS: frozenset[str] = frozenset(
        {".pqr", ".pdb", ".pdbqt", ".mol2"}
    )
    _ALLOWED_BOOLEAN_PARAMS = {"preserve_charges"}

    # NOTE: prepare_receptor4.py uses short flags so these do not match key names
    _FLAG_MAP = {
        "receptor": "-r",
        "outfile": "-o",
        "repair": "-A",
        "cleanup": "-U",
        "preserve_charges": "-C",
    }

    @override
    def __init__(self, file_path: str, data_dict: dict[str, Any] | None = None) -> None:
        """Initialize with the path to the YAML file.

        Args:
            file_path : Path to the YAML configuration file.
            data_dict : Dictionary to create the YAML file if it doesn't exist.
        """
        if data_dict is None and not os.path.exists(file_path):
            logger.warning(
                "No configuration file found. Creating default mgltools config."
            )
            data_dict = self.declare_defaults()
        super().__init__(file_path, data_dict)
        self.validate()

    def declare_defaults(self) -> dict[str, dict[str, Any]]:
        """Declare default settings for mgltools config (grouped under mgltools)."""
        default_config: dict[str, dict[str, Any]] = {
            "mgltools": {
                "input_output": {
                    "receptor": None,
                    "outfile": None,
                },
                "receptor_perception": {
                    # bonds_hydrogens|bonds|hydrogens|checkhydrogens|None
                    "repair": "checkhydrogens",
                    # see _CLEANUP_ALLOWED; default strips nonpolar H, lone
                    # pairs, waters and non-standard residues (incl. ligands)
                    "cleanup": "nphs_lps_waters_nonstdres",
                    # NOTE: -C keeps pdb2pqr charges instead of recomputing
                    # Gasteiger. Only use for PQR input
                    "preserve_charges": True,
                },
            }
        }
        return default_config

    def validate(self) -> None:
        """Validate the current configuration against allowed values."""
        default = self.declare_defaults()
        mgltools_cfg = self.data_dict.setdefault("mgltools", {})
        for section, params in default["mgltools"].items():
            section_cfg = mgltools_cfg.setdefault(section, {})
            for key, value in params.items():
                if key not in section_cfg:
                    logger.warning(
                        "Missing parameter '{key}' in section '{section}'. Using default value: {value}",
                        key=key,
                        section=section,
                        value=value,
                    )
                    section_cfg[key] = value

        config = self.data_dict.get("mgltools", {})
        sections = ["input_output", "receptor_perception"]
        input_output = config.get("input_output", {})
        perception = config.get("receptor_perception", {})

        # Validate string parameters
        for section in sections:
            for key, value in config.get(section, {}).items():
                if key not in self._ALLOWED_BOOLEAN_PARAMS:
                    if not isinstance(value, str) and value is not None:
                        raise TypeError(
                            f"Parameter '{key}' in '{section}' must be a string or None."
                        )

        # Validate boolean parameters
        for section in sections:
            for key, value in config.get(section, {}).items():
                if key in self._ALLOWED_BOOLEAN_PARAMS and not isinstance(value, bool):
                    raise ValueError(
                        f"Parameter '{key}' in '{section}' must be a boolean."
                    )

        repair = perception.get("repair")
        if repair not in self._REPAIR_ALLOWED:
            raise ValueError(
                f"Invalid repair mode '{repair}'. Allowed: {sorted(self._REPAIR_ALLOWED)}"
            )

        cleanup = perception.get("cleanup")
        if cleanup not in self._CLEANUP_ALLOWED:
            raise ValueError(
                f"Invalid cleanup mode '{cleanup}'. Allowed: {sorted(self._CLEANUP_ALLOWED)}"
            )

        receptor = input_output.get("receptor")
        if receptor is None:
            logger.warning(
                "No receptor set. 'receptor' must be provided before running "
                "prepare_receptor4.py."
            )
        else:
            if not os.path.isfile(receptor):
                raise FileNotFoundError(
                    f"Specified receptor file '{receptor}' does not exist."
                )
            extension = os.path.splitext(receptor)[1].lower()
            if extension not in self._ALLOWED_INPUT_EXTENSIONS:
                raise ValueError(
                    f"Invalid receptor extension '{extension}'. "
                    f"Allowed: {sorted(self._ALLOWED_INPUT_EXTENSIONS)}"
                )

            if perception.get("preserve_charges") and extension != ".pqr":
                logger.warning(
                    "'preserve_charges' is set but input '{receptor}' is not a PQR; "
                    "PDB/PDBQT carry no charge column so Gasteiger charges will be "
                    "computed instead.",
                    receptor=receptor,
                )
                perception["preserve_charges"] = False

        if input_output.get("outfile") is None:
            logger.warning(
                "No output file set. 'outfile' must be provided before running "
                "prepare_receptor4.py."
            )
        self.write_config(self.data_dict)
        logger.info("Configuration validation passed.")

    @override
    def to_cli(
        self, data: dict[str, Any] | None = None, prefix: str = "-"
    ) -> list[str]:
        """Convert the configuration to prepare_receptor4.py arguments.

        Args:
            data : Optional dictionary to convert. If None, uses the current config.
            prefix : Prefix for command-line flags. Default is "-".
        """
        if data is None:
            data = self.data_dict
        config = data.get("mgltools", {})

        args: list[str] = []
        for section in ("input_output", "receptor_perception"):
            for key, value in config.get(section, {}).items():
                if value is None:
                    continue
                flag = self._FLAG_MAP.get(key)
                if flag is None:
                    logger.warning(
                        "No CLI flag mapped for parameter '{key}'; skipping.", key=key
                    )
                    continue
                if isinstance(value, bool):
                    if value:
                        args.append(flag)
                else:
                    args.extend([flag, str(value)])
        return args
