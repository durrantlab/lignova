# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Implementation for the yaml class to write configuration for PDB2QR."""

import os
from typing import Any, override

from loguru import logger

from .config import YamlConfig


class ProtonationConfig(YamlConfig):
    r"""Class to handle YAML configuration files for protonation contigs."""

    # ---- Schema/constraints ----
    _FF_ALLOWED = {"AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON"}
    _PROPKA_REFERENCE_ALLOWED = {"neutral", "low-pH"}
    _LOG_LEVELS = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
    _TITRATION_METHODS = {None, "propka"}
    _MUTATOR_ALLOWED = {None, "alignment", "scwrl", "jackal"}
    _ALLOWED_BOOLEAN_PARAMS = {
        "clean",
        "nodebump",
        "noopt",
        "keep-chain",
        "neutraln",
        "neutralc",
        "assign-only",
        "whitespace",
        "drop-water",
        "include-header",
        "display-coupled-residues",
        "reuse-ligand-mol2-files",
        "keep-protons",
        "quiet",
        "protonate-all",
    }
    _ALLOWED_FLOAT_PARAMS = {"with-ph", "pH"}
    _ALLOWED_LIST_PARAMS = {"window", "grid", "file"}

    @override
    def __init__(self, file_path: str, data_dict: dict[str, Any] | None = None) -> None:
        r"""Initialize with the path to the YAML file.
        Args:
            file_path (str): Path to the YAML configuration file.
            data_dict (dict| None) : Dictionary to create the YAML file if it doesn't exist.
        """
        if data_dict is None and not os.path.exists(file_path):
            logger.warning(
                "No configuration file found. Creating default protonation config."
            )
            data_dict = self.declare_defaults()
        super().__init__(file_path, data_dict)
        self.validate()

    def declare_defaults(self) -> dict[str, dict[str, Any]]:
        r"""Declare default settings for protonation config (grouped under pdb2pqr)."""
        default_config: dict[str, dict[str, Any]] = {
            "pdb2pqr": {
                "general": {
                    "ff": "PARSE",  # AMBER|CHARMM|PARSE|TYL06|PEOEPB|SWANSON
                    "userff": None,
                    "clean": False,
                    "nodebump": False,
                    "noopt": False,
                    "keep-chain": True,
                    "assign-only": False,
                    "ffout": None,
                    "usernames": None,
                    "apbs-input": None,
                    "pdb-output": None,
                    "ligand": None,
                    "whitespace": False,
                    # NOTE: PARSE-only neutrality options — update if ff changes
                    "neutraln": True,
                    "neutralc": True,
                    "drop-water": True,
                    "include-header": True,
                },
                "pka": {
                    "titration-state-method": "propka",  # "propka" or None
                    "with-ph": 7.0,
                },
                "propka": {
                    "reference": "neutral",  # neutral|low-pH
                    "chain": None,  # e.g., "A" (use " " for no-ID chains)
                    "i": None,  # e.g., "A:10,A:11"
                    "thermophile": None,
                    "alignment": None,
                    "mutation": None,  # e.g., "N25R/N181D"
                    # NOTE: Optional include parameters path if you want it in defaults
                    # "parameters": "/path/to/propka.cfg",
                    "log-level": "INFO",  # DEBUG|INFO|WARNING|ERROR|CRITICAL
                    "pH": 7.0,
                    "window": [0.0, 14.0, 1.0],
                    "grid": [0.0, 14.0, 0.1],
                    "mutator": None,  # alignment|scwrl|jackal
                    "mutator-option": None,  # e.g., 'type="side-chain"'
                    "display-coupled-residues": False,
                    "reuse-ligand-mol2-files": False,
                    "keep-protons": False,
                    "quiet": False,
                    "protonate-all": True,
                },
            }
        }
        return default_config

    def validate(self) -> None:
        r"""Validate the current configuration against allowed values."""
        # compare values in data_dict to default allowed values if any keys are missing log a warning and use default
        default = self.declare_defaults()
        pdb2pqr_cfg = self.data_dict.setdefault("pdb2pqr", {})
        for section, params in default["pdb2pqr"].items():
            section_cfg = pdb2pqr_cfg.setdefault(section, {})
            for key, value in params.items():
                if key not in section_cfg:
                    logger.warning(
                        f"Missing parameter '{key}' in section '{section}'. Using default value: {value}"
                    )
                    section_cfg[key] = value

        self.write_config(self.data_dict)
        config = self.data_dict.get("pdb2pqr", {})
        general = config.get("general", {})
        pka = config.get("pka", {})
        propka = config.get("propka", {})
        # Validate string parameters
        for section in ["general", "propka"]:
            for key, value in config.get(section, {}).items():
                if (
                    key not in self._ALLOWED_BOOLEAN_PARAMS
                    and key not in self._ALLOWED_FLOAT_PARAMS
                    and key not in self._ALLOWED_LIST_PARAMS
                ):
                    if not isinstance(value, str) and value is not None:
                        raise TypeError(
                            f"Parameter '{key}' in section '{section}' must be a string or None."
                        )
        # Validate boolean parameters
        for section in ["general", "propka"]:
            for key, value in config.get(section, {}).items():
                if key in self._ALLOWED_BOOLEAN_PARAMS and not isinstance(value, bool):
                    raise ValueError(
                        f"Parameter '{key}' in section '{section}' must be a boolean."
                    )
        # Validate float parameters
        for section in ["pka", "propka"]:
            for key, value in config.get(section, {}).items():
                if key in self._ALLOWED_FLOAT_PARAMS and not isinstance(
                    value, (float, int)
                ):
                    raise ValueError(
                        f"Parameter '{key}' in section '{section}' must be a float."
                    )

        # Validate list parameters
        for key, value in propka.items():
            if key in self._ALLOWED_LIST_PARAMS:
                if not isinstance(value, list):
                    raise TypeError(f"Parameter '{key}' in propka must be a list.")
                if key in {"window", "grid"}:
                    if not all(isinstance(item, (float, int)) for item in value):
                        raise ValueError(
                            f"All elements of parameter '{key}' in section propka must be floats or ints"
                        )
                elif key == "file":
                    # check if the list is empty then skip
                    if len(value) == 0:
                        # drop the key from the dictionary
                        continue
                    else:
                        if not all(isinstance(item, str) for item in value):
                            raise ValueError(
                                "All elements of parameter 'file' in section propka must be strings."
                            )
                        # check if the files exist
                        for file in value:
                            if not os.path.isfile(file):
                                raise FileNotFoundError(
                                    f"PropKa file '{file}' does not exist."
                                )

        ff = general.get("ff")
        if ff not in self._FF_ALLOWED:
            raise ValueError(f"Invalid force field '{ff}'. Allowed: {self._FF_ALLOWED}")

        ffout = general.get("ffout")
        if ffout is not None and ffout not in self._FF_ALLOWED:
            raise ValueError(
                f"Invalid output force field '{ffout}'. Allowed: {self._FF_ALLOWED}"
            )
        titration_method = pka.get("titration-state-method")
        if titration_method not in self._TITRATION_METHODS:
            raise ValueError(
                f"Invalid titration method '{titration_method}'. Allowed: {self._TITRATION_METHODS}"
            )

        reference = propka.get("reference")
        if reference not in self._PROPKA_REFERENCE_ALLOWED:
            raise ValueError(
                f"Invalid PropKa reference '{reference}'. Allowed: {self._PROPKA_REFERENCE_ALLOWED}"
            )

        log_level = propka.get("log-level")
        if log_level not in self._LOG_LEVELS:
            raise ValueError(
                f"Invalid log level '{log_level}'. Allowed: {self._LOG_LEVELS}"
            )

        mutator = propka.get("mutator")
        if mutator not in self._MUTATOR_ALLOWED:
            raise ValueError(
                f"Invalid mutator '{mutator}'. Allowed: {self._MUTATOR_ALLOWED}"
            )
        ligand_file = general.get("ligand")
        if ligand_file is not None:
            if not os.path.isfile(ligand_file):
                raise FileNotFoundError(
                    f"Specified ligand file '{ligand_file}' does not exist."
                )

        # check that pH is between 0 and 14
        pH = pka.get("with-ph")
        if not (0.0 <= pH <= 14.0):
            raise ValueError("pH must be between 0 and 14.")
        pH_propka = propka.get("pH")
        if not (0.0 <= pH_propka <= 14.0):
            raise ValueError("PropKa pH must be between 0 and 14.")

        # check that window values are valid
        window = propka.get("window")
        if len(window) != 3:
            raise ValueError("Invalid PropKa window values.")

        # make sure that neutraln and neutralc are only true if ff is PARSE
        neutraln = general.get("neutraln")
        neutralc = general.get("neutralc")
        if ff != "PARSE" and (neutraln or neutralc):
            raise ValueError("neutraln and neutralc can only be True if ff is 'PARSE'.")

        logger.info("Configuration validation passed.")
