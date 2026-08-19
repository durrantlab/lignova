# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Implementation for the yaml class to write configuration for Meeko."""

import importlib.util
import os
from typing import Any, override

from loguru import logger

from .config import YamlConfig


class MeekoConfig(YamlConfig):
    """Class to handle YAML configuration files for Meeko receptor preparation."""

    # Schema for mk_prepare_receptor.py configuration
    _CHARGE_MODEL_ALLOWED = {"gasteiger", "espaloma", "zero", "read"}
    _INPUT_PARAMS = {"read_pdb", "read_pqr", "read_with_prody"}
    _OUTPUT_PARAMS = {
        "write_pdbqt",
        "write_json",
        "write_pdb",
        "write_gpf",
        "write_vina_box",
    }
    _BOX_ENVELOPING_EXTENSIONS = {".sdf", ".mol", ".mol2", ".pdb", ".pdbqt"}
    _ALLOWED_BOOLEAN_PARAMS = {
        "allow_bad_res",
        "box_center_off_reactive_res",
    }

    _EXTENSION_READERS = {
        ".pqr": "read_pqr",
        ".pdb": "read_pdb",
        ".cif": "read_with_prody",
        ".mmcif": "read_with_prody",
    }
    # NOTE: mk_prepare_receptor takes these with an optional filename, so they are
    # either True (bare flag, filename derived from output_basename) or a path string
    _ALLOWED_FLAG_OR_PATH_PARAMS = {
        "write_pdbqt",
        "write_json",
        "write_gpf",
        "write_vina_box",
        "cache_templates",
    }
    _ALLOWED_FLOAT_PARAMS = {
        "padding",
        "r_eq_12",
        "eps_12",
        "r_eq_13_scaling",
        "r_eq_14_scaling",
    }
    _ALLOWED_LIST_PARAMS = {"box_size", "box_center"}
    _MAX_REACTIVE_FLEXRES = 8

    @override
    def __init__(self, file_path: str, data_dict: dict[str, Any] | None = None) -> None:
        """Initialize with the path to the YAML file.
        Args:
            file_path : Path to the YAML configuration file.
            data_dict : Dictionary to create the YAML file if it doesn't exist.
        """
        if data_dict is None and not os.path.exists(file_path):
            logger.warning(
                "No configuration file found. Creating default meeko config."
            )
            data_dict = self.declare_defaults()
        super().__init__(file_path, data_dict)
        self.validate()

    def declare_defaults(self) -> dict[str, dict[str, Any]]:
        """Declare default settings for meeko config (grouped under meeko)."""
        default_config: dict[str, dict[str, Any]] = {
            "meeko": {
                "input_output": {
                    "read_pqr": None,
                    "read_pdb": None,
                    "read_with_prody": None,
                    "output_basename": "receptor",
                    "write_pdbqt": True,  # True so it uses output_basename, or give a path
                    "write_json": False,
                    "write_pdb": None,
                    "write_gpf": False,
                    # NOTE: gnina autoboxes off the ligand, so no vina box file is needed
                    "write_vina_box": False,
                    "debug_fn": None,
                },
                "receptor_perception": {
                    "charge_model": "read",  # gasteiger|espaloma|zero|read
                    "set_template": None,  # e.g., "A:5,7=CYX,B:17=HID"
                    "delete_residues": None,  # e.g., "A:350,B:15,16,17"
                    "blunt_ends": None,  # e.g., "A:123,200=2,A:1=0"
                    "add_templates": None,  # additional residue templates json file or resname:file.sdf
                    "cache_templates": False,  # True or a JSON cache path
                    "mk_config": None,
                    "allow_bad_res": False,
                    "default_altloc": None,  # default alternative location, overridden by wanted_altloc
                    "wanted_altloc": None,  # e.g., ":5=B,B:17=A"
                    "flexres": None,  # e.g., ":42,B:23"
                    "rot_terminal_group": None,  # e.g., ":42,B:23"
                },
                "box": {
                    "box_size": None,  # [x, y, z] in Angstrom
                    "box_center": None,  # [x, y, z] in Angstrom
                    "box_center_off_reactive_res": False,
                    "box_enveloping": None,  # [.sdf .mol .mol2 .pdb .pdbqt]
                    "padding": None,  # only used with box_enveloping [A]
                },
                "reactive": {
                    "reactive_flexres": None,  # same syntax as flexres (max 8)
                    "reactive_name": None,  # e.g., "TRP:NE1"
                    "reactive_name_specific": None,  # e.g., "A:42=NE2"
                    # NOTE: left as None so meeko applies its own defaults and the
                    # flags stay off the command line when no residue is reactive
                    "r_eq_12": None,  # equilibrium distance for 1-2 interactions, meeko default 1.8 [A]
                    "eps_12": None,  # epsilon for 1-2 interactions, meeko default 2.5
                    "r_eq_13_scaling": None,  # scaling factor for 1-3 interactions, meeko default 0.5
                    "r_eq_14_scaling": None,  # scaling factor for 1-4 interactions, meeko default 0.5
                },
            }
        }
        return default_config

    def validate(self) -> None:
        """Validate the current configuration against allowed values."""
        default = self.declare_defaults()
        meeko_cfg = self.data_dict.setdefault("meeko", {})
        for section, params in default["meeko"].items():
            section_cfg = meeko_cfg.setdefault(section, {})
            for key, value in params.items():
                if key not in section_cfg:
                    logger.warning(
                        f"Missing parameter '{key}' in section '{section}'. Using default value: {value}"
                    )
                    section_cfg[key] = value
        self.write_config(self.data_dict)
        config = self.data_dict.get("meeko", {})
        sections = ["input_output", "receptor_perception", "box", "reactive"]
        input_output = config.get("input_output", {})
        perception = config.get("receptor_perception", {})
        box = config.get("box", {})
        reactive = config.get("reactive", {})
        # Validate string parameters
        for section in sections:
            for key, value in config.get(section, {}).items():
                if (
                    key not in self._ALLOWED_BOOLEAN_PARAMS
                    and key not in self._ALLOWED_FLAG_OR_PATH_PARAMS
                    and key not in self._ALLOWED_FLOAT_PARAMS
                    and key not in self._ALLOWED_LIST_PARAMS
                ):
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
        # Validate flag parameters that optionally carry a filename
        for section in sections:
            for key, value in config.get(section, {}).items():
                if key in self._ALLOWED_FLAG_OR_PATH_PARAMS and not isinstance(
                    value, (bool, str)
                ):
                    raise TypeError(
                        f"Parameter '{key}' in '{section}' must be a boolean or a string."
                    )
        # Validate float parameters
        for section in sections:
            for key, value in config.get(section, {}).items():
                if (
                    key in self._ALLOWED_FLOAT_PARAMS
                    and value is not None
                    and not isinstance(value, (float, int))
                ):
                    raise ValueError(
                        f"Parameter '{key}' in '{section}' must be a float or None."
                    )

        # Validate list parameters
        for key, value in box.items():
            if key in self._ALLOWED_LIST_PARAMS:
                if value is None:
                    continue
                if not isinstance(value, list):
                    raise TypeError(f"Parameter '{key}' in box must be a list or None.")
                if len(value) != 3:
                    raise ValueError(
                        f"Parameter '{key}' in section box must have exactly 3 values (x, y, z)."
                    )
                if not all(isinstance(item, (float, int)) for item in value):
                    raise ValueError(
                        f"All elements of parameter '{key}' in section box must be floats or ints"
                    )

        # exactly one input source is accepted by mk_prepare_receptor
        given_inputs = [
            key for key in self._INPUT_PARAMS if input_output.get(key) is not None
        ]
        if len(given_inputs) > 1:
            raise ValueError(
                f"Only one input may be given. Got: {given_inputs}. Allowed: {self._INPUT_PARAMS}"
            )
        if len(given_inputs) == 0:
            logger.warning(
                "No input file set. One of 'read_pqr', 'read_pdb' or 'read_with_prody' must be provided before running mk_prepare_receptor."
            )
        for key in given_inputs:
            input_file = input_output.get(key)
            if not os.path.isfile(input_file):
                raise FileNotFoundError(
                    f"Specified input file '{input_file}' does not exist."
                )

        if input_output.get("read_with_prody") is not None:
            if importlib.util.find_spec("prody") is None:
                raise ImportError(
                    "ProDy is required for 'read_with_prody' but is not installed."
                )

        if not any(input_output.get(key) for key in self._OUTPUT_PARAMS):
            logger.warning(
                "No output requested. Enable one of 'write_pdbqt', 'write_json', "
                "'write_pdb', 'write_gpf' or 'write_vina_box'."
            )
        charge_model = perception.get("charge_model")
        if charge_model not in self._CHARGE_MODEL_ALLOWED:
            raise ValueError(
                f"Invalid charge model '{charge_model}'. Allowed: {self._CHARGE_MODEL_ALLOWED}"
            )
        if charge_model == "read":
            if input_output.get("read_pqr") is None and len(given_inputs) > 0:
                raise ValueError(
                    "Charge model 'read' requires 'read_pqr' as the input file."
                )
            if len(given_inputs) == 0:
                logger.warning(
                    "Charge model 'read' requires 'read_pqr' to be set before running."
                )

        mk_config = perception.get("mk_config")
        if mk_config is not None:
            if not mk_config.endswith(".json"):
                raise ValueError(f"Meeko config '{mk_config}' must be a .json file.")
            if not os.path.isfile(mk_config):
                raise FileNotFoundError(
                    f"Specified meeko config '{mk_config}' does not exist."
                )

        add_templates = perception.get("add_templates")
        if add_templates is not None and ":" not in add_templates:
            if not os.path.isfile(add_templates):
                raise FileNotFoundError(
                    f"Specified templates file '{add_templates}' does not exist."
                )

        box_size = box.get("box_size")
        box_center = box.get("box_center")
        box_enveloping = box.get("box_enveloping")
        off_reactive = box.get("box_center_off_reactive_res")
        if box_enveloping is not None and (
            box_size is not None or box_center is not None
        ):
            raise ValueError(
                "'box_enveloping' cannot be combined with 'box_size' or 'box_center'."
            )
        if off_reactive and box_center is not None:
            raise ValueError(
                "'box_center_off_reactive_res' cannot be combined with 'box_center'."
            )
        if box_enveloping is not None:
            extension = os.path.splitext(box_enveloping)[1].lower()
            if extension not in self._BOX_ENVELOPING_EXTENSIONS:
                raise ValueError(
                    f"Invalid box_enveloping extension '{extension}'. Allowed: {self._BOX_ENVELOPING_EXTENSIONS}"
                )
            if not os.path.isfile(box_enveloping):
                raise FileNotFoundError(
                    f"Specified box_enveloping file '{box_enveloping}' does not exist."
                )
        padding = box.get("padding")
        if padding is not None:
            if box_enveloping is None:
                logger.warning(
                    "Parameter 'padding' is only applied together with 'box_enveloping'."
                )
            if padding < 0.0:
                raise ValueError("padding must be non-negative.")
        if box_size is not None and not all(value > 0.0 for value in box_size):
            raise ValueError("All box_size values must be positive.")
        # a box needs a center: explicit, off a reactive residue, or from enveloping
        if box_size is not None and box_center is None and not off_reactive:
            raise ValueError(
                "'box_size' requires 'box_center' or 'box_center_off_reactive_res'."
            )
        if box_size is None and box_enveloping is None:
            # these outputs cannot be written without a defined grid box
            if input_output.get("write_vina_box") or input_output.get("write_gpf"):
                raise ValueError(
                    "'write_vina_box' and 'write_gpf' require 'box_size' or 'box_enveloping'."
                )
            logger.warning(
                "No grid box defined. Set 'box_size' or 'box_enveloping' before docking."
            )

        # max 8 reactive flexible residues
        reactive_flexres = reactive.get("reactive_flexres")
        n_reactive = (
            len([res for res in reactive_flexres.split(",") if res.strip()])
            if reactive_flexres is not None
            else 0
        )
        if n_reactive > self._MAX_REACTIVE_FLEXRES:
            raise ValueError(
                f"At most {self._MAX_REACTIVE_FLEXRES} reactive flexible residues are allowed, got {n_reactive}."
            )
        if off_reactive and n_reactive != 1:
            raise ValueError(
                "'box_center_off_reactive_res' requires exactly one reactive flexible residue."
            )
        if n_reactive == 0:
            for key in ["reactive_name", "reactive_name_specific"] + sorted(
                self._ALLOWED_FLOAT_PARAMS - {"padding"}
            ):
                if reactive.get(key) is not None:
                    logger.warning(
                        f"Parameter '{key}' in section 'reactive' is ignored without 'reactive_flexres'."
                    )

        logger.info("Configuration validation passed.")
