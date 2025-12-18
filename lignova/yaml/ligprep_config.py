r"""Implementation for the yaml class to write configuration for Gypsum-DL."""

from typing import Any, override

import os

from loguru import logger

from .config import YamlConfig


class GypsumDLConfig(YamlConfig):
    r"""Class to handle YAML configuration files for Gypsum-DL."""

    # ---- Schema/constraints ----
    _JOB_MANAGERS = {"mpi", "multiprocessing", "serial"}

    _ALLOWED_BOOLEAN_PARAMS = {
        "separate_output_files",
        "add_pdb_output",
        "add_html_output",
        "skip_optimize_geometry",
        "skip_alternate_ring_conformations",
        "skip_adding_hydrogen",
        "skip_making_tautomers",
        "skip_enumerate_chiral_mol",
        "skip_enumerate_double_bonds",
        "let_tautomers_change_chirality",
        "use_durrant_lab_filters",
        "2d_output_only",
        "cache_prerun",
    }

    _ALLOWED_INT_PARAMS = {
        "num_processors",
        "max_variants_per_compound",
        "thoroughness",
    }

    _ALLOWED_FLOAT_PARAMS = {
        "min_ph",
        "max_ph",
        "pka_precision",
    }

    @override
    def __init__(self, file_path: str, data_dict: dict[str, Any] | None = None) -> None:
        r"""Initialize with the path to the YAML file.
        Args:
            file_path (str): Path to the YAML configuration file.
            data_dict (dict| None) : Dictionary to create the YAML file if it doesn't exist.
        """
        if data_dict is None and not os.path.exists(file_path):
            logger.warning(
                "No Gypsum-DL configuration file found. Creating default Gypsum-DL config."
            )
            data_dict = self.declare_defaults()
        super().__init__(file_path, data_dict)
        self.validate()

    def declare_defaults(self) -> dict[str, dict[str, Any]]:
        r"""Declare default settings for Gypsum-DL (grouped under gypsum_dl)."""
        default_config: dict[str, dict[str, Any]] = {
            "gypsum_dl": {
                "job_specs": {
                    # IO / basic control
                    "job_manager": "multiprocessing",  # mpi|multiprocessing|serial
                    "num_processors": 4,  # -1 to use all, 1 for serial
                    # NOTE: LEAVE THIS AS THE DEFAULT VALUE FOR NOW TO AVOID MASSIVE ENUMERATIONS
                    "max_variants_per_compound": 5,
                    "thoroughness": 3,  # higher => more exhaustive / slower
                    "tasks_per_processor": 1,  # this is for mpi only
                },
                "format": {
                    # Output format
                    "separate_output_files": False,
                    "add_pdb_output": False,
                    "add_html_output": False,
                    "2d_output_only": False,
                    # Enumeration toggles
                    "skip_optimize_geometry": False,
                    "skip_alternate_ring_conformations": False,
                    "skip_adding_hydrogen": False,
                    "skip_making_tautomers": False,
                    "skip_enumerate_chiral_mol": False,
                    "skip_enumerate_double_bonds": False,
                    "let_tautomers_change_chirality": False,
                    # Filters / misc
                    "use_durrant_lab_filters": True,  # keep
                    "cache_prerun": False,
                },
                "ph": {
                    # pH / ionization (Dimorphite-DL)
                    "min_ph": 6.4,
                    "max_ph": 8.4,
                    "pka_precision": 1,
                },
            }
        }
        return default_config

    def validate(self) -> None:
        r"""Validate the current configuration against allowed values."""
        # Fill in any missing keys from defaults
        default = self.declare_defaults()
        gypsum_cfg = self.data_dict.setdefault("gypsum_dl", {})
        for section, params in default["gypsum_dl"].items():
            section_cfg = gypsum_cfg.setdefault(section, {})
            for key, value in params.items():
                if key not in section_cfg:
                    logger.warning(
                        f"Missing parameter '{key}' in Gypsum-DL section '{section}'. "
                        f"Using default value: {value}"
                    )
                    section_cfg[key] = value
        # persist any filled-in defaults
        self.write_config(self.data_dict)

        # Handles
        config = self.data_dict.get("gypsum_dl", {})
        job_specs = config.get("job_specs", {})
        fmt = config.get("format", {})
        ph = config.get("ph", {})

        # Strings for any parameter that is not explicitly typed as bool/int/float
        for section_name in ("job_specs", "format", "ph"):
            section_cfg = config.get(section_name, {})
            for key, value in section_cfg.items():
                if (
                    key not in self._ALLOWED_BOOLEAN_PARAMS
                    and key not in self._ALLOWED_INT_PARAMS
                    and key not in self._ALLOWED_FLOAT_PARAMS
                    and key != "tasks_per_processor"
                ):
                    if not isinstance(value, str) and value is not None:
                        raise TypeError(
                            f"Parameter '{key}' in section '{section_name}' "
                            f"must be a string or None."
                        )

        for key in self._ALLOWED_BOOLEAN_PARAMS:
            value = fmt.get(key)
            if not isinstance(value, bool):
                raise ValueError(
                    f"Parameter '{key}' in section 'format' must be a boolean."
                )

        for key in self._ALLOWED_INT_PARAMS:
            value = job_specs.get(key)
            if not isinstance(value, int):
                raise ValueError(
                    f"Parameter '{key}' in section 'job_specs' must be an integer."
                )

        for key in self._ALLOWED_FLOAT_PARAMS:
            value = ph.get(key)
            if not isinstance(value, (float, int)):
                raise ValueError(f"Parameter '{key}' in section 'ph' must be a float.")

        # --- Semantic validation ---

        # job_manager
        job_manager = job_specs.get("job_manager")
        if job_manager not in self._JOB_MANAGERS:
            raise ValueError(
                f"Invalid job_manager '{job_manager}'. "
                f"Allowed: {self._JOB_MANAGERS}"
            )

        # num_processors
        num_proc = job_specs.get("num_processors")
        if num_proc == 0:
            raise ValueError("num_processors must be non-zero (use -1 for all cores).")

        # thoroughness
        thoroughness = job_specs.get("thoroughness")
        if thoroughness < 1:
            raise ValueError("thoroughness must be >= 1.")

        # max_variants_per_compound
        max_variants = job_specs.get("max_variants_per_compound")
        if max_variants < 1:
            raise ValueError("max_variants_per_compound must be >= 1.")

        # pH bounds
        min_ph = float(ph.get("min_ph"))
        max_ph = float(ph.get("max_ph"))
        if not (0.0 <= min_ph <= 14.0):
            raise ValueError("min_ph must be between 0 and 14.")
        if not (0.0 <= max_ph <= 14.0):
            raise ValueError("max_ph must be between 0 and 14.")
        if min_ph > max_ph:
            raise ValueError("min_ph must be <= max_ph.")

        pka_precision = float(ph.get("pka_precision"))
        if pka_precision <= 0:
            raise ValueError("pka_precision must be > 0.")

        logger.info("Gypsum-DL configuration validation passed.")

        # ensure that tasks_per_processor is an int when job_manager is mpi
        job_manager = job_specs.get("job_manager")
        tasks_per_processor = job_specs.get("tasks_per_processor")
        num_processors = job_specs.get("num_processors")
        if job_manager == "mpi":
            # check that it is not none
            if tasks_per_processor is None:
                raise ValueError("tasks_per_processor must be specified for mpi job.")
            elif not isinstance(tasks_per_processor, int):
                raise TypeError("tasks_per_processor must be an integer for mpi job.")
            if tasks_per_processor < 1:
                raise ValueError(
                    "tasks_per_processor must be a positive integer for mpi job"
                )
            if num_processors != -1:
                if num_processors <= 1:
                    raise ValueError(
                        "num_processors must be a positive integer more than 1 for mpi job OR -1 to use all available processor."
                    )
