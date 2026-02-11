r"""Implementation for the yaml class to write configuration for GNINA."""

from typing import Any, override

import os
from collections.abc import Sequence

from loguru import logger

from .config import YamlConfig


class GninaConfig(YamlConfig):
    r"""Class to handle YAML configuration files for GNINA."""

    # ---- Allowed values / constraints ----
    _SCORING_ALLOWED = {
        "ad4_scoring",
        "dkoes_fast",
        "dkoes_scoring",
        "dkoes_scoring_old",
        "vina",
        "vinardo",
        "default",
    }

    __VECTOR_PARAMS = {
        "center": ("center_x", "center_y", "center_z"),
        "size": ("size_x", "size_y", "size_z"),
    }

    _CNN_SCORING_ALLOWED = {
        "none",
        "rescore",
        "refinement",
        "metrorescore",
        "metrorefine",
        "all",
    }
    _POSE_SORT_ORDER_ALLOWED = {"CNNscore", "CNNaffinity", "Energy"}
    # ---- CNN Model Presets ----
    # Maps preset names to their expanded model lists
    _CNN_PRESETS = {
        "fast": ["all_default_to_default_1_3_1"],
        "ensemble_v1": [
            "dense",
            "general_default2018_3",
            "dense_3",
            "crossdock_default2018",
            "redock_default2018_2",
        ],
        "ensemble_v2": [
            "dense_1_3",
            "dense_1_3_PT_KD_3",
            "crossdock_default2018_KD_4",
        ],
    }

    # Default preset when nothing is specified
    _CNN_DEFAULT_PRESET = "ensemble_v2"

    # All valid built-in CNN model names from GNINA
    _CNN_MODELS_ALLOWED = {
        "all_default_to_default_1_3_1",
        "all_default_to_default_1_3_2",
        "all_default_to_default_1_3_3",
        "crossdock_default2018",
        "crossdock_default2018_1",
        "crossdock_default2018_1_3",
        "crossdock_default2018_1_3_1",
        "crossdock_default2018_1_3_2",
        "crossdock_default2018_1_3_3",
        "crossdock_default2018_1_3_4",
        "crossdock_default2018_2",
        "crossdock_default2018_3",
        "crossdock_default2018_4",
        "crossdock_default2018_KD_1",
        "crossdock_default2018_KD_2",
        "crossdock_default2018_KD_3",
        "crossdock_default2018_KD_4",
        "crossdock_default2018_KD_5",
        "default1.0",
        "default2017",
        "dense",
        "dense_1",
        "dense_1_3",
        "dense_1_3_1",
        "dense_1_3_2",
        "dense_1_3_3",
        "dense_1_3_4",
        "dense_1_3_PT_KD",
        "dense_1_3_PT_KD_1",
        "dense_1_3_PT_KD_2",
        "dense_1_3_PT_KD_3",
        "dense_1_3_PT_KD_4",
        "dense_1_3_PT_KD_def2018",
        "dense_1_3_PT_KD_def2018_1",
        "dense_1_3_PT_KD_def2018_2",
        "dense_1_3_PT_KD_def2018_3",
        "dense_1_3_PT_KD_def2018_4",
        "dense_2",
        "dense_3",
        "dense_4",
        "general_default2018",
        "general_default2018_1",
        "general_default2018_2",
        "general_default2018_3",
        "general_default2018_4",
        "general_default2018_KD_1",
        "general_default2018_KD_2",
        "general_default2018_KD_3",
        "general_default2018_KD_4",
        "general_default2018_KD_5",
        "redock_default2018",
        "redock_default2018_1",
        "redock_default2018_1_3",
        "redock_default2018_1_3_1",
        "redock_default2018_1_3_2",
        "redock_default2018_1_3_3",
        "redock_default2018_1_3_4",
        "redock_default2018_2",
        "redock_default2018_3",
        "redock_default2018_4",
        "redock_default2018_KD_1",
        "redock_default2018_KD_2",
        "redock_default2018_KD_3",
        "redock_default2018_KD_4",
        "redock_default2018_KD_5",
    }

    _ALLOWED_BOOLEAN_PARAMS = {
        "no_lig",
        "score_only",
        "local_only",
        "minimize",
        "randomize_only",
        "accurate_line",
        "simple_ascent",
        "minimize_early_term",
        "minimize_single_full",
        "print_terms",
        "print_atom_types",
        # cnn
        "cnn_mix_emp_force",
        "cnn_mix_emp_energy",
        "cnn_verbose",
        # job
        "quiet",
        "stripH",
        "no_gpu",
        # covalent
        "covalent_fix_lig_atom_position",
        "covalent_optimize_lig",
        # output
        "atom_term_data",
        "full_flex_output",
    }

    # Options that are numeric (int/float)
    _NUMERIC = {
        # search space
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "autobox_add",
        "autobox_extend",
        # scoring/minimization
        "num_mc_steps",
        "max_mc_steps",
        "num_mc_saved",
        "temperature",
        "minimize_iters",
        "factor",
        "force_cap",
        "user_grid_lambda",
        # cnn
        "cnn_rotation",
        "cnn_empirical_weight",
        "cnn_center_x",
        "cnn_center_y",
        "cnn_center_z",
        # job
        "cpu",
        "seed",
        "exhaustiveness",
        "num_modes",
        "min_rmsd_filter",
        "device",
        # covalent
        "covalent_bond_order",
    }

    # File path options (existence checked if provided)
    _FILE_PATH = {
        "receptor",
        "ligand",
        "flex",
        "flexdist_ligand",
        "autobox_ligand",
        "custom_scoring",
        "custom_atoms",
        "user_grid",
        "cnn_model",
        "out_flex",
        "log",
        "out",
        "atom_terms",
    }

    @override
    def __init__(self, file_path: str, data_dict: dict[str, Any] | None = None) -> None:
        if data_dict is None and not os.path.exists(file_path):
            logger.warning("No GNINA config found. Creating default docking config.")
            data_dict = self.declare_defaults()
        super().__init__(file_path, data_dict)
        self.validate()

    def _expand_vectors(self) -> None:
        """expland the 3-item vector parameters into individual x/y/z entries."""
        config = self.data_dict.setdefault("gnina", {})
        region = config.setdefault("docking_region", {})

        for vec_key, (kx, ky, kz) in self.__VECTOR_PARAMS.items():
            vec = region.get(vec_key, None)
            if vec is None:
                continue

            # Accept list/tuple; reject strings and other scalars
            if not isinstance(vec, Sequence) or isinstance(vec, (str, bytes)):
                raise TypeError(
                    f"'docking_region.{vec_key}' must be a list of x, y, z or None."
                )
            if len(vec) != 3:
                raise ValueError(
                    f"'docking_region.{vec_key}' must have exactly 3 values; got {len(vec)}."
                )

            x, y, z = vec
            # allow numeric or None, but normally you'd want numeric
            for name, v in zip((kx, ky, kz), (x, y, z)):
                if v is not None and not isinstance(v, (int, float)):
                    raise TypeError(
                        f"'{name}' must be numeric or None; got {type(v).__name__}."
                    )
                region[name] = v

            del region[vec_key]

    def _expand_cnn_preset(self) -> None:
        """Expand CNN preset names (fast, ensemble_v1, ensemble_v2) into model lists.
        Modifies self.data_dict in place.
        1. If cnn is None or missing, set to default preset (ensemble_v2).
        2. If cnn is a string and matches a preset name, expand to that preset's model list.

        """
        config = self.data_dict.setdefault("gnina", {})
        cnn_section = config.setdefault("cnn", {})

        cnn_value = cnn_section.get("cnn")

        # Default case: None => default preset
        if cnn_value is None:
            logger.info(
                f"No CNN models specified. Using default preset '{self._CNN_DEFAULT_PRESET}': "
                f"{self._CNN_PRESETS[self._CNN_DEFAULT_PRESET]}"
            )
            cnn_section["cnn"] = self._CNN_PRESETS[self._CNN_DEFAULT_PRESET].copy()
            return

        # Case 2: String value - check if it's a preset name
        if isinstance(cnn_value, str):
            if cnn_value in self._CNN_PRESETS:
                logger.info(
                    f"Expanding CNN preset '{cnn_value}' to: {self._CNN_PRESETS[cnn_value]}"
                )
                cnn_section["cnn"] = self._CNN_PRESETS[cnn_value].copy()
                return
            # Single model name (not a preset) - convert to list
            cnn_section["cnn"] = [cnn_value]
            return

    def _validate_cnn_models(self) -> None:
        """Validate that all specified CNN model names are valid GNINA built-in models."""
        config = self.data_dict.get("gnina", {})
        cnn_section = config.get("cnn", {})
        cnn_value = cnn_section.get("cnn")

        if cnn_value is None:
            return

        # At this point, cnn should be a list (after _expand_cnn_preset)
        if not isinstance(cnn_value, list):
            raise TypeError(
                f"CNN models must be a list after expansion; got {type(cnn_value).__name__}."
            )

        # Validate each model name
        invalid_models = []
        for model in cnn_value:
            if not isinstance(model, str):
                raise TypeError(
                    f"CNN model name must be a string; got {type(model).__name__}."
                )
            # Check if it's a valid model or an ensemble pattern (ends with _ensemble)
            if model not in self._CNN_MODELS_ALLOWED:
                # Check for ensemble pattern (PREFIX_ensemble expands to all models starting with PREFIX)
                if model.endswith("_ensemble"):
                    prefix = model[:-9]  # Remove "_ensemble" suffix
                    matching = [
                        m for m in self._CNN_MODELS_ALLOWED if m.startswith(prefix)
                    ]
                    if not matching:
                        invalid_models.append(
                            f"{model} (no models match prefix '{prefix}')"
                        )
                else:
                    invalid_models.append(model)

        if invalid_models:
            raise ValueError(
                f"Invalid CNN model name(s): {invalid_models}. "
                f"Valid presets: {list(self._CNN_PRESETS.keys())}. "
                f"Valid models: see GNINA --help for full list."
            )

    def declare_defaults(self) -> dict[str, dict[str, Any]]:
        r"""Default GNINA config grouped by CLI sections (annotated)."""
        return {
            "gnina": {
                "input": {
                    "receptor": None,  # Receptor file path (e.g., PDBQT)
                    "ligand": None,  # Ligand file path (e.g., PDBQT)
                    "autobox_ligand": None,  # Ligand file for autoboxing (if different from docking ligand)
                },
                # TODO:REVIEW THESES
                "flexibility": {
                    "flex": None,  # Flexible receptor file / definition (enables receptor flexibility)
                    "flexres": None,  # Flexible residues selection (e.g., "A:TYR123,GLU45")
                    "flexdist_ligand": None,  # Flex residues within distance of ligand (requires ligand)
                    "flexdist": None,  # Flex residues within distance of search center / box
                    "flex_limit": None,  # Limit number of flexible residues/atoms considered
                    "flex_max": None,  # Maximum flexible residues/atoms allowed
                    "out_flex": None,  # Output flexible receptor coordinates (if flex used)
                    "full_flex_output": False,  # If True: more complete flex receptor output (larger files)
                },
                "docking_region": {
                    # NOTE: Choose EITHER (center_* + size_*) OR autobox_ligand (+ autobox_add/extend)
                    "center_x": None,  # Docking box center X (Angstrom)
                    "center_y": None,  # Docking box center Y (Angstrom)
                    "center_z": None,  # Docking box center Z (Angstrom)
                    "size_x": None,  # Docking box size X (Angstrom)
                    "size_y": None,  # Docking box size Y (Angstrom)
                    "size_z": None,  # Docking box size Z (Angstrom)
                    "autobox_add": 4.0,  # Extra padding to the autobox (Angstrom)
                    "autobox_extend": 1,  # Extend autobox to include nearby receptor space (integer behavior varies by build)
                    "no_lig": False,  # If True: run without ligand (receptor/flex-only sampling/minimization; disables ligand docking/scoring)
                },
                # ----------------------------
                # Covalent docking controls
                # ----------------------------
                # TODO:REVIEW THESES
                "covalent": {
                    "covalent_rec_atom": None,  # Receptor atom to bond to (e.g., "A:SER123:OG")
                    "covalent_lig_atom_pattern": None,  # Pattern/SMARTS-like rule to identify ligand reactive atom(s)
                    "covalent_lig_atom_position": None,  # Explicit ligand reactive atom position / index (depends on input format)
                    "covalent_fix_lig_atom_position": False,  # If True: keep ligand reactive atom fixed during docking
                    "covalent_bond_order": 1,  # Covalent bond order (1 single, 2 double, etc.)
                    "covalent_optimize_lig": False,  # If True: allow ligand geometry optimization for covalent setup
                },
                "scoring": {
                    "scoring": "vina",  # Scoring function: ad4_scoring|vina|vinardo|dkoes_*|default
                    "custom_scoring": None,  # Path to custom scoring file (if supported by build)
                    "custom_atoms": None,  # Custom atom types mapping (if using custom scoring)
                    "score_only": False,  # If True: score input pose(s) only (no docking)
                    "local_only": False,  # If True: restrict search locally around initial pose
                    "minimize": False,  # If True: energy-minimize pose(s) (often used with score_only / rescoring)
                    "randomize_only": False,  # If True: randomize poses only (no full search)
                    "num_mc_steps": None,  # Monte Carlo steps for search/minimization (None : automatically based on mobile  atoms and degrees  of freedom  within the ligand.)
                    "max_mc_steps": None,  # Hard cap on MC steps
                    "num_mc_saved": 50,  # How many MC states/poses to keep
                    "temperature": 0,  # Temperature for MC acceptance (if used)
                    "minimize_iters": 0,  # Iterations for minimizer (0 may mean default/no extra)
                    "accurate_line": False,  # Use more accurate line search (slower, can be more stable)
                    "simple_ascent": False,  # Simpler ascent/descent behavior (advanced tuning)
                    "minimize_early_term": False,  # Early terminate minimization when converged
                    "minimize_single_full": False,  # Single full minimization pass (advanced)
                    "approximation": None,  # Force field approximation: linear|spline|exact (if supported)
                    "factor": None,  # Scale factor (typically for forces/energies; build-dependent)
                    "force_cap": None,  # Cap on maximum force magnitude (stability / outliers)
                    "user_grid": None,  # External grid file (bypass internal grid generation)
                    "user_grid_lambda": -1,  # Grid mixing/weighting factor (build-dependent default -1)
                    "print_terms": False,  # Print energy term breakdown
                    "print_atom_types": False,  # Print atom typing info
                },
                "cnn": {
                    "cnn_scoring": "rescore",  # CNN scoring mode: none|rescore|refinement|metrorescore|metrorefine|all
                    # CNN model selection. Accepts:
                    #   - None: uses default ensemble (ensemble_v2)
                    #   - "fast": single fast model (all_default_to_default_1_3_1)
                    #   - "ensemble_v1": GNINA 1.0 default (dense, general_default2018_3, dense_3, crossdock_default2018, redock_default2018_2)
                    #   - "ensemble_v2": current GNINA default (dense_1_3, dense_1_3_PT_KD_3, crossdock_default2018_KD_4)
                    #   - List of specific model names for custom ensemble
                    "cnn": "ensemble_v2",  # CNN model(s) to use (preset name or list of model names; None defaults to ensemble_v2)
                    "cnn_model": None,  # Path to custom torch model file (overrides built-in models)
                    "cnn_rotation": 0,  # Number of rotations/ensembling steps for CNN scoring
                    "cnn_mix_emp_force": False,  # Mix empirical force with CNN guidance (advanced)
                    "cnn_mix_emp_energy": False,  # Mix empirical energy with CNN (advanced)
                    "cnn_empirical_weight": 1.0,  # Weight of empirical component in mixing
                    "cnn_center_x": None,  # Optional CNN grid center X (override)
                    "cnn_center_y": None,  # Optional CNN grid center Y (override)
                    "cnn_center_z": None,  # Optional CNN grid center Z (override)
                    "cnn_verbose": False,  # More CNN debug output
                },
                "output": {
                    "out": "docked.sdf.gz",  # Output poses file
                    "log": None,  # Log file path
                    "atom_terms": None,  # Per-atom term output file (if enabled)
                    "atom_term_data": False,  # If True: include atom-level scoring term data
                    "pose_sort_order": "CNNscore",  # Pose ranking: CNNscore (default), CNNaffinity, Vina energy
                },
                "misc": {
                    "cpu": None,  # CPU threads (None = auto-detect)
                    "seed": 0,  # RNG seed (0 often means seeded deterministically; varies by build)
                    "exhaustiveness": 8,  # Search thoroughness (higher = slower, more thorough)
                    "num_modes": 9,  # Number of output poses to write
                    "min_rmsd_filter": 1.0,  # Minimum RMSD between output poses (Angstrom)
                    "quiet": False,  # Reduce console output
                    "addH": None,  # Add hydrogens (behavior depends on input format/build)
                    "stripH": False,  # Remove hydrogens from inputs
                    "device": 0,  # GPU device index
                    "no_gpu": True,  # Force CPU only (disable GPU)
                },
            }
        }

    def validate(self) -> None:
        # fill in defaults
        self._expand_vectors()
        self._expand_cnn_preset()  # Expand CNN presets before other validation

        default = self.declare_defaults()
        gnina_cfg = self.data_dict.setdefault("gnina", {})

        for section, params in default["gnina"].items():
            section_cfg = gnina_cfg.setdefault(section, {})
            for key, value in params.items():
                if key not in section_cfg:
                    logger.warning(
                        f"Missing parameter '{key}' in GNINA section '{section}'. "
                        f"Using default value: {value}"
                    )
                    section_cfg[key] = value

        # ensure any filled-in defaults are saved
        self.write_config(self.data_dict)

        config = self.data_dict.get("gnina", {})
        flexibility = config.get("flexibility", {})
        region = config.get("docking_region", {})
        covalent = config.get("covalent", {})
        scoring = config.get("scoring", {})
        cnn = config.get("cnn", {})
        output = config.get("output", {})
        misc = config.get("misc", {})

        # type validation
        flat: dict[str, Any] = {}
        for section_name in (
            "flexibility",
            "docking_region",
            "covalent",
            "scoring",
            "cnn",
            "output",
            "misc",
        ):
            section_cfg = config.get(section_name, {})
            for k, v in section_cfg.items():
                flat[k] = v

        # Booleans
        for key in self._ALLOWED_BOOLEAN_PARAMS:
            if key in flat:
                value = flat.get(key)
                if not isinstance(value, bool):
                    raise ValueError(f"Parameter '{key}' must be a boolean.")

        # Numeric: allow int/float, but some fields should be int-only (validated later)
        for key in self._NUMERIC:
            if key in flat:
                value = flat.get(key)
                if value is None:
                    continue
                if not isinstance(value, (int, float)):
                    raise ValueError(
                        f"Parameter '{key}' must be numeric (int/float) or None."
                    )

        # Strings/None for anything not explicitly typed (bool/numeric), except where we allow special types
        # (e.g. addH may be bool/str in some wrappers, and cnn may be None/bool/list depending on how you expose it)
        _allow_any = {"addH", "cnn"}  # keep permissive; validate lightly below
        for section_name in (
            "flexibility",
            "docking_region",
            "covalent",
            "scoring",
            "cnn",
            "output",
            "misc",
        ):
            section_cfg = config.get(section_name, {})
            for key, value in section_cfg.items():
                if key in _allow_any:
                    continue
                if key not in self._ALLOWED_BOOLEAN_PARAMS and key not in self._NUMERIC:
                    # treat as string-ish fields (paths, enums, selectors), allow None
                    if value is not None and not isinstance(value, str):
                        raise TypeError(
                            f"Parameter '{key}' in section '{section_name}' must be a string or None."
                        )

        # addH: allow None/bool/str (some pipelines pass "true"/"false"/"polar"/etc.)
        if "addH" in misc:
            v = misc.get("addH")
            if v is not None and not isinstance(v, (bool, str)):
                raise TypeError(
                    "Parameter 'addH' in section 'misc' must be bool, string, or None."
                )

        # cnn: allow None/str/list (preset name, single model, or list of models)
        # After _expand_cnn_preset, it should always be a list
        if "cnn" in cnn:
            v = cnn.get("cnn")
            if v is not None and not isinstance(v, list):
                raise TypeError(
                    "Parameter 'cnn' in section 'cnn' must be a list of model names after expansion."
                )

        # Validate CNN model names
        self._validate_cnn_models()

        # file path validation
        # NOTE: output paths should NOT be required to exist; treat them differently.
        _output_like = {"out", "out_flex", "log", "atom_terms"}

        for key in self._FILE_PATH:
            if key not in flat:
                continue
            value = flat.get(key)
            if value is None:
                continue
            if not isinstance(value, str):
                raise TypeError(f"Parameter '{key}' must be a path string or None.")

            # Only enforce existence for input-like files; output files may not exist yet.
            if key in _output_like:
                continue

            if not os.path.exists(value):
                raise FileNotFoundError(f"Path for '{key}' does not exist: {value}")

        # --- Enums / allowed sets ---
        scoring_name = scoring.get("scoring")
        if (
            scoring_name not in (None, "vina")
            and scoring_name not in self._SCORING_ALLOWED
        ):
            raise ValueError(
                f"Invalid scoring '{scoring_name}'. Allowed: {sorted(self._SCORING_ALLOWED)} or 'vina'/None."
            )

        cnn_scoring = cnn.get("cnn_scoring")
        if cnn_scoring is not None and cnn_scoring not in self._CNN_SCORING_ALLOWED:
            raise ValueError(
                f"Invalid cnn_scoring '{cnn_scoring}'. Allowed: {sorted(self._CNN_SCORING_ALLOWED)}."
            )

        approximation = scoring.get("approximation")
        if approximation is not None and approximation not in {
            "linear",
            "spline",
            "exact",
        }:
            raise ValueError(
                "approximation must be one of: linear|spline|exact (or None)."
            )

        # --- Integer-only fields ---
        _int_only = {
            "autobox_extend",
            "num_mc_steps",
            "max_mc_steps",
            "num_mc_saved",
            "minimize_iters",
            "cnn_rotation",
            "cpu",
            "seed",
            "exhaustiveness",
            "num_modes",
            "device",
            "covalent_bond_order",
        }
        for key in _int_only:
            if key in flat and flat[key] is not None and not isinstance(flat[key], int):
                raise TypeError(f"Parameter '{key}' must be an integer (or None).")

        # --- Basic numeric ranges ---
        # docking region sizes must be positive if provided
        for k in ("size_x", "size_y", "size_z"):
            v = region.get(k)
            if v is not None and float(v) <= 0.0:
                raise ValueError(f"{k} must be > 0.")

        # exhaustiveness / num_modes / cpu
        ex = misc.get("exhaustiveness")
        if ex is not None and int(ex) < 1:
            raise ValueError("exhaustiveness must be >= 1.")

        nm = misc.get("num_modes")
        if nm is not None and int(nm) < 1:
            raise ValueError("num_modes must be >= 1.")

        cpu = misc.get("cpu")
        if cpu is not None and int(cpu) == 0:
            raise ValueError("cpu must be non-zero (use None for auto-detect).")
        # check cpu vs exhaustiveness they should be cpu>=exhaustiveness gives a warning
        if cpu is not None and ex is not None and int(cpu) < int(ex):
            logger.warning(
                "cpu is less than exhaustiveness; this may lead to suboptimal performance."
            )

        mrf = misc.get("min_rmsd_filter")
        if mrf is not None and float(mrf) < 0.0:
            raise ValueError("min_rmsd_filter must be >= 0.")

        pso = output.get("pose_sort_order")
        if pso is not None and pso not in self._POSE_SORT_ORDER_ALLOWED:
            raise ValueError(
                f"Invalid pose_sort_order '{pso}'. Allowed values: {sorted(self._POSE_SORT_ORDER_ALLOWED)}."
            )

        # user_grid_lambda: allow -1 default; otherwise expect >= 0 typically
        ugl = scoring.get("user_grid_lambda")
        if ugl is not None and float(ugl) < -1.0:
            raise ValueError("user_grid_lambda must be >= -1.")

        # covalent_bond_order: common safe constraint
        cbo = covalent.get("covalent_bond_order")
        if cbo is not None and int(cbo) < 1:
            raise ValueError("covalent_bond_order must be >= 1.")

        # --- Cross-field constraints ---
        center_keys = ("center_x", "center_y", "center_z")
        size_keys = ("size_x", "size_y", "size_z")

        has_any_center = any(region.get(k) is not None for k in center_keys)
        has_any_size = any(region.get(k) is not None for k in size_keys)
        has_all_center = all(region.get(k) is not None for k in center_keys)
        has_all_size = all(region.get(k) is not None for k in size_keys)

        autobox_lig = region.get("autobox_ligand")
        using_autobox = autobox_lig is not None

        if using_autobox and (has_any_center or has_any_size):
            raise ValueError(
                "Invalid docking_region: choose either (center_* + size_*) OR autobox_ligand, not both."
            )

        if not using_autobox:
            # If they specify any of center/size, require all of them
            if (has_any_center or has_any_size) and not (
                has_all_center and has_all_size
            ):
                raise ValueError(
                    "Invalid docking_region: when using explicit box, you must set all center_* and all size_*."
                )

        # random logic arguments to ensure the jobs makes sense
        no_lig = bool(region.get("no_lig"))
        if no_lig:
            # no ligand present => autobox by ligand doesn't make sense
            if region.get("autobox_ligand") is not None:
                raise ValueError(
                    "no_lig=True is incompatible with autobox_ligand (requires a ligand)."
                )
            # flexdist_ligand requires a ligand
            if flexibility.get("flexdist_ligand") is not None:
                raise ValueError(
                    "no_lig=True is incompatible with flexdist_ligand (requires a ligand)."
                )

        # scoring logic
        if scoring.get("score_only") and scoring.get("randomize_only"):
            raise ValueError("score_only and randomize_only cannot both be True.")

        # local_only usually implies there is an initial pose / starting structure;
        if scoring.get("local_only") and scoring.get("score_only"):
            raise ValueError(
                "local_only and score_only cannot both be True (no docking/search)."
            )

        if scoring.get("local_only") and not scoring.get("minimize"):
            logger.warning(
                "local_only=True without minimize; docking may be very limited."
            )

        # cnn_model only makes sense if cnn_scoring isn't "none"
        if cnn.get("cnn_scoring") == "none" and cnn.get("cnn_model") is not None:
            logger.warning(
                "cnn_model is set but cnn_scoring='none' (model will likely be ignored)."
            )

        # cnn models only make sense if cnn_scoring isn't "none"
        if cnn.get("cnn_scoring") == "none" and cnn.get("cnn"):
            logger.warning(
                "cnn models are set but cnn_scoring='none' (models will likely be ignored)."
            )

        # no_gpu/device consistency
        if misc.get("no_gpu") and misc.get("device") not in (None, 0):
            logger.warning(
                "no_gpu=True but device is set; device will be ignored when GPU is disabled."
            )

        logger.info("GNINA configuration validation passed.")
