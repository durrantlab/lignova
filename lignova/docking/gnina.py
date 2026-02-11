r"""Implementation for docking  GNINA docking."""

from typing import override

import os
from collections.abc import Iterable

from loguru import logger

from lignova.docking.docking import Docking
from lignova.structure.ligand import DockedLigand, PreparedLigand
from lignova.structure.protein import PreparedProtein
from lignova.yaml.docking_config import GninaConfig


class GNINA(Docking):
    """Class to dock ligands in determined protein pocket using GNINA."""

    def __init__(
        self,
        autobox: bool = True,
        box_ligand: str | None = None,
    ) -> None:
        r"""Initialize the GNINA docking class.

        Args:
            autobox
                Whether to automatically determine the docking box based on the ligand. Default is True.
            box_ligand
                Path to the ligand file to use for determining the docking box. Required if autobox is True.
        """
        self.autobox = autobox
        self.box_ligand = box_ligand

        # Validate autobox settings
        if autobox and box_ligand is None:
            raise ValueError("box_ligand is required when autobox is True.")
        if autobox and box_ligand is not None:
            if not os.path.exists(box_ligand):
                raise FileNotFoundError(f"Box ligand file {box_ligand} does not exist.")

    def _validate_target(self, target: PreparedProtein | str) -> PreparedProtein:
        """Convert target to PreparedProtein if needed."""
        if isinstance(target, PreparedProtein):
            return target

        if not os.path.exists(target):
            raise FileNotFoundError(f"Receptor file {target} does not exist.")

        if not target.endswith((".pdbqt", ".pdb")):
            raise ValueError(f"Receptor file {target} must be in PDBQT or PDB format.")

        return PreparedProtein(target)

    def _validate_ligand(
        self,
        ligand: Iterable[PreparedLigand | str] | PreparedLigand | str,
    ) -> tuple[list[PreparedLigand] | PreparedLigand, bool]:
        """Convert ligand(s) to PreparedLigand if needed.

        Returns:
            Tuple of (prepared ligand(s), is_single).
        """
        if isinstance(ligand, (PreparedLigand, str)):
            return self._validate_single_ligand(ligand), True

        prepared = [self._validate_single_ligand(lig) for lig in ligand]
        return prepared, False

    def _validate_single_ligand(self, ligand: PreparedLigand | str) -> PreparedLigand:
        """Convert a single ligand to PreparedLigand if needed."""
        if isinstance(ligand, PreparedLigand):
            return ligand

        if not os.path.exists(ligand):
            raise FileNotFoundError(f"Ligand file {ligand} does not exist.")

        if not ligand.endswith((".smi", ".mol", ".sdf")):
            raise ValueError(
                f"Ligand file {ligand} must be in smi, mol, or sdf format."
            )

        return PreparedLigand(ligand)

    def _validate_context(self, context: GninaConfig | str) -> GninaConfig:
        """Convert context to GninaConfig if needed."""
        if isinstance(context, GninaConfig):
            return context

        return GninaConfig(context)

    def _update_config_paths(
        self,
        target: PreparedProtein,
        ligand: PreparedLigand,
        context: GninaConfig,
    ) -> None:
        """Update the GNINA config with receptor, ligand, and autobox paths.

        Warns if existing config values differ from the provided inputs.
        """
        gnina_cfg = context.data_dict.get("gnina", {})
        input_cfg = gnina_cfg.get("input", {})

        # Get new paths
        new_receptor = target.file_path
        if isinstance(ligand, PreparedLigand):
            new_ligand = ligand.file_path
        else:
            new_ligand = ligand[0].file_path if ligand else None

        # Check for mismatches and warn
        existing_receptor = input_cfg.get("receptor")
        if existing_receptor is not None and existing_receptor != new_receptor:
            logger.warning(
                f"Config receptor '{existing_receptor}' differs from provided target '{new_receptor}'. "
                "Overwriting with provided value."
            )

        existing_ligand = input_cfg.get("ligand")
        if existing_ligand is not None and existing_ligand != new_ligand:
            logger.warning(
                f"Config ligand '{existing_ligand}' differs from provided ligand '{new_ligand}'. "
                "Overwriting with provided value."
            )

        if self.autobox:
            existing_autobox = input_cfg.get("autobox_ligand")
            if existing_autobox is not None and existing_autobox != self.box_ligand:
                logger.warning(
                    f"Config autobox_ligand '{existing_autobox}' differs from provided box_ligand '{self.box_ligand}'. "
                    "Overwriting with provided value."
                )

        # Update values
        input_cfg["receptor"] = new_receptor
        input_cfg["ligand"] = new_ligand

        if self.autobox:
            input_cfg["autobox_ligand"] = self.box_ligand

        # Update output paths based on ligand name if they are still defaults
        output_cfg = gnina_cfg.get("output", {})
        ligand_dir = os.path.dirname(new_ligand) if new_ligand else "."
        ligand_name = (
            os.path.splitext(os.path.basename(new_ligand))[0]
            if new_ligand
            else "docked"
        )

        # Only update if using default values
        if output_cfg.get("out") == "docked.sdf.gz":
            output_cfg["out"] = os.path.join(ligand_dir, f"{ligand_name}_docked.sdf.gz")

        # Persist changes
        context.write_config(context.data_dict)

    def _build_command(self, context: GninaConfig) -> str:
        """Build the GNINA command string from the config.

        Args:
            context: Validated GninaConfig object.

        Returns:
            Complete GNINA command string ready for SLURM script.
        """
        args = context.to_cli(context.data_dict.get("gnina", {}))
        return "gnina " + " ".join(args)

    @override
    def run(
        self,
        target: PreparedProtein | str,
        ligand: Iterable[PreparedLigand | str] | PreparedLigand | str,
        context: GninaConfig | str,
    ) -> Iterable[str] | str:
        """Dock one or multiple ligands to a single target.

        Args:
            ta: PreparedProteinrget
                Prepared protein or path to receptor file (PDBQT or PDB format).
            ligand
                Prepared ligand(s) or path(s) to ligand file(s) (smi, mol, or sdf format).
            context
                Configuration object or path to GNINA config YAML file.

        Returns:
            Docked ligand(s) with poses and scores.
        """
        # Convert inputs to objects
        validated_target = self._validate_target(target)
        validated_ligands, is_single = self._validate_ligand(ligand)
        validated_context = self._validate_context(context)

        if is_single:
            self._update_config_paths(
                validated_target, validated_ligands, validated_context
            )
            return self._build_command(validated_context)

        # Multiple ligands - generate a command for each
        commands = []
        for lig in validated_ligands:
            self._update_config_paths(validated_target, lig, validated_context)
            commands.append(self._build_command(validated_context))

        return commands
