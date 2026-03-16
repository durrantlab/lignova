r"""Implementation for docking  GNINA docking."""

import os
import re
import shutil
import subprocess
from collections.abc import Iterable
from typing import override

from loguru import logger

from lignova.docking.docking import Docking
from lignova.structure.ligand import DockedLigand, PreparedLigand
from lignova.structure.protein import PreparedProtein
from lignova.yaml.docking_config import GninaConfig


class GNINA(Docking):
    """Class to dock ligands in determined protein pocket using GNINA."""

    _VALID_REPAIRS = {
        "bonds_hydrogens",
        "bonds",
        "hydrogens",
        "checkhydrogens",
        "None",
    }

    _VALID_CLEANUPS = {
        "nphs",
        "lps",
        "waters",
        "nonstdres",
        "deleteAltB",
        "nphs_lps",
        "nphs_lps_waters",
        "nphs_lps_waters_nonstdres",
    }

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

    def _fix_pqr_spacing(
        self, pqr_path: str, cleanup: str = "nphs_lps_waters_nonstdres"
    ) -> str:
        r"""PQR files from PDB2PQR can have columns that
        run together when coordinates or charges are large negative numbers
        (e.g. ``-18.713-100.168``).MGLTools uses a strict fixed-width
        PDB parser that fails on these lines. This rewrites
        every ATOM/HETATM line with proper column widths.

        Additionally, HETATM lines that would be removed by
        prepare_receptor4.py's -U cleanup flag are stripped early,
        since MGLTools often cannot parse them from PQR files.

        The original file is backed up to ``*_org.pqr`` before overwriting.

        Args:
            pqr_path
                Path to the PQR file to fix.
            cleanup
                The cleanup mode passed to prepare_receptor4.py (-U flag).
                Used to determine which HETATM lines to strip.

        Returns:
            Path to the (overwritten) fixed PQR file."""
        backup_path = pqr_path.replace(".pqr", "_org.pqr")
        if not os.path.exists(backup_path):
            shutil.copy2(pqr_path, backup_path)
            logger.info(f"Backed up original PQR to {backup_path}")

        # Determine what to strip based on the cleanup flag
        strip_waters = "waters" in cleanup
        strip_nonstdres = "nonstdres" in cleanup

        _WATER_RESNAMES = {"HOH", "WAT", "TIP", "TIP3", "SOL"}

        float_pattern = re.compile(r"-?\d+\.\d+")

        with open(pqr_path) as fh:
            lines = fh.read().splitlines()

        removed_waters = 0
        removed_hetatm = 0
        fixed = []
        for line in lines:
            if line.startswith("HETATM"):
                tokens = line.split()
                is_water = any(t in _WATER_RESNAMES for t in tokens[1:6])

                if is_water and strip_waters:
                    removed_waters += 1
                    continue
                if not is_water and strip_nonstdres:
                    removed_hetatm += 1
                    continue

            if not line.startswith(("ATOM", "HETATM")):
                fixed.append(line)
                continue

            floats = float_pattern.findall(line)
            if len(floats) < 5:
                fixed.append(line)
                continue

            x, y, z = float(floats[-5]), float(floats[-4]), float(floats[-3])
            charge, radius = float(floats[-2]), float(floats[-1])

            # Locate where the numeric tail begins to preserve the prefix
            tail_parts = [re.escape(f) for f in floats[-5:]]
            tail_match = re.search(r"\s*".join(tail_parts), line)
            if not tail_match:
                fixed.append(line)
                continue

            prefix = line[: tail_match.start()].rstrip()

            if line.startswith("HETATM"):
                record = "HETATM"
                rest = prefix[6:]
            else:
                record = "ATOM  "
                rest = prefix[4:]

            tokens = rest.split()
            if len(tokens) < 4:
                fixed.append(line)
                continue

            serial = int(tokens[0])
            atom_name = tokens[1]
            resname = tokens[2]

            if len(tokens) >= 5:
                chain = tokens[3]
                resseq = tokens[4]
            elif len(tokens) == 4:
                field = tokens[3]
                if (
                    len(field) >= 2
                    and field[0].isalpha()
                    and field[1:].lstrip("-").isdigit()
                ):
                    chain = field[0]
                    resseq = field[1:]
                else:
                    chain = " "
                    resseq = field
            else:
                chain = " "
                resseq = "0"

            resseq_int = int(resseq)

            if len(atom_name) >= 4 or atom_name[0].isdigit():
                atom_name_fmt = f"{atom_name:<4s}"
            else:
                atom_name_fmt = f" {atom_name:<3s}"

            fixed.append(
                f"{record}{serial:>5d} "
                f"{atom_name_fmt}"
                f" "
                f"{resname:>3s}"
                f" "
                f"{chain:1s}"
                f"{resseq_int:>4d}"
                f"    "
                f"{x:8.3f}"
                f"{y:8.3f}"
                f"{z:8.3f}"
                f"{charge:8.4f}"
                f"{radius:7.4f}"
            )

        with open(pqr_path, "w") as fh:
            fh.write("\n".join(fixed) + "\n")

        if removed_waters or removed_hetatm:
            logger.info(
                f"Removed {removed_waters} water and {removed_hetatm} other "
                f"HETATM lines from {pqr_path}"
            )
        logger.info(f"Fixed PQR spacing written to {pqr_path}")
        return pqr_path

    def _prepare_protein(
        self,
        pqr_path: str,
        repair: str = "checkhydrogens",
        cleanup: str = "nphs_lps_waters_nonstdres",
        preserve_charges: bool = True,
    ) -> str:
        r"""Convert a PQR file to PDBQT using prepare_receptor4.py.

        Args:
            pqr_path
                Path to the input PQR file.
            repair
                Repair mode (-A flag). Default is checkhydrogens
                to add missing hydrogens.
            cleanup
                Cleanup mode (-U flag). Default is nphs_lps_waters_nonstdres to remove
                non-polar hydrogens, lone pairs, waters, and non-standard residues.
            preserve_charges
                Pass -C to keep the charges assigned by pdb2pqr rather than
                recalculating Gasteiger charges. Default is True.

        Returns:
            Path to the generated PDBQT file.
        """
        if not os.path.exists(pqr_path):
            raise FileNotFoundError(f"PQR file not found: {pqr_path}")
        if not pqr_path.endswith(".pqr"):
            raise ValueError(f"Expected a .pqr file, got '{pqr_path}'.")
        if repair not in self._VALID_REPAIRS:
            raise ValueError(
                f"Invalid repair mode '{repair}'. Must be one of {sorted(self._VALID_REPAIRS)}."
            )
        if cleanup not in self._VALID_CLEANUPS:
            raise ValueError(
                f"Invalid cleanup mode '{cleanup}'. Must be one of {sorted(self._VALID_CLEANUPS)}."
            )
        if not isinstance(preserve_charges, bool):
            raise TypeError(
                f"preserve_charges must be a boolean, got {type(preserve_charges)}"
            )

        out_dir = os.path.dirname(os.path.abspath(pqr_path))
        basename = os.path.splitext(os.path.basename(pqr_path))[0]
        pdbqt_path = os.path.join(out_dir, f"{basename}.pdbqt")

        # pythonsh is mgltools' own Python 2.7 wrapper.
        pythonsh = shutil.which("pythonsh")
        if pythonsh is None:
            raise RuntimeError("pythonsh not found in PATH. ")

        script = shutil.which("prepare_receptor4.py")
        if script is None:
            raise RuntimeError("prepare_receptor4.py not found in PATH.")

        cmd = [
            pythonsh,
            script,
            "-r",
            pqr_path,
            "-o",
            pdbqt_path,
            "-A",
            repair,
            "-U",
            cleanup,
        ]
        if preserve_charges:
            cmd.append("-C")

        logger.warning(
            f"Receptor is not in PDBQT format. Converting {pqr_path} to '{pdbqt_path}'"
        )
        logger.info(f"Running: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(
                f"prepare_receptor4.py failed on first attempt:\n{result.stderr}"
            )
            logger.info("Attempting to fix PQR column spacing and retrying...")
            self._fix_pqr_spacing(pqr_path, cleanup=cleanup)
            if os.path.exists(pdbqt_path):
                os.remove(pdbqt_path)
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"prepare_receptor4.py stderr:\n{result.stderr}")
                raise RuntimeError(
                    f"Even after spacing fix, prepare_receptor4.py failed with exit code {result.returncode}"
                )

        if result.stdout:
            logger.debug(f"prepare_receptor4.py stdout:\n{result.stdout}")

        if not os.path.exists(pdbqt_path):
            raise RuntimeError(
                f"prepare_receptor4.py completed but output not found at {pdbqt_path}"
            )

        logger.info(f"PDBQT receptor written to {pdbqt_path}")
        return pdbqt_path

    def _validate_target(self, target: PreparedProtein | str) -> PreparedProtein:
        """Convert target to PreparedProtein if needed."""
        if isinstance(target, PreparedProtein):
            path = target.file_path
        else:
            path = target
        if not os.path.exists(path):
            raise FileNotFoundError(f"Receptor file {path} does not exist.")

        if not path.endswith((".pdbqt", ".pdb", ".pqr")):
            raise ValueError(
                f"Receptor file {path} must be in PDBQT or PDB or PQR format."
            )

        if path.endswith(".pqr"):
            target = self._prepare_protein(path)
            return PreparedProtein(target)
        return target if isinstance(target, PreparedProtein) else PreparedProtein(path)

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
            target: PreparedProtein
                Prepared protein or path to receptor file (PDBQT or PDB format).
            ligand:
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
