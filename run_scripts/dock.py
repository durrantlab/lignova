#!/usr/bin/env python
"""
Run GNINA docking to dock ligands into protein pockets.
"""

import argparse
import sys
import os
from loguru import logger

from lignova.docking.gnina import GNINA
from lignova.yaml.docking_config import GninaConfig

def clean_sdf(input_path: str, output_path: str | None = None) -> str:
    r"""Remove the empty molecule that Gypsum-dl adds 
        at the start of SDF files. Additionally it  removes any 
        molecules with 0 atoms.
    Args:
        input_path: Path to the input SDF file.
        output_path: Path to write the cleaned SDF. If None, overwrites input.
    Returns:
        Path to the cleaned SDF file.
    """
    if output_path is None:
        output_path = input_path
    
    with open(input_path, "r") as f:
        content = f.read()
    
    # Split by molecule delimiter
    molecules = content.split("$$$$\n")
    
    valid_molecules = []
    removed = 0
    
    for mol in molecules:
        if not mol.strip():
            continue
        
        # Check for empty molecule marker: "  0  0  0  0  0  0  0  0  0  0999 V2000"
        if "  0  0  0  0  0  0  0  0  0  0999 V2000" in mol:
            logger.warning("Found empty gypsum metadata molecule, removing.")
            removed += 1
            continue
        
        valid_molecules.append(mol)
    
    if removed > 0:
        logger.info(f"Removed {removed} empty molecule(s) from {input_path}")
        
        # Write back with delimiters
        with open(output_path, "w") as f:
            f.write("$$$$\n".join(valid_molecules))
            if valid_molecules:
                f.write("$$$$\n")
        
        logger.info(f"Wrote {len(valid_molecules)} molecules to {output_path}")
    
    return output_path

def dock(
    receptor: str,
    ligand: str | list[str],
    config_path: str,
    autobox: bool = True,
    box_ligand: str | None = None,
    output_file: str | None = "gnina_commands.txt",
    clean_ligands: bool = True,
) -> None:
    """
    Dock ligand(s) to a receptor using GNINA.

    Args:
        receptor: Path to the receptor file (PDBQT or PDB format).
        ligand: Path to ligand file(s) (smi, mol, or sdf format).
        config_path: Path to the GNINA configuration file. If not exists, default config will be used.
        autobox: Whether to use autobox for docking region. Default is True.
        box_ligand: Path to the ligand file for autobox. Required if autobox is True.
        output_file: Path to write the GNINA command(s). If not provided, prints to stdout. Default is "gnina_commands.txt".
        clean_ligands: Whether to clean the input SDF files by removing empty molecules. Default is True.
    """
    gnina = GNINA(autobox=autobox, box_ligand=box_ligand)
    if clean_ligands:
        if isinstance(ligand, str):
            if ligand.endswith(".sdf"):
                ligand = clean_sdf(ligand)
        else:
            ligand = [clean_sdf(lig) if lig.endswith(".sdf") else lig for lig in ligand]

    commands = gnina.run(
        target=receptor,
        ligand=ligand,
        context=config_path,
    )

    # Handle single command or list of commands
    if isinstance(commands, str):
        commands = [commands]

    if output_file:
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            logger.warning(f"The directory {output_dir} does not exist. Creating it.")
            os.makedirs(output_dir)

        with open(output_file, "w") as f:
            for cmd in commands:
                f.write(f"{cmd}\n")
        logger.info(f"GNINA commands written to {output_file} ({len(commands)} commands)")
    else:
        for cmd in commands:
            print(cmd)


def run_cli():
    """Command line interface for GNINA docking."""
    parser = argparse.ArgumentParser(
        description="Dock ligands into protein pockets using GNINA."
    )
    parser.add_argument(
        "-r", "--receptor",
        type=str,
        required=True,
        help="Path to the receptor file (PDBQT or PDB format).",
    )
    parser.add_argument(
        "-l", "--ligand",
        type=str,
        nargs="+",
        required=True,
        help="Path(s) to the ligand file(s) (smi, mol, or sdf format). Multiple files can be specified.",
    )
    parser.add_argument(
        "-c", "--config",
        type=str,
        required=False,
        default="gnina_config.yaml",
        help="Path to the GNINA configuration file. Default: gnina_config.yaml",
    )
    parser.add_argument(
        "--autobox",
        action="store_true",
        default=True,
        help="Use autobox for docking region (default: True).",
    )
    parser.add_argument(
        "--no-autobox",
        action="store_true",
        help="Disable autobox and use explicit box coordinates from config.",
    )
    parser.add_argument(
        "-b", "--box-ligand",
        type=str,
        required=False,
        help="Path to the ligand file for autobox. Required if autobox is enabled.",
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=False,
        help="Path to write the GNINA command(s). If not provided, prints to stdout.",
    )
    parser.add_argument(
        "--no-clean",
        action="store_true",
        help="Disable automatic removal of empty molecules from SDF files.",
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Handle autobox logic
    autobox = not args.no_autobox

    # Validate box_ligand requirement
    if autobox and args.box_ligand is None:
        parser.error("--box-ligand is required when autobox is enabled. Use --no-autobox to disable.")

    # Handle single vs multiple ligands
    ligand = args.ligand[0] if len(args.ligand) == 1 else args.ligand

    dock(
        receptor=args.receptor,
        ligand=ligand,
        config_path=args.config,
        autobox=autobox,
        box_ligand=args.box_ligand,
        output_file=args.output,
        clean_ligands=not args.no_clean,

    )


if __name__ == "__main__":
    run_cli()