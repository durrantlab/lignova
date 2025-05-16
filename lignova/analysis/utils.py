"""Implementation of utility functions for the analysis module."""

from typing import TextIO

import os
import subprocess
import pandas as pd
from collections.abc import Iterable
from loguru import logger
from rdkit.Chem import rdchem, rdmolfiles, rdmolops

from ..docking.contexts import GlideContext


def interconvert_mae_sdf(
    test_file: str | TextIO,
    output_filename: str,
    ntruct: int | str | None = None,
    context: GlideContext = GlideContext.get_current(),
) -> None:
    r"""Convert ligand(s) to SDF format.

    Args:
        test_file : Test file name or file object.
        output_filename : Output file name.
        ntruct : Number of structures to convert. Default 1:5 i.e the first 5 structures.
        context : Docking context to run the command.
    """
    # GET THE path of the file extension using the os.path.splitext() function
    # if the file extension is .sdf, then the file is in SDF format
    # if the file extension is .mae, then the file is in MAE format
    # if the file extension is neither .sdf nor .mae, then the file format is not supported
    if not os.path.exists(test_file):
        logger.error(f"File {test_file} does not exist.")
        return
    logger.info(os.path.splitext(test_file)[1])
    if os.path.splitext(test_file)[1] == ".sdf":
        logger.info("Input file is in SDF format.Converting to MAE format.")
        fileformat = "-isdf"
        outformat = "-omae"
    elif (
        os.path.splitext(test_file)[1] == ".maegz"
        or os.path.splitext(test_file)[1] == ".mae"
    ):
        logger.info("Input file is in MAE format.Converting to SDF format.")
        fileformat = "-imae"
        outformat = "-osdf"
    else:
        logger.error(
            "Input file format not supported. Please provide a file in MAE or SDF format."
        )
        return

    command = [
        context.command + "/utilities/sdconvert",
        fileformat,
        test_file,
        outformat,
        output_filename,
    ]
    if ntruct is not None:
        command.extend(["-n", str(ntruct)])
    else:
        command.extend(["-all"])
    logger.info(f"Running command: {' '.join(command)}")
    try:
        with subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        ) as process:
            stdout, stderr = process.communicate()
        if process.returncode == 0:
            logger.info("File format conversion completed")
            logger.info(f"Output:\n{stdout}")
        else:
            logger.error("File format conversion failed ")
            logger.error(f"Error Output:\n{stderr}")
            raise subprocess.CalledProcessError(process.returncode, " ".join(command))
    except Exception as e:
        logger.error(f"An error occurred during rmsd calculation: {str(e)}")
        raise e


def obabel_convert(
    test_file: str | TextIO, output_filename: str, hydrogen: bool = False
) -> None:
    """Convert ligand(s) from MAE format to SDF format using obabel.

    Args:
        test_file : Test file name or file object.
        output_filename : Output file name.
        hydrogen: Whether to add hydrogens to the molecule. Default is False.
    """

    # Construct the command
    command: list[str] = ["obabel", test_file, "-O", output_filename]
    if hydrogen:
        command.append("-h")
    try:
        with subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        ) as process:
            stdout, stderr = process.communicate()
        if process.returncode == 0:
            logger.info("File format conversion completed")
            logger.info(f"Output:\n{stdout}")
        else:
            logger.error("File format conversion failed ")
            logger.error(f"Error Output:\n{stderr}")
            raise subprocess.CalledProcessError(process.returncode, " ".join(command))
    except Exception as e:
        logger.error(f"An error occurred during file format conversion: {str(e)}")
        raise e


def obabel_result_parser(output):
    """
    Parses the output from the obabel command and returns the numeric values found per line.

    Args:
        output : The output from the obabel command.

    Returns:
        A dictionary where the keys are arbitrary numbers (1, 2, 3, ...)
            and the values are lists of numeric values found per line.
    """
    # Split the output into lines
    lines = output.strip().split("\n")

    # Initialize a dictionary to store the numeric values per line
    values = {}

    # Iterate through each line and extract the numeric values
    for i, line in enumerate(lines, start=1):
        parts = line.split(",")
        line_values = []
        for part in parts:
            part = part.strip()
            if part != "inf":
                try:
                    value = float(part)
                    line_values.append(value)
                except ValueError:
                    pass  # Ignore parts that cannot be converted to float
        # Store the numeric values for the current line in the dictionary
        values[i] = line_values

    return values


# pylint: disable=no-member,c-extension-no-member
def fix_valence_and_aromaticity_issues(mol: rdchem.Mol) -> rdchem.Mol:
    """
    Fix valence issues (e.g. 4-bonded neutral N) and unmark non-ring aromatics.

    Args:
        Molecule to fix.

    Returns:
        Fixed molecule (sanitized)
    """
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() and not atom.IsInRing():
            atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        if bond.GetIsAromatic() and not bond.IsInRing():
            bond.SetIsAromatic(False)
            bond.SetBondType(rdchem.BondType.SINGLE)

    mol.UpdatePropertyCache(strict=False)

    problems = rdmolops.DetectChemistryProblems(mol)
    for problem in problems:
        if problem.GetType() == "AtomValenceException":
            atom = mol.GetAtomWithIdx(problem.GetAtomIdx())
            if (
                atom.GetAtomicNum() == 7
                and atom.GetFormalCharge() == 0
                and atom.GetExplicitValence() == 4
            ):
                atom.SetFormalCharge(1)
                atom.SetNoImplicit(True)

    rdmolops.SanitizeMol(mol)
    return mol


def load_mol(path: str, add_hs: bool = False) -> rdchem.Mol | None:
    r"""
    Load a molecule from a file using RDKit.
    Args:
        path : Path to the file.
        add_hs : Whether to add hydrogens to the molecule. Default is False.
    Returns:
        Loaded molecule or None if loading failed.
    """
    try:
        if path.endswith(".pdb"):
            return rdmolfiles.MolFromPDBFile(path, removeHs=not add_hs, sanitize=False)
        elif path.endswith(".mol2"):
            return rdmolfiles.MolFromMol2File(path, removeHs=not add_hs, sanitize=False)
    except Exception as e:
        logger.warning(f"Failed to load molecule from {path}: {e}")
    return None


# pylint: disable=no-member,c-extension-no-member
def clean_and_standardize_file(file_path: str, add_hs: bool = False) -> str:
    """
    Clean and standardize a ligand structure file (.pdb or .mol2) using RDKit.
    Attempts RDKit load, Open Babel fix, nitrogen valence fix, and atom type fix.
    Args:
        file_path : Path to the input file.
        add_hs : Whether to add hydrogens to the molecule. Default is False.
    Returns:
        Path to the cleaned .pdb file.
    """
    ext = os.path.splitext(file_path)[1].lower()
    if ext not in [".mol2", ".pdb"]:
        raise ValueError(f"Unsupported input format: {ext}")
    mol = None
    current_path = file_path
    fallback_files = []

    logger.info(f"Attempting to load file: {current_path}")
    mol = load_mol(current_path)

    # Step 2: Open Babel fallback
    if mol is None:
        logger.warning("Attempting Open Babel fallback.")
        babel_path = os.path.splitext(file_path)[0] + "_babel_fixed.mol2"
        try:
            obabel_convert(file_path, babel_path, hydrogen=add_hs)
            if not os.path.exists(babel_path):
                raise FileNotFoundError(
                    f"Open Babel did not create fallback file: {babel_path}"
                )
            current_path = babel_path
            fallback_files.append(babel_path)
            mol = load_mol(current_path)
        except Exception as e:
            logger.error(f"Open Babel fallback failed: {e}")

    # Step 3: Nitrogen valence and aromaticity fix
    if mol is not None:
        try:
            mol.UpdatePropertyCache(strict=False)
            mol = fix_valence_and_aromaticity_issues(mol)
        except Exception as e:
            logger.warning(f"Valence/aromaticity fix failed: {e}")
            mol = None  # trigger next fallback

    # Step 4: Atom type fixing
    if mol is None and ext == ".mol2":
        logger.warning("Attempting to fix atom types in mol2 file.")
        fixed_path = os.path.splitext(file_path)[0] + "_fixed_types.mol2"
        try:
            fix_mol2_atom_types(file_path, fixed_path)
            current_path = fixed_path
            fallback_files.append(fixed_path)
            mol = load_mol(current_path)
        except Exception as e:
            logger.error(f"Fixing atom types failed: {e}")

    # Final validation
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"Molecule is empty or invalid after fallbacks: {file_path}")

    # Sanitize molecule
    try:
        rdmolops.SanitizeMol(mol)
    except Exception as e:
        logger.error(f"Sanitization failed for {current_path}: {e}")
        raise ValueError(f"Failed to sanitize molecule: {file_path}") from e

    # Write final cleaned output
    output_path = os.path.splitext(file_path)[0] + "_cleaned.pdb"
    rdmolfiles.MolToPDBFile(mol, output_path)

    # Cleanup temp files
    for temp_file in fallback_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)

    return output_path


def split_mol_file(file_path: str, lig_name: str) -> str:
    r"""
    Split a .mol2 file to get the ligand and protein separately.

    Args:
        file_path : Path to the input .mol2 file.
        lig_name : Name of the ligand to extract (residue name in the mol2 file).

    Returns:
        Path to the split ligand file (PDB format).
    """
    mol = rdmolfiles.MolFromMol2File(file_path, sanitize=False)
    if mol is None:
        raise ValueError(f"Failed to load molecule from {file_path}")

    try:
        rdmolops.SanitizeMol(mol)
    except Exception as e:
        raise ValueError(f"Sanitization failed: {e}")

    fragments = rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=False)

    ligand = None
    fallback_ligand = None
    protein_frags = []

    for frag in fragments:
        resnames = set()

        for atom in frag.GetAtoms():
            for key in ("_TriposResidueName", "_TriposSubstName"):
                if atom.HasProp(key):
                    resname = atom.GetProp(key).strip()
                    resnames.add(resname)

        if lig_name in resnames:
            ligand = frag
        elif "UNK0" in resnames:
            fallback_ligand = frag
        else:
            protein_frags.append(frag)

    if ligand is None and fallback_ligand is not None:
        ligand = fallback_ligand
        lig_name = "UNK0"

    if ligand is None:
        raise ValueError(
            f"Ligand with resname '{lig_name}' not found in {file_path} "
            "and fallback 'UNK0' not found either."
        )

    # Write ligand to file
    base = os.path.splitext(os.path.basename(file_path))[0]
    out_dir = os.path.dirname(file_path)
    ligand_path = os.path.join(out_dir, f"{base}_{lig_name}.pdb")
    rdmolfiles.MolToPDBFile(ligand, ligand_path)

    return ligand_path


def guess_atom_type(atom_name: str) -> str:
    r"""Guess the atom type based on the atom name.
        This is a simple heuristic function that maps common atom names in mol2 files
    Args:
        atom_name : Atom name to guess the type for.
    Returns:
        Guessed atom type.
    """
    name = atom_name.strip().upper()
    if name.startswith("C"):
        return "C.3"
    elif name.startswith("O"):
        return "O.2"
    elif name.startswith("N"):
        return "N.3"
    elif name.startswith("S"):
        return "S.2"
    elif name.startswith("H"):
        return "H"


def fix_mol2_atom_types(input_file: str, output_file: str) -> None:
    r"""Fix atom types in a mol2 file.
    Args:
        input_file : Path to the input mol2 file.
        output_file : Path to the output mol2 file.
    """
    with open(input_file, encoding="utf-8") as f:
        lines = f.readlines()

    fixed_lines = []
    inside_atoms = False
    for line in lines:
        if line.strip() == "@<TRIPOS>ATOM":
            inside_atoms = True
            fixed_lines.append(line)
            continue
        elif line.startswith("@<TRIPOS>"):
            inside_atoms = False

        if inside_atoms and len(line.strip()) > 0:
            parts = line.split()
            if len(parts) >= 6:
                atom_type = parts[5]
                if not atom_type[0].isalpha():  # suspect it's invalid
                    parts[5] = guess_atom_type(parts[1])  # Use atom name to guess
                    line = " ".join(parts) + "\n"
        fixed_lines.append(line)

    with open(output_file, "w", encoding="utf-8") as f:
        f.writelines(fixed_lines)
