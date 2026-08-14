# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Implementation of utility functions for the analysis module."""

import math
import os
import subprocess
from typing import TextIO

import numpy as np
import pandas as pd
from loguru import logger
from rdkit import Chem
from rdkit.Chem import rdchem, rdFMCS, rdmolfiles, rdmolops

from lignova.io import decompress, get_file_ext

STANDARD_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "HIE",
    "HID",
    "HIP",
    "CYX",
    "ASH",
    "GLH",
}


def _is_protein(mol: rdchem.Mol) -> bool:
    r"""Check whether a molecule is a protein based on its PDB residue names.

    Returns True if the majority of unique residue names are standard
    amino acids (or common variants like HIE/HID/CYX).
    """
    residue_names = set()
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is not None:
            residue_names.add(info.GetResidueName().strip())
    if not residue_names:
        return False
    return len(residue_names & STANDARD_RESIDUES) / len(residue_names) > 0.5


def _set_pdb_record_type(mol: rdchem.Mol) -> None:
    r"""Set ATOM/HETATM record types based on residue identity.

    Standard amino acid residues are written as ATOM records,
    everything else (ligands, waters, ions, cofactors) as HETATM.
    Args:
        mol: RDKit molecule with PDB residue info.
    """
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is not None:
            resname = info.GetResidueName().strip()
            info.SetIsHeteroAtom(resname not in STANDARD_RESIDUES)


def mae_convert(
    input_file: str,
    output_file: str,
    concatenate: bool = True,
    remove_hs: bool = True,
    sanitize: bool = False,
    protein: bool = False,
) -> list[str]:
    r"""Convert a .mae file to SDF or PDB using maeparser from RDKit.
    Args:
        input_file: Path to the input .mae or .maegz file.
        output_file: Path to the output file (.sdf or .pdb).
        concatenate: Whether to concatenate all molecules into a single file
            (only applies to SDF). Default is True.
        remove_hs: Whether to remove hydrogens. Default is True.
        sanitize: Whether to sanitize molecules. Default is False.
            Set to False for docked poses to preserve coordinates.
        protein: Whether to output the protein if exists to a separate file with
            a _protein suffix. Default is False.
    Returns:
        List of output file paths written.
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    ext_in = get_file_ext(input_file).lower()
    if ext_in not in [".mae", ".maegz"]:
        raise ValueError(f"Invalid input format: {ext_in}. Expected .mae or .maegz")

    ext_out = get_file_ext(output_file).lower()
    if ext_out not in [".sdf", ".pdb"]:
        raise ValueError(f"Invalid output format: {ext_out}. Expected .sdf or .pdb")

    base, ext = os.path.splitext(output_file)
    output_paths = []

    with decompress(input_file) as mae_path:
        suppl = rdmolfiles.MaeMolSupplier(mae_path, removeHs=False, sanitize=sanitize)

        protein_mol = None
        ligands = []

        for i, mol in enumerate(suppl):
            if mol is None:
                logger.warning(
                    f"Molecule at index {i} in {input_file} could not be loaded."
                )
                continue

            if remove_hs:
                mol = rdmolops.RemoveAllHs(mol, sanitize=False)

            if protein_mol is None and _is_protein(mol):
                protein_mol = mol
                logger.info(
                    f"Detected protein at index {i} with ({mol.GetNumAtoms()} atoms).Removing it."
                )
                continue
            ligands.append(mol)

        if not ligands:
            raise ValueError(f"No valid ligand molecules found in {input_file}")

        logger.info(f"Loaded {len(ligands)} ligands from {input_file}")

        if protein_mol is not None and protein:
            protein_path = f"{base}_protein.pdb"
            _set_pdb_record_type(protein_mol)
            rdmolfiles.MolToPDBFile(protein_mol, protein_path)
            output_paths.append(protein_path)
            logger.info(f"Wrote protein to {protein_path}")

        if ext_out == ".sdf":
            if concatenate:
                writer = rdmolfiles.SDWriter(output_file)
                for mol in ligands:
                    writer.write(mol)
                writer.close()
                output_paths.append(output_file)
                logger.info(f"Wrote {len(ligands)} ligand(s) to {output_file}")
            else:
                for i, mol in enumerate(ligands):
                    path = f"{base}_{i}{ext}"
                    rdmolfiles.MolToMolFile(mol, path)
                    output_paths.append(path)
                logger.info(f"Wrote {len(ligands)} SDF file(s)")

        elif ext_out == ".pdb":
            for mol in ligands:
                _set_pdb_record_type(mol)
            if len(ligands) == 1:
                rdmolfiles.MolToPDBFile(ligands[0], output_file)
                output_paths.append(output_file)
                logger.info(f"Wrote 1 ligand to {output_file}")
            else:
                for i, mol in enumerate(ligands):
                    path = f"{base}_{i}{ext}"
                    rdmolfiles.MolToPDBFile(mol, path)
                    output_paths.append(path)
                logger.info(f"Wrote {len(ligands)} PDB files")
    return output_paths


def obabel_convert(
    test_file: str | TextIO, output_filename: str, hydrogen: bool = False
) -> None:
    r"""Interconvert ligand(s) file formats using obabel.

    Args:
        test_file : Test file name or file object.
        output_filename : Output file name.
        hydrogen: Whether to add hydrogens to the molecule. Default is False.
    """
    # check if the file exists
    if isinstance(test_file, str):
        if not os.path.exists(test_file):
            raise FileNotFoundError(f"Input file not found: {test_file}")
    if os.path.exists(output_filename):
        logger.warning(
            f"Output file {output_filename} already exists. Skipping conversion."
        )
        return
    # check if the file extension is .mae or .maegz and if so use the mae_convert function instead
    ext = os.path.splitext(str(test_file))[-1].lower()
    if ext in [".mae", ".maegz"]:
        logger.info(
            f"Input file {test_file} is a .mae or .maegz file. Using mae_convert instead of obabel."
        )
        mae_convert(str(test_file), output_filename, remove_hs=not hydrogen)
        return
    # Construct the command
    command: list[str] = ["obabel", test_file, "-O", output_filename]
    logger.info(f"Running command: {' '.join(command)}")
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
        elif path.endswith(".sdf"):
            suppl = rdmolfiles.SDMolSupplier(path, removeHs=not add_hs, sanitize=False)
            if len(suppl) > 0:
                logger.warning(
                    f"Multiple molecules found in this sdf file {path}.Loading only the first one."
                )
                return suppl[0]
    except Exception as e:
        logger.warning(f"Failed to load molecule from {path}: {e}")
    return None


def write_sdf(path: str | list[str], remove_hs: bool = True) -> str | list[str]:
    r"""Write a molecule or list of molecules to an SDF file, optionally removing hydrogens.
    Args:
        path: Path to the imput file (PDB or MOL2) or list of paths.
        remove_hs: Whether to remove hydrogens from the molecule(s) before writing. Default is True.
    Returns:
        Path(s) to the output SDF file(s).
    """
    if isinstance(path, str):
        path = [path]
    output_paths = []
    for p in path:
        mol = load_mol(p, add_hs=not remove_hs)
        if mol is None:
            logger.warning(f"Failed to load molecule from {p}. Skipping.")
            continue
        if remove_hs:
            mol = rdmolops.RemoveAllHs(mol, sanitize=False)
        output_path = os.path.splitext(p)[0] + ".sdf"
        rdmolfiles.MolToMolFile(mol, output_path)
        output_paths.append(output_path)
    return output_paths if len(output_paths) > 1 else output_paths[0]


def calc_mcs(
    ref_file: str,
    target_file: str,
    timeout: int = 10,
    add_hs: bool = False,
) -> tuple[list[int], list[int]]:
    """Find the MCS between two molecule files. Returns aligned atom indices.
    Args:
        ref_file : Path to the reference molecule file (PDB, MOL2, or SDF).
        target_file : Path to the target molecule file (PDB, MOL2, or SDF).
        timeout : Maximum time (in seconds) to spend on MCS calculation. Default is 10 seconds.
        add_hs : Whether to add hydrogens when loading molecules. Default is False.
    Returns:
        A tuple of two lists: (ref_indices, target_indices) corresponding to the MCS atom indices in the reference and target molecules, respectively.
    """
    ref = load_mol(ref_file, add_hs=add_hs)
    target = load_mol(target_file, add_hs=add_hs)
    if ref is None or target is None:
        raise ValueError("Failed to parse one or both input files.")
    result = rdFMCS.FindMCS(
        [ref, target],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        timeout=timeout,
    )
    if result.numAtoms == 0:
        raise ValueError("No common substructure found.")

    query = Chem.MolFromSmarts(result.smartsString)
    return (
        list(ref.GetSubstructMatch(query)),
        list(target.GetSubstructMatch(query)),
    )


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
        if os.path.exists(temp_file) and temp_file != current_path:
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


def eval_pose(
    docked_pose: str,
    reference_pose: str | None,
    protein: str,
    outfmt: str = "csv",
    multiple: bool = True,
    output_file: str | None = str,
    top_n: int = None,
    max_workers: int = None,
) -> str | pd.DataFrame:
    """
    Run posebuster (bust) command line tool to evaluate predicted pose
        whether the poses are physically valid
    Args:
        docked_pose : Path to the docked pose file.
        reference_pose : Path to the reference pose file.
        protein : Path to the protein file.
        outfmt : Output format for bust. Default is "short". Ir can be "short", "long", or "csv".
        multiple : Whether to run bust on multiple poses. Default is True.
            If True, docked_pose should be a csv file and reference_pose should be None
        output_file : Path to save the output file. If None, will return
        top_n : Number of top results to return. Default is None (return all).
        max_workers : Maximum number of workers to use for parallel processing. Default is None (use all available).
    Returns:
        Path to the output file or a pandas DataFrame with the results.
    """
    if multiple:
        # check that docked_pose is a csv file and reference_pose is None if true cool if not error
        if not isinstance(docked_pose, str) or not docked_pose.endswith(".csv"):
            raise ValueError("When multiple is True, docked_pose must be a csv file.")
        if reference_pose is not None:
            logger.warning("When multiple is True, reference_pose would be ignored. ")
        if protein is not None:
            logger.warning("When multiple is True, protein would be ignored. ")
        command = [
            "bust",
            "-t",
            docked_pose,
        ]
    else:
        # check if any of the three inputs don't exist or are empty error out
        if (
            reference_pose is None
            or not os.path.exists(reference_pose)
            or not os.path.exists(docked_pose)
            or not os.path.exists(protein)
        ):
            raise ValueError(
                "When multiple is False, all three inputs must be valid file paths."
            )
        command = [
            "bust",
            docked_pose,
            "-l",
            reference_pose,
            "-p",
            protein,
        ]
    command.extend(["-outfmt", outfmt])
    if output_file is not None:
        if not output_file.endswith(".csv"):
            output_file += ".csv"
        command.extend(["-o", output_file])
    if top_n is not None:
        command.extend(["--top-n", str(top_n)])
    if max_workers is not None:
        command.extend(["--max-workers", str(max_workers)])
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
            logger.info("Pose evaluation completed")
            logger.info(f"Output:\n{stdout}")
        else:
            logger.error("Pose evaluation failed ")
            logger.error(f"Error Output:\n{stderr}")
            raise subprocess.CalledProcessError(process.returncode, " ".join(command))
    except Exception as e:
        logger.error(f"An error occurred during pose evaluation: {str(e)}")
        raise e


_R_KCAL = 1.98721e-3

_KD_UNITS = {"M": 1e0, "mM": 1e3, "uM": 1e6, "nM": 1e9, "pM": 1e12}


def to_delta_g(
    affinity: float | np.ndarray, temperature: float = 300.0
) -> float | np.ndarray:
    """Convert CNNaffinity (pK) to ΔG in kcal/mol.

    Args:
        affinity: CNNaffinity value(s) in pK units.
        temperature: Temperature in Kelvin. Defaults to 300K (gnina default).
    Returns:
        A numpy array or float with the ΔG value(s) in kcal/mol.
    """
    result = (-_R_KCAL) * temperature * math.log(10.0) * np.asarray(affinity)
    if np.ndim(result) == 0:
        return float(result)
    return result


def to_kd(affinity: float | np.ndarray, unit: str = "uM") -> float | np.ndarray:
    """Convert CNNaffinity (pK) to Kd.

    Args:
        affinity: CNNaffinity value(s) in pK units.
        unit: Concentration unit. One of M, mM, uM, nM, pM. Defaults to uM.
    Returns:
        A numpy array or float with the Kd value(s) in the specified unit.
    """
    if unit not in _KD_UNITS:
        raise ValueError(f"Unknown unit {unit!r}. Choose from: {list(_KD_UNITS)}")
    result = np.power(10.0, -np.asarray(affinity)) * _KD_UNITS[unit]
    if np.ndim(result) == 0:
        return float(result)
    return result


def to_pActivity(value: float | np.ndarray, unit: str = "uM") -> float | np.ndarray:
    """Convert experimental affinity to pActivity (-log10(Molar)).

    Args:
        value: Kd value(s) in the specified unit.
        unit: One of M, mM, uM, nM, pM. Defaults to uM.
    Returns:
        A numpy array or float with the pActivity value(s) of the experimental data.
    """
    if unit not in _KD_UNITS:
        raise ValueError(f"Unknown unit {unit!r}. Choose from: {list(_KD_UNITS)}")
    array = np.asarray(value, np.float64)
    molar_value = array / _KD_UNITS[unit]
    pAct = -np.log10(molar_value)
    if np.ndim(pAct) == 0:
        return float(pAct)
    return pAct
