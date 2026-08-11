# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Class for calculating RMSD using Spyrmsd tool."""

import os
from typing import override

import numpy as np
from loguru import logger
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from spyrmsd.rmsd import rmsd as rmsd_fn
from spyrmsd.rmsd import symmrmsd

from lignova.structure.ligand import DockedLigand, Ligand

from ..utils import calc_mcs, load_mol, obabel_convert
from .base import RMSDBase


class spyrmsdRMSD(RMSDBase):
    r"""Calculate RMSD using spyrmsd, accounting for molecular symmetry.

    Note: This is only suitable for small molecules, not proteins.
    """

    def __init__(self, target: DockedLigand, reference: Ligand):
        r"""Initialize the spyrmsdRMSD class.

        Validates that input files exist and are in a supported format.
        Converts to SDF automatically if needed.
        """
        super().__init__(target, reference)

        if not self.reference.file_path or not os.path.exists(self.reference.file_path):
            raise FileNotFoundError(
                f"Reference file {self.reference.file_path} not found."
            )
        if not self.target.file_path or not os.path.exists(self.target.file_path):
            raise FileNotFoundError(f"Target file {self.target.file_path} not found.")

        self._ensure_compatible_format()

    def _ensure_compatible_format(self, target_format: str = "sdf") -> None:
        r"""Convert reference and target to a spyrmsd-compatible format if needed.

        Args:
            target_format: Desired format. Must be sdf or pdb.
        """
        supported = ["sdf", "pdb"]

        if self.reference.file_ext not in supported:
            ref_out = self.reference.file_path.replace(
                self.reference.file_ext, f".{target_format}"
            )
            obabel_convert(self.reference.file_path, ref_out)
            self.reference = self.reference.__class__(ref_out)

        if self.target.file_ext not in supported:
            tgt_out = str(self.target.file_path).replace(
                self.target.file_ext, f".{target_format}"
            )
            obabel_convert(self.target.file_path, tgt_out)
            self.target = self.target.__class__(tgt_out)

    def _load_all_poses(self, path: str, add_hs: bool = False) -> list:
        r"""Load all poses from an SDF, or single mol from PDB/mol2.
        Args:
          path: file path to load
          add_hs: Whether to add hydrogens during loading. Default is False.
        Returns:
          A list of RDKit Mol objects. For SDF input, this will include all poses. For PDB/mol2 input, this will be a single-molecule list.
        """
        if path.endswith(".sdf"):
            suppl = rdmolfiles.SDMolSupplier(path, removeHs=not add_hs, sanitize=False)
            return [m for m in suppl if m is not None]
        mol = load_mol(path, add_hs=add_hs)
        return [mol] if mol is not None else []

    def _run_api(
        self,
        symmetry: bool = True,
        hydrogens: bool = False,
        superimpose: bool = False,
        mcs: bool = False,
    ) -> list[float]:
        r"""In-process spyrmsd via symmrmsd / rmsd."""

        ref_mol = load_mol(self.reference.file_path, add_hs=hydrogens)
        target_poses = self._load_all_poses(self.target.file_path, add_hs=hydrogens)
        if ref_mol is None or not target_poses:
            raise ValueError("Failed to load reference or target.")
        if mcs:
            ref_idx, tgt_idx = calc_mcs(
                self.reference.file_path, self.target.file_path, add_hs=hydrogens
            )
        else:
            ref_n = ref_mol.GetNumAtoms()
            tgt_n = target_poses[0].GetNumAtoms()
            if ref_n != tgt_n:
                raise ValueError(
                    f"Atom count mismatch (ref={ref_n}, target={tgt_n}).Consider using MCS."
                )
            ref_idx = list(range(ref_n))
            tgt_idx = list(range(tgt_n))
        ref_znums = np.array([a.GetAtomicNum() for a in ref_mol.GetAtoms()])[ref_idx]
        ref_coords = ref_mol.GetConformer().GetPositions()[ref_idx]
        ref_adj = Chem.GetAdjacencyMatrix(ref_mol)[np.ix_(ref_idx, ref_idx)]
        first_pose = target_poses[0]
        pose_znums = np.array([a.GetAtomicNum() for a in first_pose.GetAtoms()])[
            tgt_idx
        ]
        pose_adj = Chem.GetAdjacencyMatrix(first_pose)[np.ix_(tgt_idx, tgt_idx)]
        pose_coords_list = [
            pose.GetConformer().GetPositions()[tgt_idx] for pose in target_poses
        ]
        if symmetry:
            values = symmrmsd(
                coordsref=ref_coords,
                coords=pose_coords_list,
                apropsref=ref_znums,
                aprops=pose_znums,
                amref=ref_adj,
                am=pose_adj,
                minimize=superimpose,
            )
        else:
            values = [
                rmsd_fn(ref_coords, pc, ref_znums, pose_znums, minimize=superimpose)
                for pc in pose_coords_list
            ]
        return list(map(float, np.atleast_1d(values)))

    @override
    def calculate(
        self,
        symmetry: bool = True,
        hydrogens: bool = False,
        superimpose: bool = False,
        mcs: bool = False,
        backend: str = "api",
        save: bool = False,
        output_filename: str | None = None,
    ) -> list[float]:
        r"""Calculate RMSD between reference and target using spyrmsd.

        Args:
            symmetry: Use symmetry-corrected RMSD. Default is True.
            hydrogens: Include hydrogens. Default is False.
            superimpose: Superimpose before calculation. Default is False.
                Set to False for in-place RMSD of docked poses.
            mcs: Apply MCS so ref and target are compared over their
                maximum common substructure. Default is False. Only
                supported with backend='api'; combining mcs=True with
                backend='cli' falls back to the API path.
            backend: api (in-process via spyrmsd library) or cli
                (subprocess to python -m spyrmsd). Default is api.
            save: If True, write the result to a text file. Default is False.
            output_filename: Output file path (without extension) if save
                is True.
        Returns:
            List of RMSD values.
        """
        if backend == "cli" and mcs:
            logger.warning(
                "MCS is not supported in CLI mode; falling back to API backend."
            )
            backend = "api"
        if backend == "api":
            values = self._run_api(symmetry, hydrogens, superimpose, mcs)
            if save:
                self._save_result(values, output_filename)
            return values
        command = ["python3", "-m", "spyrmsd"]
        if not symmetry:
            command.append("-n")
        if hydrogens:
            command.append("--hydrogens")
        if superimpose:
            command.append("-m")
        command.extend([self.reference.file_path, self.target.file_path])
        result = self._run_command(command)
        values = [
            float(line.strip()) for line in result.strip().splitlines() if line.strip()
        ]
        if save:
            self._save_result(values, output_filename)
        return values
