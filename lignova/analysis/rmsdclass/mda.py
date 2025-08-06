r"""Implementation of rmsd analysis using MDAnalysis."""

from typing import override

import MDAnalysis as mda
from loguru import logger
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

from ...structure.editing import filter_hetatoms, get_mda_universe
from .base import RMSDBase


class mdaRMSD(RMSDBase):
    r"""Class to calculate ligand rmsd using MDAnalysis."""

    @override
    def calculate(
        self, selection: str = "(not resname HOH) and (not name H*)"
    ) -> list[float]:
        r"""Calculate rmsd between docked ligand and reference ligand using MDAnalysis.

        Args:
            selection : Atom selection language for rmsd calculations.
                Default is (not resname HOH) and (not name H*)".

        Returns:
            List of rmsd values for each frame in the trajectory.
        """
        # check if the extensions are not pdb then raise an error using .file_ext
        if self.target.file_ext != "pdb":
            raise ValueError(
                f"Target file {self.target.file_path} must be in PDB format."
            )
        if self.reference.file_ext != "pdb":
            raise ValueError(
                f"Reference file {self.reference.file_path} must be in PDB format."
            )
        if not self.target.file_path:
            raise ValueError("Target file path is None.")
        docked_traj: mda.Universe = get_mda_universe(self.target.file_path)
        reference: mda.Universe = get_mda_universe(self.reference.file_path)

        # Select the ligands
        ref_ligand_het: mda.Universe = filter_hetatoms(reference)
        dock_ligand_het: mda.Universe = filter_hetatoms(docked_traj)

        ref_ligand: mda.AtomGroup = ref_ligand_het.select_atoms(selection)
        dock_ligand: mda.AtomGroup = dock_ligand_het.select_atoms(selection)
        # Check if the length of the ligands is the same
        if len(ref_ligand.atoms) != len(dock_ligand.atoms):
            logger.error(len(ref_ligand.atoms))
            logger.error(len(dock_ligand.atoms))
            raise ValueError(
                f"Reference ligand {self.reference.file_path} and docked ligand "
                + f"{self.target.file_path} have different number of atoms."
            )
        # Calculate the rmsd between the ligands
        rmsd_values: list[float] = []
        for _ in docked_traj.trajectory:
            rmsd_value: float = rmsd(
                dock_ligand.positions,
                ref_ligand.positions,
                center=True,
            )
            rmsd_values.append(float(rmsd_value))
        return rmsd_values
