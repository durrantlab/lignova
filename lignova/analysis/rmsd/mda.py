r"""Implementation of RMSD analysis using MDAnalysis."""

from typing import override

from MDAnalysis.analysis.rms import rmsd

from ...structure.editing import filter_hetatoms, get_mda_universe
from .base import RMSDBase


class mdaRMSD(RMSDBase):
    r"""Calculate ligand RMSD using MDAnalysis."""

    @override
    def calculate(
        self,
        selection: str = "(not resname HOH) and (not name H*)",
        superimpose: bool = False,
        save: bool = False,
        output_filename: str | None = None,
    ) -> list[float]:
        r"""Calculate RMSD between docked and reference ligand using MDAnalysis.

        Args:
            selection: MDAnalysis atom selection string.
                Default excludes waters and hydrogens.
            superimpose: If True, optimally superimpose the ligand onto the
                reference before computing RMSD (minimized RMSD). If False,
                compute in-place RMSD using the docked coordinates as-is.
                Default is False.
            save: If True, write results to a text file. Default is False.
            output_filename: Output file path (without extension) if save
                is True.

        Returns:
            List of RMSD values, one per frame in the trajectory.
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

        docked_traj = get_mda_universe(self.target.file_path)
        reference = get_mda_universe(self.reference.file_path)

        ref_ligand = filter_hetatoms(reference).select_atoms(selection)
        dock_ligand = filter_hetatoms(docked_traj).select_atoms(selection)

        if len(ref_ligand.atoms) != len(dock_ligand.atoms):
            raise ValueError(
                f"Atom count mismatch: reference has {len(ref_ligand.atoms)}, "
                f"target has {len(dock_ligand.atoms)}. "
                f"Files: {self.reference.file_path}, {self.target.file_path}"
            )

        rmsd_values = []
        for _ in docked_traj.trajectory:
            rmsd_value = rmsd(
                dock_ligand.positions,
                ref_ligand.positions,
                center=superimpose,
                superposition=superimpose,
            )
            rmsd_values.append(float(rmsd_value))

        if save:
            self._save_result(rmsd_values, output_filename)

        return rmsd_values
