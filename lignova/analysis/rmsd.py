r"""Implementation for RMSD analysis."""

# pylint: disable=R0801
from typing import Iterable, TextIO

import os
import subprocess

import MDAnalysis as mda
from loguru import logger
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

from ..docking.contexts import GlideContext
from ..docking.utils import convert_to_pdb, manipulate_complexes
from ..structure.editing import filter_hetatoms, find_common_atoms, select_common_atoms
from ..structure.ligand import DockedLigand
from ..structure.protein import Protein


class RMSD:
    r"""Implementation for RMSD analysis."""

    def __init__(
        self,
        ligand: Iterable[DockedLigand] | DockedLigand,
        reference: Iterable[Protein] | Protein,
        context: GlideContext.get_current(),
    ):
        r"""Initialize RMSD class.

        Parameters
        ----------
        ligand: Iterable[DockedLigand] | DockedLigand
            Docked ligand(s) Object that will be analyzed.
        reference: Iterable[Protein] | Protein
            Reference ligand(s) in a Protein object that will be used for comparison.
        context: GlideContext object
            Docking context.
        """
        self.ligand = ligand
        self.reference = reference
        self.context = context

    def rmsd_mda(self, selection: str = "not resname HOH") -> Iterable[float]:
        r"""Calculate RMSD between docked ligand and reference ligand using MDAnalysis.
        Parameters
        ----------
        selection : str
            Atom selection language for RMSD calculations.
            Default is "record_type HETATM and not resname HOH".

        Returns
        -------
        List[float]
            List of RMSD values for each frame in the trajectory.
        """
        # check if the ligand is pdb or maegz
        logger.debug(self.ligand.file_ext)
        logger.debug(self.reference.file_ext)
        if self.ligand.file_ext != "pdb":
            logger.debug(self.ligand.file_name)
            manipulate_complexes(
                self.ligand.file_path,
                context=self.context,
                mode="merge",
                outfile_name=os.path.splitext(self.ligand.file_name)[0]
                + "_complexes.maegz",
            )
            file_name = os.path.splitext(self.ligand.file_name)[0]
            convert_to_pdb(
                os.path.join(self.context.write_dir, file_name + "_complexes.maegz")
            )
            self.ligand = DockedLigand(
                os.path.join(self.context.write_dir, file_name + "_complexes.pdb")
            )
        if self.reference.file_ext != "pdb":
            manipulate_complexes(
                self.reference.file_path,
                context=self.context,
                mode="merge",
                outfile_name=os.path.splitext(self.reference.file_name)[0]
                + "_complexes.maegz",
            )
            convert_to_pdb(
                os.path.join(
                    self.context.write_dir,
                    self.reference.file_name + "_complexes.maegz",
                )
            )
            self.reference = Protein(
                os.path.join(
                    self.context.write_dir, self.reference.file_name + "_complexes.pdb"
                )
            )
        docked_traj = mda.Universe(self.ligand.file_path, format="pdb")
        reference = mda.Universe(self.reference.file_path, format="pdb")

        # Align the trajectories based on the protein
        protein1 = reference.select_atoms("name CA")
        protein2 = docked_traj.select_atoms("name CA")
        align.alignto(protein1, protein2)

        # Select the ligands
        ref_ligand = filter_hetatoms(reference).select_atoms(selection)
        dock_ligand = filter_hetatoms(docked_traj).select_atoms(selection)

        # Check if the length of the ligands is the same
        if len(ref_ligand.atoms) != len(dock_ligand.atoms):
            logger.warning(
                f"Atoms number in the reference ligand ({len(ref_ligand.atoms)}) "
                f"and docked ligand ({len(dock_ligand.atoms)}) different."
            )
            # Find the common atoms
            common_atoms = find_common_atoms(ref_ligand, dock_ligand)
            # Select the common atoms from each ligand
            ref_ligand = select_common_atoms(ref_ligand, common_atoms)
            dock_ligand = select_common_atoms(dock_ligand, common_atoms)

        # Calculate the RMSD between the ligands
        rmsd_values = []
        if len(docked_traj.trajectory) > 1:
            for _ in docked_traj.trajectory:
                rmsd_value = rmsd(dock_ligand.positions, ref_ligand.positions)
                rmsd_values.append(rmsd_value)
        else:
            rmsd_value = rmsd(dock_ligand.positions, ref_ligand.positions)
            rmsd_values.append(rmsd_value)
        logger.debug(rmsd_values)
        return rmsd_values

    def rmsd_obabel(
        self,
        firstonly: bool = True,
        save: bool = False,
        minimize: bool = False,
        output_filename: str | TextIO | None = None,
    ):
        """Calculate RMSD between reference and target file using OpenBabel,
        taking into account the symmetry of the molecules.

        Parameters
        ----------
        firstonly : bool
            Only calculate the RMSD for the first structure in the reference file. Default is True.
        save : bool
            Write the RMSD to a txt file. Default is False.
        minimize : bool
            Compute minimum RMSD. Default is False.
        output_filename : str | TextIO | None
            Output file name if csv is true. Default is None.
        """
        command = ["obrms", "-x"]
        if firstonly:
            command.append("-f")
        if minimize:
            command.append("-m")
        command.extend([self.reference.file_path, self.ligand.file_path])
        try:
            with subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            ) as process:
                stdout, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(f"RMSD calculation completed for {self.ligand.file_id}")
                logger.info(f"Output:\n{stdout}")
            else:
                logger.error(f"RMSD calculation failed for {self.ligand.file_id}")
                logger.error(f"Error Output:\n{stderr}")
                raise subprocess.CalledProcessError(
                    process.returncode, " ".join(command)
                )
        except subprocess.CalledProcessError as e:
            logger.error(f"An error occurred during RMSD calculation: {str(e)}")
            raise e

        # Parse the RMSD from the output
        rmsd_result = stdout.strip()

        # If an output filename is provided, write the RMSD to the file
        if save:
            if output_filename is not None:
                with open(output_filename + ".txt", "w", encoding="utf-8") as file:
                    file.write(f"RMSD: {rmsd_result}\n")
            else:
                raise ValueError("Output filename is required if csv is True")
        return rmsd_result

    def symmetry_rmsd(
        self,
        symmetry: bool = True,
        hydrogens: bool = False,
        superimpose: bool = False,
        save: bool = False,
        output_filename: str | TextIO | None = None,
    ):
        r"""Calculate RMSD between reference and target file using Spyrmsd,
        taking into account the symmetry of the molecules.
        Parameters
        ----------
        symmetry : bool
            Use symmetry information. Default is True.
        hydrogens : bool
            Include hydrogens in the calculation. Default is False.
        superimpose : bool
            Superimpose the molecules. Default is False. (i.e perform in-place RMSD)
        """
        # python -m spyrmsd -m -n [self.reference.file_path, self.ligand.file_path
        command = ["python", "-m", "spyrmsd"]
        if not symmetry:
            command.append("-n")
        if hydrogens:
            command.append("--h")
        if superimpose:
            command.append("-m")
        command.extend([self.reference.file_path, self.ligand.file_path])
        try:
            with subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            ) as process:
                stdout, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(f"RMSD calculation completed for {self.ligand.file_id}")
                logger.info(f"Output:\n{stdout}")
            else:
                logger.error(f"RMSD calculation failed for {self.ligand.file_id}")
                logger.error(f"Error Output:\n{stderr}")
                raise subprocess.CalledProcessError(
                    process.returncode, " ".join(command)
                )
        except subprocess.CalledProcessError as e:
            logger.error(f"An error occurred during RMSD calculation: {str(e)}")
            raise e

        # Parse the RMSD from the output
        rmsd_value = stdout.strip()
        # If an output filename is provided, write the RMSD to the file
        if save:
            if output_filename is not None:
                with open(output_filename + ".txt", "w", encoding="utf-8") as file:
                    file.write(f"RMSD: {rmsd_value}\n")
            else:
                raise ValueError("Output filename is required if save is True")
        return rmsd_value
