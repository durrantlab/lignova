r"""Implementation for RMSD analysis."""
# pylint: disable=R0801
from typing import Iterable, TextIO, Union

import os
import subprocess

import MDAnalysis as mda
from loguru import logger
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

from ..docking.contexts import GlideContext
from ..docking.utils import convert_to_pdb, get_complexes
from ..structure.editing import filter_hetatoms, find_common_atoms, select_common_atoms
from ..structure.ligand import DockedLigand
from ..structure.protein import Protein


class RMSD:
    r"""Implementation for RMSD analysis."""

    def __init__(
        self,
        ligand: Union[Iterable[DockedLigand], DockedLigand],
        reference: Union[Iterable[Protein], Protein],
        context: GlideContext.get_current(),
    ):
        r"""Initialize RMSD class.

        Parameters
        ----------
        ligand: Union[Iterable[DockedLigand], DockedLigand]
            Docked ligand(s) Object that will be analyzed.
        reference: Union[Iterable[Protein], Protein]
            Reference ligand(s) in a Protein object that will be used for comparison.
        context: GlideContext object
            Docking context.
        """
        self.ligand = ligand
        self.reference = reference
        self.context = context

    def rmsd_schrodinger(
        self,
        output_filename: Union[str, TextIO],
        cuttoff: Union[int, float, None] = None,
        asl: str = "ligand",
        align: bool = False,
        neutral_sccafold: bool = False,
    ):
        r"""Calculate RMSD between docked ligand and reference ligand using Schrodinger.
        Parameters
        ----------
        output_filename : str
            Output file name.
        cuttoff : Union[int, float, None]
            Cuttoff for RMSD. Default is None.
        asl : str
            Atom selection language for RMSD calculations. Default is "ligand".
        align : bool
            Align the docked ligand to the reference ligand. Default is True.
        neutral_sccafold : bool
            Neutralize the ligand. Default is False.
        """
        command = [self.context.command + "/run", "rmsd.py", "-a", asl]
        if neutral_sccafold:
            command.append("-use_neutral_scaffold")
        if cuttoff is not None:
            command.append("-r")
            command.append(str(cuttoff))
        if align:
            command.append("-m")
        command.extend(
            ["-c", output_filename, self.reference.file_path, self.ligand.file_path]
        )

        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
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
        except Exception as e:
            logger.error(f"An error occurred during rmsd calculation: {str(e)}")
            raise e

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
            get_complexes(self.ligand.file_path)
            file_name = os.path.splitext(self.ligand.file_name)[0]
            convert_to_pdb(
                os.path.join(self.context.write_dir, file_name + "_complexes.maegz")
            )
            self.ligand = DockedLigand(
                os.path.join(self.context.write_dir, file_name + "_complexes.pdb")
            )
        if self.reference.file_ext != "pdb":
            get_complexes(self.reference.file_path)
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
                f"Number of atoms in the reference ligand ({len(ref_ligand.atoms)}) and docked ligand ({len(dock_ligand.atoms)}) are not the same."
            )
            # Find the common atoms
            common_atoms = find_common_atoms(ref_ligand, dock_ligand)
            # Select the common atoms from each ligand
            ref_ligand = select_common_atoms(ref_ligand, common_atoms)
            dock_ligand = select_common_atoms(dock_ligand, common_atoms)

        # Calculate the RMSD between the ligands
        rmsd_values = []
        if len(docked_traj.trajectory) > 1:
            for traj in docked_traj.trajectory:
                rmsd_value = rmsd(dock_ligand.positions, ref_ligand.positions)
                rmsd_values.append(rmsd_value)
        else:
            rmsd_value = rmsd(dock_ligand.positions, ref_ligand.positions)
            rmsd_values.append(rmsd_value)
        logger.debug(rmsd_values)
        return rmsd_values

    def rmsd_obabel(
        self,
        firstonly: bool = False,
        save: bool = True,
        minimize: bool = False,
        output_filename: Union[str, TextIO, None] = None,
    ):
        """Calculate RMSD between reference and target file using OpenBabel,
        taking into account the symmetry of the molecules.

        Parameters
        ----------
        firstonly : bool
            Only calculate the RMSD for the first structure in the reference file. Default is True.
        save : bool
            Write the RMSD to a txt file. Default is True.
        minimize : bool
            Compute minimum RMSD. Default is False.
        output_filename : Union[str, TextIO,None]
            Output file name if csv is true. Default is None.
        """
        command = ["obrms", "-x"]
        if firstonly:
            command.append("-f")
        if minimize:
            command.append("-m")
        command.extend([self.reference.file_path, self.ligand.file_path])
        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            stdout, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(f"RMSD calculation completed for {self.ligand.file_id}")
                logger.info(f"Output:\n{stdout}")
                logger.debug(f"Error Output:\n{stderr}")
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
        rmsd = stdout.strip()

        # If an output filename is provided, write the RMSD to the file
        if save:
            if output_filename is not None:
                with open(output_filename + ".txt", "w") as file:
                    file.write(f"RMSD: {rmsd}\n")
            else:
                raise ValueError("Output filename is required if csv is True")
        return rmsd
