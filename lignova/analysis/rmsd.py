r"""Implementation for RMSD analysis."""
from typing import Iterable, TextIO, Union

import subprocess

import MDAnalysis as mda
from loguru import logger
from MDAnalysis.analysis import rms

from ..docking.contexts import GlideContext
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
        self.context = context.command

    def rmsd_schrodinger(
        self,
        output_filename: Union[str, TextIO],
        cuttoff: Union[int, float, None] = None,
        asl: str = "ligand",
        align: bool = True,
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
        command = [self.context + "/run", "rmsd.py", "-a", asl, "-verbose"]
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

    def rmsd_mda(self):
        r"""Calculate RMSD between docked ligand and reference ligand using MDAnalysis."""
        # Load docked ligand and reference PDB structures
        docked = mda.Universe(self.ligand.file_path)
        reference = mda.Universe(self.reference.file_path)

        # Select the atoms for RMSD calculation
        selection = "protein and name CA"  # Modify this according to your needs
