"""Implementation of ligand class."""
from typing import Union

from .base import Prepared, Structure


class Ligand(Structure):
    """Class for ligands to be docked to proteins."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    # Define a property to access _ligand_text
    @property
    def ligand_text(self):
        r"""Return the ligand text."""
        self.load(self.file_path)
        return self._ligand_text

    def load(
        self,
        file_path: Union[str, None] = None,
        write: bool = False,
        write_path: Union[None, str] = None,
        pdb_id: Union[str, None] = None,
    ) -> None:
        r"""Load structural information for a ligand.

        Parameters
          ----------
             file_path : str, optional
                 Path to ligand file.
             write_path : str, optional
                 Path to write the ligand file to disk.
             pdb_id : str, optional
                 PDB ID of ligand to download.
        """
        if file_path is not None:
            # read in file
            with open(file_path, encoding="utf-8") as f:
                self._ligand_text = f.read()
        if write_path:
            # write to file
            with open(write_path, "w") as f:
                f.write(self._ligand_text)


class PreparedLigand(Ligand, Prepared):
    r"""Class for prepared ligands to be docked to proteins."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class DockedLigand(Ligand):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
