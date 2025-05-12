"""Implementation of ligand class."""

from .base import Prepared, Structure


class Ligand(Structure):
    """Class for ligands to be docked to proteins."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._ligand_text = None

    # Define a property to access _ligand_text
    @property
    def ligand_text(self):
        r"""Return the ligand text."""
        if self._ligand_text is None:
            self.load(self.file_path)
        return self._ligand_text

    def load(
        self,
        file_path: str | None = None,
        write: bool = False,
        write_path: str | None = None,
        pdb_id: str | None = None,
    ) -> None:
        r"""Load structural information for a ligand.

        Args:
            file_path : Path to ligand file.
            write : if true write the ligand to disk. Default is False.
            write_path : Path to write ligand to disk.
            pdb_id : PDB ID of ligand to download.
        Returns:
            None
        """
        if file_path is not None:
            # read in file
            with open(file_path, encoding="utf-8") as file:
                self._ligand_text = file.read()
        if write_path:
            # write to file
            with open(write_path, "w", encoding="utf-8") as file:
                file.write(self._ligand_text)


class PreparedLigand(Ligand, Prepared):
    r"""Class for prepared ligands to be docked to proteins."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class DockedLigand(Ligand):
    r"""Class for docked ligands to be docked to proteins."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
