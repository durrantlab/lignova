"""Implementation of ligand class."""
from typing import Union

from .base import Prepared, Structure


class Ligand(Structure):
    """Class for ligands to be docked to proteins."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    # TODO _MA: figure out how best to use load in ligand context
    def load(
        self,
        file_path: Union[str, None] = None,
        write: bool = False,
        write_path: Union[None, str] = None,
        pdb_id: Union[str, None] = None,
    ) -> None:
        pass


class PreparedLigand(Ligand, Prepared):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class DockedLigand(Ligand):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
