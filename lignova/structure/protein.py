"""Implementation of protein class."""
from typing import Union

from collections.abc import Iterable

import requests
from loguru import logger

from ..io import write_text
from .base import Prepared, Structure


class Protein(Structure):
    """Proteins."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def get_pdb_from_rcsb(pdb_id: str) -> str:
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        try:
            response = requests.get(pdb_url, timeout=30)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to download PDB {pdb_id}: {e}")
            raise e
        return response.text

    def _load_from_pdb_id(
        self, pdb_id: str, write: bool = False, write_path: Union[None, str] = None
    ) -> None:
        pdb_text = Protein.get_pdb_from_rcsb(pdb_id)
        if write:
            if write_path is None:
                raise ValueError("Must provide write_path if write is True.")
            else:
                self._pdb_file_path = write_text(pdb_text, write_path, file_ext="pdb")
        else:
            self._pdb_text = pdb_text

    @property
    def pdb(self) -> Union[str, None]:
        if hasattr(self, "_pdb_text"):
            return self._pdb_text
        if hasattr(self, "_pdb_file_path"):
            with open(self._pdb_file_path, encoding="utf-8") as f:
                return f.read()
        return None

    def load(
        self,
        file_path: Union[str, None] = None,
        write: bool = False,
        write_path: Union[None, str] = None,
        pdb_id: Union[str, None] = None,
    ) -> None:
        r"""Load structural information for a protein.

        Parameters
        ----------
        file_path
            Path to file to load.
        write
            Keep structure in file and load when requested. If ``False``, this will
            keep the structure in memory.
        write_path
            Path to write to file. If ``None``, then a ``NamedTemporaryFile`` will
            be created instead.
        pdb_id
            Four-letter code to load structure from RCSB.
        """
        if file_path is not None:
            self._pdb_file_path = file_path
        if pdb_id is not None: 
            self._load_from_pdb_id(pdb_id, write, write_path)


class PreparedProtein(Protein, Prepared):
    r"""A protein that has been prepared for some downstream application."""
