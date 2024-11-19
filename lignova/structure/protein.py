"""Implementation of protein class."""

from typing import Union

import requests
from loguru import logger

from ..io import write_text
from .base import Prepared, Structure


class Protein(Structure):
    """Protein class that contains functions for loading proteins
    and preparing them for docking class."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def get_pdb_from_rcsb(pdb_id: str) -> str:
        r"""Download a PDB file from the RCSB PDB database.
        Parameters
        ----------
        pdb_id
            PDB ID of the protein to download.
        Returns
        -------
        str of PDB file
            The PDB file as a string.
        """
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        try:
            response = requests.get(pdb_url, timeout=30)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            if response.status_code == 404:
                logger.warning(
                    f"PDB file not found for {pdb_id} Trying the PDBx/mmCIF."
                )
                pdbx_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
                try:
                    response = requests.get(pdbx_url, timeout=30)
                    response.raise_for_status()
                except requests.exceptions.RequestException as exp:
                    logger.error(f"PDB file not found for {pdb_id}.")
                    raise exp
            else:
                logger.error(f"PDB file not found for {pdb_id}.")
                raise e
        return response.text

    def _load_from_pdb_id(
        self, pdb_id: str, write: bool = False, write_path: Union[None, str] = None
    ) -> None:
        r"""Load structural information for a protein from RCSB.
        Parameters
        ----------
        pdb_id
            PDB ID to load structure from RCSB.
        write
            Keep structure in file and load when requested. If ``False``, this will
            keep the structure in memory.
        write_path
            Path to write to file
        """
        pdb_text = Protein.get_pdb_from_rcsb(pdb_id)
        if write:
            if write_path is None:
                raise ValueError("Must provide write_path if write is True.")
            file_ext = "pdb" if pdb_text.startswith("HEADER") else "cif"
            self._pdb_file_path = write_text(pdb_text, write_path, file_ext=file_ext)
        else:
            self._pdb_text = pdb_text

    @property
    def pdb(self) -> Union[str, None]:
        r"""Return the PDB text."""
        if hasattr(self, "_pdb_text"):
            return self._pdb_text
        if hasattr(self, "_pdb_file_path"):
            with open(self._pdb_file_path, encoding="utf-8") as file:
                return file.read()
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
