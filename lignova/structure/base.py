r"""Implementation of structure class."""

import os
from abc import ABC, abstractmethod


class Structure(ABC):
    r"""Base class for any physical structure."""

    def __init__(self, file_path: str | None = None, file_id: str | None = None):
        self.file_path = file_path
        if file_id is None and file_path is not None:
            self.file_name = os.path.basename(file_path)
            self.file_id, self.file_ext = os.path.splitext(self.file_name)
            self.file_ext = self.file_ext.lstrip(".")
        elif file_id is not None and self.file_path is not None:
            self.file_id = file_id
            self.file_ext = os.path.splitext(self.file_path)[1].lstrip(".")

    @abstractmethod
    def load(
        self,
        file_path: str | None = None,
        write: bool = False,
        write_path: None | str = None,
        pdb_id: str | None = None,
    ) -> None:
        r"""Load structural data from a variety of sources.

        Args:
            file_path : Path to file to load.
            write : Keep structure in file and load when requested. If ``False``, this will
                keep the structure in memory.
            write_path : Path to write to file.
            pdb_id : Four-letter code to load structure from RCSB.
        """
        raise NotImplementedError()


class Prepared:
    r"""Base class for prepared structures."""

    def get_info(self, context):
        r"""This function gives information about the prepared structure.

        Args:
            context : Context object with information about the preparation.

        Returns:
            Dictionary with information about the preparation depending
                on the type of structure.
        """
        prot_info = {
            "forcefield": context.forcefield,
            "fillsidechains": context.fillsidechains,
            "disulfides": context.disulfides,
            "epik_pH": context.epik_ph,
            "epik_pHt": context.epik_pht,
            "propka_pH": context.propka_ph,
            "rmsd": context.prot_rmsd,
            "water_distance": context.prot_watdist,
            "grid_innerbox": context.grid_innerbox,
            "sample_water": context.samplewater,
            "minimize_adj_h": context.minimize_adj_h,
            "hydrogen_addition": context.rehtreat,
        }
        lig_info = {
            "forcefield": context.lig_forcefield,
            "pH": context.lig_ph,
            "max_atoms": context.lig_max_mw,
            "epik": context.lig_epik,
            "pHt": context.lig_pht,
            "stereoisomers": context.lig_stereoisomers,
        }
        if "ligand" in str(type(self)):
            return lig_info
        if "protein" in str(type(self)):
            return prot_info
        if "ligand" not in str(type(self)) and "protein" not in str(type(self)):
            return None
