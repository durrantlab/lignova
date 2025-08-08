r"""Implements the Docking class."""

from abc import ABC, abstractmethod
from collections.abc import Iterable

from ..structure.ligand import DockedLigand, PreparedLigand
from ..structure.protein import PreparedProtein
from .contexts import GlideContext


class Docking(ABC):
    r"""Abstract class for docking ."""

    @abstractmethod
    def run(
        self,
        target: PreparedProtein,
        ligand: Iterable[PreparedLigand] | PreparedLigand,
        context: GlideContext,
    ) -> Iterable[DockedLigand] | DockedLigand:
        r"""Dock one or multiple ligands to a single target.


        Args:
            target
                Prepared protein that ligands will be docked to.
            ligand
                Prepared ligand(s) that will be docked into `target`
            context
                Docking configuration context for specific program.

        """
        raise NotImplementedError()
