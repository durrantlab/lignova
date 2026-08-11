# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Implements the Docking class."""

from abc import ABC, abstractmethod
from collections.abc import Iterable

from ..structure.ligand import DockedLigand, PreparedLigand
from ..structure.protein import PreparedProtein

# from .contexts import GlideContext #NOTE:removed the type hint for the context


class Docking(ABC):
    r"""Abstract class for docking ."""

    @abstractmethod
    def run(
        self,
        target: PreparedProtein,
        ligand: Iterable[PreparedLigand] | PreparedLigand,
        context,
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
