"""Base class for RMSD calculators."""

from abc import ABC, abstractmethod

from ...docking.contexts import GlideContext
from ...structure.ligand import DockedLigand
from ...structure.protein import Protein


class RMSDBase(ABC):
    """
    Base class for RMSD calculators.
    """

    def __init__(
        self,
        target: DockedLigand,
        reference: Protein,
        context: GlideContext,
    ):
        r"""Initialize RMSD class.

        Args:
            target:
                Docked ligand(s) Object that will be analyzed.
            reference:
                Reference ligand(s) in a Protein object that will be used for comparison.
            context:
                Docking context. Default is GlideContext.get_current().
        """
        assert isinstance(target, DockedLigand), "Ligand must be a DockedLigand object."
        assert isinstance(reference, Protein), "Reference must be a Protein object."
        assert isinstance(
            context, GlideContext
        ), "Context must be an instance of GlideContext."
        self.target: DockedLigand = target
        self.reference: Protein = reference
        self.context: GlideContext = context

    @abstractmethod
    def calculate(self) -> float | list[float]:
        """Calculate RMSD."""
        raise NotImplementedError()
