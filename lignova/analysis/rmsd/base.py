"""Base class for RMSD calculators."""

from abc import ABC, abstractmethod

from ...structure.ligand import DockedLigand, Ligand


class RMSDBase(ABC):
    """
    Base class for RMSD calculators.
    """

    def __init__(
        self,
        target: DockedLigand,
        reference: Ligand,
    ):
        r"""Initialize RMSD class.

        Args:
            target:
                Docked ligand(s) Object that will be analyzed.
            reference:
                Reference ligand(s) in a Ligand object that will be used for comparison.
        """
        assert isinstance(target, DockedLigand), "Ligand must be a DockedLigand object."
        assert isinstance(reference, Ligand), "Reference must be a Ligand object."
        self.target: DockedLigand = target
        self.reference: Ligand = reference

    @abstractmethod
    def calculate(self) -> float | list[float]:
        """Calculate RMSD."""
        raise NotImplementedError()
