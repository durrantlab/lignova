"""Class to convert SMILES to Morgan fingerprints"""

from dataclasses import dataclass, field

from loguru import logger
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs.cDataStructs import ExplicitBitVect


@dataclass(frozen=True, slots=True)
class FeaturizeResult:
    """Dataclass holding the results of a featurization pass."""

    items: dict[str, ExplicitBitVect]
    """Dictionary of compound id as key and their Morgan fingerprint as value """
    skipped: dict[str, str] = field(default_factory=dict)
    """Dictionary of compound id as key and the reason for skipping as value."""

    @property
    def n_compounds(self) -> int:
        """Total compounds seen (featurized + skipped)."""
        return len(self.items) + len(self.skipped)


class MorganFeaturizer:
    """Class to convert SMILES to Morgan fingerprints.

    Args:
        radius: Morgan radius to use . Default is 2 which is ECFP4.
        fp_size: Fingerprint length in bits.
    """

    def __init__(self, radius: int = 2, fp_size: int = 2048) -> None:
        self.radius = radius
        self.fp_size = fp_size
        self._generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=fp_size
        )

    def featurize(self, raw: dict[str, str]) -> FeaturizeResult:
        """Featurize a batch of SMILES keyed by their compound id.

        Args:
            raw: dictionary of compound id as key and their SMILES as value.

        Returns:
            A `FeaturizeResult` object containing the fingerprints and any skipped compounds.
        """
        items: dict[str, ExplicitBitVect] = {}
        skipped: dict[str, str] = {}

        for cid, smiles in raw.items():
            try:
                mol = Chem.MolFromSmiles(smiles)
            except Exception as exc:
                skipped[cid] = f"SMILES parse error: {exc}"
                continue
            if mol is None:
                skipped[cid] = "invalid SMILES"
                continue
            items[cid] = self._generator.GetFingerprint(mol)

        if skipped:
            logger.warning(
                "Featurized {n_items}/{n_raw} compounds; skipped {n_skipped} (see FeaturizeResult.skipped)",
                n_items=len(items),
                n_raw=len(raw),
                n_skipped=len(skipped),
            )
        else:
            logger.info(f"Featurized {len(items)} compounds")

        return FeaturizeResult(items=items, skipped=skipped)
