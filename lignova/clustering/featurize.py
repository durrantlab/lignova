# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Class to convert SMILES to Morgan fingerprints"""

from dataclasses import dataclass, field

from loguru import logger
from rdkit import Chem
from rdkit.Chem import Mol, rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.DataStructs.cDataStructs import ExplicitBitVect


def _standardize_mol(mol) -> Mol | None:
    """Convert a molecule to its neutral parent form by stripping salts, normalizing, and uncharging.

    Args:
        mol: RDKit molecule object.
    Returns:
        The standardized RDKit molecule object, or None if standardization fails.
    """
    try:
        mol = rdMolStandardize.FragmentParent(mol)
        mol = rdMolStandardize.Normalizer().normalize(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    except Exception:
        return None
    return mol


def resolve_smiles(variants: set[str], *, standardize: bool = True) -> tuple[str, str]:
    """Resolve a set of SMILES mapping to the same InChIKey to a single deterministic representative SMILES string.

    Args:
        variants: A set of SMILES strings that map to the same InChIKey.
        standardize: Whether to standardize the molecules before comparison. Default is True.

    Returns:
        A tuple containing the chosen SMILES string and a method string indicating how the choice was made one of the following:
            - "single"              Only one SMILES string is available for these variants.
            - "same_structure"      variants that were parsed canonicalize to one structure so any of them can be returned. Thus, the minimum lexicographically ordered canonical SMILES of the variants is returned.
            - "standardized"        when standardize is True, when compared in their neutral parent forms, the neutral parents agree. Thus, the canonical SMILES of the neutral parent is returned.
            - "standardized_picked" when standardize is True, but the neutral parents still differ. Thus, the minimum lexicographically ordered canonical SMILES of the neutral parents is returned
            - "picked_raw"          when standardize is False or nothing standardized, the minimum lexicographically ordered canonical SMILES of the raw variants is returned.
    """
    if len(variants) == 1:
        return next(iter(variants)), "single"

    raw_canonical = set()
    for s in variants:
        mol = Chem.MolFromSmiles(s)
        if mol is not None:
            raw_canonical.add(Chem.MolToSmiles(mol))
    if len(raw_canonical) == 1:
        return next(iter(raw_canonical)), "same_structure"

    if standardize:
        output = set()
        for s in variants:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                continue
            mol = _standardize_mol(mol)
            if mol is None:
                continue
            output.add(Chem.MolToSmiles(mol))
        if output:
            method = "standardized" if len(output) == 1 else "standardized_picked"
            return min(output), method

    canonical = []
    for s in variants:
        mol = Chem.MolFromSmiles(s)
        canonical.append(Chem.MolToSmiles(mol) if mol is not None else s)
    return min(canonical), "picked_raw"


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
