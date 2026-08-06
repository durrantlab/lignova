# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Pydantic models for PubChem assay and compound-property payloads."""

from collections.abc import Iterable
from enum import StrEnum
from typing import Any, ClassVar

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    computed_field,
    field_validator,
)

_UNSPECIFIED = "unspecified"

DEFAULT_PROPERTIES: list[str] = [
    "SMILES",
    "XLogP",
    "ExactMass",
    "MolecularWeight",
    "Charge",
    "Complexity",
    "TPSA",
    "InChIKey",
]


def _blank_to_none(v: Any) -> Any:
    """Treat blank/whitespace cells as missing.
    Arg:
        v: The value to check.
    Returns:
        None when v is a blank string; otherwise v unchanged."""
    if isinstance(v, str) and v.strip() == "":
        return None
    return v


class AffinityStrategy(StrEnum):
    """An enum to specify how to aggregate multiple activity values for the same CID within an assay and same activity type."""

    ALL = "all"
    MIN = "min"
    AVERAGE = "average"


def _aggregate(
    values: list[tuple[float | None, str, str]], strategy: AffinityStrategy
) -> list[tuple[float | None, str, str]]:
    """Implement the correct aggregation strategy for multiple activity values of the same type for a single CID taking into account the activity outcome.
    Args:
        values: A list of tuples, each containing an activity value (float or None) and its corresponding activity name (str).
        strategy: The aggregation strategy to use. Refer to the `[AffinityStrategy][APIs.pubchem.model.AffinityStrategy]` enum for available options.
    Returns:
        A list of tuples, each containing the aggregated activity value (float or None) and its corresponding activity name (str). The order of the activity names is preserved from the input list.
    """
    groups: dict[tuple[str, str], list[float]] = {}
    order: list[tuple[str, str]] = []
    for value, affinity_type, outcome in values:
        key = (affinity_type, outcome)
        if key not in groups:
            groups[key] = []
            order.append(key)
        if value is not None:
            groups[key].append(value)

    result: list[tuple[float | None, str, str]] = []
    for affinity_type, outcome in order:
        nums = groups[(affinity_type, outcome)]
        if not nums:
            result.append((None, affinity_type, outcome))
            continue
        agg = min(nums) if strategy is AffinityStrategy.MIN else sum(nums) / len(nums)
        result.append((agg, affinity_type, outcome))
    return result


class CompoundProperties(BaseModel):
    """Compound properties for a single CID from the property endpoint."""

    model_config: ClassVar[ConfigDict] = ConfigDict(
        populate_by_name=True, extra="allow"
    )

    smiles: str = Field(default="", alias="SMILES")
    """SMILES string including stereochemistry and isotopes."""

    xlogp: float | None = Field(default=None, alias="XLogP")
    """Computed octanol-water partition coefficient (hydrophobicity)."""

    exact_mass: float | None = Field(default=None, alias="ExactMass")
    """Mass of the most likely isotopic composition (mass-spec oriented)."""

    molecular_weight: float | None = Field(default=None, alias="MolecularWeight")
    """A sum of all atomic weights of the atoms in the compound, measured in g/mol."""

    charge: int | None = Field(default=None, alias="Charge")
    """Net formal charge of the molecule."""

    complexity: float | None = Field(default=None, alias="Complexity")
    """Bertz/Hendrickson/Ihlenfeldt molecular complexity rating. Smaller values indicate simpler molecules."""

    tpsa: float | None = Field(default=None, alias="TPSA")
    """Topological polar surface area."""

    inchikey: str = Field(default="", alias="InChIKey")
    """Hashed InChI; canonical 27-character compound identifier."""

    @field_validator(
        "xlogp",
        "exact_mass",
        "molecular_weight",
        "charge",
        "complexity",
        "tpsa",
        mode="before",
    )
    @classmethod
    def _numeric_blank_to_none(cls, v: Any) -> Any:
        """To ensure  empty strings are treated as None for numeric fields."""
        return _blank_to_none(v)

    @field_validator("smiles", "inchikey", mode="before")
    @classmethod
    def _str_none_to_blank(cls, v: Any) -> Any:
        return "" if v is None else v


class _AssayRecord(BaseModel):
    """A single row from a PubChem concise assay table.

    Each row has a CID, an activity outcome, and an activity value (in uM)"""

    model_config: ClassVar[ConfigDict] = ConfigDict(populate_by_name=True)

    cid: int | None = Field(default=None, alias="CID")
    """PubChem Compound ID for the row, or None if the cell was blank."""

    outcome: str = Field(default="", alias="Activity Outcome")
    """Activity outcome (e.g. Active, Inactive, Inconclusive, Unspecified)."""

    activity: str = Field(default="", alias="Activity Value [uM]")
    """Raw activity value cell in uM, kept as text; parsed on demand."""

    assay_type: str = Field(default="", alias="Assay Type")
    """Type of assay (e.g. Screening, Confirmatory, Functional, ADME )."""

    activity_name: str = Field(default="", alias="Activity Name")
    """Per-row activity type (e.g. IC50, Ki, EC50)."""

    target_accession: str = Field(default="", alias="Target Accession")
    """Target protein accession (e.g. UniProt) for the assay row or None if not reported."""

    target_geneid: int | None = Field(default=None, alias="Target GeneID")
    """NCBI Gene ID of the assay target, or None if not reported."""

    pubmed_id: int | None = Field(default=None, alias="PMID")
    """PubMed ID associated with the CID , or None if not reported."""
    
    activity_qualifier: str = Field(default="", alias="Activity Qualifier")
    """Relation for the activity value (e.g. '<', '=', '>') for Bulk dump only."""

    activity_unit: str = Field(default="", alias="Activity Unit")
    """Unit of the activity value (e.g. 'uM', 'nM') for Bulk dump only. REST API defaults to uM """

    properties: CompoundProperties | None = None
    """Compound properties for this row's CID, filled by enrich_properties."""


    @field_validator("cid", "target_geneid", "pubmed_id", mode="before")    
    @classmethod
    def _id_blank_to_none(cls, v: Any) -> Any:
        return _blank_to_none(v)

    @field_validator(
        "outcome",
        "activity",
        "assay_type",
        "activity_name",
        "target_accession",
        "activity_qualifier",
        "activity_unit",
        mode="before",
    )
    @classmethod
    def _str_none_to_blank(cls, v: Any) -> Any:
        return "" if v is None else v

    @computed_field
    @property
    def activity_value(self) -> float | None:
        """Activity value parsed to float in uM, or None if blank/unparseable."""
        text = self.activity.strip()
        if not text:
            return None
        try:
            return float(text)
        except ValueError:
            return None


class AssayInfo(BaseModel):
    """A parsed PubChem concise assay record."""

    model_config: ClassVar[ConfigDict] = ConfigDict(populate_by_name=True)

    aid: int
    """The PubChem Assay ID."""

    records: list[_AssayRecord] = Field(default_factory=list)
    """One record per row of the concise table."""

    @computed_field
    @property
    def pubmed_id(self) -> int | None:
        """Return the PMID shared by all records, or None if they disagree or none is set. Per-row PMIDs live on record.pubmed_id."""
        seen = {r.pubmed_id for r in self.records if r.pubmed_id is not None}
        return next(iter(seen)) if len(seen) == 1 else None

    @classmethod
    def from_concise(cls, aid: int, data: dict[str, Any]) -> "AssayInfo":
        """Parse a concise assay payload into an AssayInfo.

        Args:
            aid: The PubChem Assay ID the payload belongs to.
            data: The raw concise JSON response (already confirmed to be a dict
                without a Fault by the client).

        Returns:
            An AssayInfo. Records/pubmed_id are empty/None if the table is absent.
        """
        table = data.get("Table") or {}
        columns = table.get("Columns", {}).get("Column", [])
        rows = table.get("Row", [])
        if not columns or not rows:
            return cls(aid=aid)

        pm_index = columns.index("PubMed ID") if "PubMed ID" in columns else None

        records: list[_AssayRecord] = []
        for row in rows:
            cells = row.get("Cell", [])
            row_map = {col: cells[i] for i, col in enumerate(columns) if i < len(cells)}
            if "PubMed ID" in row_map:
                row_map["PMID"] = row_map["PubMed ID"]
            records.append(_AssayRecord.model_validate(row_map))

        return cls(aid=aid, records=records)

        return cls(aid=aid, records=records, pubmed_id=pubmed_id)

    def _cids_where(self, outcome: str) -> list[int]:
        """Deduplicated CIDs whose outcome matches in a case-insensitive manner."""
        seen: set[int] = set()
        result: list[int] = []
        target = outcome.strip().lower()
        for r in self.records:
            if r.cid is not None and r.outcome.strip().lower() == target:
                if r.cid not in seen:
                    seen.add(r.cid)
                    result.append(r.cid)
        return result

    @property
    def active_cids(self) -> list[int]:
        """CIDs marked Active in the assay."""
        return self._cids_where("active")

    @property
    def inactive_cids(self) -> list[int]:
        """CIDs marked Inactive in the assay."""
        return self._cids_where("inactive")

    def cids(self, outcome: str) -> list[int]:
        """CIDs for an arbitrary outcome (e.g. 'inconclusive', 'probe')."""
        return self._cids_where(outcome)

    @property
    def unique_cids(self) -> list[int]:
        """All distinct, non-blank CIDs in the assay (any outcome)."""
        seen: set[int] = set()
        result: list[int] = []
        for r in self.records:
            if r.cid is not None and r.cid not in seen:
                seen.add(r.cid)
                result.append(r.cid)
        return result

    def binding_affinity(
        self,
        cids: Iterable[int | str] | None = None,
        strategy: AffinityStrategy | str = AffinityStrategy.ALL,
        outcome: str | None = None,
    ) -> dict[int, list[tuple[float | None, str, str]]]:
        """Map CID to a list of (activity_value, activity_name, outcome) entries.

        Args:
            cids: Optional iterable of CIDs to keep. If None, all are kept.
            strategy: "all" keeps every value; "min"/"average" collapse to one
                entry per activity_name using non-None values.
            outcome: Optional filter for activity outcome (e.g. "active"). If
                None, all outcomes excluding unspecified are kept.

        Returns:
            dict mapping CID to a list of affinity value_or_None, activity_name and outcome tuples
        """
        strategy = AffinityStrategy(strategy)
        wanted = None if cids is None else {int(c) for c in cids}
        target = outcome.strip().lower() if outcome else None
        collected: dict[int, list[tuple[float | None, str, str]]] = {}
        for r in self.records:
            oc = r.outcome.strip().lower()
            if target is None:
                if oc == _UNSPECIFIED:
                    continue
            elif oc != target:
                continue
            if r.cid is None:
                continue
            if wanted is not None and r.cid not in wanted:
                continue
            if not r.activity.strip():
                continue
            collected.setdefault(r.cid, []).append(
                (r.activity_value, r.activity_name, r.outcome)
            )
        if strategy is AffinityStrategy.ALL:
            return collected
        return {cid: _aggregate(vals, strategy) for cid, vals in collected.items()}

    def to_rows(self, outcome: str | None = None) -> list[dict[str, Any]]:
        """Flatten to one dict per record. Optionally filter by outcome."""
        target = outcome.strip().lower() if outcome else None
        rows: list[dict[str, Any]] = []
        for r in self.records:
            if target is not None and r.outcome.strip().lower() != target:
                continue
            row: dict[str, Any] = {
                "aid": self.aid,
                "cid": r.cid,
                "outcome": r.outcome,
                "activity_value": r.activity_value,
                "activity_name": r.activity_name,
                "assay_type": r.assay_type,
                "target_accession": r.target_accession,
                "target_geneid": r.target_geneid,
                "pubmed_id": r.pubmed_id,
                "activity_qualifier": r.activity_qualifier,
                "activity_unit": r.activity_unit,
            }
            if r.properties is not None:
                row.update(r.properties.model_dump(by_alias=True))
            rows.append(row)
        return rows

    def to_dataframe(self, outcome: str | None = None) -> "Any":
        """Return the flattened rows as a pandas DataFrame, optionally filtered."""
        import pandas as pd

        return pd.DataFrame(self.to_rows(outcome))
