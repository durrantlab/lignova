# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Aggregate raw compound and target data into one CompoundActivity object per compound."""

import math
import statistics
from collections.abc import Iterable
from dataclasses import dataclass

from loguru import logger

from lignova.activity.models import CompoundActivity
from lignova.analysis.utils import to_pActivity

CONVERTIBLE_TYPES = frozenset(
    {"Kd", "Ki", "Kb", "Kieq", "Kic", "IC50", "EC50", "fIC50", "fEC50"}
)
"""The assay types that convert to a comparable pActivity. Anything outside this set is ignored."""

_TYPE_RANK: dict[str, int] = {
    "Kd": 0,
    "Ki": 1,
    "Kieq": 2,
    "Kic": 2,
    "Kb": 3,
    "IC50": 4,
    "fIC50": 4,
    "EC50": 5,
    "fEC50": 5,
}
"""Preference order for winning-type selection: Kd > Ki > {Kieq, Kic} > Kb > {IC50, fIC50} > {EC50, fEC50}."""

DEFAULT_SPREAD_GATE = 1.0
"""Maximum allowed spread in log scale among the winning-type values before the quality gate fails."""


@dataclass(frozen=True, slots=True)
class Measurement:
    """A single raw, exact measurement row for a compound-target pair."""

    type: str
    """The assay type e.g kd, Ki, IC50, EC50, etc. It is used to select the winning type for aggregation."""

    value: float
    """The raw activity value, expressed in `unit`."""

    unit: str
    """The unit of `value` which is used in the pActivity conversion to a log scale."""

    cid: str = ""
    """The comppound id for the activity measurement"""

    assay_id: str = ""
    """The assay id (`aid`) attached to the measurement."""

    is_active: bool = True
    """True if the compound is active from the source assay outcome."""


def group_measurements(
    rows: Iterable[dict],
    key_col: str = "InChIKey",
    cid_col: str = "cid",
    assay_col: str = "aid",
    outcome_col: str = "outcome",
    type_col: str = "activity_name",
    value_col: str = "activity_value",
    unit_col: str = "activity_unit",
    qualifier_col: str = "activity_qualifier",
) -> dict[str, list[Measurement]]:
    """Group parquet rows for a single target into a list of exact measurements per compound.

    Args:
        rows: An iterable of row dictionaries for a single target_geneid.
        key_col: The column holding the compound key. Default is "InChIKey" which gives one structural
            node per standardized parent; pass "cid" only to deliberately key per-CID.
        cid_col: The column holding the source CID, retained on each Measurement for the crosswalk.
        assay_col: The column holding the source assay id (`aid`).
        outcome_col: The column holding the active/inactive outcome label.
        type_col: The column holding the assay type (e.g. `activity_name`).
        value_col: The column holding the raw activity value.
        unit_col: The column holding the value unit.
        qualifier_col: The column holding the relation qualifier ('=', '>', '<'); only literal '='
            rows are kept. Blank/missing qualifiers are treated as relative and dropped.

    Returns:
        A dictionary with the compound keys as the keys and the list of their exact `Measurement`
        records as the values.
    """
    grouped: dict[str, list[Measurement]] = {}
    for row in rows:
        key = row.get(key_col)
        value = row.get(value_col)
        if not key or value in (None, ""):
            continue
        qualifier = str(row.get(qualifier_col) or "").strip()
        if qualifier != "=":
            continue
        try:
            value = float(value)
        except (TypeError, ValueError):
            continue
        outcome = str(row.get(outcome_col) or "").strip().lower()
        grouped.setdefault(str(key), []).append(
            Measurement(
                type=str(row.get(type_col) or ""),
                value=value,
                unit=str(row.get(unit_col) or ""),
                cid=str(row.get(cid_col) or ""),
                assay_id=str(row.get(assay_col) or ""),
                is_active=outcome != "inactive",
            )
        )
    return grouped


def cid_crosswalk(measurements: dict[str, list[Measurement]]) -> dict[str, set[str]]:
    """Build the compound key to source CIDs map, for provenance and per-CID output expansion.

    Args:
        measurements: The output of `group_measurements`, keyed on the compound InChIKey.

    Returns:
        A dictionary with the compound InChIKey as the keys and the set of CIDs that resolved to each
        as the values. Use it to expand InChIKey-keyed cluster or cliff outputs back to per-CID rows.
    """
    original: dict[str, set[str]] = {}
    for key, rows in measurements.items():
        cids = {m.cid for m in rows if m.cid}
        if cids:
            original[key] = cids
    return original


def _type_rank(t: str) -> int:
    """Return the winning-type preference rank of an assay type.

    Args:
        t: The assay type string.

    Returns:
        The integer rank; lower is preferred. Types absent from `_TYPE_RANK` get the largest rank.
    """
    return _TYPE_RANK.get(t, len(_TYPE_RANK))


def _pactivity(value: float, unit: str) -> float | None:
    """Convert one value to pActivity via `to_pActivity` to guard against exceptions and non-finite results.

    Args:
        value: The raw activity value.
        unit: The unit of `value`.

    Returns:
        The pActivity as a float, or None if the unit is unknown or the result is not finite.
    """
    try:
        pact = to_pActivity(value, unit)
    except (ValueError, TypeError):
        return None
    if not isinstance(pact, float):
        return None
    return pact if math.isfinite(pact) else None


def aggregate_activity(
    measurements: dict[str, list[Measurement]],
    spread_gate: float = DEFAULT_SPREAD_GATE,
    convertible: frozenset[str] = CONVERTIBLE_TYPES,
) -> dict[str, CompoundActivity]:
    """Aggregate the per-compound measurements into a single `CompoundActivity` for each compound.

    Args:
        measurements: The output of `group_measurements`, keyed on the compound key.
        spread_gate: The quality gate fails when the winning-type spread exceeds this (in log units).
        convertible: The assay types that convert to pActivity.

    Returns:
        A dictionary with the compound keys as the keys and their `CompoundActivity` records as the
        values. Compounds with no convertible, exact measurement are omitted (they still cluster;
        cliff detection skips them through the missing-activity path).
    """
    activity: dict[str, CompoundActivity] = {}
    n_skipped = 0

    for key, rows in measurements.items():
        # For each convertible, positive value: its pActivity plus the row's assay id and label.
        by_type: dict[str, list[tuple[float, str, bool]]] = {}
        for m in rows:
            if m.type not in convertible or m.value <= 0:
                continue
            pact = _pactivity(m.value, m.unit)
            if pact is None:
                continue
            by_type.setdefault(m.type, []).append((pact, m.assay_id, m.is_active))

        if not by_type:
            n_skipped += 1
            continue

        winning_type = min(by_type, key=lambda t: (_type_rank(t), -len(by_type[t])))
        winners = by_type[winning_type]
        vals = [p for p, _, _ in winners]
        spread = max(vals) - min(vals)
        activity[key] = CompoundActivity(
            pActivity=statistics.median(vals),
            winning_type=winning_type,
            passes_quality_gate=spread <= spread_gate,
            is_active=any(active for _, _, active in winners),
            assay_ids=frozenset(aid for _, aid, _ in winners if aid),
        )

    logger.info(
        "Aggregated {n} compounds to activity with ({s} skipped because they had no convertible measurement).",
        n=len(activity),
        s=n_skipped,
    )
    return activity
