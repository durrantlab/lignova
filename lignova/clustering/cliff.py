# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Class for detecting activity cliffs in a set of compounds based on their similarities and activities."""

import math
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass, replace
from enum import StrEnum

from loguru import logger

from lignova.activity.models import CompoundActivity
from lignova.clustering.tanimoto import TanimotoSimilarities


class SeverityMetric(StrEnum):
    """Enum to declare which landscape index to use for severity ranking."""

    TS_SALI = "ts_sali"
    """Taylor-series Structure-activity landscape index (SALI). It is introduced by López-Pérez, K., & Miranda-Quintana, R. A. (2025).
    It is a taylor series expansion version of the classic SALI that is finite everywhere. It is the default metric."""

    SALI = "sali"
    """This is the classic Structure-activity landscape index (SALI) introduced by Guha and Van Drie in 2008. It is undefined at Tanimoto similarity (Tm) == 1, 
    so the user must choose a policy for how to handle those undefined pairs. by declaring the `sali_undefined` parameter in `CliffParams`."""


class SaliUndefined(StrEnum):
    """Enum to declare how the SALI severity metric resolves pairs that are undefined at Tanimoto similarity (Tm) == 1."""

    NEXT_LARGEST = "next_largest"
    """Replace the undefined value with the largest finite SALI in the ranked set as described by Guha & Van Drie (2008). This is the default policy."""

    MAX_SEVERITY = "max_severity"
    """Rank undefined pairs above all finite ones by assigning them to infinity. As by definition,a Tm = 1
    pair with a large activity gap is an extreme cliff where the same structure (based on the fingerprint), 
    yet different activity. Thus, it should outrank even the largest finite SALI, not merely tie it."""

    EXCLUDE = "exclude"
    """Drop undefined pairs from the ranking where their severity_label remains None. This was also used in the original Guha & Van Drie (2008) paper, 
    but it is not recommended because it discards potentially important extreme cliffs."""


class CliffSeverity(StrEnum):
    """An enum for the Severity label assigned by rank within the ranked cliff set."""

    EXTREME = "extreme"
    STRONG = "strong"
    MODERATE = "moderate"


@dataclass(frozen=True, slots=True)
class CliffParams:
    """Configuration for a single  cliff detection and severity pass."""

    min_delta: float = 2.0
    """Minimum absolute difference in pActivity for a pair to be considered a cliff. Default is 2.0 which means a 100-fold difference."""

    min_similarity: float = 0.55
    """Minimum Tanimoto similarity for a pair to be considered a cliff. Default is 0.55."""

    metric: SeverityMetric = SeverityMetric.TS_SALI
    """Which landscape index to compute and rank by. Default is TS_SALI."""

    ts_truncation: int = 3
    """TS_SALI only parameter representing the Taylor order k. Default is 3, which is the original order used by López-Pérez, K., & Miranda-Quintana, R. A. (2025)."""

    extreme_frac: float = 0.01
    """Top fraction of ranked cliffs to be labelled EXTREME."""

    strong_frac: float = 0.05
    """Cumulative top fraction of ranked cliffs to be labelled at least STRONG."""

    sali_undefined: SaliUndefined = SaliUndefined.NEXT_LARGEST
    """SALI only parameter identifying the policy for Tm == 1 pairs i.e, the undefined values."""

    def __post_init__(self) -> None:
        if self.min_delta < 0:
            raise ValueError(f"min_delta must be >= 0, got {self.min_delta}")
        if self.ts_truncation < 1:
            raise ValueError(f"ts_truncation must be >= 1, got {self.ts_truncation}")
        if not (0.0 < self.extreme_frac <= self.strong_frac <= 1.0):
            raise ValueError(
                f"require 0 < extreme_frac <= strong_frac <= 1, got "
                f"extreme_frac={self.extreme_frac}, strong_frac={self.strong_frac}"
            )
        if (
            self.metric is SeverityMetric.TS_SALI
            and self.sali_undefined is not SaliUndefined.NEXT_LARGEST
        ):
            logger.warning(
                "sali_undefined={p} is ignored on the TS_SALI branch.",
                p=self.sali_undefined.value,
            )
        if self.metric is SeverityMetric.SALI and self.ts_truncation != 3:
            logger.warning(
                "ts_truncation={k} is ignored on the SALI branch.", k=self.ts_truncation
            )


@dataclass(frozen=True, slots=True)
class ActivityCliffs:
    """Dataclass holding the results of one activity cliff detection pass over a set of compounds."""

    id_a: str
    """The compound id of the first compound in the cliff pair."""

    id_b: str
    """The compound id of the second compound in the cliff pair."""

    similarity: float
    """The Tanimoto similarity between the two compounds in the cliff pair."""

    pAct_diff: float
    """The absolute difference in pActivity between the two compounds in the cliff pair."""

    pAct_a: float
    """pActivity of compound a."""

    pAct_b: float
    """pActivity of compound b."""

    winning_type_a: str
    """Winning assay type behind pAct_a."""

    winning_type_b: str
    """Winning assay type behind pAct_b."""

    is_cross_type: bool
    """True when the two winning types differ between the two compounds in the cliff pair."""

    is_low_quality: bool
    """True if either compound failed the activity quality gate. These
    cliffs may rest on measurement noise and can be excluded downstream."""

    is_same_assay: bool
    """True when the two compounds share at least one assay id (a within-assay cliff, firmest
    ground). Cross-assay cliffs can carry a systematic offset."""

    involves_inactive: bool
    """True when at least one compound is labelled inactive by its source assay. This indicates
    an active-vs-inactive cliff."""

    landscape_index: float | None
    """The computed value of the landscape index SALI or TS_SALI (per params.metric) for this pair. 
    It is None only on the SALI branch at Tm == 1 (resolved during labelling by sali_undefined)."""

    severity_label: CliffSeverity | None = None
    """The severity label (EXTREME / STRONG / MODERATE) for this cliff pair. It is None until `label_cliff_severity` runs, or if the pair
    was not part of the ranked population."""


@dataclass(frozen=True, slots=True)
class CliffResult:
    """The standardized output of a cliff pass. It contains the pairs and the params that produced them."""

    cliffs: list[ActivityCliffs]
    """All detected cliff pairs, in detection order (all cliffs generated by `find_activity_cliffs`)."""

    params: CliffParams
    """The configuration these cliffs has used to be detected and will be used to label severity. It is stored here for reproducibility."""

    n_undefined: int = 0
    """Count of SALI-branch pairs undefined at Tm == 1. This would be 0 when using TS_SALI, or when the SALI branch had no Tm == 1 pairs."""

    @property
    def metric(self) -> SeverityMetric:
        """The landscape index that will be used to label the severity of the cliffs."""
        return self.params.metric

    @property
    def n_cliffs(self) -> int:
        """Number of cliff pairs."""
        return len(self.cliffs)

    def severity_counts(self) -> dict[CliffSeverity, int]:
        """Count of cliffs in each severity band where unlabelled pairs are ignored."""
        out: dict[CliffSeverity, int] = {}
        for c in self.cliffs:
            if c.severity_label is not None:
                out[c.severity_label] = out.get(c.severity_label, 0) + 1
        return out


def _landscape_index(
    params: CliffParams, pAct_diff: float, similarity: float
) -> float | None:
    """Compute the landscape index for a given pair of compounds using the specified metric.
    Args:
        params: The CliffParams object containing the detection and scoring configuration (min_delta, metric branch, thresholds).
        pAct_diff: The absolute difference in pActivity between the two compounds in the cliff pair.
        similarity: The Tanimoto similarity between the two compounds in the cliff pair.
    Returns:
        The computed landscape index value (float) or None if the value is undefined (e.g, for SALI at Tm == 1).
    """
    if params.metric is SeverityMetric.SALI:
        denom = 1.0 - similarity
        return None if denom <= 0.0 else pAct_diff / denom
    geom = sum(similarity**m for m in range(params.ts_truncation + 1))
    return (pAct_diff**2) * geom / (params.ts_truncation + 1)


def find_activity_cliffs(
    sims: TanimotoSimilarities,
    activity: dict[str, CompoundActivity],
    params: CliffParams = CliffParams(),
) -> CliffResult:
    """Detect activity cliffs in a set of compounds based on their similarities and activities.

    Args:
        sims: Output of `compute_pairwise` over the compounds to analyze.
        activity: A dictionary with the compound ids as the keys and their `CompoundActivity` records as the values.
        params: the `CliffParams` object containing the detection and scoring configuration (min_delta, metric branch, thresholds).
    Returns:
        A `CliffResult` holding the detected cliffs; severity labels are None until
        `label_cliff_severity` is called.
    """
    if sims.min_sim > params.min_similarity:
        raise ValueError(
            f"Edges listed in sims were built at floor {sims.min_sim}, but "
            f"params.min_similarity={params.min_similarity} is higher. "
            f"Thus activity cliff pairs are unreachable. Please rebuild sims at a "
            f"lower floor or lower params.min_similarity."
        )
    cliffs: list[ActivityCliffs] = []
    for id_a, id_b, similarity in sims.edges:
        if similarity < params.min_similarity:
            continue
        a = activity.get(id_a)
        b = activity.get(id_b)
        if a is None or b is None:
            logger.warning(
                "Missing pActivity for {id_a} or {id_b}, skipping cliff detection for this pair.",
                id_a=id_a,
                id_b=id_b,
            )
            continue

        pAct_diff = abs(a.pActivity - b.pActivity)
        if pAct_diff >= params.min_delta:
            cliffs.append(
                ActivityCliffs(
                    id_a=id_a,
                    id_b=id_b,
                    similarity=similarity,
                    pAct_diff=pAct_diff,
                    pAct_a=a.pActivity,
                    pAct_b=b.pActivity,
                    winning_type_a=a.winning_type,
                    winning_type_b=b.winning_type,
                    is_cross_type=a.winning_type != b.winning_type,
                    is_low_quality=not (
                        a.passes_quality_gate and b.passes_quality_gate
                    ),
                    is_same_assay=bool(a.assay_ids & b.assay_ids),
                    involves_inactive=not (a.is_active and b.is_active),
                    landscape_index=_landscape_index(params, pAct_diff, similarity),
                )
            )

    n_undefined = sum(1 for c in cliffs if c.landscape_index is None)
    if params.metric is SeverityMetric.SALI:
        logger.warning(
            "SALI branch: {n} of {t} cliff pairs are undefined at Tm == 1 and will be resolved "
            "by policy '{p}'. SALI rankings are NOT comparable to TS_SALI.",
            n=n_undefined,
            t=len(cliffs),
            p=params.sali_undefined.value,
        )
    logger.info(
        "Detected {n_cliffs} activity cliffs with min_delta {min_delta}.",
        n_cliffs=len(cliffs),
        min_delta=params.min_delta,
    )
    return CliffResult(cliffs=cliffs, params=params, n_undefined=n_undefined)


def high_confidence_cliffs(cliffs: list[ActivityCliffs]) -> list[ActivityCliffs]:
    """Filter a list of cliffs to only those that are high confidence, i.e. both compounds passed the quality gate and share a winning type.

    Args:
        cliffs: The full list from `find_activity_cliffs`, including flagged pairs.

    Returns:
        Only the cliffs where both compounds passed the quality gate and share a winning type.
    """
    return [c for c in cliffs if not c.is_low_quality and not c.is_cross_type]


def same_assay_cliffs(cliffs: list[ActivityCliffs]) -> list[ActivityCliffs]:
    """Filter a list of cliffs to only within-assay pairs.

    Args:
        cliffs: The full list of cliffs identified in the detection pass, including cross-assay pairs.

    Returns:
        Only the cliffs whose two compounds share at least one assay id.
    """
    return [c for c in cliffs if c.is_same_assay]


def same_assay(c: ActivityCliffs) -> bool:
    """Check if a cliff pair is within-assay, i.e. the two compounds share at least one assay id.

    Args:
        c: The cliff pair to test.

    Returns:
        True if the two compounds share at least one assay id.
    """
    return c.is_same_assay


def high_confidence(c: ActivityCliffs) -> bool:
    """Check whether a cliff is high-confidence, i.e. both compounds passed the quality gate and share a winning type.

    Args:
        c: The cliff pair to test.

    Returns:
        True if both compounds passed the quality gate and share a winning type.
    """
    return not c.is_low_quality and not c.is_cross_type


def _ranking_scores(
    result: CliffResult, restrict: Callable[[ActivityCliffs], bool] | None
) -> dict[int, float]:
    """Calculate the scores to rank the cliffs on, and solve the
    undefined values according to the parameters in the `CliffResult`.

    Args:
        result: The `CliffResult` object with the cliffs and the params that drive the ranking.
        restrict: A rule to filter the cliffs to be ranked when True. If None, all cliffs will be ranked.

    Returns:
        A dictionary where the keys are cliff indices and the values are the scores to rank on.
            Indices absent from the dictionary are not ranked and stay unlabelled.
    """
    params = result.params
    cliffs = result.cliffs

    scores: dict[int, float] = {}
    undefined: list[int] = []
    for i, c in enumerate(cliffs):
        if restrict is not None and not restrict(c):
            continue
        if c.landscape_index is None:
            undefined.append(i)
        else:
            scores[i] = c.landscape_index

    if undefined and params.sali_undefined is not SaliUndefined.EXCLUDE:
        if params.sali_undefined is SaliUndefined.MAX_SEVERITY:
            fill = math.inf
        elif scores:
            fill = max(scores.values())
        else:
            fill = math.inf
            logger.warning(
                "NEXT_LARGEST requested but all {u} ranked pairs are undefined; "
                "no finite SALI to borrow, so all are filled with inf and tie at EXTREME.",
                u=len(undefined),
            )
        for i in undefined:
            scores[i] = fill
    return scores


def label_cliff_severity(
    result: CliffResult,
    restrict: Callable[[ActivityCliffs], bool] | None = None,
) -> CliffResult:
    """Assign severity level by ranking the cliffs based on their landscape index.

    Args:
        result: The `CliffResult` object holding all cliffs and their associated thresholds and the SALI undefined
            policy are read from its params.
        restrict: A rule to filter the cliffs to be ranked. If true only cliffs that meet the criteria will be labelled.
            If None, all cliffs will be ranked.
    Returns:
        A new `CliffResult` with `severity_label` filled on the ranked cliffs.
    """
    params = result.params
    if not result.cliffs:
        return result

    scores = _ranking_scores(result, restrict)
    n = len(scores)
    total = len(result.cliffs)
    if n == 0:
        logger.warning(
            "Restriction left 0 rankable cliffs (of {t}); nothing labelled.", t=total
        )
        return result
    if n < total:
        logger.info("Ranking {n} of {t} cliffs after restriction.", n=n, t=total)
    if n < 1 / params.extreme_frac:
        logger.warning(
            "Only {n} ranked cliffs: top {f:.0%} rounds to one pair, so bands are coarse.",
            n=n,
            f=params.extreme_frac,
        )

    order = sorted(
        scores,
        key=lambda i: (-scores[i], result.cliffs[i].id_a, result.cliffs[i].id_b),
    )
    n_extreme = math.ceil(params.extreme_frac * n)
    n_strong = math.ceil(params.strong_frac * n)

    labels: dict[int, CliffSeverity] = {}
    for rank, i in enumerate(order):
        if rank < n_extreme:
            labels[i] = CliffSeverity.EXTREME
        elif rank < n_strong:
            labels[i] = CliffSeverity.STRONG
        else:
            labels[i] = CliffSeverity.MODERATE

    labelled = [
        replace(c, severity_label=labels[i]) if i in labels else c
        for i, c in enumerate(result.cliffs)
    ]
    logger.info(
        "Labelled {n} cliffs (metric={m}): {c}.",
        n=n,
        m=params.metric.value,
        c={sev.value: count for sev, count in Counter(labels.values()).items()},
    )
    return CliffResult(cliffs=labelled, params=params, n_undefined=result.n_undefined)
