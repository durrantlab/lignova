"""Class for detecting activity cliffs in a set of compounds based on their similarities and activities."""

from dataclasses import dataclass

from loguru import logger

from lignova.clustering.tanimoto import TanimotoSimilarities


@dataclass(frozen=True, slots=True)
class CompoundActivity:
    """Per-compound activity handed to cliff detection."""

    pActivity: float
    """Aggregated log-scale activity calculated from the median of the winning type's values."""

    winning_type: str
    """Winning assay type behind pActivity; enables the cross-type check."""

    passes_quality_gate: bool = True
    """False if the winning-type values disagreed by more than ~1 log unit. These cliffs are then
    flagged as low quality rather than excluded."""


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


def find_activity_cliffs(
    sims: TanimotoSimilarities,
    activity: dict[str, CompoundActivity],
    min_delta: float,
) -> list[ActivityCliffs]:
    """Detect activity cliffs in a set of compounds based on their similarities and activities.

    Args:
        sims: Output of `compute_pairwise` over the compounds to analyze.
        activity: A dictionary with the compound ids as the keys and their `CompoundActivity` records as the values.
        min_delta: Minimum absolute difference in pActivity for a pair to be considered a cliff.
    Returns:
        A list of `ActivityCliffs` objects representing the detected activity cliffs.
    """
    cliffs: list[ActivityCliffs] = []
    for id_a, id_b, similarity in sims.edges:
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
        if pAct_diff >= min_delta:
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
                )
            )
    logger.info(
        "Detected {n_cliffs} activity cliffs with min_delta {min_delta}.",
        n_cliffs=len(cliffs),
        min_delta=min_delta,
    )
    return cliffs


def high_confidence_cliffs(cliffs: list[ActivityCliffs]) -> list[ActivityCliffs]:
    """Filter a list of cliffs to only those that are high confidence, i.e. both compounds passed the quality gate and share a winning type.

    Args:
        cliffs: The full list from `find_activity_cliffs`, including flagged pairs.

    Returns:
        Only the cliffs where both compounds passed the quality gate and share a winning type.
    """
    return [c for c in cliffs if not c.is_low_quality and not c.is_cross_type]
