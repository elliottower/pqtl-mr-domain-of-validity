"""Indication rule, specificity rank and therapeutic-area flags (PREREG §Study design,
"Indication rule", "Specificity rank"; §Measured variables, "Flags", "Covariates").

Ported from `count_v4.py` (`indication_rule`); the specificity rank replaces its ontology depth.
"""
import pandas as pd
from pydantic import BaseModel, ConfigDict

ALLOWED_PREFIXES = frozenset({"MONDO", "EFO", "Orphanet"})
NON_DISEASE_TAS = frozenset({
    "GO_0008150",     # biological process
    "EFO_0000651",    # phenotype
    "EFO_0001444",    # measurement
    "EFO_0002571",    # medical procedure
    "MONDO_0005583",  # non-human animal disease
})
NEURO_TA = "MONDO_0005071"     # nervous system disorder
PSYCH_TA = "MONDO_0002025"     # psychiatric disorder
ONCOLOGY_TA = "MONDO_0045024"  # cancer or benign tumor

PASS = "pass"
NOT_IN_INDEX = "not_in_ot_disease_index"


class TherapeuticAreaFlags(BaseModel):
    model_config = ConfigDict(frozen=True)

    neuro_psych: bool
    neuro_only: bool
    psych_only: bool
    oncology: bool


def _as_set(x) -> set[str]:
    return set(x) if x is not None else set()


def indication_status(disease_id: str, is_therapeutic_area: bool, therapeutic_areas: set[str]) -> str:
    """'pass' or the first failing clause of the indication rule."""
    if disease_id.split("_")[0] not in ALLOWED_PREFIXES:
        return "prefix_not_mondo_efo_orphanet"
    if is_therapeutic_area:
        return "therapeutic_area_root"
    if not (therapeutic_areas - NON_DISEASE_TAS):
        return "no_disease_therapeutic_area"
    return PASS


def indication_rule(disease: pd.DataFrame) -> dict[str, str]:
    """disease id -> 'pass' or failing reason, over the Open Targets disease index.
    An indication absent from the index fails with NOT_IN_INDEX (callers map missing ids)."""
    return {
        did: indication_status(did, bool(ont["isTherapeuticArea"]), _as_set(tas))
        for did, ont, tas in zip(disease["id"], disease["ontology"], disease["therapeuticAreas"])
    }


def descendants_map(disease: pd.DataFrame) -> dict[str, set[str]]:
    return {i: _as_set(d) for i, d in zip(disease["id"], disease["descendants"])}


def specificity_rank(disease: pd.DataFrame) -> dict[str, int]:
    """Number of Open Targets descendants; fewer is more specific. Ties are broken by the lower
    indication ID wherever the rank orders indications (see `specificity_order_key`)."""
    return {i: len(_as_set(d)) for i, d in zip(disease["id"], disease["descendants"])}


def specificity_order_key(indication_id: str, rank: dict[str, int]) -> tuple[int, str]:
    return rank[indication_id], indication_id


def therapeutic_area_flags(therapeutic_areas: set[str]) -> TherapeuticAreaFlags:
    """`neuro_only` is the plan's "neurological (MONDO_0005071)" category of table 13 (it may
    also carry MONDO_0002025); `psych_only` is MONDO_0002025 without MONDO_0005071. The two
    partition `neuro_psych`."""
    neuro = NEURO_TA in therapeutic_areas
    psych = PSYCH_TA in therapeutic_areas
    return TherapeuticAreaFlags(neuro_psych=neuro or psych, neuro_only=neuro,
                                psych_only=psych and not neuro,
                                oncology=ONCOLOGY_TA in therapeutic_areas)


def therapeutic_area_map(disease: pd.DataFrame) -> dict[str, TherapeuticAreaFlags]:
    return {i: therapeutic_area_flags(_as_set(t)) for i, t in zip(disease["id"], disease["therapeuticAreas"])}
