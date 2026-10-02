import pandas as pd
import pytest

from stage_a.indications import PASS, indication_rule, specificity_rank, therapeutic_area_flags
from tests.synth import disease_row


def _rule(*rows: dict) -> dict[str, str]:
    return indication_rule(pd.DataFrame(list(rows)))


@pytest.mark.parametrize("row, expected", [
    (disease_row("MONDO_1", "d", tas=["EFO_0000319"]), PASS),
    (disease_row("Orphanet_1", "d", tas=["EFO_0000319"]), PASS),
    (disease_row("EFO_1", "d", tas=["EFO_0000319", "EFO_0001444"]), PASS),
    (disease_row("HP_1", "d", tas=["EFO_0000319"]), "prefix_not_mondo_efo_orphanet"),
    (disease_row("OTAR_1", "d", tas=["EFO_0000319"]), "prefix_not_mondo_efo_orphanet"),
    (disease_row("GO_1", "d", tas=["EFO_0000319"]), "prefix_not_mondo_efo_orphanet"),
    (disease_row("EFO_0000319", "d", tas=["EFO_0000319"], is_ta=True), "therapeutic_area_root"),
    (disease_row("EFO_2", "d", tas=["EFO_0001444", "EFO_0000651"]), "no_disease_therapeutic_area"),
    (disease_row("MONDO_2", "d", tas=["GO_0008150", "EFO_0002571", "MONDO_0005583"]), "no_disease_therapeutic_area"),
    (disease_row("MONDO_3", "d", tas=[]), "no_disease_therapeutic_area"),
])
def test_indication_rule_each_clause(row, expected):
    assert _rule(row)[row["id"]] == expected


def test_indication_rule_missing_therapeutic_areas_fails():
    row = disease_row("MONDO_4", "d")
    row["therapeuticAreas"] = None
    assert _rule(row)["MONDO_4"] == "no_disease_therapeutic_area"


def test_specificity_rank_counts_descendants():
    rank = specificity_rank(pd.DataFrame([disease_row("MONDO_1", "a", desc=["MONDO_2", "MONDO_3"]),
                                          disease_row("MONDO_2", "b"), disease_row("MONDO_3", "c", desc=[])]))
    assert rank == {"MONDO_1": 2, "MONDO_2": 0, "MONDO_3": 0}


@pytest.mark.parametrize("tas, neuro_psych, neuro_only, psych_only, oncology", [
    ({"MONDO_0005071"}, True, True, False, False),
    ({"MONDO_0002025"}, True, False, True, False),
    ({"MONDO_0005071", "MONDO_0002025"}, True, True, False, False),
    ({"MONDO_0045024"}, False, False, False, True),
    ({"EFO_0000319"}, False, False, False, False),
])
def test_therapeutic_area_flags_partition_neuro_psych(tas, neuro_psych, neuro_only, psych_only, oncology):
    f = therapeutic_area_flags(tas)
    assert (f.neuro_psych, f.neuro_only, f.psych_only, f.oncology) == (neuro_psych, neuro_only, psych_only, oncology)
    assert f.neuro_psych == (f.neuro_only or f.psych_only)
    assert not (f.neuro_only and f.psych_only)
