import pandas as pd
import pytest

from stage_a.mechanism import OTHER_CONFLICT, direction_of, label_row, resolve_program_targets, token_direction


@pytest.mark.parametrize("token, expected", [
    ("INVERSE AGONIST", "BLOCKING"),
    ("  inverse   agonist ", "BLOCKING"),
    ("NEGATIVE ALLOSTERIC MODULATOR", "BLOCKING"),
    ("POSITIVE ALLOSTERIC MODULATOR", "ACTIVATING"),
    ("PARTIAL AGONIST", "ACTIVATING"),
    ("ALLOSTERIC ANTAGONIST", "BLOCKING"),
    ("DISRUPTING AGENT", "BLOCKING"),
    ("ANTISENSE INHIBITOR", "BLOCKING"),
    ("DEGRADER", "BLOCKING"),
    ("SUPPRESSOR", "BLOCKING"),
    ("AGONIST", "ACTIVATING"),
    ("OPENER", "ACTIVATING"),
    ("STABILISER", "ACTIVATING"),
    ("RELEASING AGENT", "ACTIVATING"),
    ("NEGATIVE MODULATOR", "AMBIGUOUS"),
    ("MODULATOR", "AMBIGUOUS"),
    ("BINDING AGENT", "AMBIGUOUS"),
    ("OTHER", "AMBIGUOUS"),
    ("", "AMBIGUOUS"),
    (None, "AMBIGUOUS"),
])
def test_token_direction_tiers_in_order(token, expected):
    assert token_direction(token) == expected


def test_direction_maps_blocking_to_decrease_and_activating_to_increase():
    assert direction_of("INHIBITOR") == "decrease"
    assert direction_of("AGONIST") == "increase"
    assert direction_of("NEGATIVE MODULATOR") == "ambiguous"


@pytest.mark.parametrize("mol, action, blood, cls", [
    ("Antibody", "INHIBITOR", True, "aligned"),
    ("Protein", "ANTAGONIST", True, "aligned"),
    ("Enzyme", "NEGATIVE ALLOSTERIC MODULATOR", True, "aligned"),
    ("Antibody", "INHIBITOR", False, "other"),
    ("Antibody", "BLOCKER", True, "other"),
    ("Antibody", "AGONIST", True, "other"),
    ("Oligonucleotide", "ANTISENSE INHIBITOR", True, "aligned"),
    ("Oligonucleotide", "RNAI INHIBITOR", True, "aligned"),
    ("Oligonucleotide", "RNAI INHIBITOR", False, "other"),
    ("Oligonucleotide", "INHIBITOR", True, "other"),
    ("Small molecule", "DEGRADER", True, "aligned"),
    ("Antibody", "DEGRADER", True, "aligned"),
    ("Small molecule", "DEGRADER", False, "other"),
    ("Small molecule", "INHIBITOR", False, "blocking"),
    ("Small molecule", "INHIBITOR", True, "blocking"),
    ("Small molecule", "ANTAGONIST", False, "blocking"),
    ("Small molecule", "BLOCKER", False, "blocking"),
    ("Small molecule", "NEGATIVE ALLOSTERIC MODULATOR", False, "blocking"),
    ("Small molecule", "INVERSE AGONIST", False, "blocking"),
    ("Small molecule", "AGONIST", False, "other"),
    ("Small molecule", "ALLOSTERIC ANTAGONIST", False, "other"),
    ("Small molecule", "DISRUPTING AGENT", False, "other"),
    ("Small molecule", "NEGATIVE MODULATOR", False, "other"),
    ("Antibody drug conjugate", "INHIBITOR", True, "other"),
    ("Unknown", "INHIBITOR", True, "other"),
    ("Protein", "AGONIST", True, "other"),
])
def test_mechanism_class_table(mol, action, blood, cls):
    assert label_row(mol, action, "text", blood).mechanism_class == cls


def test_aligned_and_blocking_rows_always_decrease():
    for mol in ("Antibody", "Protein", "Enzyme", "Oligonucleotide", "Small molecule", "Unknown"):
        for action in ("INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "INVERSE AGONIST",
                       "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "DEGRADER", "AGONIST", "BINDING AGENT"):
            for blood in (True, False):
                lab = label_row(mol, action, "x", blood)
                if lab.mechanism_class in ("aligned", "blocking"):
                    assert lab.direction == "decrease"


@pytest.mark.parametrize("text", ["IL6 neutralising antibody", "Neutralizing anti-X", "X SEQUESTRANT",
                                  "complement C5 inhibitor", "TNF antagonist", "channel blocker"])
def test_neutralizing_binding_agent_is_inhibitor(text):
    lab = label_row("Antibody", "BINDING AGENT", text, True)
    assert (lab.mechanism_class, lab.direction, lab.neutralizing_binding_agent) == ("aligned", "decrease", True)
    sm = label_row("Small molecule", "BINDING AGENT", text, False)
    assert (sm.mechanism_class, sm.direction) == ("blocking", "decrease")


@pytest.mark.parametrize("text", ["TNF binding agent", "anti-X antibody", "", None])
def test_other_binding_agent_is_other_and_ambiguous(text):
    lab = label_row("Antibody", "BINDING AGENT", text, True)
    assert (lab.mechanism_class, lab.direction, lab.neutralizing_binding_agent) == ("other", "ambiguous", False)


def test_non_binding_agent_action_never_neutralizing_by_text():
    lab = label_row("Antibody", "AGONIST", "neutralizing inhibitor", True)
    assert (lab.mechanism_class, lab.direction, lab.neutralizing_binding_agent) == ("other", "increase", False)


def _rows(*triples) -> pd.DataFrame:
    return pd.DataFrame([{"program": "P", "gene": "G", "row_class": c, "row_direction": d, "row_reason": r}
                         for c, d, r in triples])


def test_conflicting_informative_rows_make_program_target_other():
    r = resolve_program_targets(_rows(("blocking", "decrease", "row4"), ("other", "decrease", "x"))).iloc[0]
    assert (r["pt_class"], r["pt_direction"], r["pt_reason"], r["conflict"]) == ("other", "ambiguous", OTHER_CONFLICT, True)


def test_direction_conflict_makes_program_target_other():
    r = resolve_program_targets(_rows(("blocking", "decrease", "row4"), ("other", "increase", "x"))).iloc[0]
    assert (r["pt_class"], r["pt_direction"], r["conflict"]) == ("other", "ambiguous", True)


def test_ambiguous_rows_are_not_informative():
    r = resolve_program_targets(_rows(("blocking", "decrease", "row4"), ("other", "ambiguous", "ba"))).iloc[0]
    assert (r["pt_class"], r["pt_direction"], r["pt_reason"], r["conflict"]) == ("blocking", "decrease", "row4", False)


def test_only_ambiguous_rows_give_other_without_conflict():
    r = resolve_program_targets(_rows(("other", "ambiguous", "ba"), ("other", "ambiguous", "mod"))).iloc[0]
    assert (r["pt_class"], r["pt_direction"], r["conflict"]) == ("other", "ambiguous", False)


def test_agreeing_rows_keep_their_class():
    r = resolve_program_targets(_rows(("aligned", "decrease", "row1"), ("aligned", "decrease", "row2"))).iloc[0]
    assert (r["pt_class"], r["pt_direction"], r["pt_reason"], r["conflict"]) == ("aligned", "decrease", "row1;row2", False)


# ---- the action type is read as Open Targets gives it (count_v4.py lines 144-162) -------------------

@pytest.mark.parametrize("action", ["inhibitor", " INHIBITOR", "INHIBITOR ", "Inhibitor"])
def test_the_class_table_matches_the_action_type_exactly_while_the_direction_rule_normalizes_it(action):
    lab = label_row("Small molecule", action, "x inhibitor", False)
    assert (lab.mechanism_class, lab.direction) == ("other", "decrease")
    assert label_row("Small molecule", "INHIBITOR", "x inhibitor", False).mechanism_class == "blocking"


@pytest.mark.parametrize("action", ["binding agent", "BINDING  AGENT", " BINDING AGENT"])
def test_a_binding_agent_is_neutralizing_only_under_the_exact_action_type(action):
    lab = label_row("Antibody", action, "X neutralizing antibody", True)
    assert (lab.mechanism_class, lab.direction, lab.neutralizing_binding_agent) == ("other", "ambiguous", False)


@pytest.mark.parametrize("action", [None, float("nan"), ""])
def test_a_missing_action_type_is_other_and_ambiguous(action):
    lab = label_row("Small molecule", action, None, True)
    assert (lab.mechanism_class, lab.direction, lab.neutralizing_binding_agent) == ("other", "ambiguous", False)
