import pandas as pd

from stage_a.related import related_indication_rule, s12_one_per_program

DESC = {"A": {"A1", "A2", "A11"}, "A1": {"A11"}, "A2": set(), "A11": set(), "B": set()}
RANK = {k: len(v) for k, v in DESC.items()}


def _hyp(*rows) -> pd.DataFrame:
    return pd.DataFrame([{"gene": g, "direction": d, "mechanism_class": c, "indication_id": i, "programs": frozenset(p)}
                         for g, d, c, i, p in rows])


def _kept(h: pd.DataFrame) -> set[str]:
    return set(h.loc[related_indication_rule(h, DESC, RANK), "indication_id"])


def test_ancestor_dropped_when_every_program_on_retained_descendants():
    h = _hyp(("G", "decrease", "blocking", "A", {"P1", "P2"}), ("G", "decrease", "blocking", "A1", {"P1"}),
             ("G", "decrease", "blocking", "A2", {"P2"}))
    assert _kept(h) == {"A1", "A2"}


def test_ancestor_kept_when_one_program_is_not_on_a_descendant():
    h = _hyp(("G", "decrease", "blocking", "A", {"P1", "P3"}), ("G", "decrease", "blocking", "A1", {"P1"}),
             ("G", "decrease", "blocking", "A2", {"P2"}))
    assert _kept(h) == {"A", "A1", "A2"}


def test_non_descendant_does_not_cover():
    h = _hyp(("G", "decrease", "blocking", "A", {"P1"}), ("G", "decrease", "blocking", "B", {"P1"}))
    assert _kept(h) == {"A", "B"}


def test_descendants_in_another_class_direction_or_gene_do_not_cover():
    h = _hyp(("G", "decrease", "blocking", "A", {"P1"}), ("G", "decrease", "other", "A1", {"P1"}),
             ("G", "increase", "blocking", "A2", {"P1"}), ("H", "decrease", "blocking", "A11", {"P1"}))
    assert _kept(h) == {"A", "A1", "A2", "A11"}


def test_dropped_descendant_program_still_counts_through_its_retained_descendant():
    h = _hyp(("G", "decrease", "blocking", "A", {"P1"}), ("G", "decrease", "blocking", "A1", {"P1"}),
             ("G", "decrease", "blocking", "A11", {"P1"}))
    assert _kept(h) == {"A11"}


def test_program_counted_only_at_retained_descendants():
    h = _hyp(("G", "decrease", "blocking", "A", {"P1", "P2"}), ("G", "decrease", "blocking", "A1", {"P1", "P2"}),
             ("G", "decrease", "blocking", "A11", {"P1"}))
    assert _kept(h) == {"A1", "A11"}


def test_s12_orders_by_specificity_then_lower_indication_id_then_program():
    s1 = pd.DataFrame([
        {"gene": "G", "indication_id": "MONDO_2", "programs": frozenset({"P1"})},
        {"gene": "G", "indication_id": "MONDO_1", "programs": frozenset({"P1", "P2"})},
        {"gene": "G", "indication_id": "MONDO_0", "programs": frozenset({"P2"})},
        {"gene": "G", "indication_id": "MONDO_3", "programs": frozenset({"P3"})},
        {"gene": "H", "indication_id": "MONDO_2", "programs": frozenset({"P1"})},
    ])
    rank = {"MONDO_0": 5, "MONDO_1": 0, "MONDO_2": 0, "MONDO_3": 9}
    keep = s12_one_per_program(s1, rank)
    assert list(keep) == [False, True, False, True, True]
