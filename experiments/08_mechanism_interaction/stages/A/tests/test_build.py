import hashlib
from datetime import date

import pandas as pd
import pytest

from stage_a.build import run_stage_a
from stage_a.flags import PUBLICATION_DATES, earliest_publication, require_publication_dates
from stage_a.schemas import HYPOTHESIS_COLUMNS
from tests.synth import DATES, World, gcat_study, moa, world_inputs


@pytest.fixture(scope="module")
def out():
    return run_stage_a(world_inputs(World()))


@pytest.fixture(scope="module")
def hyp(out) -> pd.DataFrame:
    return out.hypotheses.set_index(["gene_symbol", "indication_id", "mechanism_class", "direction"], drop=False)


def _one(hyp: pd.DataFrame, gene: str, ind: str) -> pd.Series:
    rows = hyp[(hyp["gene_symbol"] == gene) & (hyp["indication_id"] == ind)]
    assert len(rows) == 1, rows
    return rows.iloc[0]


def test_hypotheses_columns_match_the_interface(out):
    assert list(out.hypotheses.columns) == HYPOTHESIS_COLUMNS


def test_no_output_column_or_value_carries_phase(out):
    stage_words = ("PHASE_", "APPROVAL", "PREAPPROVAL", "EARLY_PHASE", "PRECLINICAL")
    for name in ("hypotheses", "funnel", "outcome_gwas_selection", "mechanism_crosstab"):
        df = getattr(out, name)
        assert not [c for c in df.columns if "phase" in c.lower() or "stage" in c.lower()], name
        text = df.astype(str).to_csv(index=False)
        assert not [w for w in stage_words if w in text], name
    assert out.hypotheses["eligible"].map(type).eq(bool).all()


def test_outputs_do_not_depend_on_phase_beyond_phase2_or_later():
    base = World()
    shifted = World(clinical=[(d, i, {"PHASE_2": "APPROVAL", "PHASE_3": "PHASE_2", "APPROVAL": "PHASE_2_3",
                                      "PHASE_2_3": "PREAPPROVAL", "PREAPPROVAL": "PHASE_3"}.get(s, s))
                              for d, i, s in base.clinical])
    a, b = run_stage_a(world_inputs(base)), run_stage_a(world_inputs(shifted))
    for name in ("hypotheses", "funnel", "outcome_gwas_selection", "mechanism_crosstab"):
        pd.testing.assert_frame_equal(getattr(a, name), getattr(b, name))


def test_s1_membership(out):
    h = out.hypotheses
    s1 = h[h["in_s1"] & h["heldout"]]
    pooled_only = h[h["in_s1"] & ~h["heldout"]]
    assert set(zip(pooled_only["gene_symbol"], pooled_only["indication_id"])) == {("GA", "MONDO_5000")}
    assert h.loc[~h["heldout"], "in_s1"].all()
    assert set(zip(s1["gene_symbol"], s1["indication_id"], s1["mechanism_class"])) == {
        ("GA", "MONDO_1001", "aligned"), ("GA", "MONDO_3000", "aligned"), ("GB", "MONDO_1001", "blocking"),
        ("GB", "MONDO_1002", "blocking"), ("GB", "MONDO_2000", "other"), ("GB", "MONDO_1001", "other"),
        ("GG", "MONDO_2100", "blocking")}
    assert s1["heldout"].all() and (s1["outcome_tier"] == 1).all() and (s1["overlap"] == "no").all()
    assert s1["strict_instrument"].all() and not s1["conflicted_row_restored"].any()
    assert s1["in_s4"].all() and s1["in_s10"].all() and s1["eligible"].all()


def test_hypothesis_id_is_sha256_prefix_of_key(out):
    for r in out.hypotheses.itertuples():
        key = f"{r.gene_ensembl}|{r.indication_id}|{r.direction}|{r.mechanism_class}"
        key = key if r.variant == "primary" else f"{key}|{r.variant}"
        assert r.hypothesis_id == hashlib.sha256(key.encode()).hexdigest()[:16]
    assert out.hypotheses["hypothesis_id"].is_unique


def test_salt_counts_as_parent_program(hyp):
    assert _one(hyp, "GB", "MONDO_1002")["drug_program_ids"] == "SM2"


def test_ancestor_dropped_by_related_rule_and_counted(out, hyp):
    assert not ((hyp["gene_symbol"] == "GB") & (hyp["indication_id"] == "MONDO_1000")).any()
    f = out.funnel.set_index("step")
    assert f.loc["related_indication_rule", "excluded"] == 1


def test_phase1_only_hypothesis_is_not_eligible_and_counted(out, hyp):
    assert not ((hyp["gene_symbol"] == "GB") & (hyp["direction"] == "increase")).any()
    assert out.funnel.set_index("step").loc["eligible_phase2_or_later", "excluded"] == 1


def test_indication_rule_exclusions_counted(out):
    r = out.funnel.set_index("step").loc["indication_rule"]
    assert r["excluded"] == 2
    assert r["reason"] == "no_disease_therapeutic_area=1;prefix_not_mondo_efo_orphanet=1"


def test_conflicted_program_target_never_in_both_h1_classes(out):
    h = out.hypotheses[~out.hypotheses["conflicted_row_restored"]]
    for (g, i), grp in h.groupby(["gene_symbol", "indication_id"]):
        aligned = {p for s in grp.loc[grp["mechanism_class"] == "aligned", "drug_program_ids"] for p in s.split(";")}
        blocking = {p for s in grp.loc[grp["mechanism_class"] == "blocking", "drug_program_ids"] for p in s.split(";")}
        assert not aligned & blocking, (g, i)
    conflicted = h[h["drug_program_ids"].str.split(";").map(lambda s: "SM3" in s)]
    assert list(zip(conflicted["mechanism_class"], conflicted["direction"], conflicted["class_reason"])) == [
        ("other", "ambiguous", "other:conflicting_mechanism_rows")]


def test_conflicted_rows_restored_only_under_s21(out):
    r = out.hypotheses[out.hypotheses["conflicted_row_restored"]]
    assert set(zip(r["gene_symbol"], r["indication_id"], r["mechanism_class"], r["direction"])) == {
        ("GB", "MONDO_2000", "blocking", "decrease"), ("GB", "MONDO_2000", "other", "decrease")}
    assert not r["in_s1"].any() and not r["in_s4"].any()


def test_conflict_in_blocking_and_aligned_rows_keeps_program_out_of_h1():
    w = World()
    w.molecules.append({"id": "BI1", "drugType": "Antibody", "parentId": None})
    w.molecules.append({"id": "BI1S", "drugType": "Small molecule", "parentId": "BI1"})
    w.mechanisms.append(moa("BI1", "INHIBITOR", "GA inhibitor", "GA"))
    w.clinical += [("BI1", "MONDO_1002", "PHASE_2"), ("BI1S", "MONDO_1002", "PHASE_2")]
    w.gwas_catalog.append(gcat_study("GCST1002", "MONDO_1002", "EstBB", 4000, 40000))
    h = run_stage_a(world_inputs(w)).hypotheses
    ga = h[(h["gene_symbol"] == "GA") & (h["indication_id"] == "MONDO_1002") & ~h["conflicted_row_restored"]]
    assert list(zip(ga["mechanism_class"], ga["direction"], ga["drug_program_ids"])) == [("other", "ambiguous", "BI1")]


def test_pilot_key_written_but_not_held_out(hyp):
    r = _one(hyp, "GA", "MONDO_5000")
    assert (r["heldout"], r["in_s1"], r["in_s12"], r["pilot_indication"]) == (False, True, False, True)


def test_s4_adds_tier1_overlap_unknown(hyp):
    r = _one(hyp, "GD", "MONDO_7000")
    assert (r["in_s1"], r["in_s4"], r["in_s10"], r["overlap"], r["overlap_unknown"], r["sample_overlap"]) == (
        False, True, False, "unknown", True, False)
    assert (r["mechanism_class"], r["class_reason"]) == ("aligned", "row1_biologic_inhibitor_blood_secreted")


def test_s10_adds_tier2_with_subtype_flag(hyp):
    r = _one(hyp, "GD", "MONDO_6000")
    assert (r["in_s1"], r["in_s4"], r["in_s10"], r["outcome_tier"], r["subtype_restricted"], r["outcome_accession"]) == (
        False, False, True, 2, True, "GCST61")


def test_decode_instrument_skips_icelandic_outcome_in_pipeline(out, hyp):
    r = _one(hyp, "GF", "MONDO_8000")
    assert (r["instrument_source"], r["platform"], r["instrument_assay_id"]) == ("decode", "SomaScan", "2000-2")
    assert (r["outcome_accession"], r["outcome_tier"], r["overlap"]) == ("FINNGEN_R12_Y81", 2, "no")
    sel = out.outcome_gwas_selection
    ice = sel[sel["includes_iceland"] == "yes"]
    for acc in ice["outcome_accession"]:
        assert not ((out.hypotheses["instrument_source"] == "decode") & (out.hypotheses["outcome_accession"] == acc)).any()


def test_independent_candidate_beats_larger_overlapping_one(hyp):
    assert _one(hyp, "GA", "MONDO_3000")["outcome_accession"] == "FINNGEN_R12_C3"
    r = hyp[(hyp["gene_symbol"] == "GA") & (hyp["indication_id"] == "MONDO_1001")].iloc[0]
    assert (r["outcome_accession"], r["outcome_n_case"], r["outcome_n_control"]) == ("GCST1", 5000, 50000)
    assert r["outcome_neff"] == pytest.approx(4 / (1 / 5000 + 1 / 50000))


def test_instrument_source_order_and_assay(hyp):
    gb = hyp[(hyp["gene_symbol"] == "GB") & (hyp["indication_id"] == "MONDO_1001")]
    assert set(gb["instrument_source"]) == {"decode"} and set(gb["instrument_assay_id"]) == {"1000-1"}
    ga = _one(hyp, "GA", "MONDO_3000")
    assert (ga["instrument_source"], ga["platform"], ga["instrument_assay_id"]) == ("ukbppp", "Olink", "OID001")
    assert _one(hyp, "GG", "MONDO_2100")["instrument_assay_id"] == "OID007"


def test_s20_rows_only_under_inclusive_lists(out):
    h = out.hypotheses
    s20 = h[~h["strict_instrument"]]
    assert set(zip(s20["gene_symbol"], s20["instrument_source"], s20["instrument_assay_id"])) == {
        ("GC", "ukbppp", "OID003"), ("GE", "interval", "GE.1;epigraphdb:no_assay")}
    assert not s20["in_s1"].any()
    assert not bool(h.loc[h["gene_symbol"] == "GC", "single_protein_row"].iloc[0])


def test_localization_and_area_flags(hyp):
    ga = _one(hyp, "GA", "MONDO_3000")
    assert (ga["blood_secreted_hpa"], ga["secreted_uniprot"], ga["oncology"], ga["karim_launched"]) == (True, True, True, True)
    assert not bool(_one(hyp, "GA", "MONDO_5000")["karim_launched"])
    gb = hyp[(hyp["gene_symbol"] == "GB") & (hyp["indication_id"] == "MONDO_2000")].iloc[0]
    assert (gb["neuro_psych"], gb["neuro_only"], gb["psych_only"], gb["blood_secreted_hpa"]) == (True, True, False, False)
    gg = _one(hyp, "GG", "MONDO_2100")
    assert (gg["neuro_psych"], gg["neuro_only"], gg["psych_only"]) == (True, False, True)


def test_pre_pqtl_publication_date_is_earliest_listing_source(hyp):
    assert _one(hyp, "GA", "MONDO_3000")["pre_pqtl_publication_date"] == DATES["ukbppp"]
    gb = hyp[(hyp["gene_symbol"] == "GB") & (hyp["indication_id"] == "MONDO_1001")]
    assert set(gb["pre_pqtl_publication_date"]) == {DATES["interval"]}
    assert _one(hyp, "GF", "MONDO_8000")["pre_pqtl_publication_date"] == DATES["decode"]


def test_pre_pqtl_publication_date_follows_the_dates_given():
    w = World(publication_dates={"ukbppp": date(2001, 1, 1), "decode": date(2020, 2, 2), "interval": date(2010, 1, 1)})
    h = run_stage_a(world_inputs(w)).hypotheses
    r = h[(h["gene_symbol"] == "GA") & (h["indication_id"] == "MONDO_3000")].iloc[0]
    assert r["pre_pqtl_publication_date"] == date(2001, 1, 1)


def test_s12_keeps_one_hypothesis_per_program(hyp):
    s12 = hyp[hyp["in_s12"]]
    assert not (s12["in_s12"] & ~s12["in_s1"]).any()
    assert set(zip(s12["gene_symbol"], s12["indication_id"], s12["drug_program_ids"])) == {
        ("GA", "MONDO_1001", "AB1"), ("GB", "MONDO_1001", "SM2"), ("GB", "MONDO_1001", "AB11"),
        ("GB", "MONDO_2000", "SM3"), ("GG", "MONDO_2100", "SM10")}


def test_outcome_selection_table_lists_every_tier_and_marks_selected(out):
    sel = out.outcome_gwas_selection
    m1001 = sel[sel["indication_id"] == "MONDO_1001"]
    assert set(zip(m1001["outcome_accession"], m1001["tier"])) == {("GCST1", 1), ("ieu-a-1", 1), ("GCST0", 3)}
    assert dict(zip(m1001["outcome_accession"], m1001["selected"])) == {"GCST1": True, "ieu-a-1": False, "GCST0": False}
    assert "GCSTSMALL" not in set(sel["outcome_accession"]) and "GCSTEAS" not in set(sel["outcome_accession"])


def test_mechanism_crosstab_counts_s1_rows(out):
    x = out.mechanism_crosstab
    assert x["n_hypotheses_s1"].sum() >= 7
    conflicted = x[x["conflicted"]]
    assert set(zip(conflicted["action_type"], conflicted["row_class"], conflicted["program_target_class"])) == {
        ("INHIBITOR", "blocking", "other"), ("DISRUPTING AGENT", "other", "other")}
    dup = x[(x["action_type"] == "INHIBITOR") & (x["molecule_type"] == "Small molecule") & ~x["conflicted"]
            & (x["program_target_class"] == "blocking")]
    assert dup["n_mechanism_rows"].sum() == 1  # SM2's duplicated row collapses; its salt shares the parent's row


def test_the_frozen_publication_dates_are_the_journal_online_dates():
    assert require_publication_dates(PUBLICATION_DATES) == {
        "ukbppp": date(2023, 10, 4), "decode": date(2021, 12, 2), "interval": date(2018, 6, 6)}


def test_a_trial_precedes_every_source_exactly_when_it_precedes_the_earliest():
    dates = require_publication_dates(PUBLICATION_DATES)
    both = ["ukbppp", "decode"]
    assert earliest_publication(both, dates) == date(2021, 12, 2)
    for start in (date(2021, 12, 1), date(2021, 12, 2), date(2022, 6, 1), date(2023, 10, 4)):
        assert (start < earliest_publication(both, dates)) == all(start < dates[s] for s in both)
    assert earliest_publication(["ukbppp", "decode", "interval"], dates) == date(2018, 6, 6)


# ---- class and direction come from Phase II+ rows, as in count_v4.py (lines 595, 630-632) ------------

def test_a_child_molecule_below_phase2_adds_no_mechanism_row_to_its_parents_program():
    """SM2X is a child of SM2 with its own mechanism row (an agonist) and one Phase I indication.
    count_v4.py drops its clinical_indication row before any mechanism row is derived, so SM2 on GB
    stays a plain inhibitor: blocking, decrease, no conflict."""
    base = run_stage_a(world_inputs(World()))
    w = World()
    w.molecules.append({"id": "SM2X", "drugType": "Small molecule", "parentId": "SM2"})
    w.mechanisms.append(moa("SM2X", "AGONIST", "GB agonist", "GB"))
    w.clinical.append(("SM2X", "MONDO_1002", "PHASE_1"))
    out = run_stage_a(world_inputs(w))
    pd.testing.assert_frame_equal(out.hypotheses, base.hypotheses)
    pd.testing.assert_frame_equal(out.mechanism_crosstab, base.mechanism_crosstab)
    f, fb = out.funnel.set_index("step"), base.funnel.set_index("step")
    assert f.loc["eligible_phase2_or_later", "remaining"] == fb.loc["eligible_phase2_or_later", "remaining"]
    assert f.loc["related_indication_rule", "remaining"] == fb.loc["related_indication_rule", "remaining"]


def test_the_same_child_at_phase2_does_conflict_with_its_parent():
    w = World()
    w.molecules.append({"id": "SM2X", "drugType": "Small molecule", "parentId": "SM2"})
    w.mechanisms.append(moa("SM2X", "AGONIST", "GB agonist", "GB"))
    w.clinical.append(("SM2X", "MONDO_1002", "PHASE_2"))
    h = run_stage_a(world_inputs(w)).hypotheses
    sm2 = h[h["drug_program_ids"].str.split(";").map(lambda s: "SM2" in s) & ~h["conflicted_row_restored"]]
    assert set(zip(sm2["mechanism_class"], sm2["direction"], sm2["class_reason"])) == {
        ("other", "ambiguous", "other:conflicting_mechanism_rows")}


def test_a_program_with_no_phase2_row_anywhere_still_forms_its_ineligible_hypothesis():
    out = run_stage_a(world_inputs(World()))          # SM5: an agonist of GB, Phase I only
    f = out.funnel.set_index("step")
    assert f.loc["hypotheses", "remaining"] - f.loc["eligible_phase2_or_later", "remaining"] == 1
    assert not (out.hypotheses["drug_program_ids"].str.split(";").map(lambda s: "SM5" in s)).any()
