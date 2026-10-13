import json

import numpy as np
import pandas as pd
import pytest

from followup.common import FROZEN_PLAN_SHA256, FollowupError, require_run, wilson
from followup.data import b_and_c, check_ensembl, check_unit_integrity, coloc_class, max_by_assay, midrank_quintile, s05_state
from followup.q1 import design, fit, make_spec, rcs_term, run_bootstrap, standardize
from followup.q2_q3 import ba_bootstrap, ba_by_class, balanced_accuracy
from followup.q4_q5 import augment, check_evidence, g_columns, support
from followup.q7 import survival_table


def ev_frame(rows):
    cols = ["coloc_run", "pp_h0", "pp_h1", "pp_h2", "pp_h3", "pp_h4"]
    return pd.DataFrame(rows, columns=cols)


def test_coloc_class_follows_the_order_bounds_and_tie_rule():
    df = ev_frame([
        (True, 0.0, 0.0, 0.0, 0.2, 0.8),      # shared at the bound
        (True, 0.0, 0.0, 0.0, 0.8, 0.2),      # distinct
        (True, 0.0, 0.0, 0.8, 0.1, 0.1),      # outcome_only
        (True, 0.5, 0.3, 0.1, 0.05, 0.05),    # no_outcome_signal (H0 + H1 = 0.8)
        (True, 0.2, 0.2, 0.2, 0.4, 0.0),      # unresolved: H0+H1 0.4 ties H3 0.4 -> the earlier label
        (True, 0.1, 0.1, 0.1, 0.0, 0.7),      # unresolved, H4 largest
        (False, None, None, None, None, None),
    ])
    assert list(coloc_class(df)) == ["shared", "distinct", "outcome_only", "no_outcome_signal", "unresolved_h0h1",
                                     "unresolved_h4", "not_run"]


def test_a_completed_run_with_a_missing_posterior_stops():
    with pytest.raises(FollowupError, match="missing posterior"):
        coloc_class(ev_frame([(True, 0.1, 0.1, 0.1, None, 0.7)]))


def test_b_and_c_including_the_zero_h3_case_and_the_c_floor():
    df = ev_frame([(True, 0.5, 0.0, 0.0, 0.0, 0.5), (True, 0.99, 0.0, 0.0, 0.004, 0.004), (True, 0.0, 0.0, 0.0, 0.6, 0.2),
                   (False, None, None, None, None, None)])
    b, c = b_and_c(df)
    assert list(b.round(6)[:3]) == [0.5, 0.008, 0.8] and np.isnan(b[3])
    assert c[0] == 1.0 and np.isnan(c[1]) and c[2] == pytest.approx(0.25) and np.isnan(c[3])


def test_midrank_quintile_puts_ties_together_and_cut_points_down():
    q = midrank_quintile(pd.Series([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, None]))
    assert list(q[:10]) == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5] and pd.isna(q[10])
    tied = midrank_quintile(pd.Series([5.0] * 6 + [9.0] * 4))
    assert set(tied[:6]) == {2} and set(tied[6:]) == {5}      # midrank 3.5/10 = 0.35; 8.5/10 = 0.85


def test_unit_integrity_rejects_a_unit_whose_rows_disagree():
    base = {"instrument_source": "ukbppp", "instrument_assay_id": "OID1", "gene_ensembl": "ENSG00000000001",
            "outcome_accession": "GCST1", "coloc_run": True, "not_run_reason": None, "n_shared": 100,
            "lead_variant": "rs1", "genetic_direction": 1, "pp_h0": 0.1, "pp_h1": 0.1, "pp_h2": 0.1, "pp_h3": 0.1, "pp_h4": 0.6}
    ok = pd.DataFrame([base, base])
    check_unit_integrity(ok)
    for col, val in (("pp_h4", 0.6 + 1e-9), ("lead_variant", "rs2"), ("genetic_direction", -1), ("pp_h3", None)):
        bad = pd.DataFrame([base, {**base, col: val}])
        with pytest.raises(FollowupError):
            check_unit_integrity(bad)


def test_ensembl_ids_with_a_version_suffix_stop():
    check_ensembl(pd.Series(["ENSG00000117020"]), "x")
    with pytest.raises(FollowupError):
        check_ensembl(pd.Series(["ENSG00000117020.5"]), "x")


def test_max_by_assay_takes_the_maximum_and_flags_truncation_of_the_maximum_only():
    out = max_by_assay([("a", 10), ("a", "inf"), ("b", 400), ("b", 12), ("c", 5), ("c", 7), ("d", "n/a")])
    assert out["a"] == {"value": 300.0, "rows": 2, "truncated": True}
    assert out["b"] == {"value": 300.0, "rows": 2, "truncated": True}
    assert out["c"] == {"value": 7.0, "rows": 2, "truncated": False}
    assert "d" not in out


def test_s05_uses_pp_h4_at_half_and_the_v8_direction_rule():
    df = pd.DataFrame({"coloc_run": [True, True, True, True, False], "pp_h4": [0.5, 0.49, 0.7, 0.9, np.nan],
                       "genetic_direction": [1, 1, -1, 0, 1], "direction": ["decrease", "decrease", "decrease", "decrease", "decrease"]})
    assert list(s05_state(df)) == ["supportive", "inconclusive", "contradictory", "inconclusive", "inconclusive"]


def test_rcs_term_is_zero_below_the_first_knot_and_linear_beyond_the_last():
    knots = (1.0, 2.0, 4.0)
    assert rcs_term(np.array([0.0, 0.5, 1.0]), knots) == pytest.approx([0, 0, 0])
    x = np.array([5.0, 6.0, 7.0, 8.0])
    assert np.diff(rcs_term(x, knots), 2) == pytest.approx([0, 0], abs=1e-9)


def synthetic_runs(n_genes=60, per_gene=8, seed=1):
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_genes):
        source = "ukbppp" if g % 3 else "decode"
        for _ in range(per_gene):
            rows.append({"gene_ensembl": f"ENSG{g:011d}", "instrument_source": source,
                         "platform": "Olink" if source == "ukbppp" else "SomaScan",
                         "mechanism_class": rng.choice(["aligned", "blocking", "other"]),
                         "oncology": bool(rng.random() < 0.4), "neuro_psych": bool(rng.random() < 0.2),
                         "log10_neff": rng.normal(4.3, 0.6), "assay_cis_signal_strength": rng.uniform(10, 200),
                         "coloc_run": True})
    df = pd.DataFrame(rows)
    eta = -1.0 + 0.8 * (df["log10_neff"] - 4.3) + 0.5 * (df["instrument_source"] == "decode")
    df["B"] = 1 / (1 + np.exp(-(eta + rng.normal(0, 0.5, len(df)))))
    return df


def test_platform_aliased_with_source_is_dropped_and_the_fit_recovers_the_signal():
    df = synthetic_runs()
    spec = make_spec(df, "B", with_strength=True)
    assert spec.drop_platform
    res = fit(df, spec)
    std = standardize(res, df, spec)
    assert std["log10_neff@p90-p10"] > 0.1
    assert std["instrument_source=ukbppp-decode"] < 0


def test_standardized_level_mean_is_the_average_prediction_with_every_row_set_to_that_level():
    df = synthetic_runs()
    spec = make_spec(df, "B", with_strength=False)
    res = fit(df, spec)
    X, _ = design(df.assign(mechanism_class="blocking"), spec)
    assert standardize(res, df, spec)["mechanism_class=blocking"] == pytest.approx(float(res.predict(X).mean()))


def test_a_level_absent_from_a_replicate_is_a_failed_fit_not_a_stop():
    df = synthetic_runs()
    spec = make_spec(df, "B", with_strength=False)
    assert fit(df.loc[df["instrument_source"] == "ukbppp"], spec) is None


def test_bootstrap_resumes_to_the_same_replicates(tmp_path):
    df = synthetic_runs(n_genes=30, per_gene=5)
    spec = make_spec(df, "B", with_strength=False)
    run_bootstrap(df, spec, "B_without_strength", tmp_path / "a.jsonl", reps=6)
    run_bootstrap(df, spec, "B_without_strength", tmp_path / "b.jsonl", reps=3)
    run_bootstrap(df, spec, "B_without_strength", tmp_path / "b.jsonl", reps=6)
    a = [json.loads(x) for x in (tmp_path / "a.jsonl").read_text().splitlines()]
    b = [json.loads(x) for x in (tmp_path / "b.jsonl").read_text().splitlines()]
    assert a == b and [r["rep"] for r in a] == list(range(6))


def test_balanced_accuracy_hand_example_and_undefined_without_both_outcomes():
    pred, y = np.array([1, 0, 1, 0, 0]), np.array([1, 1, 0, 0, 0])
    assert balanced_accuracy(pred, y) == pytest.approx(0.5 * (0.5 + 2 / 3))
    assert balanced_accuracy(np.array([1, 0]), np.array([1, 1])) is None


def test_ba_bootstrap_counts_undefined_replicates_and_withholds_the_interval():
    f = pd.DataFrame({"gene_ensembl": [f"g{i}" for i in range(20)], "cls": ["aligned"] * 10 + ["blocking"] * 10,
                      "S": [1, 0] * 10, "y": [1.0] * 10 + [1.0, 0.0] * 5})
    out = ba_bootstrap(f, reps=50)
    assert out["undefined"] == 50 and not out["interval_available"] and out["point"]["aligned"] is None
    assert ba_by_class(f)["difference"] is None


def test_support_counts_distinct_loci_in_the_descendant_closure_and_dates_strictly_before_phase2():
    df = pd.DataFrame({"hypothesis_id": ["h1", "h2", "h3"], "gene_ensembl": ["ENSG00000000001"] * 3,
                       "indication_id": ["D1", "D2", "D1"],
                       "earliest_phase2_start": pd.to_datetime(["2015-01-01", "2015-01-01", None])})
    ev = pd.DataFrame({"studyLocusId": ["L1", "L1", "L2", "L3", "L4"], "targetId": ["ENSG00000000001"] * 5,
                       "diseaseId": ["D1c", "D1c", "D1", "D2parent", "D1"], "resourceScore": [0.6, 0.6, 0.7, 0.9, 0.4],
                       "score": [0.6, 0.6, 0.7, 0.9, 0.4],
                       "publicationDate": ["2015-01-01", "2015-01-01", "2014-12-31", "2010-01-01", "2000-01-01"]})
    closures = {"D1": {"D1", "D1c"}, "D2": {"D2"}, "D1c": {"D1c"}, "D2parent": {"D2parent", "D2"}}
    g = support(df, ev, closures).set_index("hypothesis_id")
    assert g.loc["h1", "G"] == 1 and g.loc["h1", "G_loci"] == 2      # L1 counted once; L4 below 0.5
    assert g.loc["h1", "G_pre"] == 1                                 # L2 published the day before
    assert g.loc["h2", "G"] == 0                                     # an ancestor's evidence never counts
    assert pd.isna(g.loc["h3", "G_pre"])                             # no Phase II date


def test_evidence_check_stops_on_a_score_mismatch_or_missing_l2g():
    ev = pd.DataFrame({"resourceScore": [0.5, 0.6], "score": [0.5, 0.6]})
    check_evidence(ev)
    with pytest.raises(FollowupError):
        check_evidence(ev.assign(score=[0.5, 0.61]))
    with pytest.raises(FollowupError):
        check_evidence(ev.assign(resourceScore=[0.5, None]))


def test_survival_coding_applies_each_competing_event_rule():
    s1 = pd.DataFrame({
        "hypothesis_id": list("abcdefg"), "gene_ensembl": ["ENSG00000000001"] * 7,
        "time_to_phase3_days": [100.0, 2000, 2000, 2000, 2000, 2000, None],
        "phase3_event": [1, 0, 0, 0, 0, 0, 0],
        "status_24": ["advanced", "no_observed_advancement", "business_only", "no_observed_advancement",
                      "no_observed_advancement", "active", "no_observed_advancement"],
        "earliest_phase2_start": pd.to_datetime(["2010-01-01"] * 7),
        "last_phase2_end_date": pd.to_datetime(["2011-01-01", "2012-01-01", "2013-01-01", None, "2027-01-01",
                                                "2012-01-01", "2012-01-01"]),
        "cls": ["aligned"] * 7, "S": [1.0] * 7, "platform_somascan": [0.0] * 7, "z_log10_neff": [0.0] * 7,
        "oncology_f": [0.0] * 7})
    s1.loc[5, "last_phase2_end_date"] = pd.Timestamp("2009-01-01")
    s1.loc[5, "status_24"] = "no_observed_advancement"
    t, counts = survival_table(s1)
    assert list(t["status"]) == [1, 2, 2, 0, 0, 0]
    assert list(t["time"]) == [100.0, 730.0, 1096.0, 2000.0, 2000.0, 0.0]
    assert counts["competing_end_missing"] == 1 and counts["competing_end_after_freeze"] == 1
    assert counts["competing_end_before_start"] == 1 and counts["rows"] == 6


def test_g_columns_and_augment_keep_the_focal_and_add_the_interaction():
    from stage_d.designs import build_design
    rng = np.random.default_rng(3)
    n = 80
    frame = pd.DataFrame({"hypothesis_id": [f"h{i}" for i in range(n)], "gene_ensembl": [f"g{i % 12}" for i in range(n)],
                          "indication_id": [f"d{i % 7}" for i in range(n)], "programs": [[f"p{i % 9}"] for i in range(n)],
                          "cls": rng.choice(["aligned", "blocking"], n), "state": "inconclusive",
                          "S": rng.integers(0, 2, n).astype(float), "E": 0.0, "C": 0.0, "y": rng.integers(0, 2, n).astype(float),
                          "platform_somascan": rng.integers(0, 2, n).astype(float), "z_log10_neff": rng.normal(0, 1, n),
                          "oncology_f": rng.integers(0, 2, n).astype(float)})
    d = build_design("t", "h1", frame)
    gid = pd.Series(rng.integers(0, 2, n).astype(float), index=frame["hypothesis_id"])
    cols = g_columns(d, gid, True)
    d2 = augment(d, "t_G", cols)
    assert d2.columns[-2:] == ("G", "GxA") and d2.X.shape[1] == d.X.shape[1] + 2
    assert np.array_equal(d2.focal[: len(d.focal)], d.focal) and d2.focal[-2:].sum() == 0
    assert np.array_equal(d2.X[:, -1], d2.X[:, -2] * d.X[:, d.stratum_col])


def test_require_run_checks_the_plan_hash_and_the_token(tmp_path):
    plan, log = tmp_path / "PREREG.md", tmp_path / "PREREG.log"
    plan.write_text("not the plan")
    log.write_text("2026-10-13  RUN_START followup token=0123456789abcdef  no results seen\n")
    with pytest.raises(FollowupError, match=FROZEN_PLAN_SHA256[:8]):
        require_run(plan, log, "0123456789abcdef")


def test_wilson_matches_a_known_value():
    lo, hi = wilson(5, 20)
    assert (lo, hi) == pytest.approx((0.1119, 0.4687), abs=1e-4)
