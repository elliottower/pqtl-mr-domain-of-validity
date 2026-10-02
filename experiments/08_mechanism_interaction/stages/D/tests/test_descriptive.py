import numpy as np
import pandas as pd
import pytest
from helpers import small_frame

from stage_d.descriptive import (cells_2x2, metrics_2x2, outcome_category, stratum_2x2_table, table11_per_indication,
                                 table12_platform)


def test_metrics_2x2_hand_values():
    m = metrics_2x2(np.array([8.0, 2.0, 30.0, 60.0]))
    assert m["ppv"] == pytest.approx(0.8)
    assert m["npv"] == pytest.approx(60 / 90)
    assert m["rs"] == pytest.approx(0.8 / (30 / 90))
    assert m["or"] == pytest.approx(8 * 60 / (2 * 30))
    assert m["sensitivity"] == pytest.approx(8 / 38)
    assert m["specificity"] == pytest.approx(60 / 62)
    assert m["balanced_accuracy"] == pytest.approx((8 / 38 + 60 / 62) / 2)
    assert np.isnan(metrics_2x2(np.array([0.0, 0.0, 3.0, 4.0]))["rs"])


def test_stratum_table_point_estimates_and_component_bootstrap():
    f = small_frame(np.random.default_rng(0), n=600, n_genes=60)
    rows = {r["stratum"]: r for r in stratum_2x2_table(f, n_boot=2000)}
    for cls in ("aligned", "blocking"):
        g = f.loc[f["cls"] == cls]
        a, b, c, d = cells_2x2(g["S"].to_numpy(), g["y"].to_numpy())
        r = rows[cls]
        assert (r["a_supportive_advanced"], r["b_supportive_not"], r["c_not_supportive_advanced"],
                r["d_not_supportive_not"]) == (a, b, c, d)
        assert r["ppv_q05"] < r["ppv"] < r["ppv_q95"]
        assert r["balanced_accuracy_q05"] < r["balanced_accuracy"] < r["balanced_accuracy_q95"]


def test_component_bootstrap_resamples_whole_components():
    # one component holds every row: resampling components reproduces the data exactly, so the
    # interval collapses to the point estimate
    f = small_frame(np.random.default_rng(1), n=200)
    f["programs"] = [["shared"]] * len(f)
    r = stratum_2x2_table(f, n_boot=200)[0]
    assert r["ppv_q05"] == pytest.approx(r["ppv"]) and r["ppv_q95"] == pytest.approx(r["ppv"])


def test_table11_pools_small_indications():
    f = small_frame(np.random.default_rng(2), n=120, n_ind=40)
    rows = table11_per_indication(f)
    sizes = f.groupby("indication_id").size()
    pooled = [r for r in rows if r["indication"].startswith("pooled")]
    assert pooled[0]["n"] == int(sizes[sizes < 5].sum())
    assert pooled[0]["indications_in_row"] == int((sizes < 5).sum())
    assert sum(r["n"] for r in rows) == len(f)


def test_outcome_decomposition_categories():
    f = pd.DataFrame({"status_24": ["advanced", "no_observed_advancement", "no_observed_advancement", "active",
                                    "business_only", "undated"],
                      "stop_code": [None, "efficacy", "completed_no_successor", None, "business", None]})
    assert list(outcome_category(f)) == ["advanced", "efficacy", "completed_no_successor", "active", "business", "undated"]


def test_cross_source_agreement_hand_counts():
    f = pd.DataFrame({
        "cls": ["aligned", "aligned", "blocking", "blocking", "other", "blocking"],
        "gene_ensembl": ["g1", "g1", "g2", "g3", "g4", "g5"],
        "platform_concordant": ["concordant"] * 6,
        "evidence_state_ukbppp": ["supportive", "inconclusive", "supportive", "contradictory", None, "inconclusive"],
        "evidence_state_decode": ["supportive", "supportive", "supportive", "contradictory", "supportive", None],
    })
    out = table12_platform(f)
    s = out["cross_source_state_agreement"]
    assert (s["n_both_sources"], s["genes"], s["agree"], s["supportive_either"], s["supportive_both"]) == (4, 3, 3, 3, 2)
    assert s["agreement_fraction"] == pytest.approx(0.75)
    assert {r["class"]: (r["n"], r["agree"]) for r in out["cross_source_state_by_class"]} == {
        "aligned": (2, 1), "blocking": (2, 2), "other": (0, 0)}
    xt = {r["ukbppp_state"]: r for r in out["cross_source_state_crosstab"]}
    assert xt["inconclusive"]["decode_supportive"] == 1 and xt["supportive"]["decode_supportive"] == 2


def test_cross_source_agreement_not_available_without_both_states():
    f = pd.DataFrame({"cls": ["aligned"], "gene_ensembl": ["g"], "platform_concordant": ["untested"],
                      "evidence_state_ukbppp": [None], "evidence_state_decode": [None]})
    out = table12_platform(f)
    assert "not_available" in out["cross_source_state_agreement"]
    assert out["cross_source_state_crosstab"] == [] and out["cross_source_state_by_class"] == []
