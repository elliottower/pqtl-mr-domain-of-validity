import pytest
import pqtl_validity as pv
from pqtl_validity.taxonomy import Mechanism


def test_classify_pair_enzyme_inhibitor_predicts_uninformative():
    result = pv.classify_pair("ACE", mr_p=0.72)
    assert result.prediction == "UNINFORMATIVE"
    assert result.mr_informative_expected is False
    assert result.mr_significant is False


def test_classify_pair_soluble_ligand_predicts_informative():
    result = pv.classify_pair("TNF", mr_p=0.003)
    assert result.prediction == "INFORMATIVE"
    assert result.mr_informative_expected is True
    assert result.mr_significant is True


def test_classify_pair_unknown_gene():
    result = pv.classify_pair("BRCA1")
    assert result.prediction == "UNKNOWN"
    assert result.mechanism is None


def test_classify_pair_no_p_value():
    result = pv.classify_pair("ACE")
    assert result.mr_p is None
    assert result.mr_significant is False
    assert result.prediction == "UNINFORMATIVE"


def test_score_pairs_perfect_separation():
    genes = ["TNF", "VEGFA", "ACE", "EGFR"]
    outcomes = ["SUCCESS", "SUCCESS", "SUCCESS", "SUCCESS"]
    mr_ps = [0.001, 0.002, 0.80, 0.60]
    result = pv.score_pairs(genes, outcomes, mr_ps=mr_ps, seed=0)
    assert result.n == 4
    assert result.sensitivity > 0
    assert result.specificity >= 0


def test_score_pairs_balanced_accuracy_bounded():
    genes = ["TNF", "ACE", "EGFR", "VEGFA", "DPP4", "IL6R"]
    outcomes = ["SUCCESS", "SUCCESS", "SUCCESS", "SUCCESS", "SUCCESS", "SUCCESS"]
    mr_ps = [0.001, 0.80, 0.60, 0.002, 0.70, 0.003]
    result = pv.score_pairs(genes, outcomes, mr_ps=mr_ps, seed=42)
    assert 0.0 <= result.balanced_accuracy <= 1.0
    assert result.ba_ci_low <= result.balanced_accuracy
    assert result.balanced_accuracy <= result.ba_ci_high


def test_score_pairs_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        pv.score_pairs(["TNF", "ACE"], ["SUCCESS"])


# ----------------------------------------------------------------------
# Reproduction of published values from the shipped data tables.
# These are the end-to-end guard: if the public API stops reproducing the
# paper, they fail.
# ----------------------------------------------------------------------

import csv
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def _load(rel_path):
    with open(REPO / rel_path, newline="") as fh:
        return list(csv.DictReader(fh))


def _score(rows, **kwargs):
    return pv.score_pairs(
        genes=[r["gene"] for r in rows],
        outcomes=[r["outcome"] for r in rows],
        mr_ps=[float(r["mr_p"]) for r in rows],
        n_bootstrap=200,
        seed=0,
        **kwargs,
    )


def test_public_api_reproduces_published_primary_confusion_matrix():
    rows = _load("results/v5/classification_v5.csv")
    result = _score(rows)
    # results/v5/evaluation_v5.json, primary_metrics
    assert (result.n, result.tp, result.tn, result.fp, result.fn) == (157, 11, 101, 8, 37)
    assert result.balanced_accuracy == pytest.approx(0.5779, abs=5e-4)


def test_public_api_reproduces_published_mechanism_strata():
    rows = _load("data/classification_v34.csv")
    activity = _score([r for r in rows if r["mechanism_class"] == "activity_blocking"])
    abundance = _score([r for r in rows if r["mechanism_class"] == "abundance_modulating"])

    # Manuscript, Binary stratification section.
    assert activity.n == 76
    assert activity.balanced_accuracy == pytest.approx(0.511, abs=5e-4)
    assert activity.sensitivity == pytest.approx(0.095, abs=5e-4)

    assert abundance.n == 59
    assert abundance.balanced_accuracy == pytest.approx(0.590, abs=5e-4)


def test_pairs_without_an_adjudicated_outcome_are_excluded_and_reported():
    rows = _load("results/v5/classification_v5.csv")
    result = _score(rows)
    n_unadjudicated = sum(
        1 for r in rows if r["outcome"] not in ("SUCCESS", "FAILURE")
    )
    assert n_unadjudicated > 0, "fixture should contain EXCLUDED/PENDING rows"
    assert len(result.dropped_no_outcome) == n_unadjudicated
    assert result.n == len(rows) - n_unadjudicated


def test_mechanism_gated_rule_refuses_to_silently_drop_unknown_genes():
    # The taxonomy lookup covers fewer genes than the curated tables. Under the
    # gated rule an unknown gene must raise rather than shrink the denominator.
    genes = ["TNF", "NOT_A_REAL_GENE_XYZ"]
    with pytest.raises(ValueError, match="NOT_A_REAL_GENE_XYZ"):
        pv.score_pairs(genes, ["SUCCESS", "FAILURE"], mr_ps=[0.01, 0.02],
                       rule="mechanism_gated", seed=0)


def test_mechanism_gated_rule_accepts_curated_assignments():
    # Passing the curated column lets genes absent from the taxonomy score.
    result = pv.score_pairs(
        genes=["TNF", "NOT_A_REAL_GENE_XYZ"],
        outcomes=["SUCCESS", "FAILURE"],
        mr_ps=[0.01, 0.02],
        mechanisms=["abundance_modulating", "activity_blocking"],
        rule="mechanism_gated",
        n_bootstrap=200,
        seed=0,
    )
    assert result.n == 2
    assert result.dropped_unclassified == ()


def test_skip_mode_records_what_it_dropped():
    result = pv.score_pairs(
        genes=["TNF", "VEGFA", "NOT_A_REAL_GENE_XYZ"],
        outcomes=["SUCCESS", "FAILURE", "SUCCESS"],
        mr_ps=[0.01, 0.60, 0.02],
        rule="mechanism_gated",
        on_unclassified="skip",
        n_bootstrap=200,
        seed=0,
    )
    assert "NOT_A_REAL_GENE_XYZ" in result.dropped_unclassified
    assert result.n_dropped >= 1
    assert result.n == 2


def test_default_rule_ignores_mechanism_entirely():
    # Under the paper's screen, prediction depends only on the p-value, so
    # relabelling every mechanism must not move the result.
    genes = ["TNF", "ACE", "EGFR", "VEGFA"]
    outcomes = ["SUCCESS", "FAILURE", "SUCCESS", "FAILURE"]
    mr_ps = [0.01, 0.02, 0.60, 0.70]
    a = pv.score_pairs(genes, outcomes, mr_ps=mr_ps, n_bootstrap=200, seed=0)
    b = pv.score_pairs(genes, outcomes, mr_ps=mr_ps, n_bootstrap=200, seed=0,
                       mechanisms=["activity_blocking"] * 4)
    assert a.balanced_accuracy == b.balanced_accuracy
    assert (a.tp, a.tn, a.fp, a.fn) == (b.tp, b.tn, b.fp, b.fn)
