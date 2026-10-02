import pandas as pd
import pytest

from stage_d.gates import Diagnostics, h1_gate, h3_gate, h4_gate, prior_dominated, sampler_passed


def _frame(spec):
    """spec: list of (cls, C, n_supportive_genes, supportive_per_gene, n_nonsupportive)."""
    rows = []
    for cls, c, n_genes, per_gene, n_non in spec:
        for g in range(n_genes):
            for j in range(per_gene):
                rows.append({"cls": cls, "C": c, "S": 1.0, "state": "supportive", "gene_ensembl": f"{cls}{c}g{g}",
                             "indication_id": f"i{j}", "programs": [f"{cls}{c}p{g}"]})
        for j in range(n_non):
            rows.append({"cls": cls, "C": c, "S": 0.0, "state": "inconclusive", "gene_ensembl": f"{cls}{c}n{j}",
                         "indication_id": "i0", "programs": [f"{cls}{c}q{j}"]})
    return pd.DataFrame(rows)


def test_h1_gate_passes_at_the_registered_minimum():
    g = h1_gate(_frame([("aligned", 0, 10, 2, 5), ("blocking", 0, 10, 2, 5)]))
    assert g.structurally_estimable and g.reliable and g.reasons == []
    counts = {c.stratum: c for c in g.strata}
    assert counts["aligned"].n_supportive == 20 and counts["aligned"].n_supportive_genes == 10


@pytest.mark.parametrize("spec", [
    [("aligned", 0, 9, 3, 5), ("blocking", 0, 10, 2, 5)],    # 9 supportive genes
    [("aligned", 0, 10, 2, 5), ("blocking", 0, 19, 1, 5)],   # 19 supportive hypotheses
])
def test_h1_gate_reliability_fails_below_minimum_but_stays_estimable(spec):
    g = h1_gate(_frame(spec))
    assert g.structurally_estimable and not g.reliable


@pytest.mark.parametrize("spec", [
    [("aligned", 0, 0, 0, 30), ("blocking", 0, 20, 2, 5)],   # no supportive in aligned
    [("aligned", 0, 20, 2, 0), ("blocking", 0, 20, 2, 5)],   # no non-supportive in aligned
])
def test_h1_gate_structural_non_estimability(spec):
    g = h1_gate(_frame(spec))
    assert not g.structurally_estimable and not g.reliable


def test_h1_gate_ignores_other_class():
    g = h1_gate(_frame([("aligned", 0, 10, 2, 5), ("blocking", 0, 10, 2, 5), ("other", 0, 0, 0, 50)]))
    assert g.reliable
    assert {c.stratum for c in g.strata} == {"aligned", "blocking"}


def test_h4_gate_uses_indication_strata_over_all_classes():
    f = _frame([("aligned", 1, 5, 2, 3), ("other", 1, 5, 2, 3), ("blocking", 0, 10, 2, 5)])
    g = h4_gate(f)
    counts = {c.stratum: c for c in g.strata}
    assert counts["neuro_psych"].n_supportive_genes == 10 and counts["neuro_psych"].n_supportive == 20
    assert g.reliable
    assert not h4_gate(_frame([("other", 1, 9, 3, 3), ("blocking", 0, 10, 2, 5)])).reliable


def test_h3_gate_boundary():
    def f(n_sup, n_con):
        return pd.DataFrame({"state": ["supportive"] * n_sup + ["contradictory"] * n_con + ["inconclusive"] * 5})
    assert h3_gate(f(20, 20)).passes
    assert not h3_gate(f(20, 19)).passes
    assert not h3_gate(f(19, 40)).passes


@pytest.mark.parametrize("rhat,ess,div,ok", [
    (1.009, 401, 0, True), (1.01, 401, 0, False), (1.0, 400, 0, False), (1.0, 5000, 1, False),
    (float("nan"), 5000, 0, False),
])
def test_sampler_diagnostics_rule(rhat, ess, div, ok):
    assert sampler_passed(Diagnostics(rhat_max=rhat, ess_bulk_min=ess, divergences=div)) is ok


def test_prior_dominance_threshold_is_point_eight_of_prior_sd():
    assert not prior_dominated(1.2, 1.5)
    assert prior_dominated(1.2000001, 1.5)
    assert not prior_dominated(0.3, 1.5)
