import inspect
import itertools

import pytest

from stage_d.decisions import ModelEvidence, against_h1, decide_h1, decide_h4, h1_wording
from stage_d.fitspec import ATTEMPTS, fit_outcome, next_attempt
from stage_d.gates import StratumGate
from stage_d.pipeline import model_evidence


def ev(pr=0.99, structural=True, reliable=True, outcome="full", dominated=False, opp=None):
    return ModelEvidence(structurally_estimable=structural, reliable=reliable, fit_outcome=outcome,
                         prior_dominated=None if outcome in ("failed", "not_fitted") and dominated is None else dominated,
                         pr_predicted=pr, pr_opposite=1 - pr if opp is None else opp)


def test_h1_holds_at_threshold_and_not_below():
    assert decide_h1(ev(0.95)).status == "holds"
    assert decide_h1(ev(0.9499)).status == "does_not_hold"
    assert decide_h1(ev(0.99, outcome="full_rerun")).status == "holds"


@pytest.mark.parametrize("kw,reason", [
    ({"reliable": False}, "reliability"),
    ({"outcome": "reduced"}, "sampler"),
    ({"outcome": "reduced_rerun"}, "sampler"),
    ({"outcome": "failed"}, "sampler"),
    ({"dominated": True}, "prior-dominated"),
])
def test_h1_no_confirmatory_decision_branches(kw, reason):
    d = decide_h1(ev(0.999, **kw))
    assert d.status == "no_confirmatory_decision" and not d.permitted
    assert any(reason in r for r in d.reasons)


def test_h1_structural_non_estimability():
    d = decide_h1(ModelEvidence(structurally_estimable=False, reliable=False, fit_outcome="not_fitted"))
    assert d.status == "not_estimable" and not d.permitted


H1_CASES = [ev(0.999), ev(0.5), ev(0.999, reliable=False), ev(0.999, outcome="reduced"), ev(0.999, dominated=True),
            ModelEvidence(structurally_estimable=False, reliable=False, fit_outcome="not_fitted")]
H4_CASES = [ev(0.999), ev(0.5), ev(0.999, reliable=False), ev(0.999, outcome="reduced"), ev(0.999, dominated=True),
            ModelEvidence(structurally_estimable=False, reliable=False, fit_outcome="not_fitted")]


@pytest.mark.parametrize("m1,m4", list(itertools.product(H1_CASES, H4_CASES)))
def test_fixed_sequence_h4_confirmatory_only_after_h1_holds(m1, m4):
    h1 = decide_h1(m1)
    h4 = decide_h4(h1, m4)
    if h1.status != "holds":
        assert h4.status == "secondary_estimate" and not h4.confirmatory and not h4.permitted
    else:
        assert h4.confirmatory
        assert h4.status == decide_h1(m4).status  # same gate/threshold logic once H1 holds
    assert (h4.status == "holds") <= (h1.status == "holds")


def test_h4_holds_only_with_its_own_gates():
    h1 = decide_h1(ev(0.99))
    assert decide_h4(h1, ev(0.96)).status == "holds"
    assert decide_h4(h1, ev(0.94)).status == "does_not_hold"
    assert decide_h4(h1, ev(0.99, reliable=False)).status == "no_confirmatory_decision"


def test_against_h1_row():
    assert against_h1(ev(0.02, opp=0.98))["criterion_met"]
    assert not against_h1(ev(0.10, opp=0.90))["criterion_met"]
    assert against_h1(ev(0.02, opp=0.98, reliable=False))["decision_permitted"] is False


def test_h2_changes_wording_never_decision():
    h1 = decide_h1(ev(0.97))
    assert list(inspect.signature(decide_h1).parameters) == ["m"]
    assert h1_wording(h1, 0.5) == "difference_between_strata_without_overall_association"
    assert h1_wording(decide_h1(ev(0.6)), 0.5) == "uninformative_about_moderation"
    assert h1_wording(h1, 0.99) == "moderation_of_positive_association"


def _final(attempts):
    return {"outcome": fit_outcome(attempts), "attempts": attempts}


def _att(name, passed, pr=0.99, dominated=False):
    return {"attempt": name, "passed": passed,
            "summary": {"prior_dominated": dominated, "focal": {"pr_predicted": pr, "pr_opposite": 1 - pr}}}


GATE_OK = StratumGate(model="h1", strata=[], structurally_estimable=True, reliable=True, reasons=[])


def test_decision_uses_rerun_when_primary_fails():
    final = _final([_att("full", False, 0.5), _att("full_rerun", True, 0.97)])
    assert decide_h1(model_evidence(GATE_OK, final)).status == "holds"


def test_reduced_model_never_confirmatory_even_when_it_passes():
    final = _final([_att("full", False), _att("full_rerun", False), _att("reduced", True, 0.999)])
    d = decide_h1(model_evidence(GATE_OK, final))
    assert final["outcome"] == "reduced"
    assert d.status == "no_confirmatory_decision"
    assert decide_h4(d, ev(0.999)).status == "secondary_estimate"


def test_prior_dominance_from_fit_summary_blocks_decision():
    d = decide_h1(model_evidence(GATE_OK, _final([_att("full", True, 0.999, dominated=True)])))
    assert d.status == "no_confirmatory_decision"


def test_attempt_sequence():
    assert next_attempt([]) == ATTEMPTS[0]
    assert next_attempt([{"passed": True}]) is None
    assert next_attempt([{"passed": False}]) == ("full_rerun", "rerun", True)
    assert next_attempt([{"passed": False}] * 2) == ("reduced", "primary", False)
    assert next_attempt([{"passed": False}] * 3) == ("reduced_rerun", "rerun", False)
    assert next_attempt([{"passed": False}] * 4) is None
    assert fit_outcome([{"attempt": "full", "passed": False}] * 4) == "failed"
    assert fit_outcome([]) == "not_fitted"
