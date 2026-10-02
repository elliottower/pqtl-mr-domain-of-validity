"""Confirmatory decisions (PREREG §Inference criteria).

| hypothesis | holds when |
| H1 | Pr(beta_SxA > 0) >= 0.95 in the primary model on the held-out set, and a confirmatory
|    | decision is permitted |
| H4 | H1 holds, Pr(beta_SxC < 0) >= 0.95 on the held-out set, and a confirmatory decision is
|    | permitted |
| H2 | no confirmatory criterion |
| against H1 | Pr(beta_SxA < 0) >= 0.95, reported |

A decision is permitted only when the model is structurally estimable, passes the reliability
gate, the full model (random gene slope kept) passes its sampler diagnostics at the primary
settings or on the rerun, and the focal coefficient is not prior-dominated. The reduced model is
descriptive only. H4 is confirmatory only after H1 holds. Nothing else (H2, H3, sets, secondary
models) enters these functions.
"""
from typing import Literal

from pydantic import BaseModel

from stage_d.constants import PROB_THRESHOLD

FitOutcome = Literal["full", "full_rerun", "reduced", "reduced_rerun", "failed", "not_fitted"]
Status = Literal["holds", "does_not_hold", "no_confirmatory_decision", "not_estimable", "secondary_estimate"]


class ModelEvidence(BaseModel):
    structurally_estimable: bool
    reliable: bool
    fit_outcome: FitOutcome
    prior_dominated: bool | None = None
    pr_predicted: float | None = None
    pr_opposite: float | None = None


class Decision(BaseModel):
    hypothesis: str
    status: Status
    confirmatory: bool
    permitted: bool
    reasons: list[str]
    pr_predicted: float | None = None


def permitted(m: ModelEvidence) -> tuple[bool, list[str]]:
    reasons = []
    if not m.structurally_estimable:
        reasons.append("structurally non-estimable")
    if not m.reliable:
        reasons.append("reliability gate failed (low-information estimate)")
    if m.fit_outcome not in ("full", "full_rerun"):
        reasons.append(f"full model did not pass sampler diagnostics (fit outcome: {m.fit_outcome})")
    if m.prior_dominated is None and m.fit_outcome in ("full", "full_rerun"):
        raise ValueError("prior dominance not evaluated for a passing full fit")
    if m.prior_dominated:
        reasons.append("focal coefficient prior-dominated")
    return not reasons, reasons


def decide_h1(m: ModelEvidence) -> Decision:
    if not m.structurally_estimable:
        return Decision(hypothesis="H1", status="not_estimable", confirmatory=True, permitted=False,
                        reasons=["structurally non-estimable"])
    ok, reasons = permitted(m)
    if not ok:
        return Decision(hypothesis="H1", status="no_confirmatory_decision", confirmatory=True, permitted=False,
                        reasons=reasons, pr_predicted=m.pr_predicted)
    holds = m.pr_predicted >= PROB_THRESHOLD
    return Decision(hypothesis="H1", status="holds" if holds else "does_not_hold", confirmatory=True, permitted=True,
                    reasons=[f"Pr(beta_SxA > 0) = {m.pr_predicted:.4f}"], pr_predicted=m.pr_predicted)


def decide_h4(h1: Decision, m: ModelEvidence) -> Decision:
    if h1.status != "holds":
        return Decision(hypothesis="H4", status="secondary_estimate", confirmatory=False, permitted=False,
                        reasons=[f"fixed sequence: H1 status is {h1.status}"], pr_predicted=m.pr_predicted)
    if not m.structurally_estimable:
        return Decision(hypothesis="H4", status="not_estimable", confirmatory=True, permitted=False,
                        reasons=["structurally non-estimable"])
    ok, reasons = permitted(m)
    if not ok:
        return Decision(hypothesis="H4", status="no_confirmatory_decision", confirmatory=True, permitted=False,
                        reasons=reasons, pr_predicted=m.pr_predicted)
    holds = m.pr_predicted >= PROB_THRESHOLD
    return Decision(hypothesis="H4", status="holds" if holds else "does_not_hold", confirmatory=True, permitted=True,
                    reasons=[f"Pr(beta_SxC < 0) = {m.pr_predicted:.4f}"], pr_predicted=m.pr_predicted)


def against_h1(m: ModelEvidence) -> dict:
    """Reported row: Pr(beta_SxA < 0) >= 0.95, with whether a confirmatory decision was permitted."""
    ok, _ = permitted(m) if m.structurally_estimable else (False, [])
    met = m.pr_opposite is not None and m.pr_opposite >= PROB_THRESHOLD
    return {"pr_beta_SxA_lt_0": m.pr_opposite, "criterion_met": met, "decision_permitted": ok}


def h1_wording(h1: Decision, h2_pr_gt0: float | None) -> str:
    """How H1 is described given H2 (PREREG §Inference criteria); never changes the decision."""
    if h1.pr_predicted is None or h2_pr_gt0 is None:
        return "not_applicable"
    if h2_pr_gt0 >= PROB_THRESHOLD:
        return "moderation_of_positive_association" if h1.pr_predicted >= PROB_THRESHOLD else "no_moderation_shown"
    if h1.pr_predicted >= PROB_THRESHOLD:
        return "difference_between_strata_without_overall_association"
    return "uninformative_about_moderation"
