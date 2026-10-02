"""Descriptive diagnostics of stage C. None is a registered analysis and none changes a label.

Two of them apply the registered rules (rules.code_hypothesis, unchanged) to an alternative reading
of a hypothesis's linked trials and list the quantities whose value differs:

partial_date_boundary.csv  a ClinicalTrials.gov date given as year-month or year takes the last day
                           of that month or year instead of the first (ctgov.parse_study
                           `last_day=True`); quantities: the 12/24/36/48-month status.
phase23_diagnostic.csv     a registration listing both PHASE2 and PHASE3 counts as a Phase III
                           registration and not as a Phase II registration (`phase23_as_phase3`);
                           quantities: the four statuses and the time-to-Phase-III fields. The
                           Open Targets maxClinicalStage values are as registered, so `advanced`
                           itself cannot differ.

The third, linkage_audit.csv, is ot_phase.linkage_audit: counts at each trial-linking step.
"""
from datetime import date

from models import MATURATION_WINDOWS_MONTHS, DifferenceRow, Trial
from rules import HypothesisOutcome, code_hypothesis

STATUS_QUANTITIES = tuple(f"status_{w}" for w in MATURATION_WINDOWS_MONTHS)
PHASE3_QUANTITIES = ("phase2_start_date", "phase3_start_date", "time_to_phase3_days", "phase3_event")


def _text(v) -> str:
    return "" if v is None else v.isoformat() if isinstance(v, date) else str(v)


def outcome_values(o: HypothesisOutcome) -> dict[str, str]:
    """The compared quantities of one coded hypothesis, as the text outcomes.csv would hold."""
    return {**{f"status_{w}": o.statuses[w].value for w in MATURATION_WINDOWS_MONTHS},
            "phase2_start_date": _text(o.phase2_start), "phase3_start_date": _text(o.phase3_start),
            "time_to_phase3_days": _text(o.time_to_phase3_days), "phase3_event": _text(o.phase3_event)}


def differences(hypothesis_id: str, registered: HypothesisOutcome, alternative: HypothesisOutcome,
                quantities: tuple[str, ...]) -> list[DifferenceRow]:
    """One row per quantity of `quantities` whose value differs between the two codings."""
    reg, alt = outcome_values(registered), outcome_values(alternative)
    return [DifferenceRow(hypothesis_id=hypothesis_id, quantity=q, registered=reg[q], alternative=alt[q])
            for q in quantities if reg[q] != alt[q]]


def partial_date_boundary(hypothesis_id: str, max_stages: list[str], trials: list[Trial],
                          trials_last_day: list[Trial], freeze: date) -> list[DifferenceRow]:
    """`trials` are the linked trials as registered (partial dates on the first day),
    `trials_last_day` the same registrations with partial dates on the last day."""
    return differences(hypothesis_id, code_hypothesis(max_stages, trials, freeze),
                       code_hypothesis(max_stages, trials_last_day, freeze), STATUS_QUANTITIES)


def phase23_as_phase3(t: Trial) -> Trial:
    """A registration listing PHASE2 and PHASE3 read as a Phase III registration only."""
    if "PHASE2" in t.phases and "PHASE3" in t.phases:
        return t.model_copy(update={"phases": ("PHASE3",)})
    return t


def phase23_diagnostic(hypothesis_id: str, max_stages: list[str], trials: list[Trial], freeze: date) -> list[DifferenceRow]:
    return differences(hypothesis_id, code_hypothesis(max_stages, trials, freeze),
                       code_hypothesis(max_stages, [phase23_as_phase3(t) for t in trials], freeze),
                       STATUS_QUANTITIES + PHASE3_QUANTITIES)
