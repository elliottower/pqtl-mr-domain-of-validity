"""Outcome rules of PREREG §Measured variables, as pure functions over synthetic-testable inputs.

No I/O. Each function names the plan sentence it implements. Implementation choices the plan
does not fix are marked IMPLEMENTATION CHOICE and listed in README.md §Choices needing an
amendment.
"""
import calendar
from dataclasses import dataclass
from datetime import date

from models import (ACTIVE_STATUSES, ADVANCED_STAGES, APPROVAL_STAGE, ENDED_STATUSES,
                    MATURATION_WINDOWS_MONTHS, PRIMARY_WINDOW_MONTHS, STOP_KEYWORDS,
                    STOP_PRECEDENCE, STOPPED_STATUSES, StageCError, Status, StopCode, Trial)


# ---- dates (partial ClinicalTrials.gov dates are parsed in ctgov.parse_partial_date) ---------

def months_before(d: date, months: int) -> date:
    """Calendar subtraction; the day is clamped to the target month's length."""
    total = d.year * 12 + (d.month - 1) - months
    y, m = divmod(total, 12)
    m += 1
    return date(y, m, min(d.day, calendar.monthrange(y, m)[1]))


def is_mature(last_end: date, freeze: date, window_months: int) -> bool:
    """PREREG: the last Phase II completion or termination date is at least `window` months
    before the outcome-freeze date."""
    return last_end <= months_before(freeze, window_months)


# ---- trial state at the freeze -------------------------------------------------------------

def registered_by_freeze(t: Trial, freeze: date) -> bool:
    """PREREG §Starting and stopping rules: trial records dated after the freeze do not change
    any label. IMPLEMENTATION CHOICE: a trial first submitted after the freeze is ignored."""
    return t.first_submit_date is None or t.first_submit_date <= freeze


def post_freeze_counts(trials: list[Trial], freeze: date) -> tuple[int, int]:
    """Reported diagnostic, changes no label: (linked trials first submitted after the freeze,
    linked trials whose last update was submitted after the freeze). The second count includes
    the first; their difference is the trials that inform labels but whose archived record was
    last updated after the freeze."""
    first = sum(1 for t in trials if t.first_submit_date is not None and t.first_submit_date > freeze)
    update = sum(1 for t in trials if t.last_update_submit_date is not None and t.last_update_submit_date > freeze)
    return first, update


def status_at_freeze(t: Trial, freeze: date) -> str:
    """The retrieved status, except that a trial whose completion or termination date falls
    after the freeze was still running at the freeze. IMPLEMENTATION CHOICE: such a trial is
    read as ACTIVE_NOT_RECRUITING, so neither its end nor its why_stopped counts."""
    if (t.overall_status in ENDED_STATUSES | STOPPED_STATUSES and t.completion_date is not None
            and t.completion_date > freeze):
        return "ACTIVE_NOT_RECRUITING"
    return t.overall_status


def is_phase2(t: Trial) -> bool:
    """IMPLEMENTATION CHOICE: a Phase II trial is one whose phases include PHASE2 (so Phase
    I/II and Phase II/III registrations count)."""
    return "PHASE2" in t.phases


def is_phase3(t: Trial) -> bool:
    """IMPLEMENTATION CHOICE: a Phase III trial lists PHASE3 and not PHASE2, matching Open
    Targets, where PHASE_2_3 is not advanced."""
    return "PHASE3" in t.phases and "PHASE2" not in t.phases


# ---- why_stopped ---------------------------------------------------------------------------

def classify_why_stopped(text: str) -> StopCode:
    """PREREG keyword table, case-insensitive substring match. IMPLEMENTATION CHOICE: a text
    matching keywords of more than one code takes the precedence order, as stops do."""
    low = (text or "").lower()
    hits = {code for code, kws in STOP_KEYWORDS.items() if any(k in low for k in kws)}
    if not hits:
        return StopCode.unclassified
    return combine_stop_codes([StopCode(c) for c in hits])


def combine_stop_codes(codes: list[StopCode]) -> StopCode:
    """PREREG: where stops carry different codes, efficacy > safety > unclassified > business."""
    if not codes:
        raise StageCError("combine_stop_codes needs at least one code")
    for c in STOP_PRECEDENCE:
        if StopCode(c) in codes:
            return StopCode(c)
    raise StageCError(f"codes outside the keyword table: {codes}")


# ---- hypothesis coding ---------------------------------------------------------------------

@dataclass(frozen=True)
class HypothesisOutcome:
    statuses: dict[int, Status]
    advanced: dict[int, int | None]
    stop_code: StopCode
    why_stopped_texts: str
    approved: bool
    phase3_success: int | None
    phase2_start: date | None
    phase3_start: date | None
    time_to_phase3_days: int | None
    phase3_event: int | None
    last_phase2_end: date | None


def advanced_value(status: Status) -> int | None:
    if status == Status.advanced:
        return 1
    if status == Status.no_observed_advancement:
        return 0
    return None


def stop_decomposition(trials: list[Trial], freeze: date) -> tuple[StopCode, str, bool]:
    """PREREG §Outcome decomposition. Returns (code, `||`-joined why_stopped texts,
    every-stop-is-business). Stops are the hypothesis's terminated, suspended or withdrawn
    trials of any phase. With no stop, a hypothesis with a completed trial carries
    completed_no_successor; a completed trial's missing why_stopped is never read."""
    stops = sorted((t for t in trials if status_at_freeze(t, freeze) in STOPPED_STATUSES),
                   key=lambda t: t.nct_id)
    if stops:
        codes = [classify_why_stopped(t.why_stopped) for t in stops]
        texts = "||".join(t.why_stopped for t in stops)
        return combine_stop_codes(codes), texts, all(c == StopCode.business for c in codes)
    if any(status_at_freeze(t, freeze) == "COMPLETED" for t in trials):
        return StopCode.completed_no_successor, "", False
    return StopCode.none, "", False


def last_end(trials: list[Trial], freeze: date, phase_pred) -> date | None:
    ends = [t.completion_date for t in trials
            if phase_pred(t) and status_at_freeze(t, freeze) in ENDED_STATUSES
            and t.completion_date is not None]
    return max(ends) if ends else None


def primary_status(advanced: bool, trials: list[Trial], all_business: bool, freeze: date,
                   window_months: int) -> Status:
    """PREREG primary outcome, with the business-only exclusion and the pre-stage undated
    rule. IMPLEMENTATION CHOICE on labels only (all three are excluded): an ongoing trial is
    checked first, then business-only, then undated."""
    if advanced:
        return Status.advanced
    if any(status_at_freeze(t, freeze) in ACTIVE_STATUSES for t in trials):
        return Status.active
    if all_business:
        return Status.business_only
    end = last_end(trials, freeze, is_phase2)
    if end is None:
        return Status.undated
    return Status.no_observed_advancement if is_mature(end, freeze, window_months) else Status.active


def phase3_success(advanced: bool, approved: bool, trials: list[Trial], freeze: date) -> int | None:
    """PREREG secondary: among hypotheses that reached Phase III, 1 if any drug reached
    APPROVAL; 0 if none did, no Phase III trial is active and the last Phase III completion or
    termination is at least 24 months before the freeze; otherwise pending (empty)."""
    if not advanced:
        return None
    if approved:
        return 1
    p3 = [t for t in trials if is_phase3(t)]
    if any(status_at_freeze(t, freeze) in ACTIVE_STATUSES for t in p3):
        return None
    end = last_end(p3, freeze, is_phase3)
    if end is None or not is_mature(end, freeze, PRIMARY_WINDOW_MONTHS):
        return None
    return 0


def time_to_phase3(trials: list[Trial], freeze: date) -> tuple[date | None, date | None, int | None, int | None]:
    """PREREG secondary: time from the earliest dated Phase II start to the first Phase III
    start, censored at the freeze. Returns (phase2_start, phase3_start, days, event).
    IMPLEMENTATION CHOICE: withdrawn trials (never enrolled) and starts after the freeze are
    not starts; the first Phase III start is the first on or after the Phase II origin."""
    def starts(pred):
        return [t.start_date for t in trials
                if pred(t) and t.start_date is not None and t.start_date <= freeze
                and t.overall_status != "WITHDRAWN"]

    p2 = starts(is_phase2)
    if not p2:
        return None, None, None, None
    origin = min(p2)
    p3 = [d for d in starts(is_phase3) if d >= origin]
    if p3:
        first = min(p3)
        return origin, first, (first - origin).days, 1
    return origin, None, (freeze - origin).days, 0


def code_hypothesis(max_stages: list[str], trials: list[Trial], freeze: date) -> HypothesisOutcome:
    """All stage C outcomes for one hypothesis. `max_stages` are the Open Targets 26.09
    maxClinicalStage values of every clinical_indication row of the hypothesis's drugs for the
    indication; `trials` its linked ClinicalTrials.gov records."""
    trials = [t for t in trials if registered_by_freeze(t, freeze)]
    p2_start, p3_start, days, event = time_to_phase3(trials, freeze)
    # Descriptive table 7 (follow-up from the last Phase II to the freeze): the same date the
    # maturation rule compares with the window, for every hypothesis with one.
    p2_end = last_end(trials, freeze, is_phase2)
    if not max_stages:
        # PREREG §Missing data: no retrievable phase -> excluded and counted.
        return HypothesisOutcome({w: Status.no_phase for w in MATURATION_WINDOWS_MONTHS},
                                 {w: None for w in MATURATION_WINDOWS_MONTHS}, StopCode.none, "",
                                 False, None, p2_start, p3_start, days, event, p2_end)
    advanced = any(s in ADVANCED_STAGES for s in max_stages)
    approved = APPROVAL_STAGE in max_stages
    if advanced:
        code, texts, all_business = StopCode.none, "", False
    else:
        code, texts, all_business = stop_decomposition(trials, freeze)
    statuses = {w: primary_status(advanced, trials, all_business, freeze, w)
                for w in MATURATION_WINDOWS_MONTHS}
    return HypothesisOutcome(
        statuses=statuses, advanced={w: advanced_value(s) for w, s in statuses.items()},
        stop_code=code, why_stopped_texts=texts, approved=approved,
        phase3_success=phase3_success(advanced, approved, trials, freeze),
        phase2_start=p2_start, phase3_start=p3_start, time_to_phase3_days=days, phase3_event=event,
        last_phase2_end=p2_end)
