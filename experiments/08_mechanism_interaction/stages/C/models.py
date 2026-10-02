"""Row schemas and fixed constants for stage C (outcome coding).

Every constant here is transcribed from ../../PREREG.md (frozen at b946087) or fixed by
../INTERFACES.md. Nothing in this package reads stage B output or any evidence value.
"""
from datetime import date
from enum import Enum

from pydantic import BaseModel, ConfigDict, field_validator

# PREREG §Starting and stopping rules: the outcome-freeze date is the date of the freeze commit.
OUTCOME_FREEZE_DATE = date(2026, 9, 30)

# PREREG §Measured variables: 24 months primary; 12, 36, 48 in S14.
PRIMARY_WINDOW_MONTHS = 24
MATURATION_WINDOWS_MONTHS = (12, 24, 36, 48)

# Open Targets 26.09 clinical_indication.maxClinicalStage values.
ADVANCED_STAGES = frozenset({"PHASE_3", "PREAPPROVAL", "APPROVAL"})
PHASE2_PLUS_STAGES = frozenset({"PHASE_2", "PHASE_2_3", "PHASE_3", "PREAPPROVAL", "APPROVAL"})
APPROVAL_STAGE = "APPROVAL"

# ClinicalTrials.gov API v2 OverallStatus values (AACT `studies.overall_status`, normalized by
# ctgov.parse_status).
# PREREG: "recruiting, active, not yet recruiting or enrolling by invitation".
ACTIVE_STATUSES = frozenset({"RECRUITING", "ACTIVE_NOT_RECRUITING", "NOT_YET_RECRUITING",
                             "ENROLLING_BY_INVITATION"})
# PREREG §Outcome decomposition: "terminated, suspended or withdrawn trials" carry why_stopped.
STOPPED_STATUSES = frozenset({"TERMINATED", "SUSPENDED", "WITHDRAWN"})
# PREREG: "the last Phase II completion or termination date".
ENDED_STATUSES = frozenset({"COMPLETED", "TERMINATED"})

# PREREG §Outcome decomposition keyword table, in table order.
STOP_KEYWORDS: dict[str, tuple[str, ...]] = {
    "business": ("business", "strategic", "commercial", "portfolio", "sponsor decision",
                 "funding", "financial"),
    "efficacy": ("efficacy", "futility", "futile", "lack of benefit", "did not meet"),
    "safety": ("safety", "adverse", "toxicity", "tolerability"),
}
# PREREG: "the order efficacy > safety > unclassified > business assigns one code".
STOP_PRECEDENCE = ("efficacy", "safety", "unclassified", "business")


class Status(str, Enum):
    advanced = "advanced"
    no_observed_advancement = "no_observed_advancement"
    active = "active"
    business_only = "business_only"
    undated = "undated"
    # PREREG §Missing data: "A hypothesis with no retrievable phase is excluded and counted."
    no_phase = "no_phase"


class StopCode(str, Enum):
    efficacy = "efficacy"
    safety = "safety"
    unclassified = "unclassified"
    business = "business"
    completed_no_successor = "completed_no_successor"
    none = ""


class Trial(BaseModel):
    """One ClinicalTrials.gov registration as the AACT snapshot of 2026-09-30 holds it, reduced to
    the fields stage C uses (AACT `studies` column in each comment)."""
    model_config = ConfigDict(frozen=True)

    nct_id: str
    phases: tuple[str, ...]           # phase, split on "/", e.g. ("PHASE2",) or ("PHASE2", "PHASE3")
    overall_status: str               # overall_status (API v2 enum)
    why_stopped: str = ""             # why_stopped
    start_date: date | None = None    # start_month_year, first day of a partial date
    completion_date: date | None = None  # completion_month_year, first day of a partial date
    first_submit_date: date | None = None  # study_first_submitted_date
    last_update_submit_date: date | None = None  # last_update_submitted_date (diagnostic only)

    @field_validator("nct_id")
    @classmethod
    def _nct(cls, v: str) -> str:
        if not (v.startswith("NCT") and len(v) == 11 and v[3:].isdigit()):
            raise ValueError(f"not an NCT id: {v!r}")
        return v


class HypothesisDrugs(BaseModel):
    """The three stage A fields stage C may read, plus what stage C derives from Open Targets."""
    model_config = ConfigDict(frozen=True)

    hypothesis_id: str
    indication_id: str
    drug_program_ids: tuple[str, ...]
    max_stages: tuple[str, ...]       # maxClinicalStage of every clinical_indication row linked
    nct_ids: tuple[str, ...]          # linked ClinicalTrials.gov ids (trial-linking rule)


class OutcomeRow(BaseModel):
    """One row of outcomes.csv, columns exactly as in INTERFACES.md §Stage C."""
    hypothesis_id: str
    status_24: Status
    advanced_24: int | None
    status_12: Status
    status_36: Status
    status_48: Status
    advanced_12: int | None
    advanced_36: int | None
    advanced_48: int | None
    stop_code: StopCode
    why_stopped_texts: str
    efficacy_coded: bool
    phase2_start_date: date | None
    phase3_start_date: date | None
    time_to_phase3_days: int | None
    phase3_event: int | None
    phase3_success: int | None
    approved: bool
    chembl_max_phase_for_ind: float | None
    earliest_phase2_start: date | None
    last_phase2_end_date: date | None   # descriptive table 7; the date `rules.last_end` gives for Phase II


OUTCOME_COLUMNS = list(OutcomeRow.model_fields)


class FreezeDiagnosticRow(BaseModel):
    """One row of post_freeze_updates.csv (reported diagnostic; no label reads it)."""
    hypothesis_id: str
    n_linked_trials: int
    n_not_in_archive: int
    n_first_submitted_after_freeze: int
    n_last_update_after_freeze: int


FREEZE_DIAGNOSTIC_COLUMNS = list(FreezeDiagnosticRow.model_fields)


class LinkageAuditRow(BaseModel):
    """One row of linkage_audit.csv: the count at one step of the trial-linking rule (descriptive)."""
    step: str
    count: int
    description: str


LINKAGE_AUDIT_COLUMNS = list(LinkageAuditRow.model_fields)


class DifferenceRow(BaseModel):
    """One row of partial_date_boundary.csv or phase23_diagnostic.csv: a quantity of one hypothesis
    whose value under the registered rule differs from its value under the alternative reading the
    file is about (descriptive; no label reads it)."""
    hypothesis_id: str
    quantity: str
    registered: str
    alternative: str


DIFFERENCE_COLUMNS = list(DifferenceRow.model_fields)


class StageCError(Exception):
    """Raised when an input violates a stage C contract."""


class BlindingViolation(StageCError):
    """Raised when stage C is pointed at a stage B path or an evidence column."""
