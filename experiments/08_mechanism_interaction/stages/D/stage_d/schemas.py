"""Row schemas for the three stage outputs D reads (stages/INTERFACES.md), and the loader.

Every row of every input is validated; a row that fails raises SchemaError naming the file,
the row and the field. Blank CSV cells are None.
"""
import types
from datetime import date
from pathlib import Path
from typing import Annotated, Literal, Union, get_args, get_origin

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator


class SchemaError(ValueError):
    pass


class _Row(BaseModel):
    model_config = ConfigDict(extra="ignore", frozen=True)

    @model_validator(mode="before")
    @classmethod
    def _blank_to_none(cls, data):
        return {k: (None if isinstance(v, str) and v.strip() == "" else v) for k, v in data.items()}


MechanismClass = Literal["aligned", "blocking", "other"]
EvidenceState = Literal["supportive", "contradictory", "inconclusive"]
Binary = Annotated[int, Field(ge=0, le=1)]
Tier = Annotated[int, Field(ge=1, le=3)]
Sign = Annotated[int, Field(ge=-1, le=1)]
# no_phase: PREREG §Missing data, "A hypothesis with no retrievable phase is excluded and counted."
OutcomeStatus = Literal["advanced", "no_observed_advancement", "active", "business_only", "undated", "no_phase"]
S19Arm = Literal["neutralizing_biologic", "small_molecule_blocker"]


class HypothesisRow(_Row):
    """Stage A `hypotheses.csv`."""
    hypothesis_id: str
    gene_symbol: str
    gene_ensembl: str
    indication_id: str
    indication_name: str
    direction: Literal["decrease", "increase", "ambiguous"]
    mechanism_class: MechanismClass
    class_reason: str
    drug_program_ids: str
    eligible: bool
    heldout: bool
    specificity_rank: int
    instrument_source: Literal["ukbppp", "decode", "interval"]
    instrument_assay_id: str
    platform: Literal["Olink", "SomaScan"] | None = None
    outcome_accession: str
    outcome_source: str
    outcome_tier: Tier
    overlap: Literal["no", "unknown", "yes"]
    outcome_n_case: int | None = None
    outcome_n_control: int | None = None
    outcome_neff: float | None = None
    blood_secreted_hpa: bool
    secreted_uniprot: bool
    neuro_psych: bool
    neuro_only: bool
    psych_only: bool
    oncology: bool
    pilot_indication: bool
    karim_launched: bool
    single_protein_row: bool
    conflicted_row_restored: bool
    strict_instrument: bool
    in_s1: bool
    in_s4: bool
    in_s10: bool
    in_s12: bool
    subtype_restricted: bool
    phenotype_broader: bool
    sample_overlap: bool
    overlap_unknown: bool
    pre_pqtl_publication_date: date | None = None
    # S8 and S9 (INTERFACES.md stage A): hypotheses re-formed under each set's class table.
    # A row with variant s8 / s9 exists only there. Absent -> set not formed.
    variant: Literal["primary", "s8", "s9", "s21"] | None = None
    in_s8: bool | None = None
    in_s9: bool | None = None
    in_s21: bool | None = None
    # S19 arm (INTERFACES.md stage A); empty -> None, in neither arm.
    s19_arm: S19Arm | None

    @model_validator(mode="after")
    def _check(self):
        if not self.drug_program_ids.strip(";"):
            raise ValueError("drug_program_ids is empty")
        if self.specificity_rank < 0:
            raise ValueError("specificity_rank < 0")
        if (self.neuro_only or self.psych_only) and not self.neuro_psych:
            raise ValueError("neuro_only or psych_only set without neuro_psych")
        if self.neuro_only and self.psych_only:
            raise ValueError("neuro_only and psych_only both set")
        return self


class EvidenceRow(_Row):
    """Stage B `evidence.csv`."""
    hypothesis_id: str
    coloc_run: bool
    not_run_reason: str | None = None
    pp_h0: float | None = None
    pp_h1: float | None = None
    pp_h2: float | None = None
    pp_h3: float | None = None
    pp_h4: float | None = None
    n_shared: int | None = None
    frac_pqtl_retained: float | None = None
    frac_outcome_retained: float | None = None
    sentinel_or_proxy_retained: bool | None = None
    low_coverage: bool | None = None
    lead_variant: str | None = None
    genetic_direction: Sign
    evidence_state: EvidenceState
    S: Binary
    E: float
    protein_altering: bool | None = None
    platform_concordant: Literal["concordant", "discordant", "untested"] | None = None
    splicing_candidate: bool | None = None
    s15a_pp_h4: float | None = None
    s15b_pp_h4: float | None = None
    s15c_pp_h4: float | None = None
    s15d_pp_h4: float | None = None
    s15e_pp_h4: float | None = None
    s15f_pp_h4: float | None = None
    s15g_pp_h4: float | None = None
    s15f_low_coverage_excluded: bool | None = None
    s16_evidence_state: EvidenceState | None = None
    # S17 (INTERFACES.md stage B). Absent -> set not formed.
    s17_sentinel_p: float | None = None
    # Descriptive table 12 (INTERFACES.md stage B): the primary rule with each source's instrument;
    # empty where the hypothesis names no assay in that source. Absent -> agreement not reported.
    evidence_state_ukbppp: EvidenceState | None = None
    evidence_state_decode: EvidenceState | None = None
    # Frozen-rule sensitivity set (INTERFACES.md stage B): whether the outcome file passed the frozen
    # CI-versus-p rule; empty where that check did not run on it. Absent -> set not formed.
    outcome_file_frozen_ci_p_pass: bool | None = None

    @model_validator(mode="after")
    def _check(self):
        if self.S != int(self.evidence_state == "supportive"):
            raise ValueError("S disagrees with evidence_state")
        if not -1.0 <= self.E <= 1.0:
            raise ValueError("E outside [-1, 1]")
        if not self.coloc_run and (self.evidence_state != "inconclusive" or self.E != 0.0):
            raise ValueError("coloc not run but state is not inconclusive or E != 0")
        if self.coloc_run and self.pp_h4 is None:
            raise ValueError("coloc_run without pp_h4")
        return self


class OutcomeRow(_Row):
    """Stage C `outcomes.csv`."""
    hypothesis_id: str
    status_24: OutcomeStatus
    advanced_24: Binary | None = None
    status_12: OutcomeStatus
    status_36: OutcomeStatus
    status_48: OutcomeStatus
    advanced_12: Binary | None = None
    advanced_36: Binary | None = None
    advanced_48: Binary | None = None
    stop_code: Literal["efficacy", "safety", "unclassified", "business", "completed_no_successor"] | None = None
    why_stopped_texts: str | None = None
    efficacy_coded: bool
    phase2_start_date: date | None = None
    phase3_start_date: date | None = None
    time_to_phase3_days: float | None = None
    phase3_event: Binary | None = None
    phase3_success: Binary | None = None
    approved: bool
    chembl_max_phase_for_ind: float | None = None
    earliest_phase2_start: date | None = None
    # Follow-up column of table 7 (INTERFACES.md stage C). Absent -> reported missing.
    last_phase2_end_date: date | None = None

    @model_validator(mode="after")
    def _check(self):
        for m in (12, 24, 36, 48):
            status, adv = getattr(self, f"status_{m}"), getattr(self, f"advanced_{m}")
            expected = {"advanced": 1, "no_observed_advancement": 0}.get(status)
            if adv != expected:
                raise ValueError(f"advanced_{m}={adv} disagrees with status_{m}={status}")
        if (self.time_to_phase3_days is None) != (self.phase3_event is None):
            raise ValueError("time_to_phase3_days and phase3_event must be both set or both empty")
        if self.time_to_phase3_days is not None and self.time_to_phase3_days < 0:
            raise ValueError("negative time_to_phase3_days")
        return self


def _dtype_for(annotation) -> str:
    args = get_args(annotation)
    optional = get_origin(annotation) in (Union, types.UnionType) and type(None) in args
    base = next((a for a in args if a is not type(None)), annotation) if optional else annotation
    if get_origin(base) is Annotated:
        base = get_args(base)[0]
    if get_origin(base) is Literal:
        values = get_args(base)
        base = type(values[0])
    if base is bool:
        return "boolean" if optional else "bool"
    if base is int:
        return "Int64" if optional else "int64"
    if base is float:
        return "float64"
    if base is date:
        return "datetime64[ns]"
    return "object"


def load_rows(path: Path, model: type[_Row]) -> pd.DataFrame:
    """Read a CSV as text, validate every row against `model`, return typed columns."""
    raw = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [f for f, info in model.model_fields.items() if info.is_required() and f not in raw.columns]
    if missing:
        raise SchemaError(f"{path}: missing required columns {missing}")
    records = []
    for i, rec in enumerate(raw.to_dict(orient="records")):
        try:
            records.append(model.model_validate(rec).model_dump())
        except ValidationError as err:
            raise SchemaError(f"{path}: data row {i + 1}: {err}") from err
    frame = pd.DataFrame(records, columns=list(model.model_fields))
    for field, info in model.model_fields.items():
        dtype = _dtype_for(info.annotation)
        frame[field] = pd.to_datetime(frame[field]) if dtype == "datetime64[ns]" else frame[field].astype(dtype)
    frame.attrs["present_columns"] = sorted(set(raw.columns) & set(model.model_fields))
    if frame["hypothesis_id"].duplicated().any():
        dup = frame.loc[frame["hypothesis_id"].duplicated(), "hypothesis_id"].iloc[0]
        raise SchemaError(f"{path}: duplicate hypothesis_id {dup}")
    return frame
