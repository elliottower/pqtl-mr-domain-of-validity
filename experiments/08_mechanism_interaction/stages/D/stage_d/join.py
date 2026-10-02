"""Join of stages A, B and C on hypothesis_id, and the columns D derives from the joined table.

The join is one-to-one and exhaustive: B and C each carry exactly A's hypothesis ids, or the
join raises. B's evidence state and score are recomputed from PP.H4 and the two directions; a
disagreement raises, because S15 recomputes states the same way.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import BaseModel

from stage_d.evidence import evidence_score, evidence_state
from stage_d.schemas import EvidenceRow, HypothesisRow, OutcomeRow, load_rows


class JoinError(ValueError):
    pass


class Imputation(BaseModel):
    reference_rows: int
    neff_missing_reference: int
    neff_missing_all: int
    neff_median: float
    log10_neff_mean: float
    log10_neff_sd: float
    platform_missing_reference: int
    platform_missing_all: int
    platform_mode: str


class DerivedInfo(BaseModel):
    n_hypotheses: int
    n_evidence: int
    n_outcomes: int
    n_joined: int
    pre_pqtl_publication_true: int
    pre_pqtl_publication_false: int
    pre_pqtl_publication_missing_trial_date: int
    pre_pqtl_publication_missing_publication_date: int
    imputation: Imputation
    present_optional_columns: dict[str, list[str]]


OPTIONAL_COLUMNS = {
    "A": ("variant", "in_s8", "in_s9", "in_s21"),
    "B": ("s17_sentinel_p", "evidence_state_ukbppp", "evidence_state_decode"),
    "C": ("last_phase2_end_date",),
}
TABLE12_SOURCES = ("ukbppp", "decode")


def load_stage_tables(stages_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    hyp = load_rows(stages_root / "A" / "output" / "hypotheses.csv", HypothesisRow)
    ev = load_rows(stages_root / "B" / "output" / "evidence.csv", EvidenceRow)
    out = load_rows(stages_root / "C" / "output" / "outcomes.csv", OutcomeRow)
    return hyp, ev, out


def join_stages(hyp: pd.DataFrame, ev: pd.DataFrame, out: pd.DataFrame) -> pd.DataFrame:
    ids = set(hyp["hypothesis_id"])
    for name, frame in (("B evidence", ev), ("C outcomes", out)):
        other = set(frame["hypothesis_id"])
        if other != ids:
            raise JoinError(f"{name}: {len(other - ids)} ids not in A, {len(ids - other)} A ids missing")
    joined = hyp.merge(ev, on="hypothesis_id", validate="one_to_one").merge(
        out, on="hypothesis_id", validate="one_to_one")
    if len(joined) != len(hyp):
        raise JoinError("join changed the row count")
    joined.attrs["present_columns"] = {
        "A": hyp.attrs.get("present_columns", []), "B": ev.attrs.get("present_columns", []),
        "C": out.attrs.get("present_columns", [])}
    return joined


def check_evidence_consistency(df: pd.DataFrame) -> None:
    state = evidence_state(df["pp_h4"], df["genetic_direction"], df["direction"])
    bad = state != df["evidence_state"].to_numpy()
    if bad.any():
        raise JoinError(f"{int(bad.sum())} rows: B evidence_state disagrees with PP.H4 and directions "
                        f"(first {df.loc[bad, 'hypothesis_id'].iloc[0]})")
    score = evidence_score(df["pp_h4"], df["genetic_direction"], df["direction"])
    bad_e = ~np.isclose(score, df["E"].to_numpy(dtype=float), atol=1e-9)
    if bad_e.any():
        raise JoinError(f"{int(bad_e.sum())} rows: B score E disagrees with direction x PP.H4")
    present = df.attrs.get("present_columns", {}).get("B", [])
    for src in TABLE12_SOURCES:
        col = f"evidence_state_{src}"
        if col not in present:
            continue
        own = (df["instrument_source"] == src).to_numpy()
        bad_s = own & (df[col].astype(object).to_numpy() != df["evidence_state"].to_numpy())
        if bad_s.any():
            raise JoinError(f"{int(bad_s.sum())} rows: {col} disagrees with evidence_state on the selected "
                            f"{src} instrument (first {df.loc[bad_s, 'hypothesis_id'].iloc[0]})")


def programs_of(value: str) -> list[str]:
    return [p for p in str(value).split(";") if p]


def derive(joined: pd.DataFrame) -> tuple[pd.DataFrame, DerivedInfo]:
    check_evidence_consistency(joined)
    df = joined.copy()
    df["programs"] = df["drug_program_ids"].map(programs_of)
    df["k_programs"] = df["programs"].map(len).astype(int)
    df["neurological"] = df["neuro_psych"] & ~df["psych_only"]

    trial = df["earliest_phase2_start"]
    pub = df["pre_pqtl_publication_date"]
    known = trial.notna() & pub.notna()
    df["pre_pqtl_publication"] = pd.array(np.where(known, trial < pub, pd.NA), dtype="boolean")

    reference = df["in_s1"] & df["heldout"]
    if not reference.any():
        raise JoinError("no held-out S1 hypotheses; covariate standardization has no reference set")
    neff_ref = df.loc[reference, "outcome_neff"].dropna()
    if neff_ref.empty or (neff_ref <= 0).any():
        raise JoinError("held-out effective N missing everywhere or non-positive")
    neff_median = float(neff_ref.median())
    neff = df["outcome_neff"].fillna(neff_median)
    log_neff = np.log10(neff.to_numpy(dtype=float))
    ref_log = log_neff[reference.to_numpy()]
    mean = float(ref_log.mean())
    sd = float(ref_log.std(ddof=1)) if len(ref_log) > 1 else 0.0
    df["z_log10_neff"] = (log_neff - mean) / sd if sd > 0 else 0.0
    platform_mode = str(df.loc[reference, "platform"].dropna().mode().sort_values().iloc[0])
    platform = df["platform"].fillna(platform_mode)
    df["platform_somascan"] = (platform == "SomaScan").astype(float)
    df["oncology_f"] = df["oncology"].astype(float)
    df["S_primary"] = df["S"].astype(float)

    info = DerivedInfo(
        n_hypotheses=len(df), n_evidence=len(df), n_outcomes=len(df), n_joined=len(df),
        pre_pqtl_publication_true=int((df["pre_pqtl_publication"] == True).sum()),  # noqa: E712
        pre_pqtl_publication_false=int((df["pre_pqtl_publication"] == False).sum()),  # noqa: E712
        pre_pqtl_publication_missing_trial_date=int(trial.isna().sum()),
        pre_pqtl_publication_missing_publication_date=int(pub.isna().sum()),
        imputation=Imputation(
            reference_rows=int(reference.sum()),
            neff_missing_reference=int(df.loc[reference, "outcome_neff"].isna().sum()),
            neff_missing_all=int(df["outcome_neff"].isna().sum()),
            neff_median=neff_median, log10_neff_mean=mean, log10_neff_sd=sd,
            platform_missing_reference=int(df.loc[reference, "platform"].isna().sum()),
            platform_missing_all=int(df["platform"].isna().sum()),
            platform_mode=platform_mode),
        present_optional_columns={
            s: [c for c in cols if c in joined.attrs.get("present_columns", {}).get(s, [])]
            for s, cols in OPTIONAL_COLUMNS.items()},
    )
    return df, info
