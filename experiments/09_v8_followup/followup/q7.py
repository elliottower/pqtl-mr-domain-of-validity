"""Q7 data: the population and terms of V8's registered Cox secondary (`stage_d.pipeline.cox_models`)
with the plan's event coding. Event 1 = Phase III entry at time_to_phase3_days; event 2 = the
competing ending (status_24 no_observed_advancement or business_only) at last_phase2_end_date minus
earliest_phase2_start; otherwise 0, censored at time_to_phase3_days, which stage C measures to the
outcome freeze. The R script q7.R fits the models."""
import numpy as np
import pandas as pd

from followup.common import FollowupError
from stage_d.designs import COVARIATES

FREEZE = pd.Timestamp("2026-09-30")
COMPETING = ("no_observed_advancement", "business_only")


def survival_table(s1: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    dated = s1.loc[s1["time_to_phase3_days"].notna()].copy()
    t = dated["time_to_phase3_days"].astype(float)
    ev = dated["phase3_event"].astype(float)
    start = pd.to_datetime(dated["earliest_phase2_start"])
    end = pd.to_datetime(dated["last_phase2_end_date"])
    competing = (ev == 0) & dated["status_24"].astype(str).isin(COMPETING)
    counts = {"rows": int(len(dated)), "phase3_entry": int((ev == 1).sum()), "competing_status": int(competing.sum()),
              "competing_end_missing": int((competing & end.isna()).sum()),
              "competing_end_after_freeze": int((competing & (end > FREEZE)).sum()),
              "competing_end_before_start": int((competing & (end < start)).sum())}
    usable = competing & end.notna() & (end <= FREEZE) & (end >= start)
    status = np.where(ev == 1, 1, np.where(usable, 2, 0))
    time = np.where(usable, (end - start).dt.days.astype(float), t)
    time = np.where(competing & (end < start), 0.0, time)
    if (np.asarray(time) < 0).any():
        raise FollowupError(f"{int((np.asarray(time) < 0).sum())} negative survival times")
    counts["competing_event"] = int((status == 2).sum())
    counts["censored"] = int((status == 0).sum())
    out = pd.DataFrame({"hypothesis_id": dated["hypothesis_id"].to_numpy(), "gene_ensembl": dated["gene_ensembl"].to_numpy(),
                        "time": time, "status": status, "cls": dated["cls"].to_numpy(),
                        "S": dated["S"].to_numpy(dtype=float),
                        "A": (dated["cls"] == "aligned").to_numpy(dtype=float)})
    out["SxA"] = out["S"] * out["A"]
    for col, label in COVARIATES:
        out[label] = dated[col].to_numpy(dtype=float)
    out["group"] = out["cls"].astype(str) + "|S=" + out["S"].astype(int).astype(str)
    return out, counts


def registered_terms(frame: pd.DataFrame, model: str) -> list[str]:
    """The terms cox_models uses: S, A, SxA (H1) or S (H2), plus each covariate that varies."""
    terms = ["S", "A", "SxA"] if model == "H1" else ["S"]
    return terms + [label for _, label in COVARIATES if np.ptp(frame[label].to_numpy(dtype=float)) > 0]
