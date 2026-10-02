"""Gates that withhold a fit or a confirmatory decision (PREREG §Inference criteria, §Statistical models).

Order of application for H1 (and H4 with the neurologic-or-psychiatric strata):
  1. structural non-estimability: a stratum has no supportive or no non-supportive hypothesis
     -> not fitted, reported as not estimable with the counts;
  2. reliability: a stratum has < 10 genes with a supportive hypothesis or < 20 supportive
     hypotheses -> fitted and reported as a low-information estimate, no confirmatory decision;
  3. sampler diagnostics (R-hat < 1.01, bulk ESS > 400, zero divergences), rerun at target
     acceptance 0.99 with 4,000 draws, then the reduced model (gene slope dropped), which is
     descriptive only;
  4. prior dominance: posterior SD of the focal coefficient > 0.8 x prior SD -> no decision.
Counts for 1 and 2 are taken on the pre-outcome frame (stage B states, before outcomes join).
"""
import pandas as pd
from pydantic import BaseModel

from stage_d.constants import (ESS_BULK_MIN, H3_MIN_EACH, PRIOR_DOMINANCE_FRACTION, RELIABILITY_MIN_SUPPORTIVE_GENES,
                               RELIABILITY_MIN_SUPPORTIVE_HYPOTHESES, RHAT_MAX)


class StratumCounts(BaseModel):
    stratum: str
    n: int
    n_supportive: int
    n_nonsupportive: int
    n_supportive_genes: int
    n_supportive_indications: int
    n_supportive_programs: int


class StratumGate(BaseModel):
    model: str
    strata: list[StratumCounts]
    structurally_estimable: bool
    reliable: bool
    reasons: list[str]


class H3Gate(BaseModel):
    n_supportive: int
    n_contradictory: int
    passes: bool


class Diagnostics(BaseModel):
    rhat_max: float
    ess_bulk_min: float
    divergences: int


def stratum_counts(frame: pd.DataFrame, stratum: str) -> StratumCounts:
    sup = frame.loc[frame["S"] == 1]
    return StratumCounts(
        stratum=stratum, n=len(frame), n_supportive=len(sup), n_nonsupportive=int((frame["S"] == 0).sum()),
        n_supportive_genes=int(sup["gene_ensembl"].nunique()),
        n_supportive_indications=int(sup["indication_id"].nunique()),
        n_supportive_programs=len({p for ps in sup["programs"] for p in ps}))


def _gate(model: str, parts: dict[str, pd.DataFrame]) -> StratumGate:
    strata = [stratum_counts(f, name) for name, f in parts.items()]
    reasons = []
    structural = True
    reliable = True
    for c in strata:
        if c.n_supportive == 0 or c.n_nonsupportive == 0:
            structural = False
            reasons.append(f"{c.stratum}: {c.n_supportive} supportive, {c.n_nonsupportive} non-supportive")
        if c.n_supportive_genes < RELIABILITY_MIN_SUPPORTIVE_GENES:
            reliable = False
            reasons.append(f"{c.stratum}: {c.n_supportive_genes} genes with a supportive hypothesis "
                           f"(< {RELIABILITY_MIN_SUPPORTIVE_GENES})")
        if c.n_supportive < RELIABILITY_MIN_SUPPORTIVE_HYPOTHESES:
            reliable = False
            reasons.append(f"{c.stratum}: {c.n_supportive} supportive hypotheses "
                           f"(< {RELIABILITY_MIN_SUPPORTIVE_HYPOTHESES})")
    return StratumGate(model=model, strata=strata, structurally_estimable=structural,
                       reliable=structural and reliable, reasons=reasons)


def h1_gate(pre_outcome: pd.DataFrame) -> StratumGate:
    return _gate("h1", {"aligned": pre_outcome.loc[pre_outcome["cls"] == "aligned"],
                        "blocking": pre_outcome.loc[pre_outcome["cls"] == "blocking"]})


def h4_gate(pre_outcome: pd.DataFrame) -> StratumGate:
    return _gate("h4", {"neuro_psych": pre_outcome.loc[pre_outcome["C"] == 1],
                        "other_indications": pre_outcome.loc[pre_outcome["C"] == 0]})


def h3_gate(pre_outcome: pd.DataFrame) -> H3Gate:
    n_sup = int((pre_outcome["state"] == "supportive").sum())
    n_con = int((pre_outcome["state"] == "contradictory").sum())
    return H3Gate(n_supportive=n_sup, n_contradictory=n_con, passes=n_sup >= H3_MIN_EACH and n_con >= H3_MIN_EACH)


def sampler_passed(d: Diagnostics) -> bool:
    return d.rhat_max < RHAT_MAX and d.ess_bulk_min > ESS_BULK_MIN and d.divergences == 0


def prior_dominated(posterior_sd: float, prior_sd: float) -> bool:
    return posterior_sd > PRIOR_DOMINANCE_FRACTION * prior_sd
