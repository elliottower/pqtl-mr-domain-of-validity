"""Fit specifications and the attempt sequence (PREREG §Statistical models).

Every Bayesian model runs the same sequence, one checkpoint file per attempt:
  full           gene slope kept, 4 chains x 2,000 draws after 2,000 tuning, target 0.95
  full_rerun     if `full` misses R-hat < 1.01, bulk ESS > 400 or zero divergences:
                 target 0.99, 4,000 draws
  reduced        if the rerun still fails: gene slope dropped, primary settings (descriptive only)
  reduced_rerun  if `reduced` fails its diagnostics: target 0.99, 4,000 draws
The sequence stops at the first attempt that passes. Every attempt run is reported.
"""
from pydantic import BaseModel, ConfigDict

from stage_d.constants import PRIMARY_SAMPLER, RERUN_SAMPLER, SamplerSettings

ATTEMPTS = (("full", "primary", True), ("full_rerun", "rerun", True),
            ("reduced", "primary", False), ("reduced_rerun", "rerun", False))


class FitSpec(BaseModel):
    model_config = ConfigDict(frozen=True)
    fit_id: str
    design_id: str
    kind: str
    set_id: str
    prior: str
    role: str
    primary: SamplerSettings = PRIMARY_SAMPLER
    rerun: SamplerSettings = RERUN_SAMPLER

    def sampler(self, which: str) -> SamplerSettings:
        return self.primary if which == "primary" else self.rerun


def next_attempt(done: list[dict]) -> tuple[str, str, bool] | None:
    """The attempt to run given the attempts already on disk (in order), or None when finished."""
    if any(a["passed"] for a in done):
        return None
    return ATTEMPTS[len(done)] if len(done) < len(ATTEMPTS) else None


def fit_outcome(done: list[dict]) -> str:
    passed = [a["attempt"] for a in done if a["passed"]]
    if passed:
        return passed[0]
    return "failed" if done else "not_fitted"
