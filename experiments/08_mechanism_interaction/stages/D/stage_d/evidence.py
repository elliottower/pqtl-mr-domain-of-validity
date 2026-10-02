"""Evidence state and score from colocalization output (PREREG §Measured variables).

supportive     PP.H4 >= threshold and the genetic direction matches the intervention direction
               (decrease requires a risk-increasing protein, +1; increase a protective one, -1)
contradictory  PP.H4 >= threshold and the direction is opposite
inconclusive   otherwise (PP.H4 below threshold, coloc not run, either direction undetermined)

E = d * PP.H4 with d = +1 matching, -1 opposite, 0 otherwise.
"""
import numpy as np
import pandas as pd

PRIMARY_PP_H4_THRESHOLD = 0.80


def _match_opposite(genetic_direction, intervention_direction) -> tuple[np.ndarray, np.ndarray]:
    gd = np.asarray(genetic_direction, dtype=float)
    idir = np.asarray(intervention_direction, dtype=object)
    match = ((idir == "decrease") & (gd == 1)) | ((idir == "increase") & (gd == -1))
    opposite = ((idir == "decrease") & (gd == -1)) | ((idir == "increase") & (gd == 1))
    return match, opposite


def evidence_state(pp_h4, genetic_direction, intervention_direction,
                   threshold: float = PRIMARY_PP_H4_THRESHOLD) -> np.ndarray:
    pp = pd.to_numeric(pd.Series(pp_h4), errors="raise").to_numpy(dtype=float)
    match, opposite = _match_opposite(genetic_direction, intervention_direction)
    coloc = np.nan_to_num(pp, nan=-1.0) >= threshold
    return np.where(coloc & match, "supportive", np.where(coloc & opposite, "contradictory", "inconclusive"))


def evidence_score(pp_h4, genetic_direction, intervention_direction) -> np.ndarray:
    pp = np.nan_to_num(pd.to_numeric(pd.Series(pp_h4)).to_numpy(dtype=float), nan=0.0)
    match, opposite = _match_opposite(genetic_direction, intervention_direction)
    return np.where(match, 1.0, np.where(opposite, -1.0, 0.0)) * pp
