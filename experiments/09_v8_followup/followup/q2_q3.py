"""Q2 (every hypothesis with PP.H4 >= 0.5) and Q3 ("S17 with Phase III success": the V8 S17 set
with y = phase3_success; its Bayesian fits run through stage_d, its balanced accuracy here)."""
import numpy as np
import pandas as pd

from followup.common import BOOT_FAIL_MAX, BOOT_REPS, SEED

Q2_COLUMNS = ("hypothesis_id", "gene_symbol", "gene_ensembl", "indication_id", "indication_name", "drug_program_ids",
              "mechanism_class", "direction", "instrument_source", "instrument_assay_id", "platform", "lead_variant",
              "outcome_accession", "outcome_neff", "pp_h0", "pp_h1", "pp_h2", "pp_h3", "pp_h4", "genetic_direction",
              "evidence_state", "S_05", "advanced_24", "phase3_success", "heldout", "in_s1", "unit_id")


def q2_table(df: pd.DataFrame) -> pd.DataFrame:
    t = df.loc[df["coloc_run"].astype(bool) & (df["pp_h4"].astype(float) >= 0.5), list(Q2_COLUMNS)]
    return t.sort_values(["pp_h4", "gene_ensembl", "indication_id", "hypothesis_id"],
                         ascending=[False, True, True, True], kind="mergesort").reset_index(drop=True)


def s17_phase3_frame(s17: pd.DataFrame) -> pd.DataFrame:
    """The S17 frame with y replaced by phase3_success; rows with it missing are left out by the
    design builder (non-null y only), as for V8's h1_p3s."""
    f = s17.copy()
    f["y"] = pd.to_numeric(f["phase3_success"], errors="raise").astype(float)
    return f


def balanced_accuracy(pred: np.ndarray, y: np.ndarray) -> float | None:
    if not ((y == 1).any() and (y == 0).any()):
        return None
    return 0.5 * (pred[y == 1].mean() + (1 - pred[y == 0]).mean())


def ba_by_class(frame: pd.DataFrame) -> dict[str, float | None]:
    out = {}
    for cls in ("aligned", "blocking"):
        g = frame.loc[frame["cls"] == cls]
        out[cls] = balanced_accuracy(g["S"].to_numpy(dtype=float), g["y"].to_numpy(dtype=float))
    out["difference"] = None if None in (out["aligned"], out["blocking"]) else out["aligned"] - out["blocking"]
    return out


def ba_bootstrap(frame: pd.DataFrame, reps: int = BOOT_REPS) -> dict:
    """Gene-clustered bootstrap of the aligned - blocking difference; a replicate lacking either
    outcome value in either stratum is undefined and counted."""
    pop = frame.loc[frame["cls"].isin(["aligned", "blocking"]) & frame["y"].notna()].reset_index(drop=True)
    genes = np.array(sorted(pop["gene_ensembl"].unique()))
    by_gene = pop.groupby("gene_ensembl").indices
    diffs = []
    for rep in range(reps):
        rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(10, rep)))
        idx = np.concatenate([by_gene[g] for g in genes[rng.integers(0, len(genes), len(genes))]])
        d = ba_by_class(pop.iloc[idx])["difference"]
        if d is not None:
            diffs.append(d)
    undefined = reps - len(diffs)
    available = undefined / reps <= BOOT_FAIL_MAX
    point = ba_by_class(pop)
    return {"n": int(len(pop)), "by_class_n": pop["cls"].value_counts().to_dict(), "point": point,
            "replicates": reps, "undefined": undefined, "interval_available": available,
            "difference_ci95": [float(np.quantile(diffs, 0.025)), float(np.quantile(diffs, 0.975))] if available else None}
