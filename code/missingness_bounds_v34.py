"""
V3.4 Stage C: Missingness sensitivity / bounding analysis.

Computes best-case / worst-case BA under all possible assignments
of the 57 missing pairs, plus the Perplexity-flagged base-rate
asymmetry between analyzed and missing sets.

Usage:
    cd ~/Documents/GitHub/transport-wrapper
    uv run python -m DRUGS_EXPANDED.missingness_bounds_v34
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

HERE = Path(__file__).parent
OUT_DIR = HERE / "results" / "v34"


def compute_ba(tp, tn, fp, fn):
    sens = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
    spec = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
    return (sens + spec) / 2


def main():
    adj = pd.read_csv(OUT_DIR / "adjudicated_v34.csv")
    evaluable = adj[
        (adj["outcome"].isin(["SUCCESS", "FAILURE"]))
        & (adj["excluded"] != True)
        & (adj["excluded"] != "True")
    ].copy()

    classif = pd.read_csv(OUT_DIR / "classification_v34.csv")
    has_mr = set(classif["pair_id"])

    analyzed = evaluable[evaluable["pair_id"].isin(has_mr)]
    missing = evaluable[~evaluable["pair_id"].isin(has_mr)]

    n_analyzed = len(analyzed)
    n_missing = len(missing)
    n_total = len(evaluable)

    print(f"=== V3.4 Missingness Bounding Analysis ===\n")
    print(f"Evaluable total: {n_total}")
    print(f"Analyzed (have MR): {n_analyzed}")
    print(f"Missing (no MR):    {n_missing}\n")

    # ---- Base rate asymmetry ----
    s_analyzed = (analyzed["outcome"] == "SUCCESS").sum()
    f_analyzed = (analyzed["outcome"] == "FAILURE").sum()
    s_missing = (missing["outcome"] == "SUCCESS").sum()
    f_missing = (missing["outcome"] == "FAILURE").sum()
    s_total = (evaluable["outcome"] == "SUCCESS").sum()
    f_total = (evaluable["outcome"] == "FAILURE").sum()

    print("--- Base Rate Asymmetry ---")
    print(f"Evaluable set:  {s_total}S / {f_total}F  ({s_total/n_total:.1%} SUCCESS)")
    print(f"Analyzed set:   {s_analyzed}S / {f_analyzed}F  ({s_analyzed/n_analyzed:.1%} SUCCESS)")
    print(f"Missing set:    {s_missing}S / {f_missing}F  ({s_missing/n_missing:.1%} SUCCESS)")
    print(f"Gap: {s_missing/n_missing - s_analyzed/n_analyzed:+.1%}\n")

    table = np.array([[s_analyzed, f_analyzed], [s_missing, f_missing]])
    chi2, p_chi2, _, _ = chi2_contingency(table, correction=True)
    _, p_fisher = fisher_exact(table)
    print(f"Chi-square (Yates): {chi2:.3f}, p = {p_chi2:.4f}")
    print(f"Fisher exact: p = {p_fisher:.4f}\n")

    # ---- Current results ----
    tp = ((classif["prediction"] == "SUCCESS") & (classif["outcome"] == "SUCCESS")).sum()
    tn = ((classif["prediction"] == "FAILURE") & (classif["outcome"] == "FAILURE")).sum()
    fp = ((classif["prediction"] == "SUCCESS") & (classif["outcome"] == "FAILURE")).sum()
    fn = ((classif["prediction"] == "FAILURE") & (classif["outcome"] == "SUCCESS")).sum()

    ba_current = compute_ba(tp, tn, fp, fn)
    print(f"--- Current Results (n={n_analyzed}) ---")
    print(f"TP={tp}, TN={tn}, FP={fp}, FN={fn}")
    print(f"BA = {ba_current:.4f}\n")

    # ---- Bounding: if we had all 195 pairs ----
    # Missing pairs have no MR result, so prediction = FAILURE (not significant)
    # This means: missing SUCCESS → FN, missing FAILURE → TN
    print(f"--- Bounding Analysis (all {n_total} pairs) ---")
    print(f"Missing pairs default to prediction=FAILURE (no MR signal)\n")

    # Default assignment: all missing get prediction=FAILURE
    tp_full = tp
    tn_full = tn + f_missing
    fp_full = fp
    fn_full = fn + s_missing

    ba_full_default = compute_ba(tp_full, tn_full, fp_full, fn_full)
    print(f"DEFAULT (missing → FAILURE prediction):")
    print(f"  TP={tp_full}, TN={tn_full}, FP={fp_full}, FN={fn_full}")
    print(f"  BA = {ba_full_default:.4f}")
    print(f"  (sensitivity = {tp_full/(tp_full+fn_full):.4f}, "
          f"specificity = {tn_full/(tn_full+fp_full):.4f})\n")

    # Best case: all missing SUCCESS are actually causal (prediction=SUCCESS, correct)
    # and all missing FAILURE stay as TN
    tp_best = tp + s_missing
    tn_best = tn + f_missing
    fp_best = fp
    fn_best = fn

    ba_best = compute_ba(tp_best, tn_best, fp_best, fn_best)
    print(f"BEST CASE (all missing SUCCESS would have been causal):")
    print(f"  TP={tp_best}, TN={tn_best}, FP={fp_best}, FN={fn_best}")
    print(f"  BA = {ba_best:.4f}")
    print(f"  (sensitivity = {tp_best/(tp_best+fn_best):.4f}, "
          f"specificity = {tn_best/(tn_best+fp_best):.4f})\n")

    # Worst case: all missing FAILURE would have been causal (FP)
    # and all missing SUCCESS stay as FN
    tp_worst = tp
    tn_worst = tn
    fp_worst = fp + f_missing
    fn_worst = fn + s_missing

    ba_worst = compute_ba(tp_worst, tn_worst, fp_worst, fn_worst)
    print(f"WORST CASE (all missing FAILURE would have been causal):")
    print(f"  TP={tp_worst}, TN={tn_worst}, FP={fp_worst}, FN={fn_worst}")
    print(f"  BA = {ba_worst:.4f}")
    print(f"  (sensitivity = {tp_worst/(tp_worst+fn_worst):.4f}, "
          f"specificity = {tn_worst/(tn_worst+fp_worst):.4f})\n")

    # Realistic scenario: missing pairs have same causal rate as analyzed
    # Among analyzed: 16/138 = 11.6% predicted causal
    causal_rate = (tp + fp) / n_analyzed
    expected_causal_missing = round(n_missing * causal_rate)
    print(f"PROPORTIONAL SCENARIO (same causal rate {causal_rate:.1%} in missing):")
    print(f"  Expected ~{expected_causal_missing} causal predictions among {n_missing} missing")

    # Distribute proportionally between S and F
    s_causal_prop = round(expected_causal_missing * s_missing / n_missing)
    f_causal_prop = expected_causal_missing - s_causal_prop

    tp_prop = tp + s_causal_prop
    tn_prop = tn + (f_missing - f_causal_prop)
    fp_prop = fp + f_causal_prop
    fn_prop = fn + (s_missing - s_causal_prop)

    ba_prop = compute_ba(tp_prop, tn_prop, fp_prop, fn_prop)
    print(f"  TP={tp_prop}, TN={tn_prop}, FP={fp_prop}, FN={fn_prop}")
    print(f"  BA = {ba_prop:.4f}\n")

    # ---- Mechanism-stratified missingness ----
    print("--- Mechanism-Stratified Missingness ---")
    for mech in ["abundance_modulating", "activity_blocking", "mixed"]:
        m_analyzed = analyzed[analyzed["mechanism_class"] == mech]
        m_missing = missing[missing["mechanism_class"] == mech]
        s_a = (m_analyzed["outcome"] == "SUCCESS").sum()
        f_a = (m_analyzed["outcome"] == "FAILURE").sum()
        s_m = (m_missing["outcome"] == "SUCCESS").sum()
        f_m = (m_missing["outcome"] == "FAILURE").sum()
        rate_a = s_a / len(m_analyzed) if len(m_analyzed) > 0 else float("nan")
        rate_m = s_m / len(m_missing) if len(m_missing) > 0 else float("nan")
        print(f"  {mech}:")
        print(f"    Analyzed: {len(m_analyzed)} ({s_a}S/{f_a}F, {rate_a:.1%} SUCCESS)")
        print(f"    Missing:  {len(m_missing)} ({s_m}S/{f_m}F, {rate_m:.1%} SUCCESS)")

    # ---- Disease-level missingness table ----
    print("\n--- Per-Disease Missingness ---")
    print(f"{'Disease':<35} {'Analyzed':>8} {'Missing':>7} {'%Miss':>6} {'S_miss':>6} {'F_miss':>6}")
    diseases = sorted(evaluable["disease"].unique())
    for d in diseases:
        d_eval = evaluable[evaluable["disease"] == d]
        d_analyzed = analyzed[analyzed["disease"] == d]
        d_missing = missing[missing["disease"] == d]
        pct = len(d_missing) / len(d_eval) if len(d_eval) > 0 else 0
        s_m = (d_missing["outcome"] == "SUCCESS").sum()
        f_m = (d_missing["outcome"] == "FAILURE").sum()
        if len(d_missing) > 0:
            print(f"  {d:<33} {len(d_analyzed):>8} {len(d_missing):>7} {pct:>6.0%} {s_m:>6} {f_m:>6}")

    # ---- Save results ----
    results = {
        "n_evaluable": n_total,
        "n_analyzed": n_analyzed,
        "n_missing": n_missing,
        "base_rate_evaluable": s_total / n_total,
        "base_rate_analyzed": s_analyzed / n_analyzed,
        "base_rate_missing": s_missing / n_missing,
        "base_rate_gap": s_missing / n_missing - s_analyzed / n_analyzed,
        "chi2_yates": float(chi2),
        "chi2_p": float(p_chi2),
        "fisher_p": float(p_fisher),
        "current_ba": ba_current,
        "full_default_ba": ba_full_default,
        "full_best_ba": ba_best,
        "full_worst_ba": ba_worst,
        "full_proportional_ba": ba_prop,
        "missing_success": int(s_missing),
        "missing_failure": int(f_missing),
    }

    out_path = OUT_DIR / "missingness_bounds_v34.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
