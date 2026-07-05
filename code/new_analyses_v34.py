"""
V3.4 Analyses 13-14: Effect direction concordance + Leave-one-disease-out.

Runs on existing classification_v34.csv + adjudicated_v34.csv — no API calls.

Usage:
    cd ~/Documents/GitHub/transport-wrapper
    uv run python -m DRUGS_EXPANDED.new_analyses_v34
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

HERE = Path(__file__).parent
OUT_DIR = HERE / "results" / "v34"

REDUCING_ACTIONS = {"INHIBITOR", "ANTAGONIST", "ANTISENSE INHIBITOR", "DISRUPTING AGENT"}
ENHANCING_ACTIONS = {"AGONIST", "ACTIVATOR", "EXOGENOUS PROTEIN", "STABILISER"}


def get_expected_sign(action_types_str: str) -> str | None:
    if pd.isna(action_types_str):
        return None
    actions = {a.strip() for a in action_types_str.split(";")}
    has_reducing = bool(actions & REDUCING_ACTIONS)
    has_enhancing = bool(actions & ENHANCING_ACTIONS)
    if has_reducing and not has_enhancing:
        return "positive"
    if has_enhancing and not has_reducing:
        return "negative"
    return None


def compute_ba(preds, labs):
    tp = sum(1 for p, l in zip(preds, labs) if p == "SUCCESS" and l == "SUCCESS")
    tn = sum(1 for p, l in zip(preds, labs) if p == "FAILURE" and l == "FAILURE")
    fp = sum(1 for p, l in zip(preds, labs) if p == "SUCCESS" and l == "FAILURE")
    fn = sum(1 for p, l in zip(preds, labs) if p == "FAILURE" and l == "SUCCESS")
    sens = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
    spec = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
    return (sens + spec) / 2, {"tp": tp, "tn": tn, "fp": fp, "fn": fn}


def main():
    classif = pd.read_csv(OUT_DIR / "classification_v34.csv")
    adj = pd.read_csv(OUT_DIR / "adjudicated_v34.csv")

    df = classif.merge(adj[["pair_id", "action_types"]], on="pair_id", how="left")

    ba_full, m_full = compute_ba(df["prediction"].tolist(), df["outcome"].tolist())
    print(f"Current: n={len(df)}, BA={ba_full:.4f}")
    print(f"  TP={m_full['tp']} TN={m_full['tn']} FP={m_full['fp']} FN={m_full['fn']}\n")

    # =================================================================
    # Analysis 13: Effect direction concordance
    # =================================================================
    print("=" * 70)
    print("  ANALYSIS 13: Effect Direction Concordance")
    print("=" * 70)

    df["expected_sign"] = df["action_types"].map(get_expected_sign)
    df["beta_sign"] = df["mr_beta"].map(
        lambda b: "positive" if b > 0 else "negative" if b < 0 else None
    )
    df["direction_concordant"] = df["expected_sign"] == df["beta_sign"]

    has_direction = df["expected_sign"].notna() & df["mr_beta"].notna()
    dir_pairs = df[has_direction]
    n_conc = dir_pairs["direction_concordant"].sum()
    n_dir = len(dir_pairs)
    print(f"\n  All pairs with directional prediction: {n_conc}/{n_dir} concordant ({n_conc/n_dir:.1%})")

    # Among significant pairs
    sig_dir = dir_pairs[dir_pairs["mr_causal"] == True]
    n_sig_conc = sig_dir["direction_concordant"].sum()
    n_sig = len(sig_dir)
    print(f"  Significant pairs (mr_causal=True): {n_sig_conc}/{n_sig} concordant ({n_sig_conc/n_sig:.1%})")

    # TP vs FP breakdown
    sig_tp = sig_dir[sig_dir["outcome"] == "SUCCESS"]
    sig_fp = sig_dir[sig_dir["outcome"] == "FAILURE"]
    tp_conc = sig_tp["direction_concordant"].sum()
    fp_conc = sig_fp["direction_concordant"].sum()
    print(f"    TP (drug works, MR sig): {tp_conc}/{len(sig_tp)} concordant")
    print(f"    FP (drug fails, MR sig): {fp_conc}/{len(sig_fp)} concordant")

    # By mechanism
    for mech in ["abundance_modulating", "activity_blocking"]:
        m_dir = dir_pairs[dir_pairs["mechanism_class"] == mech]
        m_conc = m_dir["direction_concordant"].sum()
        m_sig = m_dir[m_dir["mr_causal"] == True]
        m_sig_conc = m_sig["direction_concordant"].sum()
        print(f"\n  {mech}:")
        print(f"    All: {m_conc}/{len(m_dir)} concordant ({m_conc/len(m_dir):.1%})")
        if len(m_sig) > 0:
            print(f"    Significant: {m_sig_conc}/{len(m_sig)} concordant ({m_sig_conc/len(m_sig):.1%})")

    # Non-significant pairs direction
    nonsig_dir = dir_pairs[dir_pairs["mr_causal"] != True]
    nonsig_conc = nonsig_dir["direction_concordant"].sum()
    print(f"\n  Non-significant pairs: {nonsig_conc}/{len(nonsig_dir)} concordant "
          f"({nonsig_conc/len(nonsig_dir):.1%}) — expect ~50% by chance")

    # Detail table
    print(f"\n  Detail (significant pairs with direction):")
    print(f"  {'pair_id':<42} {'beta':>8} {'expect':>8} {'conc':>5} {'out':>7} {'mech':<20} {'actions'}")
    for _, r in sig_dir.sort_values("mr_p").iterrows():
        conc = "Y" if r["direction_concordant"] else "N"
        actions = str(r["action_types"])[:35] if pd.notna(r["action_types"]) else ""
        mech = r["mechanism_class"][:18] if pd.notna(r["mechanism_class"]) else ""
        print(f"  {r['pair_id']:<42} {r['mr_beta']:>8.3f} {r['expected_sign']:>8} "
              f"{conc:>5} {r['outcome']:>7} {mech:<20} {actions}")

    # Binomial tests
    if n_sig > 0:
        binom_sig = binomtest(n_sig_conc, n_sig, 0.5)
        print(f"\n  Binomial test (significant concordance vs 50%): p={binom_sig.pvalue:.4f}")
    if n_dir > 0:
        binom_all = binomtest(n_conc, n_dir, 0.5)
        print(f"  Binomial test (all-pair concordance vs 50%):   p={binom_all.pvalue:.4f}")

    # =================================================================
    # Analysis 14: Leave-one-disease-out
    # =================================================================
    print(f"\n{'='*70}")
    print("  ANALYSIS 14: Leave-One-Disease-Out BA Stability")
    print("=" * 70)

    diseases = sorted(df["disease"].unique())
    lodo = {}
    print(f"\n  {'Disease dropped':<35} {'n_rem':>5} {'BA':>7} {'delta':>7}")
    print(f"  {'-'*35} {'-'*5} {'-'*7} {'-'*7}")
    for disease in diseases:
        left = df[df["disease"] != disease]
        if len(left) < 2:
            continue
        ba_left, _ = compute_ba(left["prediction"].tolist(), left["outcome"].tolist())
        delta = ba_left - ba_full
        lodo[disease] = {"ba": ba_left, "delta": delta, "n_removed": len(df[df["disease"] == disease])}
        marker = " ***" if abs(delta) > 0.02 else ""
        print(f"  {disease:<35} {len(left):>5} {ba_left:>7.4f} {delta:>+7.4f}{marker}")

    ba_vals = [v["ba"] for v in lodo.values()]
    print(f"\n  BA range: [{min(ba_vals):.4f}, {max(ba_vals):.4f}]")
    print(f"  BA std:   {np.std(ba_vals):.4f}")
    most = max(lodo.items(), key=lambda x: abs(x[1]["delta"]))
    print(f"  Most influential: {most[0]} (n={most[1]['n_removed']}, delta={most[1]['delta']:+.4f})")

    # Check if any single disease removal crosses 0.50
    below_half = [d for d, v in lodo.items() if v["ba"] < 0.50]
    if below_half:
        print(f"\n  WARNING: BA < 0.50 when dropping: {', '.join(below_half)}")
    else:
        print(f"\n  All leave-one-out BAs remain >= 0.50")

    # Save
    results = {
        "direction_concordance": {
            "n_directional": int(n_dir),
            "n_concordant": int(n_conc),
            "rate": float(n_conc / n_dir) if n_dir > 0 else None,
            "n_significant": int(n_sig),
            "n_significant_concordant": int(n_sig_conc),
            "rate_significant": float(n_sig_conc / n_sig) if n_sig > 0 else None,
            "tp_concordant": int(tp_conc),
            "tp_total": len(sig_tp),
            "fp_concordant": int(fp_conc),
            "fp_total": len(sig_fp),
        },
        "leave_one_disease_out": {
            "ba_range": [float(min(ba_vals)), float(max(ba_vals))],
            "ba_std": float(np.std(ba_vals)),
            "most_influential": most[0],
            "most_influential_delta": float(most[1]["delta"]),
            "crosses_050": below_half,
            "per_disease": {k: v for k, v in lodo.items()},
        },
    }

    out_path = OUT_DIR / "new_analyses_v34.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n  Saved: {out_path}")


if __name__ == "__main__":
    main()
