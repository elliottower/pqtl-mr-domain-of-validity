"""
V3.3 Classifier — evaluate MR-only predictions against outcome labels.

Takes the V3.3 MR catalog (cis-pQTL Wald ratios from EpiGraphDB, 33 diseases)
and the adjudicated outcome labels, joins them, and computes:
  1. MR-only balanced accuracy (primary)
  2. Per-disease breakdown
  3. Bootstrap 90% CI
  4. Comparison with V3.2 (n=30) results

Usage:
    cd ~/Documents/GitHub/transport-wrapper
    uv run python -m DRUGS_EXPANDED.classify_v33
"""
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

OUT_DIR = Path(__file__).parent / "results"


def compute_metrics(predictions: list[str], labels: list[str]) -> dict:
    tp = sum(p == "SUCCESS" and l == "SUCCESS" for p, l in zip(predictions, labels))
    tn = sum(p == "FAILURE" and l == "FAILURE" for p, l in zip(predictions, labels))
    fp = sum(p == "SUCCESS" and l == "FAILURE" for p, l in zip(predictions, labels))
    fn = sum(p == "FAILURE" and l == "SUCCESS" for p, l in zip(predictions, labels))
    sens = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
    spec = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
    ba = (sens + spec) / 2 if not (np.isnan(sens) or np.isnan(spec)) else float("nan")
    ppv = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
    npv = tn / (tn + fn) if (tn + fn) > 0 else float("nan")
    return {
        "n": len(predictions), "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "sensitivity": sens, "specificity": spec, "balanced_accuracy": ba,
        "ppv": ppv, "npv": npv,
    }


def bootstrap_ba(predictions: list[str], labels: list[str], n_boot: int = 10000) -> dict:
    rng = np.random.default_rng()
    bas = []
    n = len(predictions)
    preds_arr = np.array(predictions)
    labs_arr = np.array(labels)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        p_boot = preds_arr[idx]
        l_boot = labs_arr[idx]
        tp = ((p_boot == "SUCCESS") & (l_boot == "SUCCESS")).sum()
        fn = ((p_boot == "FAILURE") & (l_boot == "SUCCESS")).sum()
        tn = ((p_boot == "FAILURE") & (l_boot == "FAILURE")).sum()
        fp = ((p_boot == "SUCCESS") & (l_boot == "FAILURE")).sum()
        sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
        spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
        if not (np.isnan(sens) or np.isnan(spec)):
            bas.append((sens + spec) / 2)
    bas = np.array(bas)
    return {
        "ba_mean": float(np.mean(bas)),
        "ba_ci_lo": float(np.percentile(bas, 5)),
        "ba_ci_hi": float(np.percentile(bas, 95)),
        "n_boot": n_boot,
    }


def main():
    OUT_DIR.mkdir(exist_ok=True)
    ts = datetime.now().isoformat()
    print(f"[{ts}] V3.3 Classifier — MR-only evaluation (n=50 target)")

    # Load MR catalog and outcome labels
    catalog = pd.read_csv(OUT_DIR / "v33_mr_catalog.csv")
    labels = pd.read_csv(Path(__file__).parent / "OUTCOME_LABELS_V33.csv")

    # Filter to evaluable (non-excluded, SUCCESS or FAILURE)
    evaluable_labels = labels[
        (~labels["excluded"].astype(bool)) &
        (labels["outcome"].isin(["SUCCESS", "FAILURE"]))
    ].copy()
    print(f"  Evaluable labels: {len(evaluable_labels)}")
    print(f"    SUCCESS: {(evaluable_labels['outcome'] == 'SUCCESS').sum()}")
    print(f"    FAILURE: {(evaluable_labels['outcome'] == 'FAILURE').sum()}")

    # Join with MR catalog
    merged = evaluable_labels.merge(
        catalog[["pair_id", "beta", "se", "p", "ci_lower", "ci_upper", "mr_supports_causal"]],
        on="pair_id",
        how="inner",
    )
    print(f"  Joined pairs (label + MR): {len(merged)}")

    if len(merged) < len(evaluable_labels):
        missing = set(evaluable_labels["pair_id"]) - set(merged["pair_id"])
        print(f"  WARNING: {len(missing)} labels without MR match: {missing}")

    # MR-only classification: CAUSAL → SUCCESS, NULL → FAILURE
    merged["prediction"] = merged["mr_supports_causal"].map({True: "SUCCESS", False: "FAILURE"})
    merged["correct"] = merged["prediction"] == merged["outcome"]

    # Print each pair
    print(f"\n  {'pair_id':<40} {'MR':>6} {'pred':>8} {'label':>8} {'ok':>3}")
    print(f"  {'-'*40} {'-'*6} {'-'*8} {'-'*8} {'-'*3}")
    for _, row in merged.sort_values("p").iterrows():
        tag = "CAUSAL" if row["mr_supports_causal"] else "NULL"
        ok = "Y" if row["correct"] else "N"
        print(f"  {row['pair_id']:<40} {tag:>6} {row['prediction']:>8} {row['outcome']:>8} {ok:>3}")

    # Compute metrics
    preds = merged["prediction"].tolist()
    labs = merged["outcome"].tolist()
    metrics = compute_metrics(preds, labs)
    boot = bootstrap_ba(preds, labs)

    print(f"\n  === MR-only Results (n={metrics['n']}) ===")
    print(f"  TP={metrics['tp']}  TN={metrics['tn']}  FP={metrics['fp']}  FN={metrics['fn']}")
    print(f"  Sensitivity: {metrics['sensitivity']:.3f}")
    print(f"  Specificity: {metrics['specificity']:.3f}")
    print(f"  Balanced Accuracy: {metrics['balanced_accuracy']:.3f}")
    print(f"  Bootstrap 90% CI: [{boot['ba_ci_lo']:.3f}, {boot['ba_ci_hi']:.3f}]")
    print(f"  PPV: {metrics['ppv']:.3f}")
    print(f"  NPV: {metrics['npv']:.3f}")

    # Naive baselines
    base_rate = (merged["outcome"] == "SUCCESS").mean()
    print(f"\n  Base rate (SUCCESS): {base_rate:.3f}")
    print(f"  Always-SUCCESS BA: 0.500")
    print(f"  Always-FAILURE BA: 0.500")

    # Per-disease
    print(f"\n  Per-disease:")
    per_disease = {}
    for disease in sorted(merged["canonical_disease"].unique()):
        d = merged[merged["canonical_disease"] == disease]
        d_preds = d["prediction"].tolist()
        d_labs = d["outcome"].tolist()
        d_metrics = compute_metrics(d_preds, d_labs)
        per_disease[disease] = d_metrics
        print(f"    {disease:<35} n={d_metrics['n']:>2}  BA={d_metrics['balanced_accuracy']:.2f}  "
              f"TP={d_metrics['tp']} TN={d_metrics['tn']} FP={d_metrics['fp']} FN={d_metrics['fn']}")

    # Save results
    n_excluded = labels["excluded"].astype(bool).sum()
    n_pending = labels["exclusion_reason"].str.contains("PENDING", na=False).sum()
    n_mismatch = labels["exclusion_reason"].str.contains("TARGET_MISMATCH", na=False).sum()
    results = {
        "timestamp": ts,
        "version": "V3.3",
        "analysis": "MR-only (cis-pQTL Wald, EpiGraphDB catalog, 33 diseases)",
        "n_total_pairs": len(labels),
        "n_evaluable": metrics["n"],
        "n_success": int((merged["outcome"] == "SUCCESS").sum()),
        "n_failure": int((merged["outcome"] == "FAILURE").sum()),
        "n_excluded": int(n_excluded),
        "exclusion_reasons": {
            "PENDING": int(n_pending),
            "TARGET_MISMATCH": int(n_mismatch),
        },
        "metrics": metrics,
        "bootstrap": boot,
        "base_rate": float(base_rate),
        "per_disease": per_disease,
    }

    eval_path = OUT_DIR / "evaluation_v33.json"
    with open(eval_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  Wrote evaluation to {eval_path}")

    # Save per-pair results
    merged_out = merged[["pair_id", "gene_symbol", "canonical_disease", "outcome",
                          "beta", "se", "p", "mr_supports_causal", "prediction", "correct",
                          "phase3_drugs", "notes"]].copy()
    class_path = OUT_DIR / "classification_v33.csv"
    merged_out.to_csv(class_path, index=False)
    print(f"  Wrote classification to {class_path}")

    ts2 = datetime.now().isoformat()
    print(f"\n[{ts2}] Done")


if __name__ == "__main__":
    main()
