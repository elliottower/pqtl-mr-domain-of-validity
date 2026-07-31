"""Close the four pre-declared analyses that were never reported.

Three are pre-specified secondary analyses 9, 10 and 11 in
protocol/DEVIATION_LOG.md:396-398; the fourth is the cross-platform sign-flip
check declared as a STOP CONDITION at DEVIATION_LOG.md:177-184 and :408.

All four run off data already in the repo. Nothing is downloaded.

  1. Cross-platform concordance (STOP CONDITION). For every protein with
     instruments in more than one source, compare sentinel betas. Flag sign
     flips and order-of-magnitude changes.
  2. Instrument-source stratification (analysis 9).
  3. Sample-overlap sensitivity (analysis 10): recompute the headline excluding
     pairs flagged sample_overlap=True.
  4. Oncology vs non-oncology stratification (analysis 11).

Scoring goes through the public API (pqtl_validity.core.score_pairs), which
reproduces results/v5/evaluation_v5.json exactly and the manuscript's headline
BA of 0.559 on V3.4.

Usage:
    uv run --with pandas python scripts/deviation_log_gaps.py
"""

import json
from pathlib import Path

import pandas as pd

from pqtl_validity.core import score_pairs

REPO = Path(__file__).resolve().parent.parent
# The analyses were declared against V3.4, but the manuscript now reports V5
# (SUBMISSION_STATUS.md: "V5 data (n=157), deCODE instruments"). Run both, and
# treat V5 as primary since that is what the paper claims.
DATASETS = {
    "v5_primary": REPO / "results" / "v5" / "classification_v5.csv",
    "v34_as_declared": REPO / "data" / "classification_v34.csv",
}
ADJUDICATED = REPO / "data" / "adjudicated_v34.csv"
OUT = REPO / "results" / "deviation_log_gaps.json"

SEED = 20260731
N_BOOTSTRAP = 10_000


def score(df: pd.DataFrame, label: str) -> dict:
    """Score a subset through the public API, which is the paper's screen."""
    if len(df) == 0:
        return {"stratum": label, "n_pairs": 0, "note": "empty stratum"}

    res = score_pairs(
        genes=df["gene"].tolist(),
        outcomes=df["outcome"].tolist(),
        mr_ps=df["mr_p"].tolist(),
        n_bootstrap=N_BOOTSTRAP,
        ci_level=0.90,
        seed=SEED,
    )
    return {
        "stratum": label,
        "n_pairs": int(len(df)),
        "n_scored": int(res.n),
        "tp": res.tp, "tn": res.tn, "fp": res.fp, "fn": res.fn,
        "sensitivity": round(float(res.sensitivity), 4),
        "specificity": round(float(res.specificity), 4),
        "balanced_accuracy": round(float(res.balanced_accuracy), 4),
        "ba_ci_low": round(float(res.ba_ci_low), 4),
        "ba_ci_high": round(float(res.ba_ci_high), 4),
        "ci_level": float(res.ci_level),
        "n_dropped": int(res.n_dropped),
        "dropped_no_outcome": list(res.dropped_no_outcome),
    }


def cross_platform_concordance(adj: pd.DataFrame) -> dict:
    """STOP CONDITION: any cross-platform beta sign flip in the overlap test."""
    have_beta = adj[adj["sentinel_beta"].notna() & adj["sentinel_rsid"].notna()]
    # A protein is testable only if it carries instruments from >1 source.
    per_gene = have_beta.groupby("gene")["instrument_source"].nunique()
    multi = sorted(per_gene[per_gene > 1].index.tolist())

    comparisons = []
    for gene in multi:
        rows = have_beta[have_beta["gene"] == gene]
        by_source = rows.groupby("instrument_source").first()
        sources = by_source.index.tolist()
        for a, b in [(sources[i], sources[j])
                     for i in range(len(sources)) for j in range(i + 1, len(sources))]:
            ba_, bb_ = float(by_source.loc[a, "sentinel_beta"]), float(by_source.loc[b, "sentinel_beta"])
            flip = (ba_ * bb_) < 0
            ratio = abs(bb_) / abs(ba_) if ba_ != 0 else float("inf")
            comparisons.append({
                "gene": gene,
                "source_a": a, "beta_a": ba_, "rsid_a": by_source.loc[a, "sentinel_rsid"],
                "source_b": b, "beta_b": bb_, "rsid_b": by_source.loc[b, "sentinel_rsid"],
                "sign_flip": bool(flip),
                "abs_ratio": round(float(ratio), 3),
                "order_of_magnitude_change": bool(ratio > 10 or ratio < 0.1),
            })

    n_flip = sum(c["sign_flip"] for c in comparisons)
    n_oom = sum(c["order_of_magnitude_change"] for c in comparisons)
    # Diagnose *why* the test is or is not supported: the declared check needs
    # (a) a gene present under >1 source and (b) a recorded beta under each.
    beta_by_source = {
        str(s): int(adj[adj["instrument_source"] == s]["sentinel_beta"].notna().sum())
        for s in adj["instrument_source"].dropna().unique()
    }
    return {
        "n_genes_with_multiple_sources": len(multi),
        "n_comparisons": len(comparisons),
        "n_sign_flips": n_flip,
        "n_order_of_magnitude_changes": n_oom,
        "stop_condition_triggered": bool(n_flip > 0),
        "instrument_source_counts": {k: int(v) for k, v in
                                     adj["instrument_source"].value_counts().items()},
        "sentinel_beta_present_by_source": beta_by_source,
        "max_sources_per_gene": int(adj.groupby("gene")["instrument_source"].nunique().max()),
        "comparisons": comparisons,
        "interpretation": (
            "STOP CONDITION TRIGGERED: at least one cross-platform sign flip."
            if n_flip > 0 else
            "No cross-platform sign flip among testable proteins."
            if comparisons else
            "Not testable: no protein carries instruments from more than one "
            "source, so the declared overlap test has no support in this dataset."
        ),
    }


def run_one(cls: pd.DataFrame) -> dict:
    """Analyses 9, 10 and 11 for a single classification table."""
    out: dict = {"n_rows_classification": int(len(cls))}

    out["headline_reproduction"] = score(cls, "all pairs (headline)")

    out["analysis_09_instrument_source"] = [
        score(cls[cls["instrument_source"] == s], f"instrument_source={s}")
        for s in sorted(cls["instrument_source"].dropna().unique())
    ]

    # sample_overlap is categorical ('clean' / 'overlap_flagged'), not boolean.
    overlap = cls["sample_overlap"].astype(str).str.strip() == "overlap_flagged"
    out["analysis_10_sample_overlap"] = {
        "n_flagged": int(overlap.sum()),
        "excluding_flagged": score(cls[~overlap], "sample_overlap excluded"),
        "flagged_only": score(cls[overlap], "sample_overlap only"),
    }

    # disease_area is already exactly 'oncology' / 'non_oncology'. Match exactly:
    # a substring test on "oncol" also matches "non_oncology".
    area = cls["disease_area"].astype(str).str.strip()
    out["analysis_11_oncology"] = {
        "n_oncology": int((area == "oncology").sum()),
        "oncology": score(cls[area == "oncology"], "oncology"),
        "non_oncology": score(cls[area == "non_oncology"], "non-oncology"),
        "disease_area_values": {k: int(v) for k, v in
                                cls["disease_area"].value_counts(dropna=False).items()},
    }
    return out


def show(label: str, res: dict) -> None:
    print(f"\n{'='*66}\n{label}  (n_rows={res['n_rows_classification']})\n{'='*66}")
    h = res["headline_reproduction"]
    print(f"headline BA = {h['balanced_accuracy']:.3f} "
          f"[{h['ba_ci_low']:.3f}, {h['ba_ci_high']:.3f}]  n={h['n_scored']}")
    print("\nanalysis 9 - instrument source:")
    for r in res["analysis_09_instrument_source"]:
        if r.get("n_scored"):
            print(f"  {r['stratum']:<34} BA={r['balanced_accuracy']:.3f} "
                  f"[{r['ba_ci_low']:.3f}, {r['ba_ci_high']:.3f}]  n={r['n_scored']}")
    print(f"\nanalysis 10 - sample overlap "
          f"({res['analysis_10_sample_overlap']['n_flagged']} flagged):")
    r = res["analysis_10_sample_overlap"]["excluding_flagged"]
    print(f"  {r['stratum']:<34} BA={r['balanced_accuracy']:.3f} "
          f"[{r['ba_ci_low']:.3f}, {r['ba_ci_high']:.3f}]  n={r['n_scored']}")
    print(f"\nanalysis 11 - oncology ({res['analysis_11_oncology']['n_oncology']} pairs):")
    for k in ("oncology", "non_oncology"):
        r = res["analysis_11_oncology"][k]
        if r.get("n_scored"):
            print(f"  {r['stratum']:<34} BA={r['balanced_accuracy']:.3f} "
                  f"[{r['ba_ci_low']:.3f}, {r['ba_ci_high']:.3f}]  n={r['n_scored']}")


def main() -> None:
    adj = pd.read_csv(ADJUDICATED)
    out: dict = {
        "seed": SEED,
        "n_bootstrap": N_BOOTSTRAP,
        "adjudicated_file": ADJUDICATED.name,
        "stop_condition_cross_platform": cross_platform_concordance(adj),
        "datasets": {},
    }
    for name, path in DATASETS.items():
        out["datasets"][name] = run_one(pd.read_csv(path))
        out["datasets"][name]["source_file"] = str(path.relative_to(REPO))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, default=str))

    sc = out["stop_condition_cross_platform"]
    print(f"STOP CONDITION (cross-platform sign flip):")
    print(f"  genes carrying >1 instrument source: {sc['n_genes_with_multiple_sources']}"
          f"  (max sources per gene = {sc['max_sources_per_gene']})")
    print(f"  sentinel_beta recorded by source: {sc['sentinel_beta_present_by_source']}")
    print(f"  comparisons possible: {sc['n_comparisons']}, sign flips: {sc['n_sign_flips']}")
    print(f"  {sc['interpretation']}")
    for name, res in out["datasets"].items():
        show(name, res)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
