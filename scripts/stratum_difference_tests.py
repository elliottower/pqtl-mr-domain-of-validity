"""Test whether balanced accuracy differs between strata.

The manuscript's central claim is a difference between mechanism classes, but no
test of that difference exists anywhere in the repo -- each stratum is reported
separately and the contrast is read off two overlapping intervals.

This applies the same test to two groupings at the same standard:

  mechanism_class    abundance_modulating vs activity_blocking  (the paper's claim)
  instrument_source  EpiGraphDB vs UKB-PPP                      (dismissed in S10 as
                                                                 overlapping intervals)

Applying it to both is the point. A gap can only be dismissed for overlapping
intervals if the same reasoning is allowed to touch the headline gap.

Two tests per grouping, no modelling assumptions beyond the classifier itself:

  paired bootstrap   resample pairs within each group, recompute the BA difference,
                     take the percentile interval of the difference distribution
  permutation        shuffle the group label across pairs, recompute the BA
                     difference, count how often |shuffled| >= |observed|

Usage:
    python scripts/stratum_difference_tests.py [--data PATH] [--iters N]
"""

import argparse
import csv
import pathlib
import random

REPO = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_DATA = REPO / "results" / "v5" / "classification_v5.csv"


def load(path):
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh) if r["outcome"] in ("SUCCESS", "FAILURE")]
    return rows


def balanced_accuracy(rows):
    tp = sum(1 for r in rows if r["mr_causal"] == "True" and r["outcome"] == "SUCCESS")
    fp = sum(1 for r in rows if r["mr_causal"] == "True" and r["outcome"] == "FAILURE")
    tn = sum(1 for r in rows if r["mr_causal"] == "False" and r["outcome"] == "FAILURE")
    fn = sum(1 for r in rows if r["mr_causal"] == "False" and r["outcome"] == "SUCCESS")
    if (tp + fn) == 0 or (tn + fp) == 0:
        return None
    return ((tp / (tp + fn)) + (tn / (tn + fp))) / 2


def confusion(rows):
    return (
        sum(1 for r in rows if r["mr_causal"] == "True" and r["outcome"] == "SUCCESS"),
        sum(1 for r in rows if r["mr_causal"] == "True" and r["outcome"] == "FAILURE"),
        sum(1 for r in rows if r["mr_causal"] == "False" and r["outcome"] == "FAILURE"),
        sum(1 for r in rows if r["mr_causal"] == "False" and r["outcome"] == "SUCCESS"),
    )


def paired_bootstrap(a, b, iters, seed):
    """Resample within each group; return the difference distribution."""
    rnd = random.Random(seed)
    diffs = []
    for _ in range(iters):
        ra = [a[rnd.randrange(len(a))] for _ in range(len(a))]
        rb = [b[rnd.randrange(len(b))] for _ in range(len(b))]
        ba_a, ba_b = balanced_accuracy(ra), balanced_accuracy(rb)
        if ba_a is not None and ba_b is not None:
            diffs.append(ba_a - ba_b)
    diffs.sort()
    return diffs


def permutation(a, b, observed, iters, seed):
    """Shuffle group membership; count |shuffled diff| >= |observed diff|."""
    rnd = random.Random(seed)
    pool = a + b
    n_a = len(a)
    hits = valid = 0
    for _ in range(iters):
        rnd.shuffle(pool)
        ba_a = balanced_accuracy(pool[:n_a])
        ba_b = balanced_accuracy(pool[n_a:])
        if ba_a is None or ba_b is None:
            continue
        valid += 1
        if abs(ba_a - ba_b) >= abs(observed):
            hits += 1
    return (hits + 1) / (valid + 1), valid


def report(label, name_a, rows_a, name_b, rows_b, iters, seed):
    ba_a, ba_b = balanced_accuracy(rows_a), balanced_accuracy(rows_b)
    if ba_a is None or ba_b is None:
        print(f"\n{label}: undefined (a stratum has no positives or no negatives)")
        return
    observed = ba_a - ba_b
    diffs = paired_bootstrap(rows_a, rows_b, iters, seed)
    lo, hi = diffs[int(0.05 * len(diffs))], diffs[int(0.95 * len(diffs))]
    lo95, hi95 = diffs[int(0.025 * len(diffs))], diffs[int(0.975 * len(diffs))]
    p_perm, n_valid = permutation(rows_a, rows_b, observed, iters, seed + 1)
    frac_le0 = sum(1 for d in diffs if d <= 0) / len(diffs)

    print(f"\n{label}")
    print(f"  {name_a:<24} n={len(rows_a):<4} BA={ba_a:.4f}   TP/FP/TN/FN = {confusion(rows_a)}")
    print(f"  {name_b:<24} n={len(rows_b):<4} BA={ba_b:.4f}   TP/FP/TN/FN = {confusion(rows_b)}")
    print(f"  observed difference      {observed:+.4f}")
    print(f"  paired bootstrap 90% CI  [{lo:+.4f}, {hi:+.4f}]   {'EXCLUDES 0' if lo > 0 or hi < 0 else 'includes 0'}")
    print(f"  paired bootstrap 95% CI  [{lo95:+.4f}, {hi95:+.4f}]   {'EXCLUDES 0' if lo95 > 0 or hi95 < 0 else 'includes 0'}")
    print(f"  bootstrap P(diff <= 0)   {frac_le0:.4f}")
    print(f"  permutation p (2-sided)  {p_perm:.4f}   ({n_valid} valid permutations)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=str(DEFAULT_DATA))
    ap.add_argument("--iters", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20260816)
    args = ap.parse_args()

    rows = load(args.data)
    print(f"data: {args.data}")
    print(f"n = {len(rows)} scored pairs; {args.iters} iterations")

    report(
        "MECHANISM CLASS  (the paper's central claim)",
        "abundance_modulating", [r for r in rows if r["mechanism_class"] == "abundance_modulating"],
        "activity_blocking", [r for r in rows if r["mechanism_class"] == "activity_blocking"],
        args.iters, args.seed,
    )
    report(
        "INSTRUMENT SOURCE  (dismissed in Table S10)",
        "EpiGraphDB", [r for r in rows if r["instrument_source"] == "EpiGraphDB"],
        "UKB-PPP", [r for r in rows if r["instrument_source"] == "UKB-PPP"],
        args.iters, args.seed,
    )
    report(
        "DISEASE AREA  (confound check)",
        "non_oncology", [r for r in rows if r["disease_area"] == "non_oncology"],
        "oncology", [r for r in rows if r["disease_area"] == "oncology"],
        args.iters, args.seed,
    )


if __name__ == "__main__":
    main()
