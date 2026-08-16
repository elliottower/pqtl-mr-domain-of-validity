"""Power for the V7 mechanism-interaction test. Design calculation only.

Uses published marginal structure (n=157, base success rate 0.306) and hypothetical
interaction sizes. Does not read outcomes, tier assignments, or any result file.
Its purpose is to choose 2-tier vs 3-tier before the corrected dataset exists.

Model: success ~ mr_evidence + tier + mr_evidence:tier, Firth-penalised,
one-sided, gene-clustered permutation.
"""
import math, random

N_TOTAL, BASE = 157, 0.306
MAIN = 0.45          # log-OR per unit -log10(p), assumed
SEED = 20260816


def firth_interaction(rows, iters=120, lr=0.5, pen=0.5):
    """Ridge-penalised logistic as a stable stand-in for Firth in simulation."""
    b = [0.0] * 4
    for _ in range(iters):
        g = [0.0] * 4
        for e, t, y in rows:
            x = (1.0, e, t, e * t)
            z = sum(bi * xi for bi, xi in zip(b, x))
            p = 1 / (1 + math.exp(-max(-30, min(30, z))))
            for i in range(4):
                g[i] += (y - p) * x[i]
        for i in range(4):
            b[i] = b[i] + lr * (g[i] / len(rows) - pen * b[i] / len(rows))
    return b[3]


def draw(n_per_tier, tiers, inter, rnd):
    rows = []
    for t in range(1, tiers + 1):
        for _ in range(n_per_tier):
            e = rnd.expovariate(1 / 1.1)               # -log10(p), skewed
            lo = math.log(BASE / (1 - BASE)) + e * (MAIN - inter * (t - 1))
            y = 1 if rnd.random() < 1 / (1 + math.exp(-lo)) else 0
            rows.append((e, t, y))
    return rows


def power(n_per_tier, tiers, inter, reps=300, nperm=40):
    rnd = random.Random(SEED)
    hits = 0
    for _ in range(reps):
        rows = draw(n_per_tier, tiers, inter, rnd)
        obs = firth_interaction(rows)
        null = []
        for _ in range(nperm):
            ts = [r[1] for r in rows]
            rnd.shuffle(ts)
            null.append(firth_interaction([(e, t, y) for (e, _, y), t in zip(rows, ts)]))
        # one-sided: hypothesis is a NEGATIVE interaction
        p = (sum(1 for v in null if v <= obs) + 1) / (len(null) + 1)
        if p < 0.05:
            hits += 1
    return hits / reps


if __name__ == "__main__":
    print("V7 interaction power, continuous evidence, one-sided, alpha=0.05")
    print("base success 0.306, main effect %.2f log-OR per unit -log10(p)\n" % MAIN)
    for tiers, npt in ((3, 52), (2, 78)):
        print("  %d tiers, %d pairs/tier (n=%d)" % (tiers, npt, tiers * npt))
        for inter in (0.10, 0.20, 0.30, 0.45):
            print("     interaction %.2f  ->  power %3.0f%%" % (inter, 100 * power(npt, tiers, inter)))
        print()
