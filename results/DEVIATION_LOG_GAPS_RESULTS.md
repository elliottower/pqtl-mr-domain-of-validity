# Four pre-declared analyses, now run

Closes the three pre-specified secondary analyses (nos. 9, 10, 11 at
`protocol/DEVIATION_LOG.md:396-398`) and the cross-platform stop condition
(`:177-184`, `:408`), none of which had been reported.

Script: `scripts/deviation_log_gaps.py`. Output: `results/deviation_log_gaps.json`.
Run: `uv run --with pandas python scripts/deviation_log_gaps.py`.
Seed 20260731, 10,000 bootstrap resamples, 90% CI (matching the paper).

**Validation.** Scoring reproduces `results/v5/evaluation_v5.json` exactly
(tp=11, tn=101, fp=8, fn=37, BA=0.5779) and reproduces the manuscript's headline
**BA = 0.559** on V3.4. No claim below rests on a re-implementation of the
classifier.

---

## Stop condition: not testable with this data

The protocol declares that **any cross-platform beta sign flip halts the
analysis** — comparing SomaScan (EpiGraphDB/INTERVAL) against Olink (UKB-PPP)
for every protein instrumented in both.

**No protein is instrumented in both.** Max instrument sources per gene = 1
across all 211 adjudicated pairs. Separately, `sentinel_beta` is recorded for
135 UKB-PPP rows and **zero** EpiGraphDB rows.

So the check cannot be run, and it was never runnable — the design assumed an
overlap the dataset does not contain. This is a protocol defect, not a failed
test.

**Recommended wording:** state that the declared overlap test was not
evaluable because instrument sources are disjoint by protein, and that
cross-platform concordance therefore remains unaddressed. Do not report it as
passed.

---

## Analysis 9 — instrument source: no significant difference

| stratum | n | BA | 90% CI |
|---|---|---|---|
| EpiGraphDB | 59 | 0.640 | [0.549, 0.735] |
| UKB-PPP | 98 | 0.532 | [0.471, 0.598] |

V3.4: EpiGraphDB 0.612 [0.522, 0.709], UKB-PPP 0.519 [0.454, 0.592].

The intervals overlap substantially in both versions. EpiGraphDB scores higher,
but the data do not support a claim that instrument platform drives the
headline. **Report as a null.**

---

## Analysis 10 — sample overlap: headline is robust

Excluding the 5 pairs flagged `overlap_flagged`: BA **0.578 → 0.584**
[0.530, 0.640], n 157 → 152. V3.4: 0.559 → 0.565.

Two-sample overlap bias does not account for the result.

---

## Analysis 11 — oncology vs non-oncology: a real split ★

| stratum | n | BA | 90% CI |
|---|---|---|---|
| oncology | 55 | 0.489 | [0.466, **0.500**] |
| non-oncology | 102 | 0.595 | [**0.525**, 0.664] |

**The intervals do not overlap**, and it replicates on V3.4 (0.486
[0.462, 0.500] vs 0.581 [0.506, 0.659]).

Oncology pairs sit at chance — the upper CI bound is exactly 0.500. The
headline's signal comes from non-oncology pairs. This is the one stratification
of the three that changes the interpretation, and it is consistent with the
paper's mechanism argument: the 53 RTK-flagged oncology pairs are dominated by
activity-blocking targets, where cis-pQTL MR is expected to be uninformative.

**Recommended:** report this as a supplement table. It strengthens rather than
weakens the paper's thesis, but it must be stated — the headline BA is a blend
of an at-chance stratum and an above-chance one.

---

## Fixed: the shipped API now reproduces the paper

`pqtl_validity.core.score_pairs` did not reproduce any published number. Two
causes, the first more serious than it first appeared:

**1. It implemented a different decision rule.** The paper evaluates the
standard pQTL-MR screen: predict SUCCESS when MR is significant, whatever the
mechanism, applied to all pairs for the composite and within each stratum for
the mechanism comparison. `score_pairs` instead gated on mechanism, predicting
FAILURE for every activity-blocking target regardless of its p-value. That rule
appears nowhere in the manuscript. It disagreed with the stored predictions on
7 of 161 pairs.

**2. It silently shrank the denominator.** Mechanisms were re-derived through
`classify_gene()` rather than read from the curated `mechanism_class` column.
That lookup is missing 17 genes the tables classify, so 20 of 161 pairs were
dropped without warning. Where `classify_gene` does return a class it agrees
with the curated column on every row; the defect was coverage, not correctness.

**Now fixed.** `score_pairs` takes a `rule` argument defaulting to
`"mr_significance"` (the paper's screen); the gated variant is available as
`"mechanism_gated"` and is documented as not evaluated in the paper. Under the
gated rule an unclassifiable gene raises by default and names the offending
genes, rather than quietly reducing *n*; `on_unclassified="skip"` restores the
old behaviour but records what was dropped in the result. Pairs without an
adjudicated outcome (EXCLUDED, PENDING) are excluded and reported rather than
counted.

Verified end to end against the shipped data tables:

| | reproduced | published |
|---|---|---|
| V5 composite | n=157, TP 11, TN 101, FP 8, FN 37, BA 0.5779 | identical |
| V3.4 composite | BA 0.559 | identical |
| V3.4 activity arm | n=76, BA 0.511, sens 0.095 | identical |
| V3.4 abundance arm | n=59, BA 0.590 | identical |

Locked in by `tests/test_core.py`, which loads the CSVs and asserts the
confusion matrices. Suite: **26 passed**.
