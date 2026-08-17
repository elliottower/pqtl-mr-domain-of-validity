# Test of the mechanism difference

Run 2026-08-17 on `results/v5_1/classification_v5_1.csv` after two independent simulated
reviewers identified that the manuscript asserts a difference between mechanism strata
without testing it. Both named the Gelman–Stern error: comparing whether each stratum's
interval excludes 0.50 is not a test of whether the strata differ.

## Result

| quantity | value |
|---|---|
| abundance_modulating | n = 65, 44 genes, BA 0.6165 |
| activity_blocking | n = 75, 42 genes, BA 0.5136 |
| difference | +0.1028 |
| gene-clustered bootstrap, 90% CI | [−0.0198, +0.2262] — includes 0 |
| gene-clustered bootstrap, 95% CI | [−0.0424, +0.2504] — includes 0 |
| P(difference ≤ 0) | 0.085 |
| gene-clustered permutation, two-sided | p = 0.2144 |

20,000 iterations each. Genes are the resampling and permutation unit: pairs sharing a gene
share an exposure GWAS and an instrument, so pair-level resampling treats 143 correlated
observations as independent.

One-sided reading, which V3.4's pre-declared direction would license, gives 0.085 and
approximately 0.107. Neither clears 0.05.

## What follows

Supported: pooled BA 0.5783 with a 90% CI excluding 0.50; abundance_modulating BA 0.6165
with a 90% CI excluding 0.50; activity_blocking not reaching the pre-registered criterion.

Not supported: that the two strata differ. The mechanism-dependence claim is a direction
consistent with the data at n = 140 and not an established difference.

The direction is as predicted and the point estimate is substantial, so this is a question
of power rather than of sign. Both reviewers independently converged on the same remedy:
apply the mechanism annotation to a corpus large enough to test the interaction.
