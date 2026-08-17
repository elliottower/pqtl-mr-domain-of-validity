# cis-pQTL Mendelian randomization and Phase III drug-target outcomes

Analysis archive and pre-registration chain. Deposited 2026-08-17.

Repository: https://github.com/elliottower/pqtl-mr-domain-of-validity
ORCID: 0000-0001-7004-8884

## What this is

A pre-registered evaluation of whether cis-pQTL Mendelian randomization predicts Phase III
success for drug targets, across 143 drug-target-indication pairs spanning 91 genes and 24
diseases. Outcomes were adjudicated from FDA and EMA approvals and ClinicalTrials.gov. The
classifier is a single-SNP Wald ratio with `p < 0.05` predicting SUCCESS, with no tunable
parameters, fixed in advance.

`analysis_report.pdf` states the results and their limits in full. No manuscript is
included: the current draft asserts a mechanism-dependent difference that these data do not
support, and depositing it would archive a claim the accompanying analysis contradicts.

## Results

| stratum | n | balanced accuracy | 90% CI | pre-registered criterion (lower bound > 0.50) |
|---|---|---|---|---|
| Primary | 143 | 0.5783 | [0.5223, 0.6364] | met |
| Effective N ≥ 2,000 | 134 | 0.5810 | [0.5243, 0.6415] | met |
| Abundance-modulating | 65 | 0.6165 | [0.5309, 0.7040] | met |
| Activity-blocking | 75 | 0.5136 | [0.4554, 0.5826] | not met |

Sensitivity 0.227, specificity 0.929, against a base rate of 0.308 SUCCESS. Positive
predictive value 0.588, a 1.9-fold lift over the base rate, resting on 17 positive calls of
which 10 are correct, drawn from a small number of genes.

## What the analysis does not show

The registered prediction was that predictive performance differs between
abundance-modulating and activity-blocking targets. **It does not, at this sample size.**
The observed gap is +0.1028, with a gene-clustered permutation p of 0.2144 and a bootstrap
90% interval of [−0.0198, +0.2262]. `results/MECHANISM_DIFFERENCE_TEST.md` has the detail.

The gap also does not survive three checks. It falls to +0.0400 among pairs whose outcome
GWAS can be identified, to +0.0148 outside immune-mediated disease, and reverses to −0.0349
within UKB-PPP instruments. Two other pre-specified stratifications produce differences of
the same size and the same non-significance: disease area +0.0941 and instrument source
+0.0757.

Detecting a gap of this size at 80% power requires roughly 345 to 435 analyzed pairs. The
applicability funnel puts the entire universe of Phase III drug-target pairs with a plasma
cis-pQTL instrument at 211. The question is not resolvable with data that currently exist.

## Provenance and corrections

An audit of every outcome-GWAS accession found eight errors: five pointing at an unrelated
trait, three withdrawn from the OpenGWAS index. A separate `DISEASE_GWAS` mapping in
`code/expand_v5.py` supplied two further diseases and had never been audited; one used a
216-case substitute where the deviation log pre-specified a 2,816-case source, with an
unflagged sample-overlap breach. All are corrected in `code/disease_gwas.py`, each restoring
a source named in `protocol/DEVIATION_LOG.md`. `results/OUTCOME_GWAS_TRACE.md` is the full
trace.

**The corrections did not change the conclusion.** The mechanism gap moved from +0.0955
before correction to +0.1028 after, and to +0.1029 after LD-proxy recovery was applied
symmetrically. What changed the conclusion was testing it.

## Contents

```
analysis_report.pdf   the results and their limits
protocol/             every pre-registration version V1-V7, amendments, and the V5.1 and
                      V5.2 corrections, each with its sibling SHA-256 freeze record, plus
                      the deviation log
results/              the corrected analysis, its caches, the accession audit and trace,
                      and the test of the mechanism difference
results/prior/        the uncorrected analyzed set, for comparison
data/                 frozen candidates, adjudicated outcomes, and the EpiGraphDB catalog
code/, scripts/       the analysis and the correction harness
MANIFEST.sha256       a hash per file
```

`data/epigraphdb/` holds three files that were not previously public. Forty of the 143
pairs take their estimate from the EpiGraphDB MR-EvE catalog, and without these files those
estimates cannot be regenerated. The catalog records source, method, number of SNPs and
rsID, but **not which outcome GWAS was used**, so for those 40 pairs the outcome phenotype,
its case and control counts, and any sample overlap with the exposure cohort remain
unverifiable. They reproduce exactly from the deposited file; their provenance cannot be
reconstructed.

## Verifying

Protocol freeze records are plain SHA-256 over the file:

```
shasum -a 256 protocol/v5/CORRECTIONS_V5_1.md    # against CORRECTIONS_V5_1_sha256.txt
shasum -a 256 -c MANIFEST.sha256                 # the whole deposit
```

`data/frozen_candidates_v34.csv` hashes to
`c428c733aa66cbdb45ae4761424871f4dfa93e173f56ccd56f101b0dc33ebe62`, matching the value
recorded in the deviation log before outcomes were adjudicated. That hash is what
establishes the candidate set predates the outcome labels.

Protocol V1 and V3 do not verify against their recorded hashes; both name pre-rename
filenames and entered the repository already renamed, so those records describe files that
never existed under those names. Both are superseded.

## License

Data and documents CC-BY-4.0. Code MIT.
