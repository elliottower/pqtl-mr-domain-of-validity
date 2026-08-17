# V5.1 — outcome-GWAS corrections and a power floor

**Written:** 2026-08-16, before any correction is applied and before the corrected dataset
exists. No V5.1 quantity has been computed.
**Basis:** `results/OUTCOME_GWAS_TRACE.md`, committed at `2782ab2` prior to this document.
**Scope:** data corrections and one new sensitivity analysis. No hypothesis, classifier, or
outcome label changes.

## 1. What is corrected, and why it is a bug fix

Eight outcome-GWAS accessions were wrong: five pointed at an unrelated trait, three were
withdrawn from the OpenGWAS index after they were queried. A ninth, Juvenile idiopathic
arthritis, used a substitute for the source pre-specified in `DEVIATION_LOG.md`.

Every replacement below restores a source that the July deviation log already named. None
is a new choice of instrument or outcome, and none was selected by consulting a result.

| disease | from | to | consortium | pairs |
|---|---|---|---|---|
| Autism spectrum disorder | `ieu-b-87` (oral/pharyngeal cancer) | `ieu-a-1185` | iPSYCH-PGC | 1 |
| Chronic kidney disease | `ieu-b-4874` (bladder cancer) | `ieu-a-1102` | CKDGen, Pattaro 2015 | 6 |
| Systemic lupus erythematosus | `ieu-a-1073` (serum copper) | `ebi-a-GCST003156` | Bentham 2015 | 4 |
| Anorexia nervosa | `ieu-b-61` (withdrawn) | `ieu-a-1186` | PGC-ED, Duncan 2017 | 3 |
| Pancreatic cancer | `ieu-b-4866` (withdrawn) | `ieu-a-822` | PanScan1, Amundadottir 2009 | 10 |
| Juvenile idiopathic arthritis | `ebi-a-GCST90018873` (Sakaue 2021) | `ebi-a-GCST005528` | Hinks 2013 | 6 |
| Melanoma | `ieu-a-62` (waist circumference) | corrected in both dicts | — | 0 |
| Glioma | `ieu-b-4987` (withdrawn) | `finn-b-C3_BRAIN` | FinnGen | 3 |

Glioma is the single case with no match. GICC, the pre-specified source, is absent from
OpenGWAS; the only index entry matching "glioma" is `prot-a-1217`, a protein assay.
`finn-b-C3_BRAIN` is "malignant neoplasm of brain", a broader phenotype, and rests on the
deviation log's pre-declared FinnGen fallback. Its 3 pairs fall below the floor in §3 and
are therefore reported in the sensitivity stratum regardless.

### 1.1 Juvenile idiopathic arthritis

`code/expand_v5.py` holds a second `DISEASE_GWAS` dict which supplied JIA with
`ebi-a-GCST90018873` (Sakaue 2021, 216 cases / 409,001 controls). The deviation log
pre-specifies Hinks 2013. The substitution was never logged, and Sakaue 2021 is a Biobank
Japan, UK Biobank and FinnGen meta-analysis whose controls include UK Biobank, while five
of the six JIA pairs use UKB-PPP instruments — an unflagged breach of the log's own
sample-overlap policy.

Hinks 2013 is available as `ebi-a-GCST005528` (2,816 cases / 13,056 controls, European,
pre-dating UK Biobank). Restoring it removes the overlap breach rather than flagging it,
and raises JIA's effective sample size from 864 to 9,266. No retroactive overlap flag is
required because the corrected source contains no UK Biobank participants.

The two `DISEASE_GWAS` dicts are consolidated into one shared constant so the divergence
cannot recur.

## 2. PCSK9 | Hypercholesterolemia is dropped

The outcome GWAS is an LDL cholesterol GWAS and the Phase III endpoint for PCSK9
inhibitors is LDL lowering, so the pair tests a proposition against itself.
PCSK9 | Myocardial infarction is retained.

## 3. Effective-N floor — new, pre-specified here

For a single-SNP Wald ratio the MR z-statistic equals the SNP-outcome association
z-statistic, so detectable effect size is a property of the outcome GWAS alone. Effective
sample size is `4 / (1/ncase + 1/ncontrol)`.

**The floor is effective N = 2,000.** At that size a study has 80% power at α = 0.05
two-sided to detect a per-allele odds ratio of 1.10 at MAF 0.30. Established disease loci
sit largely between OR 1.05 and 1.15, so a study below the floor cannot detect a typical
effect and its nulls carry no information about the drug target.

The floor is derived from that power argument and fixed before the corrected dataset
exists. It is not chosen by reference to which diseases it includes or excludes.

**Application.** The floor is applied uniformly to every disease. The primary result is
reported on the full corrected set. A pre-specified sensitivity analysis is reported
alongside it, restricted to diseases at or above the floor. Neither is designated the
headline on the basis of which is more favourable.

Diseases below the floor after the §1 corrections:

| disease | effective N | pairs |
|---|---|---|
| Thyroid cancer | 1,036 | 6 |
| Glioma | 1,852 | 3 |

Nine pairs, 6% of the analysed set. Thyroid cancer was never previously questioned and
sits below glioma, which is why the floor is applied as a uniform rule rather than as a
per-disease judgment.

## 4. Disclosure — 40 pairs with unrecoverable outcome-GWAS provenance

Forty of 157 pairs (25%) take their MR estimate from the EpiGraphDB/INTERVAL MR-EvE
catalog. `v34_mr_catalog.csv` records `source`, `method` (Wald), `n_snps` (1) and `rsid`,
but no outcome-GWAS accession. Which outcome GWAS EpiGraphDB used cannot be determined
from any artifact in either repository.

All 40 reproduce their stored p-values exactly from the catalog file, so the values are
not in question. Only their provenance is. This is stated in Limitations and no attempt is
made to reconstruct it.

## 5. Order of operations

1. This document is frozen and committed.
2. The §1 corrections are applied to both dicts, which are then consolidated.
3. The corrected dataset is rebuilt.
4. The analysis runs **once**. Primary on the full set, sensitivity above the floor,
   sample-overlap sensitivity as already specified in `DEVIATION_LOG.md`.
5. Whatever comes out is reported.
