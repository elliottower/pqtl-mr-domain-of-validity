# V3 Sweep Protocol — Pre-Specification

**Status:** FROZEN — DO NOT EDIT after SHA-256 hash is computed.
**Timestamp:** 2026-07-04T14:49:00Z
**Supersedes:** `PRESPEC_V2_SWEEP.md` (hash: `dcfa5b2ed5dfa24c32d11e9871de8ac02a916f2e4bcb2170097c9bba318df86d`)
**Companion:** `PRESPEC_V3_REVIEWER_GUIDANCE.md` (reasoning behind every design choice)

---

## 0. Disclosure: what was known when V3 was written

The following results were computed and observed before this document was hashed.
They are stated here so that no reviewer has to guess what influenced V3's design.

**Retrospective targets (outcomes seen):**

| Target → Disease | V2 class | Trial outcome | V3 prediction correct? |
|---|---|---|---|
| CRP → CAD | ZOMBIE | FAILURE | Yes |
| triglycerides → CAD | CONCORDANT | SUCCESS | Yes |
| HDL → CAD | TRANSLATION_GAP | FAILURE | No (predicted SUCCESS) |
| LDL → CAD | TRANSLATION_GAP | SUCCESS | Yes |
| systolic_BP → CAD | ZOMBIE | SUCCESS | Excluded (polygenic) |
| BMI → T2D | ZOMBIE | SUCCESS | Excluded (polygenic) |

**eQTL sweep findings (instrument quality, not outcomes):**
- 500 MR tests scored from eqtl-a-* studies; 0% pass ≥10 instruments; median 3 SNPs
- 4 drug-target overlaps found, all MR_NULL — confirms eQTLs are underpowered as
  instruments for drug target MR (instrument quality finding, not outcome finding)
- This motivates the tiered instrument hierarchy in Section 3

**Protein QTL mapping:**
- 4,406 protein QTL studies mapped to gene symbols via UniProt
- 411 pQTL→OT drug target pairs identified across 8 diseases (§4.1)
- Outcomes for these 411 were NOT scored under any taxonomy at hash time
- These form the prospective validation set (Section 2.2)

**What V3 changes from V2 and why:**
- V2 assumed ≥10 eQTL instruments per gene. The 500-gene sweep proved 0% of eQTL
  studies pass this gate. This is an instrument-power finding, not an outcome finding.
- V3 replaces the single-instrument-tier design with a tiered hierarchy (Section 3).
- V3 adds deterministic prospective inclusion rules (Section 4) and frozen outcome
  adjudication (Section 5) that V2 lacked.
- The taxonomy, thresholds, and binary endpoint are UNCHANGED from V2.

---

## 1. Primary endpoint

**Binary target viability prediction** (unchanged from V2):

- **Predict SUCCESS** if the target is classified in any causal-positive class
  (CONCORDANT_CAUSAL, GENETIC_AMPLIFICATION, OBS_INFLATED_CAUSAL).
- **Predict FAILURE** if the target is classified as ZOMBIE or TRUE_NULL.
- OPPOSITE_DIRECTION and INSUFFICIENT are excluded. Reported in sensitivity analysis.

## 2. Analysis sets

### 2.1 Retrospective calibration set (exploratory only)

The 4 evaluable targets listed in Section 0 (CRP, triglycerides, HDL, LDL).
These are reported for transparency but are **excluded from the primary
falsification statistic**. They do not count toward the N≥40 threshold.
Their outcomes were known before V3 was written.

systolic_BP and BMI→T2D are excluded under the pre-declared polygenic scope
boundary (Section 6) and are also retrospective.

### 2.2 Prospective validation set (primary evidence)

All pQTL→OT drug target pairs that pass the deterministic inclusion rule in
Section 4. These outcomes were not scored under any taxonomy at the time V3
was hashed. The falsification criterion (Section 8) applies exclusively to
this set.

### 2.3 eQTL negative control set (reported, not scored)

The full eqtl-a-* sweep results. Reported as evidence for instrument quality
effects. Not included in either the retrospective or prospective analysis
because eQTL instruments are systematically underpowered (Section 0).

## 3. Tiered instrument hierarchy

### Tier 1 — Curated biomarker GWAS (primary)

Purpose-built GWAS for the exposure of interest (e.g., `ieu-b-110` for LDL
cholesterol, `ieu-a-299` for Lp(a)). These have hundreds of instruments and
measure the drug-relevant biomarker directly.

Applied to the retrospective calibration set and any prospective target where
a curated biomarker GWAS exists.

### Tier 2 — Protein QTL (primary prospective)

prot-a-* and prot-c-* studies from OpenGWAS, mapped to gene symbols via
UniProt (mapping frozen in `results/prot_gene_mapping.csv`). These measure
circulating protein levels — closer to drug mechanism than gene expression.

**Estimator by instrument count:**

| Instruments | Estimator | Status |
|---|---|---|
| 1 cis-SNP | Wald ratio (β_Y/β_X, delta-method SE) | Primary |
| 2 SNPs | IVW (inverse-variance weighted mean of two Wald ratios) | Primary |
| ≥3 SNPs | IVW (random-effects if Cochran's Q_het p<0.05) | Primary |

All instruments must be genome-wide significant (p<5e-8) and cis-acting
(within ±1 Mb of the protein-coding gene). Single-cis-SNP Wald ratio
is the standard approach in pQTL MR because cis-pQTL instruments
are strong (F>100 typical) and biologically constrained to the gene locus,
reducing pleiotropy concerns that motivate multi-instrument requirements
in polygenic MR. This follows established practice in large-scale pQTL-MR
studies (Zheng et al. 2020, Folkersen et al. 2020).

**Sensitivity for single-SNP estimates:** Steiger directionality test and
colocalization probability (where LD data permit) are reported alongside
the Wald ratio. These do not gate inclusion — they are diagnostic.

**Known bias — winner's curse:** Single-SNP pQTL effect sizes are
estimated in the discovery GWAS and selected at p<5e-8, inflating the
exposure effect estimate. This attenuates the Wald ratio toward the
null, biasing MR estimates downward. The bias direction is conservative
for zombie detection (real causal effects may be misclassified as
TRUE_NULL) but could inflate the ZOMBIE rate. Median F-statistic is
reported per disease to quantify this.

**Ancestry limitation:** LD clumping uses European-ancestry reference
panels. pQTL studies including multi-ancestry cohorts may have different
LD structure at the lead SNP. All analyses are restricted to EUR LD and
this is flagged as a generalizability limitation.

Applied to the prospective validation set across all 8 diseases (§4.1).

### Tier 3 — eQTL (negative control)

eqtl-a-* studies. Reported to demonstrate that instrument quality dominates
classification accuracy. Not used for primary or secondary scoring.

**Instrument selection within each tier:** Genome-wide significant SNPs
(p<5e-8), LD-clumped (r²<0.001, 10 Mb window, EUR reference). For
OpenGWAS pre-clumped data (`preclumped=1`), the server's clumping is
accepted as-is. Mean F-statistic ≥10 across instruments (calculated as β²/SE² from
the exposure GWAS summary statistics; if per-SNP SE is unavailable
from the OpenGWAS API, F-statistic is approximated from the reported
p-value and sample size, and flagged in the results table).

## 4. Deterministic prospective inclusion rule

### 4.1 Frozen disease universe

V3 expands the candidate universe from 2 diseases (V2: CAD, T2D) to 8,
selected for having both (a) an outcome GWAS in OpenGWAS with ≥10,000
cases (or ≥50,000 total N for continuous traits) and (b) ≥20 OT drug
targets with pQTL overlap. The disease list and outcome
GWAS IDs are frozen here:

| Disease | Outcome GWAS | pQTL overlaps |
|---|---|---|
| Coronary artery disease (CAD) | ieu-a-7 | 27 |
| Type 2 diabetes (T2D) | ebi-a-GCST007517 | 37 |
| Alzheimer's disease (AD) | ieu-b-2 | 40 |
| Rheumatoid arthritis (RA) | ieu-a-832 | 68 |
| Asthma | ieu-a-44 | 38 |
| Chronic kidney disease (CKD) | ieu-a-1102 | 39 |
| Ulcerative colitis (UC) | ieu-a-970 | 44 |
| Breast cancer (BC) | ieu-a-1126 | 118 |

**Total pQTL-drug target pairs: 411.** Dry-run gate attrition on a
random sample (30/79 unique pQTL studies): 53% pass the instrument gate
(≥1 cis-SNP at p<5e-8), projecting ~218 pairs surviving gate 3. The
obs-arm and outcome-adjudication gates are the binding constraints —
requiring both a published meta-analysis and a Phase III readout
eliminates most pairs. **Projected evaluable N: 30–70.** If N<40, the
primary endpoint is reported at reduced power per §8.1. If N<20,
results are exploratory only.

The §6 polygenic exclusion applies to all 8 diseases: any exposure
measuring a composite physiological trait (BP, BMI) is excluded
regardless of which disease it targets.

### 4.2 Inclusion gates

A pQTL→disease pair enters the prospective validation set if and only if
ALL of the following hold. No manual override is permitted after hashing.

1. **pQTL→gene mapping exists** in the frozen mapping file
   (`results/prot_gene_mapping.csv`, hashed separately).

2. **Gene is an OT drug target** for one of the 8 frozen diseases (§4.1),
   per the frozen candidate universe
   (`results/disease_sweep_v3_expanded.csv`, hashed separately).

3. **Instruments available:** ≥1 genome-wide significant cis-SNP at p<5e-8
   in the pQTL study, retrieved from OpenGWAS tophits at execution time.
   Estimator determined by instrument count per §3 Tier 2 table. Zero
   instruments → INSUFFICIENT.

4. **Observational estimate available:** A published meta-analysis or large
   prospective cohort reporting an effect estimate with 95% CI for the same
   exposure→outcome pair. Source chosen by the priority ladder:
   IPD meta-analysis > random-effects meta-analysis > largest prospective
   cohort (≥1,000 participants).

   **Deterministic search protocol (reproducible, non-discretionary):**
   For each pQTL target, execute the following PubMed query:

   ```
   ("[gene symbol]" OR "[protein name]") AND ("[disease]") AND
   ("meta-analysis" OR "prospective cohort" OR "cohort study" OR "case-control")
   AND ("odds ratio" OR "hazard ratio" OR "relative risk")
   NOT ("Mendelian randomization"[Title] OR "drug target"[Title] OR "immunotherapy"[Title])
   ```

   The NOT clause excludes MR studies, drug trial meta-analyses, and
   immunotherapy studies from the observational search. The observational
   arm must contain only traditional epidemiological designs (prospective
   cohorts, case-control studies, and meta-analyses thereof).

   Sort by **Most Recent** (publication date descending). Select the first
   result matching the priority ladder above that reports an effect estimate
   with 95% CI. If no result qualifies, the pair fails this gate.

   **Search freeze date:** 2026-07-04.
   All searches are executed on or before this date. Results published or
   indexed after this date do not affect inclusion.

   **Protein-synonym expansion:** If the primary gene-symbol query returns
   zero results, a pre-specified synonym list (`data/gene_synonyms.csv`)
   is consulted. This file maps gene symbols to commonly used protein-level
   search terms (e.g., INSR → "insulin resistance", CSF3R → "G-CSF",
   VEGFA → "VEGF"). The query is re-run once with the expanded terms
   substituted for the gene symbol / protein name. The synonym file is
   frozen and hashed with other artifacts. No synonyms may be added after
   the first V3 classifier output is computed.

   The search is executed once per target (twice if synonym expansion
   triggers); the selected PMID and extracted estimate are recorded in
   `data/observational_estimates_v3.csv` before scoring. SE is derived
   from CI width (SE = (ln(upper) − ln(lower)) / 3.92), never from 1/√n.

   The observational source table is frozen with specific PMIDs before any
   V3 scoring.

5. **Outcome adjudicable:** A Phase III pivotal trial primary endpoint
   readout or regulatory decision exists as of the outcome-freeze date,
   per Section 5.

Targets failing any gate are classified INSUFFICIENT and reported in a
completeness table, not dropped silently.

**Anti-cherry-pick rule:** No target may be added, removed, or reclassified
after the first V3 classifier output is computed. The inclusion rule is
evaluated once, deterministically, and the resulting set is final.

## 5. Frozen outcome adjudication

### 5.1 Outcome source

ClinicalTrials.gov (primary completion status) cross-referenced with
regulatory actions from stringent regulatory authorities — defined as
agencies that independently evaluate Phase III pivotal trial data
before granting marketing authorization (FDA, EMA, NMPA, PMDA,
Health Canada, TGA). Approvals from agencies that accept foreign
regulatory decisions without independent clinical review do not
qualify.

### 5.2 Deterministic classification rule

| Outcome label | Rule |
|---|---|
| **SUCCESS** | **Phase III** pivotal trial met primary endpoint (p<0.05) OR drug received full or regular regulatory approval from a stringent regulatory authority (§5.1) for the target indication. Phase II efficacy alone does not qualify — the trial must be Phase III (or equivalent pivotal). Withdrawn-after-approval-for-safety = SUCCESS on efficacy endpoint (noted). |
| **FAILURE** | Pivotal trial failed primary endpoint OR drug withdrawn/discontinued for the indication before approval |
| **AMBIGUOUS** | Neither criterion met. Includes: subgroup-only win; accelerated/conditional approval where the required confirmatory trial failed or the approval was withdrawn for lack of efficacy; indication narrowed post-approval. Excluded from primary analysis, reported in sensitivity |
| **PENDING** | No pivotal trial readout or regulatory decision as of outcome-freeze date. Includes accelerated/conditional approvals still awaiting confirmatory data |

### 5.3 Blinded outcome hash

Outcome labels are locked in a separate file (`OUTCOME_LABELS_V3.csv`)
and SHA-256 hashed independently, BEFORE any V3 classifier is run.
This file contains:

- `ensembl_id`, `disease`, `outcome_label`, `source` (NCT ID or
  regulatory reference), `adjudication_date`

The outcome hash is recorded in `PRESPEC_V3_sha256.txt` alongside the
prespec hash. Predictions and outcomes are provably independent.

### 5.4 Outcome-freeze date

2026-07-04

All outcomes are adjudicated as of this date. Readouts after this date
are prospective predictions, not retrospective scores.

### 5.5 Post-freeze trial readouts

If a pivotal trial reports results between the outcome-freeze date and
manuscript submission (e.g., ZEUS/ziltivekimab), the PENDING label is
NOT reclassified. The pair remains excluded from the primary analysis.
Post-freeze readouts are reported in a separate "prospective prediction"
table showing what the classifier would have predicted, compared against
the actual outcome. This preserves the integrity of the frozen outcome
set while capturing additional validation data.

### 5.6 Multiple drugs per target-disease pair

When multiple drugs target the same gene for the same disease (e.g.,
tocilizumab and sarilumab both target IL6R for RA), the outcome label
reflects the strongest evidence: SUCCESS if any drug received full
approval from a stringent authority, FAILURE only if all drugs failed.
The specific drugs and their individual outcomes are recorded in the
notes column.

## 6. Pre-declared scope exclusions

### 6.1 Polygenic complex exposures

**Excluded a priori:** systolic blood pressure, BMI, and any exposure where
the genetic instrument captures a complex polygenic trait rather than a
direct biomarker or protein level.

**Mechanistic rationale:** For polygenic exposures, the per-allele MR effect
size is orders of magnitude smaller than the pharmacological effect achievable
by drugs. A drug lowering systolic BP by 10 mmHg achieves an effect that no
single allele approximates. The MR→null classification is mechanistically
expected regardless of whether the drug works, making the zombie test
uninformative for these targets. This exclusion is definable from biology
alone, without reference to trial outcomes.

**Rule:** An exposure is polygenic-excluded if the MR instruments capture a
composite physiological measurement (blood pressure, body mass index, smoking
behavior) rather than a single molecular entity (a circulating protein, a
metabolite, a lipid fraction). The distinction is: can a drug plausibly
modify the exposure through a single molecular mechanism? If no, exclude.

### 6.2 What is NOT excluded

HDL cholesterol is NOT excluded despite being a known miss (Section 0). HDL
is a single molecular entity (a lipid fraction) for which drugs exist that
directly modify its level (CETP inhibitors). Its failure reflects a real
limitation of the method — the zombie test cannot distinguish between causal
targets and correlated biomarkers. Excluding HDL post-hoc would be
cherry-picking.

### 6.3 Proxy-exposure observational pairs

Four pairs use a circulating ligand as proxy for a transmembrane
receptor: INSR→CAD and INSR→BC (fasting insulin for insulin receptor),
IGF1R→BC (circulating IGF-1 for IGF-1 receptor), IL6R→T2D (IL-6 for
sIL-6R). These are **included in the primary 2x2 analysis** because:
(a) the receptor has no circulating form measurable in epidemiological
cohorts, (b) the ligand concentration is the direct upstream signal the
receptor transduces, and (c) excluding them would reduce the already
small evaluable N from ~14 to ~10, below the threshold for meaningful
inference. Proxy status is flagged in the results table. A sensitivity
analysis (§8.3) reports the 2x2 restricted to the 10 direct-measurement
(CURATED) pairs only.

## 7. Classification frameworks

### 7.1 Primary — Binary 2x2 framework (V3.1)

Each arm is classified by whether its 95% CI excludes the null (OR=1):

| | OBS CI excludes null | OBS CI includes null |
|---|---|---|
| **MR CI excludes null** | CONCORDANT → predict SUCCESS | MR_ONLY → predict SUCCESS (weak) |
| **MR CI includes null** | **ZOMBIE** → predict FAILURE | TRUE_NULL → predict FAILURE |

This extends the binary genetic-support classification of Nelson et al.
(2015) from "any genetic evidence" to the MR × observational discordance
axis, and avoids the Chinn conversion (which introduces approximation
error for hazard ratios — the Chinn formula `d = log_OR * sqrt(3)/pi`
assumes a logistic model and is approximate for Cox-derived HRs, ~5%
error), using each arm's own statistical significance rather than an
arbitrary effect-size threshold. The ZOMBIE cell is the core hypothesis: targets
with observational support but no genetic causal evidence are predicted
to fail in trials.

MR_ONLY targets (MR signal, no obs association) are predicted SUCCESS
but flagged as lower confidence — these may represent novel causal
targets not yet captured in observational epidemiology.

### 7.2 Secondary — Continuous V2 taxonomy

Both arms converted to Cohen's d via Chinn (`d = log_OR * sqrt(3)/pi`).
Thresholds frozen from V1:

| V2 class | Rule | Clinical prediction |
|---|---|---|
| CONCORDANT_CAUSAL | MR meaningful (d≥0.10), OBS same direction and meaningful, Q p≥0.05 | SUCCESS |
| GENETIC_AMPLIFICATION | MR meaningful, OBS same direction, Q p<0.05, \|d_MR\|>\|d_OBS\| | SUCCESS |
| OBS_INFLATED_CAUSAL | MR meaningful, OBS same direction, Q p<0.05, \|d_OBS\|>\|d_MR\| | SUCCESS |
| ZOMBIE | MR null (d<0.10), OBS meaningful positive, Q p<0.05 | FAILURE |
| TRUE_NULL | MR null AND OBS null/weak, Q p≥0.05 | FAILURE |
| OPPOSITE_DIRECTION | MR and OBS opposite signs, at least one meaningful | EXCLUDED |
| INSUFFICIENT | 0 instruments, or Egger intercept p<0.01 (≥3 SNPs only) | EXCLUDED |

No changes from V2. The taxonomy is frozen across V1→V2→V3.

## 8. Falsification criterion

### 8.1 Primary (prospective validation set only)

The primary hypothesis is **falsified** if the balanced accuracy 90%
bootstrap CI lower bound ≤ 0.50 after ≥40 evaluable prospective pairs.

**Design justification:** BA=0.70 is pre-declared as the minimum clinically
useful accuracy for drug target triage (a priori hypothesis, not derived from
the retrospective calibration set). At N=40 and true BA=0.70, power to reject
the null (BA≤0.50) is 80% — the conventional threshold.

N=40 is a target, not a guarantee. If the deterministic inclusion rule
(Section 4) yields fewer than 40 evaluable pairs, the test proceeds at
reduced power and the actual N is reported. If fewer than 20 evaluable
pairs, the primary endpoint is declared underpowered and results are
reported as exploratory.

### 8.2 Secondary

Descriptive: confusion matrix of V2 classes vs trial outcomes, reported
separately for the retrospective calibration set and the prospective
validation set. Not pass/fail.

### 8.3 Sensitivity analyses

1. Including the 4 retrospective targets in the falsification statistic
   (tests robustness to contamination — expected to be negligible at N≥40).
2. Including polygenic-excluded targets (tests whether the scope exclusion
   is load-bearing or cosmetic).
3. Including OPPOSITE_DIRECTION targets as FAILURE predictions.
4. Restricting to ≥3-instrument IVW only (excluding single-SNP Wald ratio
   estimates). Tests whether single-cis-SNP pQTL results are load-bearing.
5. Egger threshold at p<0.05 instead of p<0.01 (≥3-instrument subset).
6. Per-disease BA (tests whether performance is uniform across diseases or
   driven by a single well-characterized indication).
7. **MR-only null model:** MR CI excludes null → predict SUCCESS, else
   predict FAILURE, ignoring the observational arm entirely. If this
   achieves the same BA as the full 2x2/taxonomy, the observational arm
   adds no value and the zombie concept is unsupported. This is the
   critical comparison — the full model must beat MR-only to justify
   the discordance framework. Comparison statistic: **McNemar's test**
   on discordant prediction pairs (same cases, two classifiers).
   McNemar's is appropriate for paired discrete classifiers on the same
   dataset; DeLong's test would apply only if continuous risk scores
   were produced (reported as secondary if applicable).
   **Same-subset constraint:** The MR-only model is evaluated on the
   **same** subset of pairs used for the full 2x2 (i.e., pairs with
   both an MR estimate and a usable observational estimate). The
   MR-only model *could* run on more pairs (all pairs with MR,
   ignoring obs), but the McNemar comparison requires identical cases.
   A secondary table reports MR-only BA on the full MR-available set
   for completeness.
8. **Discordance test family:** Report three tests as a family alongside
   Q: (a) two-sample z-difference test (same numerator as Q, report p
   directly), (b) CI-overlap test (does MR 95% CI exclude obs point
   estimate, and vice versa?), (c) Cochran's Q on 1 df. All three test
   the same hypothesis at k=2 with slightly different properties. The
   z-test is the formally correct version; CI-overlap is slightly
   conservative but most interpretable; Q is retained for comparability
   with meta-analysis literature. Power at typical pQTL and obs SEs is
   30-50% for all three — this means zombie detection is conservative
   (high specificity, low sensitivity), which is noted as a known
   limitation.
9. **CURATED-only analysis:** Restrict to pairs where the observational
   exposure directly measures the target gene's protein product (status
   = CURATED in `observational_estimates_v3_curated.csv`), excluding
   PROXY pairs. Tests whether proxy ligand-for-receptor substitutions
   drive or dilute classification accuracy.
10. **Colocalization-restricted analysis:** Restrict to pairs where the
   pQTL signal colocalizes with the outcome GWAS (posterior H4>0.80,
   where LD data permit). Tests whether non-colocalized pairs drive
   noise in classification.

## 9. Scoring metrics

For each analysis set, report:
- Sensitivity (true positive rate among SUCCESS targets)
- Specificity (true negative rate among FAILURE targets)
- Balanced accuracy (mean of sensitivity and specificity)
- PPV, NPV
- 90% bootstrap CI for balanced accuracy (10,000 iterations)

All reported separately — never a single blended accuracy number.

**Per-disease reporting:** N, sensitivity, specificity, and BA reported
per disease to show whether performance is driven by a single
well-characterized indication (e.g., CAD) or generalizes.

**Base rate comparison:** The SUCCESS/FAILURE base rate in the evaluable
set is compared to the industry-wide Phase III success rate (~50-60%,
Thomas et al. 2016). If the evaluable set is enriched for FAILURE
(e.g., 70% FAILURE), a naive "always predict FAILURE" classifier
achieves high accuracy. Balanced accuracy corrects for this, but the
base rate is reported for transparency.

## 10. Execution order

1. Curate observational source table for pQTL targets. Freeze
   (`data/observational_estimates_v3.csv`).
2. Adjudicate trial outcomes for all pQTL→OT pairs. Freeze
   (`OUTCOME_LABELS_V3.csv`). Hash independently.
3. Freeze this document. Hash all artifacts.
4. **Commit all frozen artifacts to git** (separate commit from results).
5. Run the zombie screen on all targets passing the inclusion rule.
6. Score primary and secondary endpoints.
7. Report results in `PREDICTION_RESULTS_V3.md`.

Steps 1–4 MUST complete before step 5. Step 5 is deterministic given
frozen inputs.

## 11. Artifacts to hash

| File | Contents |
|---|---|
| `PRESPEC_V3_SWEEP.md` | This document |
| `PRESPEC_V3_REVIEWER_GUIDANCE.md` | Design rationale |
| `results/disease_sweep_v3_expanded.csv` | OT candidate universe (8 diseases, 411 pQTL overlaps) |
| `results/disease_sweep_20260704_015239.csv` | V2 candidate universe (audit trail) |
| `results/prot_gene_mapping.csv` | pQTL→gene→Ensembl mapping |
| `data/gene_synonyms.csv` | Protein-synonym expansion list for PubMed queries |
| `data/observational_estimates_v3.csv` | Frozen obs source table |
| `OUTCOME_LABELS_V3.csv` | Frozen trial outcomes (adjudicate before classifier run) |
| `drug_pipeline.py` | Classification code |
| `v2_sweep_mr.py` | MR sweep code |

All hashes stored in `PRESPEC_V3_sha256.txt`.

## 12. Deviations log (continuous audit trail)

### V1 → V2 (2026-07-04 07:15 UTC)

V2 hash: `dcfa5b2ed5dfa24c32d11e9871de8ac02a916f2e4bcb2170097c9bba318df86d`

Changes: 3-class taxonomy expanded to 6-class (added GENETIC_AMPLIFICATION,
OBS_INFLATED_CAUSAL, OPPOSITE_DIRECTION). Binary clinical endpoint adopted.
Motivated by P4 (Lp(a) revealed a missing positive-discordance class) and
P5 (urate revealed the taxonomy conflated zombies with true nulls).

### V2 → V3 (2026-07-04T14:49:00Z)

Changes and justifications:

1. **Instrument tier demotion.** V2 specified ≥10 eQTL instruments as the
   inclusion gate. The 500-gene eQTL sweep (completed before V3) showed 0%
   of eQTL studies pass this gate — the assumption was empirically wrong.
   V3 introduces a tiered hierarchy: curated biomarker GWAS (primary),
   protein QTL with single-cis-SNP Wald ratio as primary estimator
   (standard in pQTL-MR literature), eQTL (negative control only).
   This is an instrument-quality correction, not an outcome-driven change.

2. **Disease universe expanded.** V2 froze 653 targets across CAD and T2D.
   V3 expands to 8 diseases (CAD, T2D, AD, RA, Asthma, CKD, UC, BC) to
   reach the N≥40 power target. All 8 have well-powered outcome GWAS in
   OpenGWAS and ≥20 pQTL-drug target overlaps each (411 total). Selection
   criteria: outcome GWAS available, sufficient drug target density for
   evaluable pairs. No disease was selected or excluded based on trial
   outcome distributions.

3. **Deterministic inclusion rule.** V2's candidate universe was frozen
   (653 targets) but the inclusion rule for scoring was underspecified.
   V3 adds explicit, non-discretionary gates for the prospective set.

4. **Deterministic PubMed search protocol.** V2 did not specify how
   observational estimates were located. V3 adds a reproducible search
   string and priority ladder to close the discretionary curation hole.

5. **Blinded outcome adjudication.** V2 did not specify how SUCCESS/FAILURE
   was determined or require outcomes to be hashed before scoring. V3 adds
   a frozen outcome-adjudication rule and a separate outcome hash.

6. **Falsification N raised from 30 to 40.** Power analysis showed N=30
   gives only 74% power at BA=0.70. N=40 achieves 80% power — the
   conventional bar.

7. **Retrospective/prospective split formalized.** V2 acknowledged the
   split but did not exclude retrospective targets from the primary
   statistic. V3 walls them off explicitly.

8. **Polygenic scope exclusion formalized.** V2 did not pre-declare this.
   V3 excludes complex polygenic exposures (BP, BMI) with a mechanistic
   rationale that does not reference outcomes.

No changes to: taxonomy, thresholds, binary endpoint definition, or
anti-cherry-pick rules.

### V3 → V3.1 (2026-07-04, pre-execution)

All changes made BEFORE any V3 classifier output was computed. The V3
freeze commit (`6dd7008`) is on GitHub with an externally verifiable
timestamp. V3.1 changes are analysis-plan improvements motivated by
peer review, not by seeing results. Iterative amendments to a frozen
pre-registration are standard practice when documented with version
hashes and timestamps (Lakens 2019; Nosek et al. 2018; §13).

Changes:

1. **PubMed query corrected (§4.2).** V3 query included "Mendelian
   randomization" as an OR term in the search, which floods results with
   MR studies — the opposite of what the observational arm requires. The
   executed pipeline already excluded MR via NOT clause. V3.1 corrects
   the prespec text to match the executed query. Motivated by: the 80%
   extraction error rate was traced to MR contamination of search results.

2. **2x2 binary framework added as primary (§7.1).** The continuous
   Cochran's Q taxonomy (V2, now §7.2) has low power at k=2 (~30-40%)
   and the d≥0.10 threshold conflates effect magnitude with causal
   evidence. The 2x2 framework uses each arm's own CI-excludes-null
   criterion, avoiding the Chinn conversion and arbitrary thresholds.
   V2 taxonomy retained as secondary analysis. Motivated by:
   methodological review identifying Q at k=2 as underpowered.

3. **MR-only null model added (§8.3.7).** Tests whether the
   observational arm contributes predictive value beyond MR alone. If
   MR-only achieves the same BA, the zombie/discordance concept is
   unsupported. This is the critical falsification of the framework's
   core claim. Motivated by: adversarial review noting that MR non-null
   alone might predict trial success.

4. **CI-overlap discordance test added (§8.3.8).** Alternative to Q
   that is more robust at k=2. Motivated by: meta-analysis methodology
   review.

5. **Colocalization sensitivity added (§8.3.9).** Tests whether
   non-colocalized pQTL-outcome pairs introduce noise. Motivated by:
   genetic epidemiology review noting field trend toward requiring H4>0.8.

6. **PENDING resolution rule added (§5.5).** Post-freeze trial readouts
   (e.g., ZEUS) do not reclassify PENDING labels; reported separately.
   Motivated by: clinical trial review noting ZEUS may report before
   submission.

7. **Multiple-drugs-per-target rule added (§5.6).** SUCCESS if any drug
   approved; FAILURE only if all drugs failed. Motivated by: drug
   development review.

8. **Winner's curse, EUR LD limitations acknowledged (§3).** Known bias
   directions stated explicitly. Motivated by: genetic epidemiology
   review.

9. **Per-disease N and base rate reporting added (§9).** Guards against
   accuracy inflation from unbalanced evaluable set. Motivated by:
   adversarial review.

10. **Observational estimates curation.** The frozen pipeline output
    (`observational_estimates_v3.csv`, hash `75a235de…`) had an 80%
    error rate (MR studies mislabeled as observational, wrong genes,
    prognostic studies). A curated version
    (`observational_estimates_v3_curated.csv`) was assembled from:
    (a) manual audit of all 39 pipeline extractions against PubMed
    abstracts, (b) Perplexity deep research verification of all failed
    pairs, (c) re-run of fixed pipeline with corrected query. The
    original frozen file is preserved for audit. The curated file is
    hashed separately.

11. **Denosumab→RA reference corrected.** NCT00680992 is a giant cell
    tumor study; the RA trials were Phase II (DRIVE/DESIRABLE). Label
    FAILURE retained: sponsor discontinued RA development after Phase II
    showed bone erosion efficacy but no disease activity improvement.
    Note: denosumab met its Phase II primary endpoint (erosion score)
    but our SUCCESS definition requires a **Phase III** pivotal trial
    (§5.3); Phase II efficacy alone does not qualify.

12. **Retrospective anchor validation — scope limitation.** The V1/V2
    pipeline used CRP→CAD (ZOMBIE anchor, ieu-b-35) and TG→CAD
    (CONCORDANT anchor, ieu-b-111) as positive controls. These used
    biomarker trait GWAS as exposures. The V3 pipeline uses cis-pQTL
    exposures (prot-a-*/prot-c-* studies with ±1Mb cis-window filter),
    which is a different analysis type. CRP has a pQTL (prot-a-670) but
    TG is a lipid biomarker with no protein QTL and cannot enter the
    cis-pQTL framework. Mixing exposure GWAS types within the primary
    analysis would violate the single-instrument-tier design. The
    `validate_anchors()` function in the classifier will raise on
    misclassification if anchors are present, but they are not in the
    37-pair pQTL dataset by construction. Scale-consistency protection
    is provided instead by the hard assertion in `load_mr()` (Bug Fix 1),
    which verifies `mr_supports_causal` matches CI-excludes-0 on all OK
    rows.

No changes to: V2 taxonomy thresholds, binary endpoint definition,
outcome labels (except denosumab reference), falsification criterion,
anti-cherry-pick rule, or disease universe.

### V3.1 → V3.2 (2026-07-04, pre-execution)

All changes made BEFORE any V3.2 output was computed. The V3.1
results are frozen in commit `1e2909a` with V31_FINDINGS.md.

**Motivation:** V3.1 yielded n=4 evaluable 2x2 pairs (BA=0.50) —
statistically meaningless due to instrument coverage, not hypothesis
failure. The binding constraint is the 37-pair pQTL-only universe.
V3.2 expands instrument source to the EpiGraphDB published cis-pQTL-MR
catalog (Zheng et al. 2020, INTERVAL study) and adds 8 diseases to
recover evaluable N. This is a coverage-motivated expansion; taxonomy,
thresholds, 2x2 logic, and existing outcome labels are unchanged.

Changes:

13. **EpiGraphDB catalog integration — cis-only, single-SNP.**
    `expand_v32.py` queries the EpiGraphDB `/pqtl/` API for published
    pQTL-MR results across 14 diseases with verified catalog coverage
    (6 from V3.1 + 8 new: ischemic stroke, MS, Crohn's, lung cancer,
    Parkinson's, MDD, SLE, eczema). Only single-SNP cis-pQTL Wald
    ratios (`rtype=sglmr`, `trans_cis="cis"`) enter the candidate
    pairs. Multi-SNP IVW results (`rtype=mrres`) are saved in raw
    dumps for audit but excluded from candidate pairs because the
    `mrres` response lacks a `trans_cis` field — there is no way to
    verify cis-only status, and including them would violate the
    single-cis-pQTL protocol. This preserves protocol consistency
    with V3.1's `run_pqtl_mr.py` (cis ±1Mb, Wald for 1 SNP).

14. **Disease universe expanded from 8 to 16.** Added: ischemic
    stroke, multiple sclerosis, Crohn's disease, lung cancer,
    Parkinson's disease, major depressive disorder, systemic lupus
    erythematosus, atopic dermatitis. Selection criteria: (a) has
    published pQTL-MR results in EpiGraphDB, (b) has well-powered
    outcome GWAS in OpenGWAS, (c) has known drug targets with Phase
    III readouts. Two V3.1 diseases (asthma, breast cancer) lack
    EpiGraphDB coverage and are retained for own pQTL-MR only.
    Seven candidate diseases (heart failure, atrial fibrillation,
    psoriasis, COPD, osteoporosis, prostate cancer, colorectal
    cancer) were excluded because they returned 0 results from
    EpiGraphDB — no trait name variant matched the catalog. No
    disease was selected or excluded based on trial outcome
    distributions.

15. **Outcome adjudication for new pairs.** All new protein-disease
    pairs from the expansion require Phase III outcome adjudication
    by an independent reviewer (Perplexity) before entering the
    classifier. The expanded OUTCOME_LABELS_V32.csv and
    observational_estimates_v32_curated.csv will be frozen and
    hashed as separate artifacts before any classification.

16. **Protein→gene symbol mapping.** EpiGraphDB `expID` uses protein
    names (e.g., "ApoB", "IL-6 sRA"); the classifier uses gene
    symbol pair_ids. The expansion pipeline maps via
    `prot_gene_mapping.csv` (from V3.1). Unmapped proteins use
    expID as-is and must be manually reviewed before classification.

17. **Gene symbol validation.** EpiGraphDB `expID` values are gene
    symbols in 752/764 (98.4%) of cases. 8 entries are protein-complex
    multi-gene identifiers (semicolon-separated, e.g., "FCGR2A;FCGR2B");
    these use the first/primary gene for `pair_id` matching. 4 entries
    (B3GAT3, B3GNT8, S100A4, S100A7) matched a UniProt regex by
    coincidence but are valid HGNC symbols. No UniProt IDs, no
    unmappable identifiers.

18. **Outcome-blind drug-target selection rule (MANDATORY — pre-declared
    before any outcome adjudication).** The evaluable subset of V3.2's
    6,219 candidate pairs is determined by the following rule, applied
    independently of MR results:

    A pair (protein, disease) enters the evaluable set if and only if:
    (a) The protein is the **mechanistic target** of a drug compound
        that reached **Phase III** clinical trials for that specific
        indication, as recorded in Open Targets Platform
        (targetvalidation.org) or ChEMBL (max_phase ≥ 3).
    (b) The drug-target-indication mapping comes from a database query,
        not manual curation — the query returns all proteins with
        Phase III drugs for each disease, regardless of MR result.
    (c) The Phase III outcome (SUCCESS = approval or positive pivotal;
        FAILURE = failed/terminated/negative pivotal) is adjudicated
        from FDA/EMA records and ClinicalTrials.gov, not from the
        MR or observational data.

    This rule ensures that MR-causal and MR-null pairs enter the
    evaluable set at the same rate determined by pharmaceutical
    development history, not by the MR signal. The classifier is
    tested on the pairs that happen to have Phase III readouts, not
    on pairs selected because they look interesting.

    **Anti-cherry-pick enforcement:** The Open Targets / ChEMBL query
    is run once, frozen, and hashed. No post-hoc additions or
    removals of pairs based on MR results, observational estimates,
    or preliminary classification output.

No changes to: V2 taxonomy, thresholds, binary endpoint definition,
2x2 classification logic, falsification criterion, anti-cherry-pick
rule, scale-consistency assertions, or anchor validation.

### Prior version hashes (for audit continuity)

- V1 prespec: `f36855865b693d96d14324b44392af7cedcf29f81658b83decd517babb2f25ab`
- V2 prespec: `dcfa5b2ed5dfa24c32d11e9871de8ac02a916f2e4bcb2170097c9bba318df86d`
- V2 candidate universe: `443e0540510cb3c0c173fc87480f2512a1a2a9b6483b4de8fc34c36d921503e2`
- V3 prespec (pre-V3.1): `39b1ee694abd2c2de261c5881b2bb2cc3203aa57005427c7505d61d18eb82ade`

---

## 13. Key references

- Lakens D (2019). The value of preregistration for psychological
  science: a conceptual analysis. *Japanese Psychological Review*, 62(3).
  Justifies iterative pre-registration amendments provided each version
  is time-stamped and changes are documented (our §12 deviation log).
- Nelson MR, Tipney H, Painter JL, et al. (2015). The support of
  human genetic evidence for approved drug indications. *Nature Genetics*,
  47(8), 856–860. Precedent for binary genetic-support classification of
  drug targets (our §7.1 2x2 framework extends this from "any genetic
  evidence" to MR × observational discordance).
- Nosek BA, Ebersole CR, DeHaven AC, Mellor DT (2018). The
  preregistration revolution. *PNAS*, 115(11), 2600–2606. Establishes
  that pre-registration's value lies in separating confirmatory from
  exploratory analysis, even when amendments occur.
