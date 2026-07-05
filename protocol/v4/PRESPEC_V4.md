# Pre-Specification V4: LoF Burden Instrument + Evidence Routing

**Author:** Elliot Tower  
**Date:** 2026-07-05  
**Status:** FROZEN (pending SHA-256 hash)  
**Depends on:** V3.4 (Zenodo Version 2)

## Overview

V4 adds a second genetic instrument — rare-variant loss-of-function (pLoF) burden from Genebass (UKB exome sequencing) — to the existing cis-pQTL MR instrument from V3.4. The goal is to test whether the two instruments, which have orthogonal failure modes, produce concordant or discordant evidence on the same drug targets.

This is NOT a replication of V3.4. It is an extension that asks: does a knockout-flavored instrument agree with an abundance-flavored instrument, and does agreement predict Phase III success better than either alone?

## Instrument Definitions

### Instrument 1: cis-pQTL MR (from V3.4)

- Source: EpiGraphDB + UKB-PPP cis-pQTL catalogs
- Test: Single-SNP Wald ratio MR
- Signal: Effect of lifelong variation in circulating protein abundance on disease risk
- Known failure mode: Blind to activity-blocking drugs (V3.4 primary finding)
- Results: Frozen in V3.4 classification (138 analyzed pairs)

### Instrument 2: pLoF Burden (NEW)

- Source: Genebass (Karczewski et al., UKB 450k exome sequencing)
- Test: Gene-level pLoF burden test (SKAT-O framework)
- Signal: Effect of predicted loss-of-function variants on disease phenotype
- Known failure mode: Underpowered for most genes (rare variants); tests knockout, not graded inhibition
- Phenotype mapping: One canonical phenocode per disease, chosen by scientific relevance and expected case count (see canonical_disease_mapping.json)

### Why these are orthogonal

cis-pQTL MR asks: "what happens when this protein's circulating abundance changes?" pLoF burden asks: "what happens when this gene is knocked out by rare coding variants?" These are mechanistically distinct questions. A drug that blocks enzymatic activity without changing protein level would be invisible to pQTL but potentially visible to LoF (knockout removes both abundance and activity). The failure modes are independent.

## Evaluable Set

- V3.4 analyzed pairs: 138
- Pairs with LoF burden result: 132 (96%)
- Structural missingness: 6 pairs across 4 diseases
  - Anorexia nervosa (3 pairs): absent from Genebass, low UKB prevalence
  - Autism spectrum disorder (1 pair): underascertained in adult UKB
  - Juvenile idiopathic arthritis (1 pair): pediatric disease
  - Neuroblastoma (1 pair): pediatric cancer
- Pairs with both pQTL and LoF instruments: 132

## Pre-Registered Analyses

All analyses use the 132 pairs with both instruments. Outcomes (SUCCESS/FAILURE) are carried from V3.4 adjudication — no new outcome adjudication.

### Analysis 1 (Primary): LoF Burden AUC

Continuous ROC analysis using the signed LoF burden z-score as the predictor of Phase III outcome. Permutation test for significance (10,000 permutations of outcome labels).

Rationale: With only 7/132 pairs significant at p<0.05, a thresholded classifier has near-zero sensitivity. The continuous AUC extracts signal from the full z-score distribution.

### Analysis 2 (Primary): Directional Concordance Between Instruments

For each of the 132 pairs, check whether the pQTL MR beta and the LoF burden beta agree in sign. Report:
- Overall concordance rate (and binomial test vs 50%)
- Concordance stratified by mechanism class (activity-blocking vs abundance-modulating)
- Concordance stratified by V3.4 outcome (SUCCESS vs FAILURE)

Rationale: If both instruments point the same direction, the evidence is stronger regardless of whether either crosses a significance threshold. This is the core "evidence routing" test.

### Analysis 3 (Primary): Depth-2 Licensing

A target is "depth-2 licensed" when:
1. Both pQTL and LoF burden results exist for that gene-disease pair
2. Both instruments agree in effect direction (same sign of beta)
3. At least one instrument is nominally significant (p < 0.05)

Report: BA and AUC for depth-2-licensed targets vs all targets. Test whether depth-2-licensed targets have higher predictive accuracy than the full set.

Rationale: The invariance-depth hypothesis predicts that concordance across orthogonal instruments is more predictive than any single instrument. This is the minimal test of that claim.

NOTE: Colocalization (coloc/SuSiE) is the gold standard for confirming shared causal signal. We do not have individual-level data to run formal colocalization. The depth-2 rule above uses directional agreement as a proxy. This is a known limitation — depth-2 licensing without colocalization is weaker evidence than with it.

### Analysis 4 (Secondary): Thresholded LoF Classifier

Binary classifier: LoF burden p < 0.05 predicts SUCCESS. Report BA, sensitivity, specificity, PPV.

Rationale: Expected to be near-chance due to sparse significance. Included for completeness and comparison with V3.4 pQTL classifier.

### Analysis 5 (Secondary): Mechanism-Stratified LoF AUC

Repeat Analysis 1 separately for activity-blocking and abundance-modulating targets.

Prediction: If LoF burden is orthogonal to pQTL, then LoF should show LESS mechanism-dependence than pQTL did. Specifically, LoF should not be blind to activity-blocking targets the way pQTL was, because knockout removes both abundance and activity. This is the key test of orthogonality.

### Analysis 6 (Exploratory): Combined Instrument Score

Average the signed z-scores from pQTL and LoF for each pair. Compute AUC for the combined score. Compare to individual-instrument AUCs.

Rationale: If the instruments are truly orthogonal and both carry signal, combining should improve discrimination. If they are redundant, combining should not help.

## Phenotype Mapping Caveats (Pre-Declared)

1. **Glioma mapped to C71 (all brain cancer):** ICD-10 C71 is broader than glioma, including non-glioma brain tumors. Chosen for power (more cases). The LoF signal is diluted toward any brain malignancy.

2. **Composite phenotypes (IBD, Crohn's, Depression):** Purpose-built Genebass phenotypes may aggregate sub-phenotypes differently than the specific disease a drug trial targeted.

3. **CKD mapped to serum creatinine (30700):** Continuous biomarker rather than binary diagnosis. Higher creatinine = impaired renal function. Direction convention: positive LoF beta = increased creatinine = disease-promoting.

4. **Thyroid cancer mapped to C73 (ICD-10):** Thyroid_custom may include non-cancer thyroid conditions; ICD-10 is more specific for cancer.

## What Does NOT Change

- Classifier for pQTL: MR p < 0.05 predicts SUCCESS (identical to V3.3/V3.4)
- Outcome adjudication: identical to V3.4
- Anti-cherry-pick rule: mechanical intersection of independent databases
- Falsification criterion: 90% CI lower bound > 0.50

## Deviation Rules

Any post-hoc deviation from this specification must be:
1. Documented in DEVIATION_LOG_V4.md with rationale
2. Clearly labeled as post-hoc in the paper
3. Not used to replace a pre-registered analysis that produced an unfavorable result

## Files to Hash

1. `PRESPEC_V4.md` (this document)
2. `canonical_disease_mapping.json` (one phenocode per disease)
3. `query_genebass_lof_v2.py` (instrument collection script)
4. `genebass_lof_results_v2.json` (collected instrument data, 132 pairs)
