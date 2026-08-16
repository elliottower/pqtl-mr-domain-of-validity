# Pre-Specification V5: pQTL Instrument Expansion

**Author:** Elliot Tower  
**Date:** 2026-07-13  
**Status:** FROZEN (pending SHA-256 hash)  
**Depends on:** V3.4 (Zenodo Version 2), V4 (LoF Extension)

## Motivation

V3.4 reports BA=0.590 for the abundance-modulating stratum, with 90% bootstrap CI [0.496, 0.691]. The lower bound includes 0.50 (chance), weakening the claim that pQTL MR is informative for abundance-modulating drugs. Of the 101 abundance-modulating candidate pairs identified in Stage B, only 59 (58%) obtained MR results. The remaining 42 pairs are structurally missing due to three causes:

1. **No outcome GWAS ID** for 2 diseases (juvenile idiopathic arthritis, neuroblastoma): 8 pairs (6 UKB-PPP, 2 EpiGraphDB)
2. **Non-rs variant identifiers** in UKB-PPP that OpenGWAS cannot resolve: 8 pairs (CD86 x4 diseases, PDCD1 x2, MMP13 x1, MIF x1)
3. **Instrument SNP absent from outcome GWAS**: remaining 26 pairs where the sentinel rs-number exists but returned no association from OpenGWAS

This protocol expands the analyzed set by addressing all three causes. No new drug targets are added. No outcome adjudication is changed. The same zero-parameter classifier (MR p<0.05 -> predict SUCCESS) applies.

## Expansion Steps

### Step 1: Add Outcome GWAS for JIA and Neuroblastoma

Add OpenGWAS IDs for two diseases currently mapped to `None`:

| Disease | OpenGWAS ID | Trait Label | Sample Size | SNPs |
|---------|-------------|-------------|-------------|------|
| Juvenile idiopathic arthritis | ebi-a-GCST90018873 | Juvenile rheumatoid arthritis | 409,217 | 24,190,478 |
| Neuroblastoma | ieu-a-816 | Neuroblastoma | 4,881 | 468,788 |

Expected recovery: up to 8 pairs (6 UKB-PPP + 2 EpiGraphDB with no sentinel RSIDs that might be rescued by the Step 3 deCODE lookup).

### Step 2: Resolve Non-rs Variant Identifiers

For UKB-PPP sentinel variants in chr:pos_ref_alt format, resolve to rs-numbers using the ENSEMBL Variant Recoder API (GRCh37/hg19). If no rs-number exists for the exact variant, attempt LD proxy lookup via the ENSEMBL LD endpoint (r-squared >= 0.8, EUR population, window 500kb).

Affected variants:

| Variant ID | Gene | Chr:Pos | Diseases |
|------------|------|---------|----------|
| 3:121822003_GAATGGGTTGTCTAACTCTACTGGTTTTTCCC_G | CD86 | 3:121822003 | CKD, Crohn's, RA, UC |
| 2:242801752_CG_C | PDCD1 | 2:242801752 | Lung cancer, Ovarian cancer |
| 11:102826877_GC_G | MMP13 | 11:102826877 | Alzheimer's |
| 22:24256692_AT_A | MIF | 22:24256692 | RA |

Expected recovery: up to 8 pairs if rs-numbers are found and present in outcome GWAS.

### Step 3: LD Proxy Lookup for Remaining Failures

For the 26 pairs where a valid rs-number exists but returned no association from the outcome GWAS, attempt LD proxy lookup:

1. Query ENSEMBL LD endpoint for each sentinel SNP (EUR, r-squared >= 0.8, window 500kb)
2. Select the proxy with highest r-squared that is present in the outcome GWAS
3. Compute Wald ratio using the proxy SNP's outcome association, adjusting the exposure beta/SE for LD (multiply by signed sqrt of r-squared)

This is a standard LD proxy MR approach used in TwoSampleMR.

Expected recovery: 5-15 pairs (depends on LD structure and GWAS coverage).

### Step 4: Add deCODE Genetics cis-pQTL Source

For the 8 EpiGraphDB pairs that lack sentinel RSIDs entirely (EpiGraphDB catalog has these genes but no SNP-level data), look up cis-pQTL instruments from Ferkingstad et al. 2021 (Nature Genetics):

- Platform: SomaScan (4,907 aptamers)
- Sample: 35,559 Icelanders
- Summary statistics: decode.com/summarydata/

Target genes: GHR, CSF2RB, TNC, PCSK9 (x2 diseases), SNCA, CLEC4C

For each gene, extract the top cis-pQTL (within 1Mb of the gene body, genome-wide significant p < 5e-8) as the instrument. Run standard Wald ratio MR against the appropriate outcome GWAS.

Expected recovery: 4-6 pairs (deCODE covers ~4900 proteins; PCSK9 and GHR are almost certainly present).

## Analysis Plan

All analyses from V3.4 are re-run on the expanded dataset. The primary analysis remains:

**Analysis 1 (Primary):** Balanced accuracy of the zero-parameter classifier (MR p<0.05 -> SUCCESS), stratified by mechanism class (abundance-modulating vs activity-blocking). 10,000 bootstrap resamples for 90% CI.

**Analysis 2 (Primary):** Permutation AUC test (10,000 permutations) for the abundance-modulating stratum.

**Analysis 3 (Sensitivity):** Compare expanded results to V3.4 results. Report whether the point estimate and CI shift are consistent with random sampling from the same generating process.

**Analysis 4 (Missingness):** Re-run missingness bounds from V3.4 on any remaining missing pairs.

## Decision Rules

- If expanded abundance CI lower bound > 0.50: report as primary result in paper, note V3.4 as preliminary.
- If expanded abundance CI still includes 0.50 but point estimate remains > 0.55: report both V3.4 and V5 results, note that the abundance arm trends informative but sample size limits power.
- If expanded abundance point estimate drops below 0.55: the abundance-modulating claim is weakened; discuss possible explanations.

## Code

Pipeline implemented in `code/expand_v5.py`. All intermediate results saved to `results/v5/`. No manual data entry — everything computed from APIs.

## Integrity

This document will be SHA-256 hashed before any Step 1-4 queries are executed. The hash and timestamp will be recorded in `protocol/v5/PRESPEC_V5_sha256.txt`.
