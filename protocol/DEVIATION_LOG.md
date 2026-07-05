# V3.4 Deviation Log

**Status:** STAGE B COMPLETE — evaluable N=195 (63S/132F), ready for Stage C.
**Supersedes:** V3.3 (33 diseases, N=50 evaluable, BA=0.580)
**Stage A frozen:** 2026-07-04T23:51Z
**Stage B completed:** 2026-07-05T00:44Z

## Rationale

V3.3 (N=50) is underpowered for composite significance: the 90% CI
for balanced accuracy [0.496, 0.667] includes 0.50. The mechanism-
stratified result (abundance BA=0.708, CI excludes chance) is the
paper's thesis, but a reviewer can legitimately object that the
composite fails the pre-declared falsification criterion. Expanding
to N>=80 would either confirm the composite is a genuine mixture of
signal and noise, or reveal that the mechanism-stratified result was
a small-sample artifact.

All genetic filters, classifier logic, and analysis framework are
unchanged from V3.3. The ONLY changes are: (1) wider disease
universe, (2) additional pQTL instrument sources, (3) a pre-declared
Phase II secondary tier.

## What does NOT change from V3.3

| Parameter | V3.3 value | V3.4 value | Changed? |
|-----------|-----------|-----------|----------|
| Instrument type | cis-pQTL, single-SNP Wald ratio | same | NO |
| cis window | as defined in EpiGraphDB catalog | same | NO |
| Significance threshold | genome-wide (p<5e-8) or catalog default | same | NO |
| trans-pQTL | excluded | excluded | NO |
| multi-SNP IVW | excluded | excluded | NO |
| Classifier | MR p<0.05 -> CAUSAL -> predict SUCCESS | same | NO |
| Classifier parameters | zero | zero | NO |
| Falsification criterion | BA 90% CI lower bound > 0.50 | same | NO |
| Outcome adjudication | FDA/EMA + ClinicalTrials.gov | same | NO |
| TARGET_MISMATCH rule | exclude if drug target != pQTL protein | same | NO |
| Anti-cherry-pick rule | mechanical intersection of independent DBs | same | NO |

## Change 1: Disease universe expansion (33 -> ~70)

### Pre-specified inclusion rule

A disease is included if and only if:
(a) it has >=1 mapped EFO/MONDO term in Open Targets with >=3
    distinct drugs at Phase III+, AND
(b) it has a publicly available GWAS in IEU OpenGWAS or the GWAS
    Catalog with N>=10,000.

Diseases are enumerated by querying Open Targets for all indications
meeting (a), then filtering by (b). No disease is added or removed
by hand. If the rule returns diseases not anticipated, they are kept.
If it excludes anticipated diseases (e.g., psoriasis), that is
documented but not overridden.

### Oncology deduplication

Each tumor type counts as a distinct indication, provided:
- The drug's mechanism is the same across indications
- The same pivotal NCT trial does not appear under two disease labels
- Dedup key: (gene, drug, NCT_or_regulatory_ref)

## Change 2: Instrument source expansion

### Added sources

| Source | Proteins | Approx N | Assay |
|--------|----------|----------|-------|
| EpiGraphDB/INTERVAL | ~1,740 | ~3,300 | SomaScan |
| UKB-PPP | ~2,900+ | ~54,000 | Olink |
| deCODE | ~4,900 | ~35,500 | SomaScan |

### Instrument selection rule (pre-specified)

**Sentinel cis-pQTL, single-SNP Wald ratio.** To maintain estimator
consistency with V3.1--V3.3, every protein uses one instrument: the
single sentinel cis-pQTL (lowest p among conditionally independent
cis signals in the source study). The MR estimate is a single-SNP
Wald ratio (beta_outcome / beta_protein). Multi-SNP IVW and trans
instruments remain excluded.

### Instrument source precedence (pre-specified)

When cis-pQTLs exist in multiple sources for the same protein, use
the source with the largest discovery N for that protein:
UKB-PPP (N~54,000 Olink) > INTERVAL (N~3,300 SomaScan), resolved
per-protein by actual sample size. Ties broken alphabetically by
source name. The MR p-value is NEVER used to select the instrument.

deCODE instruments are deferred to V3.5 as a cross-platform
robustness check (SomaScan aptamer vs Olink antibody).

### Gene-symbol harmonization rule

UKB-PPP protein annotations use Olink assay targets. Before
inclusion, every rescued gene must pass:
1. Map Olink assay target → HGNC gene symbol via UKB-PPP ST3
   Gene symbol column (not assay target column)
2. Reject multi-gene assays (e.g., IL12A_IL12B) unless one gene
   is unambiguously the Phase III drug target
3. Confirm gene symbol matches the Open Targets target gene
   (approvedSymbol), not a paralog or alias

### Sample overlap policy (pre-specified)

UKB-PPP pQTL instruments are derived from UK Biobank participants.
Two-sample MR requires non-overlapping samples between exposure
and outcome GWAS. For each protein-disease pair:

- **Primary (no overlap):** Use a disease GWAS with NO UK Biobank
  participants where available: FinnGen, CARDIoGRAM (CHD/MI),
  DIAGRAM (T2D), IMSGC (MS), PGC (psychiatric), IARC (cancers),
  IIBDGC (IBD/CD/UC), etc. These pairs count toward primary N.
- **Secondary (overlap-flagged):** Where no non-UKB GWAS exists,
  use UKB-derived GWAS but flag as `sample_overlap=True`. These
  pairs are reported separately and excluded from headline BA.
  A sensitivity analysis including them is pre-declared.
- **EpiGraphDB/INTERVAL pairs:** INTERVAL is a separate cohort
  (N~3,300). EpiGraphDB's outcome GWAS vary but most pre-date
  UKB inclusion. These pairs are presumed non-overlapping unless
  the specific outcome GWAS is known to include UKB.

### Outcome GWAS source mapping (computed, outcome-blind)

Non-UKB outcome GWAS selected per disease. Full mapping in
`pqtl_sources/non_ukb_gwas_sources.json`.

**Confirmed clean (pre-UKB or UKB-free consortium, 18 diseases):**
CARDIoGRAMplusC4D (CHD/MI, 2015), Okada (RA, 2014), IMSGC (MS),
Bentham (SLE, 2015), ILCCO/TRICL (lung, 2017), OCAC (ovarian,
2017), PanScan (pancreatic), PGC AN (anorexia), iPSYCH-PGC
(autism), GLGC (LDL, 2013), IIBDGC (IBD/UC, 2017), Hinks (JIA,
2013), COG (neuroblastoma), deCODE (thyroid, 2017), GICC (glioma),
GenoMEL (melanoma).

**Verified clean after checking UKB inclusion (3 diseases, 18 pairs):**
- Alzheimer's: Kunkle 2019/IGAP confirmed no UKB participants
  (UKB proxy cases only in separate Jansen meta-analysis)
- Crohn's: de Lange 2017 IIBDGC published Feb 2017, pre-UKB
  genotype release (mid-2017/2018)
- MDD: PGC2 (Wray 2018) includes UKB, but a non-UKB release is
  available: mdd-rmUKBB on PGC GitHub / Figshare

**Switched to FinnGen R10 (4 diseases, 17 pairs):**
- Bipolar: Mullins 2021 PGC BIP confirmed includes UKB cohort
- Parkinson's: Nalls 2019 includes 18,618 UKB proxy cases
- ALS: van Rheenen 2021 UKB inclusion uncertain → FinnGen fallback
- Schizophrenia: PGC3 (Trubetskoy 2022) UKB inclusion uncertain
  → FinnGen fallback

**Overlap-flagged (1 disease, 7 UKB-PPP pairs):**
CKD (CKDGen includes UKB, no clean alternative with adequate N).
These pairs are excluded from primary BA, reported in sensitivity.

### Sample overlap summary (verified)

```
clean (non-UKB outcome or INTERVAL instrument):  204
overlap_flagged (sensitivity analysis only):         7
```

FinnGen R10 is pre-declared as universal fallback for any disease
where the primary non-UKB GWAS fails verification.

### Oncology RTK caveat (pre-noted)

Several rescued genes (EGFR, ERBB2/3/4, MET, KIT, PDGFRA/B, FLT1)
encode receptor tyrosine kinases measured as shed/circulating
protein in plasma. The drug mechanism targets membrane-bound
receptor on tumor cells, not circulating protein abundance. These
pairs are expected to show MR-null (false negatives) because the
cis-pQTL proxies abundance, not tumor-cell signaling. This is the
same abundance-vs-activity blind spot documented in V3.3. These
pairs are retained but flagged as mechanism_class=activity_blocking
in adjudication.

### Required verification

For every protein with instruments in both EpiGraphDB/INTERVAL and
UKB-PPP:
1. Compute both MR estimates
2. Plot old beta vs new beta — should correlate near identity line
3. Flag any sign-flip or order-of-magnitude change
4. Record assay type (SomaScan vs Olink) per protein

### Yield funnel (computed)

```
A = Phase III drug-target genes across clinical indications:  576
B = with valid cis-instrument in INTERVAL (EpiGraphDB):        37
C_ukbppp = with valid cis-instrument in UKB-PPP:              108
C_any = with cis-instrument in EpiGraphDB OR UKB-PPP:         112
newly_rescued = C_any - B:                                      75
```

STOP CONDITION: if C - B < 10, the instrument expansion underperformed.
Result: C - B = 75. PASSED — proceed to pair enumeration.

### Evaluable pairs (computed)

```
Total unique gene-disease pairs:                               211
  EpiGraphDB/INTERVAL instrument:                               76
  UKB-PPP instrument (newly rescued):                          135
Unique genes with at least one pair:                           112
Unique clinical indications:                                    24
RTK/kinase-flagged pairs (expected false negatives):            53
Non-RTK pairs (primary abundance analysis):                    158
```

Pair list saved to `pqtl_sources/evaluable_pairs_v34.json` and
`results/v34/evaluable_pairs_v34.csv`.

## Change 3: Phase II sensitivity tier (SECONDARY analysis only)

### Specification

- Primary analysis remains Phase III only. The Phase II tier is a
  pre-declared secondary analysis, reported alongside, NEVER blended
  into the headline BA.
- Phase II inclusion rule: include a Phase II pair only if the
  outcome is definitive — a met/missed pre-registered primary
  efficacy endpoint in a randomized controlled Phase II, OR a
  program discontinuation explicitly attributed to efficacy failure.
- Exclude: safety-only, single-arm, terminated for business/funding,
  or ambiguous Phase II trials.
- Tag: evidence_tier in {phase3, phase2_definitive}

## Stage architecture (outcome-blindness enforcement)

```
STAGE A: candidate_selection.py
  -> outputs frozen_candidates_v34.csv (NO outcome column)
  -> hash immediately, commit hash, THEN STOP

STAGE B: adjudication.py
  -> reads frozen_candidates_v34.csv
  -> adds: outcome, drug, reference, mechanism_class, evidence_tier
  -> mechanism_class assigned from pharmacology BLIND to MR p-value
  -> runs ONLY after Stage A hash is logged

STAGE C: analysis.py
  -> reads adjudicated file, computes metrics
  -> anything computed for first time after seeing V3.4 outcomes
     is labeled post-hoc
```

CONTAMINATION CHECK: candidate_selection.py must NOT contain any
string matching: outcome, success, approval, NCT, phase3_result.
If it does, the script is contaminated and must be rejected.

## Change 4: Mechanism classification (pre-declared covariate)

### Classification rule

Each gene-disease pair is assigned a `mechanism_class` based on
whether the drug's therapeutic effect is mediated by changing the
target protein's circulating level (abundance) or by blocking a
specific enzymatic/receptor activity independent of protein level
(activity). The classification is determined by target biology and
drug modality, assigned BLIND to MR p-values, and frozen before
MR computation.

**abundance_modulating**: drug works by changing circulating protein
level. Includes antibodies that neutralize circulating proteins
(anti-TNF, anti-VEGF, anti-IL12/23, anti-PD1/PDL1), antisense
oligonucleotides (tofersen, inclisiran), protein replacement
(tPA, G-CSF, IFN-beta), receptor decoys (abatacept), and antibodies
that sequester secreted ligands (denosumab, romosozumab). For
secreted enzymes (MMPs, sPLA2, kallikrein), plasma level
approximates activity, so these are classified abundance.

**activity_blocking**: drug blocks enzymatic activity, receptor
signaling, or catalytic function without changing protein level.
Includes small-molecule kinase inhibitors (EGFR/ERBB/KIT/MET/PDGFR
inhibitors), enzyme active-site inhibitors (AChE, PARP, DPP4, ACE,
PDE5A, COMT, carbonic anhydrase), and receptor antagonists
(azeliragon, vorapaxar, gabapentinoids). Includes receptor tyrosine
kinases measured as shed ectodomain in plasma, where abundance does
not reflect membrane-bound signaling.

**mixed**: gene has drugs of both modalities (APP: anti-amyloid
antibodies + BACE inhibitors; ERBB2: trastuzumab + lapatinib;
IGF1R: mecasermin agonist + kinase inhibitors).

### Classification method

The classification is a **deterministic gene-level lookup** (Python
dict keyed by HGNC symbol), not per-pair judgment. The lookup was
constructed by examining each gene's Phase III drug portfolio:

- Genes whose Phase III drugs are predominantly mAbs, fusion
  proteins, ADCs, antisense oligos, or siRNAs → abundance
- Genes whose Phase III drugs are predominantly small-molecule
  enzyme/kinase inhibitors or receptor antagonists → activity
- Genes with both modalities → mixed

The `action_type` field from Open Targets (INHIBITOR, ANTAGONIST,
etc.) does NOT determine the classification. 69 pairs have
`action_type=INHIBITOR` but are classified `abundance_modulating`
because the "inhibitor" is a monoclonal antibody that neutralizes
the circulating protein (e.g., anti-TNF, anti-VEGF, anti-PD1), not
a small-molecule active-site blocker. The full cross-tab:

```
action_type          abundance  activity  mixed  total
INHIBITOR                   66        99      4    169
ANTAGONIST                  11         3      0     14
AGONIST                     12         0      2     14
ACTIVATOR                    5         0      0      5
MODULATOR                    0         4      0      4
BINDING AGENT                3         0      1      4
ANTISENSE INHIBITOR          3         0      0      3
OTHER                        2         0      1      3
EXOGENOUS PROTEIN            2         0      0      2
RNAI INHIBITOR               1         0      0      1
DISRUPTING AGENT             0         0      1      1
STABILISER                   0         0      1      1
PROTEOLYTIC ENZYME           1         0      0      1
```

### Distribution (computed, outcome-blind)

```
abundance_modulating:  101 pairs
activity_blocking:     106 pairs
mixed:                   4 pairs  (APP, ERBB2, IGF1R×2)
```

### Edge cases (documented before MR computation)

- MMP1/3/7/8/9/12/13: secreted enzymes, drugs are enzyme inhibitors,
  but plasma level ≈ activity → classified abundance (conservative
  for the drug-failure direction)
- PLA2G2A/PLA2G10: secreted phospholipases, same reasoning → abundance
- ACE: ectoenzyme, ACE inhibitors are small molecules, circulating ACE
  may not reflect tissue activity → classified activity
- F10/F2: circulating coagulation factors, drugs block active site,
  but level ≈ activity for circulating enzymes → classified activity
  (conservative)
- GHR: receptor antagonist (pegvisomant), but lower GHR level → less
  signaling → classified abundance
- APP: anti-amyloid antibodies (lecanemab, aducanumab) + BACE
  inhibitors (verubecestat) → mixed
- ERBB2: trastuzumab (antibody) + lapatinib (small molecule) → mixed
- IGF1R: mecasermin IGF-1 agonist (ALS/anorexia) + kinase inhibitors
  (oncology) → mixed

### Prediction (pre-declared)

MR predicts drug success for abundance_modulating pairs (cis-pQTL
proxies drug mechanism) and is expected to be uninformative for
activity_blocking pairs (cis-pQTL is blind to active-site
inhibition). This replicates and extends the V3.3 finding
(abundance BA=0.708, activity BA=0.480) on a larger, independently
selected set.

## Evaluable-set definition (locked before Stage C)

A pair is **evaluable** if and only if:
1. `outcome` in {SUCCESS, FAILURE} (excludes EXCLUDED and PENDING), AND
2. `excluded` = False (excludes TARGET_MISMATCH pairs where the linked
   drug does not actually target the pQTL gene)

This yields **N = 195** evaluable pairs (63 SUCCESS, 132 FAILURE).

The 4 pairs with outcome=SUCCESS AND excluded=True are:
- RET|Lung cancer: alectinib targets ALK, not RET
- RRM2B|Lung cancer: gemcitabine targets RRM1/RRM2, not RRM2B
- RRM2B|Ovarian cancer: same
- RRM2B|Pancreatic cancer: same

These TARGET_MISMATCH exclusions carry forward from V3.3 unchanged.
The base rate for the evaluable set is 63/195 = 32.3% SUCCESS.

## MR classifier threshold (pre-locked, unchanged from V3.3)

The MR classifier is: **Wald ratio p < 0.05** → predict CAUSAL →
predict SUCCESS. This threshold was pre-registered in V3.1 and has
not changed. The classifier has zero tunable parameters: no
threshold selection, no model fitting, no optimization against
outcome labels. Stage C MUST NOT tune the p-value cutoff or any
other parameter to maximize BA — doing so would reintroduce the
circularity the outcome-blind freeze was built to prevent.

## Pre-specified analyses (V3.4)

1. Primary: BA with 10,000-iteration bootstrap 90% CI, Phase III only
2. Mechanism stratification (abundance vs activity) — CO-PRIMARY
3. Gene-level deduplication — CO-PRIMARY (upgraded from V3.3 post-hoc)
4. Exact binomial tests on sensitivity and specificity
5. Continuous ROC/AUC with permutation p
6. Dipyridamole artifact removal
7. TARGET_MISMATCH inclusion
8. Phase II tier BA (secondary, reported separately)
9. Instrument source stratification (EpiGraphDB vs UKB-PPP) — SECONDARY
10. Sample overlap sensitivity (exclude overlap-flagged pairs) — SECONDARY
11. Oncology vs non-oncology stratification — SECONDARY
12. Exclude single-pair diseases (autism, hypercholesterolemia) from
    any per-disease figure — they are kept in pooled but not estimable
    per-disease

## Stop conditions

Halt and reassess if:
- Newly-rescued proteins (C - B) < 10
- Any cross-platform beta sign-flip in overlap test
- Adjudication Cohen's kappa < 0.8 on blind re-check of ~15% sample
- Known-answer trap pair is mislabeled by adjudicator
- Final Phase-III-only N < 70

## Prior version hashes

- V1: f36855865b693d96d14324b44392af7cedcf29f81658b83decd517babb2f25ab
- V2: dcfa5b2ed5dfa24c32d11e9871de8ac02a916f2e4bcb2170097c9bba318df86d
- V3 (pre-V3.1): 39b1ee694abd2c2de261c5881b2bb2cc3203aa57005427c7505d61d18eb82ade
- V3.1: 753569cbbd3613d4deb40bc34dacec229746fc795632d1262f4f36ff0ea3a93b
- V3.3 evaluation: 8cec062fd8cad1011eb23a6bcb2a42085255a602b8ed0f0887d8e52926f46804
- V3.3 classification: 305dd26265225838744a6403aef124886a9b213db27a73b1f314237647b763a6
- V3.4 frozen candidates: c428c733aa66cbdb45ae4761424871f4dfa93e173f56ccd56f101b0dc33ebe62
