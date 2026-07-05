# Pre-Specification & Registered-Prediction Protocol
## Zombie-Target Transportability Screen

**Timestamp (UTC):** 2026-07-04 04:57 UTC
**Status:** FROZEN before any prospective target is scored. This document defines the rules; predictions in Part C are locked as of the timestamp above.
**Purpose:** Remove investigator degrees of freedom (the criticism every reviewer raised) and register falsifiable predictions BEFORE trial readouts, so the method cannot be accused of post-hoc storytelling.

---

## Part A — Frozen Method Parameters

These are locked and inherited from prior work (`paper_v1.tex`). They are NOT to be retuned per target.

- **Test statistic:** Cochran's Q across exactly two strata — one pooled MR estimate (GEN) and one observational estimate (OBS). df = 1.
- **Scale:** both arms converted to Cohen's d via Chinn, `d = log_OR * sqrt(3)/pi`. Per-unit observational effects converted to per-SD using documented exposure SDs BEFORE conversion.
- **Thresholds (frozen):** MR null = |d| < 0.10; observational positive = |d| > 0.10; discordance = Q p < 0.05.
- **Classes:** ZOMBIE (MR null + obs positive, Q sig) / TRANSLATION_GAP (both positive, Q sig) / CONCORDANT (Q not sig, both positive) / TRUE_NULL (neither) / DISCORDANT_OTHER.
- **Validity gate (diagnostic, non-voting):** Egger intercept p < 0.05 -> PLEIOTROPIC_MR flag. A flagged verdict is reported but down-weighted; it does NOT change the class.

---

## Part B — Pre-Specified Selection Rules (the anti-cherry-pick rules)

### B1. Target inclusion
A target enters the screen iff ALL hold, decided before any Q is computed:
1. A continuous or binary exposure with a well-powered OpenGWAS GWAS (>= 10 genome-wide significant instruments at p < 5e-8, mean F >= 10).
2. A published observational effect for the exposure->disease pair from a meta-analysis or pooled cohort (not a single small study).
3. The exposure and outcome GWAS are on non-overlapping samples where feasible (to limit weak-instrument bias toward the confounded estimate).

### B2. Observational-source rule (ONE source per target, chosen by fixed priority — NOT by which gives the cleanest story)
Pick the single observational estimate by this priority ladder, stopping at the first available tier:
1. Largest individual-participant-data meta-analysis (e.g., ERFC, PSC).
2. Largest published random-effects meta-analysis reporting an OR/HR with a 95% CI.
3. Largest single prospective cohort with a 95% CI.
- SE is reconstructed from the reported 95% CI: `se = (log(hi) - log(lo)) / (2 * 1.96)`. Sources without a usable CI are EXCLUDED (never `1/sqrt(n)`).
- The chosen PMID is recorded before scoring. No swapping after seeing Q.

### B3. MR-instrument rule
- Instruments: all SNPs at p < 5e-8 for the exposure GWAS named in the frozen panel; F >= 10 filter; LD-clumped locally against the 1000G EUR panel (r2 < 0.001, 10 Mb) — no post-hoc instrument curation.
- Pooled MR estimate = inverse-variance-weighted (IVW). Egger run only for the pleiotropy flag.

### B4. Outcome-blinding rule (the core of the whole protocol)
- A target is a **RETROSPECTIVE calibration case** only if its pivotal trial outcome is already public.
- A target is a **PROSPECTIVE registered prediction** only if, as of this document's timestamp, NO pivotal Phase 3 readout exists. Its verdict is locked in Part C.
- No target may be moved between categories after scoring.

---

## Part C — Registered Prospective Predictions (LOCKED as of 2026-07-04 04:57 UTC)

Each row is a falsifiable bet made BEFORE the trial reads out. The method wins only if the aggregate hit rate beats chance across these, evaluated at readout.

| # | Target -> disease | Screen verdict | Trial / drug | Prediction (falsifiable) | Falsified if... |
|---|---|---|---|---|---|
| P1 | IL-6 pathway -> CAD | **CONCORDANT** (causal) | ziltivekimab, ZEUS trial | ZEUS meets primary MACE endpoint | ZEUS null / fails primary |
| P2 | CRP -> CAD | **ZOMBIE** (marker, not cause) | any direct anti-CRP agent | No direct CRP-lowering drug reduces MACE | A direct anti-CRP agent succeeds on MACE |
| P3 | LDL -> CAD | CONCORDANT/translation (causal) | ongoing LDL-lowering (e.g. new PCSK9/oral) | Continues to succeed | A well-powered LDL-lowering trial nulls |

**Aggregate falsification of the METHOD (not just a row):**
- The method is considered FALSIFIED as a prospective tool if, across >= 8 registered predictions accrued over time, its balanced accuracy (mean of sensitivity and specificity) does not exceed 0.5 at the 90% bootstrap CI lower bound.
- Reported honestly: sensitivity and specificity SEPARATELY, never a single blended "accuracy."

**The headline contrast to feature:** P1 vs P2. CRP and IL-6 look nearly identical in observational cardiology, yet the screen calls CRP ZOMBIE and IL-6 CONCORDANT. ZEUS is the natural experiment that adjudicates them. This is the single cleanest prospective test the panel can offer.

---

## Part D — What to actually run next on OpenGWAS (scaling the prospective arm)

The 12 calibration targets have known outcomes and are done. New registerable predictions come from exposures with a live/pending trial and no Phase 3 readout yet. Run the frozen pipeline on a druggable-exposure sweep and register any target that lands ZOMBIE or CONCORDANT AND has an active trial:

1. **Lp(a) -> CAD** — strong causal MR expected; olpasiran / pelacarsen Phase 3 pending. Predict CONCORDANT -> success. (Clean expected positive.)
2. **IL-6 / IL-6R -> other indications** (e.g., depression vs CAD) — tests whether the SAME exposure is CONCORDANT in one disease and ZOMBIE in another (IL6_depression was ZOMBIE; IL6_CAD CONCORDANT). A within-exposure dissociation is a very strong result.
3. **Homocysteine -> CAD** — classic suspected ZOMBIE (observational positive, MR null, B-vitamin trials failed): use as a held-out calibration check, not a prediction.
4. **Urate -> CAD / gout-CVD** — contested; register whatever the screen says.
5. **A druggable-genome eQTL sweep** (cis-eQTL exposures -> disease) — the true scale-up: hundreds of targets, register the subset with active trials.

**Operational note:** switch MR to LOCAL LD clumping (download 1000G panel + PLINK, pass `bfile`/`plink_bin`) before the sweep — remote clumping is the ~3-min-per-query bottleneck and forbids parallelism.

---

## Part E — Integrity commitments
- This file is timestamped and frozen; Part C predictions are not edited after 2026-07-04 04:57 UTC (append-only log for new predictions, each newly timestamped).
- Every verdict reports Q with a 95% CI and the pleiotropy flag.
- Every verdict is explicitly a TWO-SOURCE (1 df) comparison; this limitation is stated in the main text.
- Sensitivity and specificity reported separately; the two working-drug zombies (amyloid, serotonin) are presented as "confounded mechanism with downstream-acting drug," not hidden.

---

## Part F — Registered Predictions, Batch 2 (LOCKED as of 2026-07-04 06:30 UTC)

Appended per Part E append-only rule. Same frozen method parameters (Part A), same selection rules (Part B).

### F1. Individual target predictions

| # | Target -> disease | Expected screen verdict | Trial / drug | Prediction (falsifiable) | Falsified if... |
|---|---|---|---|---|---|
| P4 | Lp(a) -> CAD | **CONCORDANT** (causal) | olpasiran (OCEAN(a)-Outcomes), pelacarsen (Lp(a)HORIZON) | Phase III meets primary MACE endpoint | Both Lp(a)-lowering Phase III trials null on MACE |
| P5 | urate -> CAD | **ZOMBIE** (marker, not cause) | allopurinol (ALL-HEART, already null 2022) | Screen correctly classifies as ZOMBIE; no future urate-lowering drug reduces CVD events | A urate-lowering agent succeeds on CVD MACE in Phase III |

### F2. Druggable genome sweep predictions (registered BEFORE running the sweep)

The following aggregate predictions are registered before the Open Targets sweep is executed. They are falsifiable properties of the method applied to the full druggable genome, not individual target calls.

| # | Prediction | Falsified if... |
|---|---|---|
| S1 | Among targets with APPROVED drugs for the indication, the screen classifies >60% as CONCORDANT | ≤60% of approved-drug targets are CONCORDANT |
| S2 | Among targets the screen classifies as ZOMBIE, <30% have an approved drug for that indication | ≥30% of ZOMBIE-classified targets have approved drugs |
| S3 | The screen's ZOMBIE/CONCORDANT split correlates with max clinical stage: higher stage targets are more likely CONCORDANT (Spearman rho > 0.15, p < 0.05) | rho ≤ 0.15 or p ≥ 0.05 |
| S4 | Lp(a) -> CAD is classified CONCORDANT with Q p > 0.05 (non-significant discordance) | Q p < 0.05 for Lp(a) |
| S5 | Urate -> CAD is classified ZOMBIE with Q p < 0.05 (significant discordance, MR null, obs positive) | Q p > 0.05 or MR not null for urate |
