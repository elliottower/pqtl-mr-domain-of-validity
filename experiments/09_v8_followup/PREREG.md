# Where does direction-aware cis-pQTL evidence for drug-indication hypotheses disappear, and does the pilot mechanism split return under the pilot's own outcome?

Sections use the [OSF Preregistration](https://osf.io/prereg/) question titles verbatim, so
this maps onto a registration without being rewritten. A question that does not apply is
answered **N/A** with the reason, never deleted.

**Stage D results of V8 were accessed before this plan was written.** The analyses below are post
hoc, specified on 2026-10-13 after the V8 stage D outputs were opened, and are fixed here before any
of their outputs exist. They do not alter, replace or satisfy any V8 confirmatory decision rule
(OSF 9tzfk; `experiments/08_mechanism_interaction/PREREG.md`, plan sha256 6cbba084…). Every result
produced under this plan is reported as post hoc, with this date.

**Declaration (repeated in every output file).** These analyses were selected after inspection of the
frozen V8 stage D results. They are post hoc and exploratory, not independent tests of V8 H1 or H4.
Their estimates and intervals do not alter the registered no-decision status of any V8 hypothesis,
and no result here is described as validating, confirming, replicating or rescuing anything.

## Research questions or hypotheses

V8 found 18 supportive hypotheses among 11,235 (stage B), too few for its registered mechanism test
(H1: no confirmatory decision). This plan asks where the evidence is lost, what the rare strong
evidence looks like, and whether the pilot's mechanism split depends on the pilot's outcome.

**Q1 (where evidence is lost).** Among completed colocalizations, is the scarcity of PP.H4 ≥ 0.8
located at the stage of regional evidence for both traits (B = PP.H3 + PP.H4 low), or at the stage of
a shared rather than distinct signal (C = PP.H4 / (PP.H3 + PP.H4) low given B high)? And how do both
stages vary with outcome-GWAS effective N, assay cis-signal strength, instrument source and platform, mechanism
class, and indication group?

**Q2 (the strong evidence, itemized).** Which hypotheses have PP.H4 ≥ 0.5, and what are their
targets, indications, directions and outcomes?

**Q3 (the pilot's own outcome).** Under the pilot classifier (S17: outcome-association p < 0.05 at the
sentinel, direction ignored), is the evidence-by-mechanism interaction for Phase III success
different from the interaction V8 estimated for Phase II → III advancement?

**Q4 (comparison with locus-level genetic support).** On the same hypotheses, how does strict
direction-aware cis-pQTL support (S) relate to Open Targets GWAS credible-set genetic support (G),
and does S carry information about advancement beyond G?

**Q5 (one relaxed threshold).** With PP.H4 ≥ 0.5 in place of 0.8 and the frozen direction rule
unchanged, how many hypotheses are supportive, and what do the V8 H1 and H2 models estimate?

**Q6 (information needed).** How many supportive hypotheses and supportive genes per mechanism
stratum would the registered H1 rule need to pass its gates and reach a decision with 80%
probability, at odds-ratio ratios 1, 1.5, 2 and 3?

**Q7 (the registered Cox secondary, checked).** Do the V8 time-to-Phase III estimates hold when the
proportional-hazards assumption is tested and programs that end without Phase III are treated as a
competing outcome rather than as censored?

No question here is confirmatory. There are no directional predictions; each answer is an estimate
with its interval. Two possible patterns of particular substantive interest, reported on the same
footing as all alternatives, are whether Q1 locates evidence loss primarily at low outcome-GWAS
information and whether Q4 estimates information from S beyond G.

## Foreknowledge of data or evidence

Authors have observed the data, but have not performed the proposed analyses

## Explanation of foreknowledge and managing unintended influences

The following has been seen in full: every V8 stage A–D output (`experiments/08_mechanism_interaction/
stages/{A,B,C,D}/output`, sealed), including the H1 no-decision, every fit under all priors, the 31
analysis sets, the frequentist checks and the Cox secondary; the coloc.susie time-limit figure; and a
post hoc count of the evidence states (`planning/.../V_round19_after_stage_d/evidence_breakdown.json`:
10,905 completed colocalizations, 37 with PP.H4 ≥ 0.8, 35 with 0.5 ≤ PP.H4 < 0.8, 744 with PP.H3 ≥ 0.8,
4,783 with PP.H0 + PP.H1 ≥ 0.8; among the 37, 18 supportive, 12 contradictory, 7 ambiguous). An
external review of these results (Perplexity, 2026-10-13) proposed the analyses below.

Not seen: S17 with Phase III success as the outcome (stage D fitted Phase III success only with the
strict evidence variable, interaction near zero under all priors); any Open Targets GWAS
credible-set evidence for these hypotheses; the PP.H4 ≥ 0.5 refits; the decomposition models; the
information simulation; any proportional-hazards diagnostic or competing-risk model.

The questions were chosen after seeing that strict evidence is scarce and that S17 shows no
interaction, so the choice of questions is not blind and the paper says so. What this plan fixes is
everything after that choice: populations, variables, thresholds (one relaxed threshold, 0.5, chosen
because PP.H4 then exceeds all other hypotheses combined, not by its result), models, priors,
bootstrap settings and every reported table. No threshold grid, no added subgroup and no change of
denominator follows the outputs.

## Study type

Non-randomized study

## Intention for causal interpretation

No causal relationship inferred

## Blinding of experimental treatments

No blinding is involved

## Additional blinding during research or analysis

N/A — the inputs are sealed V8 outputs that have already been read; outputs of this plan are produced once, by committed code, after the freeze.

## Study design

Secondary analysis of the sealed V8 tables, joined by `hypothesis_id`, plus one external table
(Open Targets 26.09 GWAS credible-set evidence) for Q4. The unit is the V8 therapeutic hypothesis.
For Q1 and Q2 a **unique genetic unit** is also defined: one colocalization computation, keyed by the
`hypotheses.csv` columns (`instrument_source`, `instrument_assay_id`, `gene_ensembl`,
`outcome_accession`); the colocalization window (±500 kb) and priors are fixed in V8, so they do not
vary. Hypotheses with the same key share one colocalization result. Every row sharing a key must agree,
in `evidence.csv`, on `coloc_run`, `not_run_reason`, `pp_h0`–`pp_h4` (absolute tolerance 1e-12),
`n_shared`, `lead_variant` and `genetic_direction`; any disagreement is a data error and stops the run. V8 records no sentinel
per hypothesis, so independent signals of one assay at one gene cannot be distinguished; the unit is
the narrowest recoverable assay-level unit. Every Q1/Q2 count is reported at both levels.

**Resampling cluster.** The gene (`gene_ensembl`). Every unique genetic unit belongs to exactly one gene,
so a gene carries all of its units and all their hypotheses; resampling genes is coarser than
resampling units and keeps the dependence among hypotheses that share a pQTL locus.

## Randomization

N/A — no assignment.

## Data collection procedures

No new data are collected except for Q4: the three parquet parts of Open Targets 26.09
`output/evidence_gwas_credible_sets/` (`part-0000{0,1,2}-c9f0fb93-8f5e-4da3-93c3-703223d7d0ae-c000.zstd.parquet`)
and the 26.09 `output/disease/00000000.parquet`, downloaded to the Modal
volume, and `evidence.md` of `opentargets/platform-docs` at commit 02978a52, with each file's sha256,
and those of the two pinned pQTL spreadsheets, written to `results/inputs_sha256.json` and logged
below the line before any join with V8 outcomes. Third-party files stay in the gitignored `inputs/`.

## Data collection procedures - File upload

N/A — the files are named and hashed above.

## Sample size

Fixed by V8: 11,235 hypotheses (stage B), 10,905 completed colocalizations, the S1 held-out set
(H1 model: 3,625 analysed aligned and blocking hypotheses; H2 model: all S1 hypotheses analysed),
and for Q3 the S1 hypotheses with a non-missing Phase III success outcome. Q2 has 72 rows. These
numbers cannot resolve small interactions; Q3 in particular is expected to be imprecise, because only
hypotheses that reached Phase III enter it.

## Sample size rationale

N/A — the data are fixed; no sample is chosen.

## Starting and stopping rules

N/A — each analysis runs once on the fixed data. A run that fails for an operational reason is rerun unchanged and logged.

## Manipulated variables

N/A — nothing is manipulated.

## Measured variables

From V8 (sealed): `evidence.csv` (stage B: coloc_run, not_run_reason, pp_h0–pp_h4, genetic_direction,
evidence_state, E, s17_sentinel_p, lead_variant); `hypotheses.csv` (stage A: gene_ensembl,
indication_id, direction, mechanism_class, instrument_source, instrument_assay_id, platform,
outcome_accession, outcome_neff, oncology, neuro_psych, heldout, in_s1); `outcomes.csv`
(stage C: advanced_24, status_24, phase3_success, approved, earliest_phase2_start,
last_phase2_end_date, time_to_phase3_days, phase3_event).

**assay_cis_signal_strength.** Not the instrument's strength (the sentinel is not recorded per
hypothesis) and not an F statistic: the largest finite −log10 p among the cis associations of the
hypothesis's exact assay (`instrument_assay_id`) in the pinned table of its source, computed without
any disease outcome. UKB-PPP: Sun et al. 2023 Supplementary Table 9 (`inputs/ukbppp/
sun2023_MOESM3_ESM.xlsx`, sheet ST9, header on row 4; rows with `cis/trans` = cis; column
`log10(p) (discovery)`), matched on the `OID…` field of `UKBPPP ProteinID` = `instrument_assay_id`.
deCODE: Ferkingstad et al. 2021 Supplementary Table 8, "Information on cis pQTL" (`inputs/decode/
ferkingstad2021_MOESM4_ESM.xlsx`, sheet ST08, header on row 3; column `-Log10(P) (unadj.)`), matched
on `SeqId` = `instrument_assay_id`. The maximum is taken over all rows of that exact assay (several rows per assay are expected:
several cis signals). A reported p of 0 (or an infinite −log10 p) is set to 300 and flagged. Reported
by source: assays matched, unmatched, with several rows, truncated at 300, and missing. Never taken from another assay or platform. INTERVAL instruments,
and any assay not found, are missing; the count by source is reported. Both spreadsheets are hashed into `results/inputs_sha256.json` with the sheet and column names used;
if a named sheet or column is absent the run stops.

**Locus-level genetic support G (Q4).** Open Targets 26.09 documents this evidence as "Evidence
scoring: Locus-to-Gene score, filtered to use scores above 0.05" (`opentargets/platform-docs`
`evidence.md` at commit 02978a52a708d8105a04c63cd89f9abe10a4d6cb, saved and hashed with the inputs). In
the 1,032,808 rows of part 2 of `evidence_gwas_credible_sets`, read on 2026-10-13, `resourceScore` equals
`score` in every row, both ranging 0.050–0.997; `resourceScore` is used, and the run stops if any row of
the three parts has `resourceScore` ≠ `score` or a missing `resourceScore`. Disease closure: D(i) =
{indication i} ∪ the `descendants` of i in the 26.09 `disease` table (`output/disease/00000000.parquet`);
ancestors are never used.

    G = 1 if |{ studyLocusId : targetId = gene_ensembl, diseaseId ∈ D(indication_id),
                               resourceScore >= 0.5 }| >= 1, else 0
    G_loci = that count (descriptive)
    G_pre = 1 if at least one such studyLocusId has a row with non-missing publicationDate strictly
            earlier than earliest_phase2_start, else 0; rows with a missing publicationDate do not
            count toward G_pre; hypotheses with a missing earliest_phase2_start have G_pre missing
            and are excluded from the G_pre analyses (counted)

Every `gene_ensembl` and every `targetId` used must match `ENSG` followed by 11 digits, with no version
suffix (all 11,235 V8 values do); any other format stops the run, and no ID is stripped or remapped.
Rows sharing (`studyLocusId`, `targetId`, `diseaseId`) count once. 0.5 is a prespecified operational L2G
threshold, not a validated natural cutoff. G is contemporaneous with the release, not with the
trials, and is reported as a retrospective association, not a target-selection validation; G_pre
measures documented publication before Phase II, not whether the evidence was used in target
selection.

## Measured variables - File upload

N/A — the variables are columns of the sealed V8 tables and the named Open Targets files.

## Indices

- B = PP.H3 + PP.H4 (regional evidence that both traits associate).
- C = PP.H4 / (PP.H3 + PP.H4), undefined when B < 0.01.
- Colocalization class. `not_run` when `coloc_run` is false (by `not_run_reason`). Otherwise, in this
  order, the first that holds: `shared` (PP.H4 ≥ 0.8), `distinct` (PP.H3 ≥ 0.8), `outcome_only` (PP.H2 ≥
  0.8), `no_outcome_signal` (PP.H0 + PP.H1 ≥ 0.8), else `unresolved`. Because the five posteriors sum to 1,
  at most one of the first four can hold, so the order only fixes the reading; all bounds are ≥.
  `unresolved` is subdivided by its largest posterior (H0 + H1 together, H2, H3, H4; ties to the
  earlier in that list). A completed run with any missing posterior stops the run (none is expected).
- C is reported only where B ≥ 0.01; where PP.H3 = 0 < PP.H4, C = 1 exactly.
- S_0.5, per hypothesis, from the sealed `evidence.csv` and `hypotheses.csv` rows:

      if not coloc_run or pp_h4 is missing or pp_h4 < 0.5:  S_0.5 = inconclusive
      else: d = direction_match(direction, genetic_direction)   # stages/B/stage_b/evidence.py
            S_0.5 = supportive if d = +1, contradictory if d = -1, inconclusive if d = 0

  This is `evidence_state(..., threshold=0.5)` exactly; no aggregation is needed because V8 has one
  colocalization per hypothesis. S_0.5 involves PP.H4 only, never the Open Targets L2G score.

## Indices - File upload

N/A — defined above.

## Statistical models

**Q1, descriptive (primary).** The colocalization classes and the distributions (median, quartiles)
of B and C, by mechanism class, instrument source, platform, oncology, neuro_psych, quintile of
outcome_neff and quintile of assay_cis_signal_strength, at the hypothesis and unique-unit levels.
Quintiles: among non-missing values, ranks with ties given the same midrank, divided by n, cut at
0.2, 0.4, 0.6 and 0.8 (a rank exactly at a cut falls in the lower quintile); missing values form a
separately counted category in the tables and are handled in the models as below. One figure: the
stacked share of classes by outcome_neff quintile.

**Q1, models (secondary).** Two fractional logistic regressions (quasi-binomial GLM, logit link,
continuous response in [0, 1]) at the hypothesis level on completed runs: (a) B, on all completed
runs; (b) C, on runs where C is defined (B ≥ 0.01). The binary counts B ≥ 0.8 and PP.H4 > PP.H3 are
descriptive only, in the tables above. Predictors: mechanism class (aligned, blocking, other), instrument source, platform, oncology,
neuro_psych, and restricted cubic splines with three knots (10th, 50th, 90th percentiles) for
log10 outcome_neff and assay_cis_signal_strength. Source and platform are partly aliased and are
reported, never interpreted causally. No variable selection. The models have no random effects.
**Rank rule:** before every fit, including every bootstrap replicate, the design-matrix rank is
computed; if source and platform make it rank deficient, platform is omitted and source kept (platform
stays in the descriptive tables); any other rank deficiency stops the run.

**Marginal standardization.** For a categorical predictor X with levels x1…xm: copy the full analysis
population once per level, set X to that level in every row of the copy, predict every row, and
average; the reported value for each level is that average, and contrasts are differences between
levels. For a spline term, the same at the 10th, 25th, 50th, 75th and 90th percentiles of the
observed variable. Never predictions at mean covariates.

**Bootstrap.** (1) Build the gene table. (2) Draw genes with replacement, as many as there are genes,
without strata. (3) Take every hypothesis of each drawn gene, as many times as the gene is drawn.
(4) Refit both models in full. (5) Recompute every standardized value and contrast. (6) 10,000
replicates, seed 20261013 (numpy `default_rng`). (7) Percentile 95% intervals. A replicate whose fit
fails or does not converge is counted and left out; if more than 5% of replicates fail for a model,
its intervals are reported as unavailable. If fewer than 50 completed runs satisfy B ≥ 0.01, model (b) is not fitted and C is reported
descriptively only.

**Q2.** One table of every hypothesis with PP.H4 ≥ 0.5: target, indication, drug program IDs, mechanism
class, direction, instrument source, assay, platform, lead variant, outcome accession and effective N,
PP.H0–PP.H4, genetic direction, evidence state, S_0.5, advanced_24, phase3_success, held-out and S1
membership, and unique-unit ID; sorted by PP.H4 descending, then gene, indication, hypothesis_id. No
row is singled out in the text except by a rule written here: none.

**Q3.** The population is the V8 S17 set exactly: `stages/D/stage_d/sets.py` `form_sets(df)["S17"]`
on the stage D join (held-out S1 hypotheses, aligned and blocking classes for the h1 model, S17
supportive when `s17_sentinel_p` < 0.05, otherwise inconclusive, so hypotheses without a usable
sentinel p stay in the denominator with S = 0). The outcome is replaced by the registered
`phase3_success` (stage C: 1 if any drug reached approval for the indication; 0 if none did, no Phase
III trial is active and the last Phase III completion or termination is ≥ 24 months before the
outcome freeze; otherwise pending and missing); rows with missing `phase3_success` (did not reach
Phase III, or pending) are left out, which is the rule V8's registered `h1_p3s` model already applies.
A hypothesis is one row whatever its number of drugs, as in V8. This set is called "S17 with Phase III
success", never S17. The model is the V8 H1 model (`stages/D/stage_d`, `build_design` kind `h1` with
y = `phase3_success`, the registered priors normal15, normal1 and student_t3, the registered sampler
and rerun sequence). Reported: the S×A posterior (median, 90% interval, Pr(β > 0), Pr(β > log 1.5), mass
within ±log 1.5), stratum odds ratios, the V8 reliability counts and prior-dominance diagnostic
(descriptive only). Beside it, at pilot scale: balanced accuracy of the S17 classifier (S17 supportive
→ predict success) for phase3_success within each class, their difference, and a 95% interval from
10,000 gene-clustered bootstrap resamples (steps 1–3 and 6–7 above, seed 20261013). In each replicate
balanced accuracy is recomputed within the aligned and blocking strata and differenced; a replicate
lacking either outcome value in either stratum gives an undefined difference and is counted, not
imputed; percentiles use the defined replicates, and if more than 5% are undefined the interval is
reported as unavailable. No pooling with the
pilot pairs. Posterior marginal quantities are those of V8's `summarize_fit` (new levels of the random
effects drawn from their fitted distribution, as registered).

**Q4.** Cross-tabulation of S × G, and of each with advanced_24, by class. The V8 H2 model with G added
(fixed effects S, G), and the V8 H1 model with G and G×A added (S, A, S×A, G, G×A), prior normal15, on
S1 with y = advanced_24. Reported: posterior of each focal coefficient, and the advancement rate in
each S × G cell with its 95% Wilson interval. The same with G_pre in place of G.

**Q5.** The V8 H1 and H2 models, prior normal15, on S1 with S_0.5 in place of S; the V8 gates and
prior-dominance rule reported descriptively. Beside them, the 0.8 results unchanged. One descriptive
figure: advancement rate by PP.H4 bin [0, 0.1), [0.1, 0.5), [0.5, 0.8), [0.8, 1] × direction state
(concordant, discordant, ambiguous), counts printed in each cell.

**Q6 (surrogate decision calibration, not power).** The registered planning simulation (`power/power_v9/power_v9_logic.py`: `simulate`, the
surrogate `fit_focal` in mode `h1`) on the structure of the realized S1 H1 set (its genes,
indications, drug programs and classes), with supportive shares set so that the expected number of
supportive hypotheses is k in each stratum: blocking share k / n_blocking and `aligned_support_mult`
= n_blocking / n_aligned, for k ∈ {14, 20, 30, 40, 60, 80, 120, 160, 240}. `simulate` draws the
evidence and the outcomes for every dataset, so realized counts vary around k; the data-generating
truth is the scenario's settings below, not a fitted model; gene, indication and drug-program
random effects are generated as in V8's planning simulation, so dependence among hypotheses of one
gene is kept. Odds-ratio ratios 1, 1.5, 2 and 3; other
settings at the V8 reference cell (random-effect SD 0.7, gene-slope SD 0.7, blocking odds ratio 1,
base advancement 0.30); 1,000 datasets per cell. Second pass, for each odds-ratio ratio: k* is the smallest grid value whose
estimated decision probability is at least 0.8; 1,000 additional datasets are run at k* and at the two
grid values immediately below and the two immediately above it (truncated at the ends of the grid),
giving 2,000 at those cells; if no k reaches 0.8, 1,000 additional datasets are run at the two largest
k. Master seed 20261013, with one non-overlapping numpy `SeedSequence` child stream per cell and pass.
Every probability is reported with its Monte Carlo 95% Wilson interval. Per dataset: the V8 reliability gate
from the simulated counts (≥ 20 supportive hypotheses and ≥ 10 supportive genes per stratum), prior
dominance (surrogate posterior SD > 0.8 × the registered prior SD), and a positive decision (gate
passed, not prior-dominated, Pr(β > 0) ≥ 0.95). A fit that fails (`fit_focal` returns None) is a non-decision and is counted. Reported per cell:
the three probabilities, the failure share, and the median 90% interval width (3.29 × posterior SD),
all called estimated decision probabilities under the surrogate; and, for each ratio, the smallest k on the grid at
which the decision probability reaches 0.8. The surrogate is the V8 planning approximation, not the
exact model. These are scenarios, not forecasts.

**Q7.** On the population of the registered Cox secondary (stage D `pipeline.cox_models`: S1
hypotheses with a dated Phase II start; H1 on aligned and blocking, H2 on all), with the registered
terms and covariates: (a) the proportional-hazards test of Grambsch and Therneau (scaled Schoenfeld
residuals against rank of time), per term and global, for the registered H1 and H2 models; (b)
cumulative incidence (Aalen–Johansen) of Phase III entry by S × class, with the competing event
defined below, plotted to 15 years with numbers at risk; (c) the registered H1 and H2 fixed effects
in a Fine–Gray subdistribution-hazard model for Phase III entry. **Implementation:** R with the
`survival` package (Debian bookworm `r-base-core` and `r-cran-survival`; versions written to
`results/environment.json`): (a) `cox.zph` on `coxph` fits of the registered terms (Efron ties),
which tests the proportional-hazards assumption for the fitted covariate effects; (b) `survfit` on a multi-state `Surv` (Aalen–
Johansen); (c) `finegray` then `coxph` with `cluster(gene_ensembl)`, Efron ties, so its variance is
**gene-clustered robust**, which differs from the registered two-way clustered Cox and is labelled so. The
Fine–Gray fit is `coxph(Surv(fgstart, fgstop, fgstatus) ~ registered terms + cluster(gene_ensembl),
weights = fgwt, ties = "efron")` on the expanded `finegray` data; omission of `fgwt` is a run-stopping
error.
**Event coding,** in this order: Phase III entry (`phase3_event` = 1) at `time_to_phase3_days`;
otherwise the competing event, for `status_24` ∈ {`no_observed_advancement`, `business_only`}, at
`last_phase2_end_date` minus `earliest_phase2_start`; otherwise censored at the outcome freeze
(2026-09-30). A competing-event status with a missing `last_phase2_end_date`, or with that date after
the freeze, is censored at the freeze; with that date before the Phase II start, it is censored at
time 0. Each of the three is counted. The competing event is an operational composite of business-only
and no-observed-advancement endings, not a biological-failure endpoint. Curves are shown to 15 years,
with numbers at risk (hypotheses not yet absorbed by either event and not censored) at 0, 2, 4, 6, 8,
10, 12 and 15 years.

## Statistical models - File upload

N/A — the models are the V8 models in `stages/D/stage_d`, named above, plus the two fractional logistic regressions and the survival models defined here.

## Transformations

log10 for outcome_neff; −log10 p capped at 300 for assay_cis_signal_strength; splines as above.
Nothing else.

## Inference criteria

| question | reported | decision rule |
|---|---|---|
| Q1 | class shares, B and C distributions, standardized probabilities with 95% bootstrap intervals | none (descriptive) |
| Q2 | the full table | none |
| Q3 | S×A posterior summaries and balanced-accuracy difference with 95% interval | none; no claim of replication or non-replication rests on a threshold |
| Q4 | S × G cells, posterior of S given G | none |
| Q5 | counts and posteriors at 0.5 beside 0.8 | none |
| Q6 | probabilities per cell, smallest k for 0.8 | none |
| Q7 | PH test statistics and p-values, cumulative-incidence curves, Fine–Gray estimates beside the registered Cox estimates | none |

**Nothing in this plan can satisfy or fail V8's H1 or H4.** The paper describes every result here as a
post hoc estimate and gives its interval; it does not describe any as significant or as confirming.

## Data inclusion and exclusion

As in V8 for each population (S1 for Q3–Q5; all stage B hypotheses for Q1–Q2). Q1 models use
completed colocalizations only; the 330 not-run hypotheses are counted by reason. Q3 uses S1
hypotheses with non-missing phase3_success in the aligned and blocking classes.

## Missing data

Missing assay_cis_signal_strength: the Q1 models are fitted on complete cases and again without assay_cis_signal_strength;
both are reported with the counts. Missing covariates otherwise as in V8 (median or modal
imputation, counted). Missing G cannot occur (absence of a row is G = 0).

## Other planned analysis

None. Anything not specified above, if done, is reported as unplanned and dated.

## Context and additional information

**Maximum claim under this plan.** Descriptive and estimation statements about where strict cis-pQTL
evidence is lost (Q1), what the strongest cases are (Q2), how the pilot classifier's mechanism
interaction looks under Phase III success (Q3), how strict cis-pQTL support relates to GWAS
credible-set support (Q4), the effect of one relaxed threshold (Q5), the information a decision
would need (Q6) and whether the registered Cox estimates survive a proportional-hazards check and a
competing-outcome model (Q7), each labelled post hoc. No confirmatory claim of any kind.

**Code and environment.** One Modal app, `pqtl-v9-followup`: the V8 stage D image definition for the
Bayesian fits (Q3–Q5) and the Python analyses (Q1, Q2, Q4 cross-tabulations), the V8 power image for
Q6, and a Debian bookworm image with `r-base-core` and `r-cran-survival` from apt for Q7. Code is committed before the freeze's
log entry `RUN_START followup token=<token>`; outputs are written to `results/` and recorded with
`results run` before they are opened.
