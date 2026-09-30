# Does drug mechanism change the association between cis-pQTL evidence and clinical advancement?

**Status:** FROZEN at `b946087070e4`
**Plan sha256:** `6cbba0840b213f586ef455b33e569a84a94fc676103b74adef6c46e2e0e70040`
**Log:** 12 entries, head `d1aa4a36`
**Frozen:** 2026-09-30

Sections use the [OSF Preregistration](https://osf.io/prereg/) question titles verbatim, so
this maps onto a registration without being rewritten. A question that does not apply is
answered **N/A** with the reason, never deleted.

## Research questions or hypotheses

A cis-pQTL instruments lifelong variation in the circulating amount of one protein. Some drugs
act on that same quantity, by neutralizing the circulating protein or lowering its production;
others leave the amount unchanged and block what the protein does. The question is whether
colocalized, direction-concordant cis-pQTL evidence is associated with clinical advancement
more strongly for the first kind of drug than for the second.

The unit is a **therapeutic hypothesis**: target gene + indication + intervention direction
(decrease or increase target function) + mechanism class (§Measured variables). Drugs that
share all four are one hypothesis.

**H1.** The association between supportive cis-pQTL evidence and Phase II → Phase III
advancement is stronger for abundance-aligned hypotheses than for function-blocking hypotheses.

**H2 (key secondary).** Across all therapeutic hypotheses, supportive cis-pQTL evidence is
associated with Phase II → Phase III advancement.

**H4.** The association between supportive cis-pQTL evidence and advancement is weaker for
neurologic or psychiatric indications (Open Targets `therapeuticAreas` containing
MONDO_0005071, nervous system disorder, or MONDO_0002025, psychiatric disorder) than for other
indications, across all held-out hypotheses.

H1 carries the design and is the sole primary confirmatory test. H4 tests a second prediction
of the same account: a plasma measurement should say less about targets that act in the
brain. H4 is confirmatory only through the fixed sequence H1 → H4. H2 is a key secondary
claim: its model is fitted and reported in full, it enters no confirmatory family, and its
failure does not void H1 (§Inference criteria).

**Confirmatory tests.** H1, and H4 in fixed sequence after H1, are the only confirmatory tests.
H2, H3 and everything else in this document are registered and reported, and none of them
can satisfy or fail H1 or H4.

**The H4 grouping.** The union of MONDO_0005071 and MONDO_0002025 is not a Mondo subclass
relation: Mondo does not classify psychiatric disorder under nervous system disorder. The
union is an indication-level proxy for likely central target action, not evidence that every
program in it requires penetration of the blood–brain barrier.

**What H1 compares.** Mechanism class is almost entirely a property of the target gene: 11
genes carry hypotheses in both classes (§Sample size). H1 therefore compares the
evidence–advancement association across two largely different sets of target genes, not
mechanisms within the same target. A within-gene diagnostic on the crossover genes is
reported (§Statistical models); it has little power and does not replace H1.

**Unit of inference.** Inference is about therapeutic hypotheses in the pinned Open Targets
26.09 universe under the registered rules, not about future programs or other databases.

**Key secondary (not confirmatory).** Among hypotheses whose protein and disease signals
colocalize, advancement is compared between genetic effects pointing in the direction the
drug acts and genetic effects pointing the opposite way (§Statistical models, H3). This tests
the sign requirement built into the evidence definition.

The confirmatory set is the held-out set (§Study design). The pairs analyzed in V5/V5.1 are a
disclosed pilot and enter only the pooled secondary analysis.

## Foreknowledge of data or evidence

Authors have observed the data, but have not performed the proposed analyses

## Explanation of foreknowledge and managing unintended influences

What was seen before this document was written, stated as fact:

1. **Pilot results.** On 143 Phase III pairs (V5.1/V5.2, `results/v5_1/`), balanced accuracy
   of the p < 0.05 classifier was 0.617 for abundance-modulating (n = 65) and 0.514 for
   activity-blocking (n = 75); the difference, +0.103, had a gene-clustered 95% interval of
   [−0.042, +0.250] and permutation p = 0.214 (`results/v5_1/MECHANISM_DIFFERENCE_TEST.md`).
   Earlier versions (V3.3, V3.4, V5) on overlapping sets were also seen. Phase III outcomes
   of all 161 pairs in `results/v5/classification_v5.csv` are known.
2. **Direction-concordance run.** A direction-concordance variant of the classifier was run
   unblinded on 138 pilot pairs on 2026-08-16 and displayed in a working session
   (`experiments/PREREG_direction_concordance.md` §0). The author states he did not read the
   numbers. It touched pilot pairs only.
3. **Phase marginals of the universe.** The all-indication count
   (`feasibility/v2_all_indications/coverage_all_indications.json`, 2026-09-30) is known in
   full: 6,853 held-out gene–indication pairs across 389 genes and 692 indications, of which
   4,078 have a most advanced drug at Phase II and 2,775 at Phase III or later. This is the
   marginal of the variable from which the primary outcome is built, with no evidence state,
   mechanism label or pair identity attached. The count's output files carry no per-pair phase.
   The strict count under the registered rules
   (`feasibility/v3_strict/coverage_strict.json`, 2026-09-30) is also known: S1 has 5,134
   hypotheses (385 genes, 496 indications), of which 584 are abundance-aligned (69 genes),
   2,734 function-blocking (170 genes) and 1,816 other; S4 has 5,422. It reports no phase.
   The round-3 recount (`feasibility/v4_round3/coverage_v4.json`, 2026-09-30), which applies
   the rules in this document, is also known; its counts are in §Sample size. It reports no
   phase; its membership files carry opaque integer IDs only.
4. **Per-indication phase splits for the pilot indications.**
   `feasibility/coverage_counts.json` (2026-09-29, sha256 `eb1d0339…7e3acf`) contains, for the
   24 pilot indications, per-indication counts of pairs whose most advanced drug is Phase III
   versus Phase II only. Feasibility agents read part of it on 2026-09-29 and 2026-09-30. The
   author and the main drafting session have not read those breakdowns. Set S13 removes every
   hypothesis on those 24 indications.
5. **Published results.** Karim et al. 2026 `@karim2026proteogenomic`: pQTL-supported
   target–indication pairs launch at 4.73 times the rate of unsupported pairs; enzymes show
   the largest family-level enrichment (7/9 launched); the 23 launched pQTL-supported pairs
   come from 12 targets (APOB, CSF3, CSF3R, EGLN1, F2, FLT3, IL12B, IL4R, KIT, PCSK9,
   SERPINA1, VWF). This list was read on 2026-09-29 and is partial outcome information for
   any matching hypothesis; the workbook is pinned by sha256 in `inputs/karim2026/`. Ravarani
   et al. 2026 `@ravarani2026retrospective`: MR p < 0.05 does not enrich Phase II success at
   11,482 pairs. Minikel et al. 2024 `@minikel2024refining`: genetic support RS ≈ 2.6.
   Ferolito et al. 2026 `@ferolito2026biobanks`: among approved pairs the MR direction matches
   the mechanism of action 84% of the time.
6. **Mechanism boundary.** The abundance / function-blocking distinction, and the rule in
   §Measured variables, were written after item 1. The rule codes F2 and EGLN1 inhibitors as
   function-blocking knowing from item 5 that both were launched with pQTL support.
7. **Power simulations.** Perplexity executed `power/power_v8.py` (the 24-indication design)
   before the script was reviewed; its summary is in
   `power/power_v8_perplexity_run_2026-09-30.json` and was read. `power/power_v8b.py` and
   `power/power_v8c.py` were run and read on 2026-09-30; both use no study data. The
   power_v9 grid and the exact-model calibration (`power/power_v9/`) were run on Modal and
   read on 2026-09-30; they simulate outcomes on the opaque v4 membership structure and use no
   outcome, evidence state or phase (§Sample size).
8. **Pilot outcome lines seen by an agent.** On 2026-09-30 a feasibility agent printed lines
   ~355–400 of `protocol/DEVIATION_LOG.md`, which carry V3.3/V3.4 pilot results (balanced
   accuracies, the 195-pair outcome split, four excluded pilot pairs). These concern pilot
   pairs only and are covered by item 1.
9. **Reviews.** Five reviewers' reports on a related paper (Frontiers MS 1933481) and four
   Perplexity reviews of the V8 plan and drafts
   (`planning/perplexity_v8_design/01_ADJUDICATION.md`) shaped this design.

No evidence state, colocalization result, mechanism label under the V8 rule, or outcome for
any held-out hypothesis has been computed.

The pilot supplied the hypothesis; it does not enter the confirmatory test. The held-out set
is every eligible hypothesis whose gene+indication key is absent from
`results/v5/classification_v5.csv` and `results/v5_1/classification_v5_1.csv`.

Mechanism labels come from a written rule over database fields (Open Targets / ChEMBL molecule
type, action type and mechanism text; Human Protein Atlas secretome location), applied by
script. No person or model assigns a label. The two cases item 6 names are coded against H1.

Work runs in four stages, each committed with its output hashed before the next starts:
(A) enumeration, eligibility, mechanism labels, outcome-GWAS selection and instrument
selection; (B) colocalization, harmonization and evidence states, by a script that does not
read trial data; (C) outcome coding, by a script that does not read stage B output; (D) the
join and every analysis below. The plan is frozen first; each stage's script and output are
committed afterwards and logged below the line, and no stage D script runs before stages A–C
are logged.

**Phase blinding in stage A.** Eligibility requires a drug at Phase II or later, so stage A
reads the phase field. The stage A script writes only a boolean `eligible` per hypothesis and
never writes or prints phase; stage C recomputes phase from the same pinned table.

Items 3 and 4 bear on the base rate of the primary outcome, not on evidence or mechanism. Set
S13 removes every held-out hypothesis matching a Karim launched pair (item 5) and every
hypothesis on the 24 pilot indications (item 4).

## Study type

Non-randomized study

## Intention for causal interpretation

No causal relationship inferred

## Blinding of experimental treatments

No blinding is involved

## Additional blinding during research or analysis

Stage separation (§Explanation of foreknowledge): the evidence script (B) never reads trial
records and the outcome script (C) never reads stage B output. The author does not open stage
B or C output before both are committed and hashed. Mechanism labels are assigned in stage A,
before either exists.

**Never read, by any stage:** Supplementary Tables ST16, ST17, ST18, ST38 and ST39 of Eldjarn
et al. 2023 `@eldjarn2023proteomics` (protein–phenotype, protein–IBD and pQTL–disease-GWAS
associations). Stage A reads ST20 and ST29 of that workbook and no other sheet.

## Study design

**Universe.** Open Targets Platform release 26.09 `@ochoa2021opentargets`, bulk tables
`clinical_indication` (drug × disease, `maxClinicalStage`), `drug_mechanism_of_action`,
`drug_molecule` and `disease`, as downloaded and hashed in
`feasibility/v2_all_indications/coverage_all_indications.json`. A drug's targets are its own
mechanism rows, or its parent molecule's where it has none. Every indication in the table is
in scope, subject to the indication rule below.

**Indication rule.** An indication enters if its ID is in the Open Targets 26.09 disease
index, its prefix is MONDO, EFO or Orphanet, it is not itself a therapeutic-area root
(`ontology.isTherapeuticArea` false), and at least one of its `therapeuticAreas` is a disease
area, that is, not one of GO_0008150
(biological process), EFO_0000651 (phenotype), EFO_0001444 (measurement), EFO_0002571
(medical procedure) or MONDO_0005583 (non-human animal disease). An indication with no
`therapeuticAreas` fails. HP, MP, GO, OBA and OTAR terms are excluded.

**Eligible hypothesis.** A therapeutic hypothesis enters when (1) at least one of its drugs
reached Phase II or later for the indication (`maxClinicalStage` ∈ PHASE_2, PHASE_2_3,
PHASE_3, PREAPPROVAL, APPROVAL); (2) the target gene is on a source-reported cis-pQTL list
(below); and (3) the outcome-GWAS rule (§Measured variables) selects a tier-1 GWAS with no
known sample overlap with the selected instrument source.

**Instrument lists.** A gene has a cis-pQTL in a source when it is on that source's published
cis list, at the source's own significance threshold, and the assay maps to exactly one gene
and one protein accession, as implemented in `feasibility/v4_round3/count_v4.py`: UKB-PPP
`@sun2023ukbppp`, the cis-pQTL gene list (`inputs/ukbppp/ukbppp_cis_pqtl_genes.json`, pinned
by sha256) restricted to Olink assays with one HGNC symbol and one UniProt accession (isoform
suffix stripped; the `UniProt2` accession where `UniProt` holds a label); deCODE
`@ferkingstad2021decode`, the cis rows of Supplementary Table ST02 on SOMAmers mapping to one
gene; INTERVAL `@sun2018interval`, the cis rows of Supplementary Table ST4 on SOMAmers mapping
to one UniProt accession and one HGNC symbol. The V3.4 EpiGraphDB gene list carries no assay
or accession and is not used in S1. Source-native inclusive lists (every gene each source
names, multi-gene assays split) form sensitivity set S20. The sentinel and regional
statistics used in stage B come from the same source.

**Drug programs.** A drug program is the canonical parent molecule: Open Targets
`drug_molecule.parentId` where present, otherwise the `drugId`. Salts, formulations and other
child molecules count as their parent's program.

**Related indications.** For one gene, intervention direction and mechanism class, an
indication X is dropped when every Phase II+ program of the X hypothesis is also a Phase II+
program of some retained hypothesis on an ontology descendant of X. Terms are processed in
`specificity_rank` order, most specific first, so a drug program is counted at its most
specific indication. The rule is applied after the overlap step and before the held-out
filter.

**Specificity rank.** An indication's `specificity_rank` is its number of Open Targets 26.09
descendants; fewer is more specific, and ties go to the lower indication ID. This one rank is
used wherever this document orders indications by specificity (the related-indication rule
and S12).

**Sets.** The held-out set (confirmatory) is every eligible hypothesis whose gene+indication
key is absent from `results/v5/classification_v5.csv` and `results/v5_1/classification_v5_1.csv`.
The pooled set adds every key of `classification_v5.csv` that survives the V8 rules,
re-derived under them. S4 is S1 plus the held-out hypotheses whose selected candidate is tier 1
with overlap *unknown* and that survive the related-indication rule run on the tier-1
*no* + *unknown* pool. Every hypothesis carries a `heldout` flag fixed in stage A; the
remaining sets are in §Other planned analysis.

**Outcome aggregation.** A hypothesis *advanced* if any of its drugs reached Phase III or
later for the indication. Pairs with drugs in two mechanism classes form two hypotheses, each
with its own outcome. Mechanism rows on protein complexes are kept; S18 restricts to
hypotheses with a single-protein row.

**Conflicting mechanism rows.** Per drug program and target gene, exact duplicate
mechanism-of-action rows are collapsed and a class and direction are derived for each
remaining row. Rows whose direction is not ambiguous are informative. If the informative rows
agree, that class and direction are used; if they conflict, the program–target is *other*
with ambiguous direction and leaves H1. The same drug outcome never enters both H1 classes.
S21 restores the conflicted rows under the earlier rule (each row forms its own hypothesis).

**Funnel.** Stage A writes a flow table: universe rows → gene–indication pairs → each
indication-rule, instrument, outcome-GWAS-tier, overlap and related-indication exclusion with
its count and reason → hypotheses. It is committed before stage B runs.

## Randomization

N/A — observational; nothing is assigned.

## Data collection procedures

| source | use | pin |
|---|---|---|
| Open Targets Platform 26.09 bulk tables | universe, phase, drug type, action type, mechanism text, disease ontology, FinnGen study index | sha256 in `coverage_all_indications.json` |
| ChEMBL (release current at freeze) | `max_phase_for_ind`, cross-check of phase | release number and dump sha256 in stage A |
| Human Protein Atlas (version in `coverage_all_indications.json`) | secretome location | sha256 in `coverage_all_indications.json` |
| UniProt 2026_03 | KW-0964 / SL-0243 Secreted (sensitivity) | sha256 in `coverage_all_indications.json` |
| UKB-PPP, Synapse `syn51365303` (European discovery) | cis-pQTL regional summary statistics, one tar per protein | per-protein file sha256, stage B |
| UKB-PPP metadata, Synapse `syn51396728`, `syn51396727` | Olink assay → gene / UniProt map; GRCh38 position → rsID maps | sha256, stage A |
| deCODE, Ferkingstad 2021, non-normalized per-protein files | cis-pQTL regional summary statistics | per-protein file sha256, stage B |
| deCODE, Ferkingstad 2021 `assocvariants.annotated.txt.gz`, `assocvariants.excluded.txt.gz`, readme | effect-allele frequency; multiallelic correction; excluded variants | sha256 in `inputs/decode/` |
| INTERVAL, OpenGWAS `prot-a-*` | cis-pQTL regional summary statistics | retrieval log, stage B |
| Eldjarn 2023 supplementary workbook, sheets ST20 and ST29 only | platform concordance; protein-altering-variant flags | workbook sha256, stage A |
| GWAS Catalog r2026-09-13 | outcome GWAS candidates (mapped trait URI, ancestry, cases, controls) and full summary statistics | sha256 in `coverage_all_indications.json`; per-study file sha256, stage B |
| FinnGen R12 | outcome GWAS candidates and summary statistics | per-endpoint file sha256, stage B |
| IEU OpenGWAS | outcome GWAS candidates, regional associations, `gwasinfo` | listing sha256 in `coverage_all_indications.json`; retrieval log, stage B |
| deCODE outcome GWAS | outcome GWAS candidates, only where the rule selects them | sha256, stage A |
| eQTL Catalogue API (GTEx v8 eQTL and leafcutter sQTL) | splicing flag | response cache, stage B |
| 1000 Genomes phase 3 EUR, GRCh38 | LD for the protein-altering flag and `coloc.susie` | file sha256, stage B |
| ClinicalTrials.gov API v2 | trial status, dates, `why_stopped` | snapshot of every retrieved record, stage C |
| ENSEMBL VEP (source build) | protein-altering consequence of lead variants | response cache, stage B |

Regional pQTL and outcome statistics are streamed on Modal, filtered to the ±1 Mb window
around each instrument (the ±500 kb primary window is a subset), and written per protein with
a checkpoint after each; the full files are not stored.

**Data manifest.** Stage A writes `manifest.tsv` with one row per input file: source,
release or version, URL or accession, download date and sha256. Stage B appends one row per
regional extract. The manifest is committed; raw and extracted data files are not.

## Data collection procedures - File upload

N/A — the scripts and snapshots are committed in the repository at the logged stage commits.

## Sample size

Every eligible hypothesis enters; there is no target n.

The outcome-blind round-3 count (`feasibility/v4_round3/coverage_v4.json`) applies the rules
in this document. S1 has 5,064 held-out hypotheses (4,778 gene–indication pairs, 384 genes,
493 indications). The H1 model has 3,253: 582 abundance-aligned in 69 genes and 2,671
function-blocking in 170 genes, with 11 genes in both classes; 1,811 are *other*. H4 uses all
of S1: 698 hypotheses in 168 genes on neurologic or psychiatric indications and 4,366 on
other indications. S4 has 5,355. Restricted to blood-secreted targets on both sides (S19), H1 has
582 aligned hypotheses in 69 genes and 408 blocking in 30 genes. Active programs are removed
at stage C, so the analysed counts will be lower.

**Planning power (size-corrected surrogate test).** The figures below are size-corrected
surrogate-test planning power: they describe a two-way cluster-robust test, not the
registered Bayesian decision. `power/power_v9/` simulates outcomes on the realized membership
structure (the v4
membership files: gene, indication and every drug program of each hypothesis, opaque IDs).
The data-generating model has a gene random intercept, a gene random slope on S, an
indication random intercept and a multiple-membership drug-program intercept, with 324
scenarios: random-effect SD 0.3, 0.7, 1.0; gene-slope SD 0.3, 0.7; supportive share 3%, 6%,
10%, equal or 1.5 times higher in the aligned class; supportive-evidence odds ratio 0.8, 1.0,
1.25 in the blocking stratum; H1 odds-ratio ratio 1.5, 2, 3 (H4: neurologic-or-psychiatric to other ratio 1/1.5,
1/2, 1/3); base advancement 0.30. Each scenario has 1,000 simulated datasets, analysed by a
two-way cluster-robust logistic model (gene–drug connected components × indication). That
test rejects too often under no effect (up to 8.9% for H1, 19.4% for H2), so power is
size-corrected: the critical value is the 95th percentile of the statistic in a matched
no-effect scenario with the same variance structure. Results are in
`power/power_v9/results/power_v9_summary.json`.

| supportive share | H1, OR ratio 2 | H1, OR ratio 3 | H4, ratio 1/2 | H4, ratio 1/3 |
|---|---|---|---|---|
| 3% | 0.22 (0.17–0.36) | 0.42 (0.28–0.66) | 0.27 (0.18–0.32) | 0.47 (0.30–0.57) |
| 6% | 0.35 (0.25–0.57) | 0.60 (0.49–0.89) | 0.30 (0.28–0.54) | 0.56 (0.54–0.83) |
| 10% | 0.46 (0.27–0.75) | 0.74 (0.56–0.97) | 0.44 (0.37–0.66) | 0.80 (0.67–0.93) |

Each cell is size-corrected power in the reference scenario (random-effect SD 0.7, slope SD
0.7, equal supportive shares, blocking odds ratio 1), with the range over the 36 scenarios
of random-effect SD, slope SD, supportive-share ratio and blocking odds ratio in parentheses. At an odds-ratio ratio of 1.5, H1 power is 0.13–0.23 in
the reference scenario and at most 0.41 in any scenario. H2 power depends on whether
supportive evidence is associated with advancement in the function-blocking majority: with a
blocking-stratum odds ratio of 1.25 it is 0.25–0.96 (median 0.53–0.62 by slope SD); with 1.0,
0.05–0.48; with 0.8, at most 0.04.

**Limited exact-model validation.** On datasets simulated with no effect, the exact
registered Bayesian model (PyMC, NUTS, the model in §Statistical models without covariates;
4 chains × 1,000 draws after 1,000 tuning steps, half the registered length) declared the focal effect in 8 of 200 datasets for H1 (4.0%, 95% interval 1.7–7.7%), 9 of 200
for H2 (4.5%, 2.1–8.4%) and 9 of 200 for H4 (4.5%, 2.1–8.4%), with no divergent transitions;
the two-way cluster-robust test declared it in 7.5%, 13.0% and 6.5% of the same datasets. On
30 datasets each at H1 odds-ratio ratios of 2 and 3, the exact model declared the effect in
40% and 73% of datasets and agreed with the fast test's decision in 90% and 87% of them
(`power/power_v9/results/calibrate_v9_summary.json`, `calibrate_v9_null_size_all.json`).
Earlier design-stage grids, `power/power_v8b.py` and `power/power_v8c.py`, used a gene random
intercept only and are superseded.

**Representative exact-model cells.** Before the freeze, the exact model (NUTS, as above) was
run on the same simulated datasets as the surrogate test in these cells, all at base
advancement 0.30 and a blocking-stratum odds ratio of 1: the reference settings (random-effect
SD 0.7, gene-slope SD 0.7, supportive share 6%) at H1 odds-ratio ratios of 1 (200 datasets),
1.5, 2 and 3 (100 each); odds-ratio ratio 3 at supportive shares of 3% and 10% (100 each);
odds-ratio ratio 3 at low heterogeneity (SDs 0.3 and 0.3) and high heterogeneity (SDs 1.0 and
0.7) (100 each); and no effect at a supportive share of 3% and at high heterogeneity (200
each) (`power/power_v9/power_v9_calibrate.py`). The table covers 1,300 simulated datasets
(600 with no effect, 700 with an effect), one exact-model fit per dataset and no reruns: 260 of
them (the 200 reference no-effect datasets and the first 30 at odds-ratio ratios 2 and 3) are
those of the validation above, and 1,040 were added for this table.

| cell | OR ratio | supportive share | RE SD / slope SD | datasets | exact model declares H1 (95% interval) | surrogate, unadjusted | same decision | fits with divergences |
|---|---|---|---|---|---|---|---|---|
| reference, no effect | 1 | 6% | 0.7 / 0.7 | 200 | 4.0% (1.7–7.7%) | 7.5% | 95% | 0 |
| no effect, low share | 1 | 3% | 0.7 / 0.7 | 200 | 2.0% (0.5–5.0%) | 4.5% | 97% | 0 |
| no effect, high heterogeneity | 1 | 6% | 1.0 / 0.7 | 200 | 3.0% (1.1–6.4%) | 8.5% | 94% | 0 |
| reference | 1.5 | 6% | 0.7 / 0.7 | 100 | 18% (11–27%) | 22% | 92% | 1 |
| reference | 2 | 6% | 0.7 / 0.7 | 100 | 37% (28–47%) | 43% | 90% | 0 |
| reference | 3 | 6% | 0.7 / 0.7 | 100 | 74% (64–82%) | 74% | 88% | 0 |
| low share | 3 | 3% | 0.7 / 0.7 | 100 | 43% (33–53%) | 52% | 85% | 0 |
| high share | 3 | 10% | 0.7 / 0.7 | 100 | 83% (74–90%) | 83% | 88% | 0 |
| low heterogeneity | 3 | 6% | 0.3 / 0.3 | 100 | 80% (71–87%) | 86% | 88% | 2 |
| high heterogeneity | 3 | 6% | 1.0 / 0.7 | 100 | 55% (45–65%) | 57% | 86% | 0 |

Figures from `power/power_v9/results/calibrate_v9_round4_cells.json`; 95% intervals are
Clopper–Pearson over datasets. In the three no-effect cells the exact decision rule declared
H1 in 2.0–4.0% of datasets. At a threefold difference the exact model declared H1 in 74% of
datasets at the reference settings, against 0.60 for the size-corrected surrogate test in the
same scenario; the surrogate grid therefore understates exact-model power in the cells checked.
Three of the 700 non-null fits had divergent transitions, the case the registered rerun rule
covers (§Statistical models); none of the 600 no-effect fits did.

These runs check the null behavior of the exact Bayesian decision and its agreement with the
surrogate test in the listed configurations. They do not establish 5% frequentist error
control for the Bayesian decision across the grid, nor exact-model power over the grid; the
grid figures remain planning power for the surrogate test.

Planning power is 0.60–0.74 for a threefold H1 difference at a supportive share of 6–10% in the
reference scenario (exact model in the checked cells: 0.74 and 0.83), and lower for smaller
differences or a smaller supportive share. A
twofold difference is detected about half the time or less. The supportive share is not known
before stage B.

The Phase III success outcome (§Measured variables) is reported as an estimate without a
powered test.

## Sample size rationale

N/A — the sample is every eligible hypothesis in the pinned universe.

## Starting and stopping rules

Data are collected once. The outcome-freeze date is the date of the freeze commit; trial
records dated after it do not change any label. There are no interim analyses.

## Manipulated variables

N/A — nothing is manipulated.

## Measured variables

**Primary outcome: advancement.** 1 (*advanced*) if any drug in the hypothesis reached Phase
III or later for the indication (Open Targets 26.09 `maxClinicalStage` ∈ PHASE_3,
PREAPPROVAL, APPROVAL); 0 (*no observed advancement*) if none did, no trial of any of its
drugs for the indication is recruiting, active, not yet recruiting or enrolling by
invitation, and the last Phase II completion or termination date is at least 24 months
before the outcome-freeze date; otherwise *active*, and excluded. The maturation window is
24 months in the primary analysis and 12, 36 and 48 months in S14. ChEMBL
`max_phase_for_ind` is recorded beside each hypothesis and the agreement is reported.

**Outcome decomposition (descriptive).** Each hypothesis with no observed advancement is
coded from `why_stopped` of its terminated, suspended or withdrawn trials by the keyword
table below, applied by script; the original text is preserved beside the code. Trials that
completed with no successor carry *completed, no successor*; a completed trial has no
`why_stopped`, and the code does not read it as a failure of efficacy.

| code | `why_stopped` contains (case-insensitive) |
|---|---|
| business | business, strategic, commercial, portfolio, sponsor decision, funding, financial |
| efficacy | efficacy, futility, futile, lack of benefit, did not meet |
| safety | safety, adverse, toxicity, tolerability |
| unclassified | anything else, or empty |

A hypothesis whose every stop is *business* is excluded from the primary outcome. Where stops
carry different codes, the order efficacy > safety > unclassified > business assigns one code.
Regulatory approval (`maxClinicalStage` = APPROVAL) is recorded and descriptive.

**Secondary outcome: time to Phase III.** Time from the earliest dated Phase II start of the
hypothesis on ClinicalTrials.gov to the first Phase III start, censored at the outcome-freeze
date; hypotheses without a dated Phase II trial are excluded from this outcome and counted.

**Secondary outcome: Phase III success.** Among hypotheses that reached Phase III: 1 if any
drug reached APPROVAL for the indication; 0 if none did, no Phase III trial is active, and
the last Phase III completion or termination is at least 24 months before the outcome-freeze
date; otherwise pending, and excluded.

**Target localization.** *Blood-secreted* means Human Protein Atlas secretome location
"Secreted to blood" for the gene's Ensembl ID. UniProt Secreted (KW-0964 or SL-0243 on the
reviewed entry) is the sensitivity definition (S9). No per-target exceptions are made.

**Mechanism class.** Assigned per drug by the first matching row; a hypothesis takes the class
its drugs share.

| molecule type | action type | target localization | class |
|---|---|---|---|
| Antibody, Protein, Enzyme | INHIBITOR, ANTAGONIST, NEGATIVE ALLOSTERIC MODULATOR, or a neutralizing BINDING AGENT | blood-secreted | abundance-aligned |
| Oligonucleotide | ANTISENSE INHIBITOR, RNAI INHIBITOR | blood-secreted | abundance-aligned |
| any | DEGRADER | blood-secreted | abundance-aligned |
| Small molecule | INHIBITOR, ANTAGONIST, BLOCKER, NEGATIVE ALLOSTERIC MODULATOR, INVERSE AGONIST | any | function-blocking |
| anything else | anything else | any | other |

A **neutralizing BINDING AGENT** is a BINDING AGENT row whose `mechanismOfAction` text
contains, case-insensitively, "inhibitor", "antagonist", "blocker", "neutraliz",
"neutralis" or "sequestr"; it is treated as INHIBITOR for both mechanism class and
intervention direction. Any other BINDING AGENT row is *other* and its direction ambiguous.
Antibodies against targets that are not blood-secreted, agonists, activators and exogenous
proteins are *other*. Molecule types and action types not named in the table are *other*,
including antibody–drug conjugates, molecules of unknown type, ALLOSTERIC ANTAGONIST,
DISRUPTING AGENT and NEGATIVE MODULATOR. Under this rule PCSK9 antibodies and siRNA are abundance-aligned; F2
direct thrombin inhibitors and EGLN1 HIF-prolyl hydroxylase inhibitors are function-blocking.
The realized cross-tabulation of molecule type × action type × localization × class is
reported.

**Intervention direction.** From action type by the tiered token rule of
`experiments/PREREG_direction_concordance.md` §3.2–§3.4, adopted verbatim, with the
neutralizing-BINDING-AGENT rule above: *decrease* for BLOCKING tokens, *increase* for
ACTIVATING tokens, *ambiguous* otherwise. The rule is applied to each mechanism row's action
type; ALLOSTERIC ANTAGONIST and DISRUPTING AGENT give *decrease*, NEGATIVE MODULATOR
*ambiguous*.

**deCODE release.** Every deCODE protein uses the non-normalized per-protein summary
statistics of Ferkingstad et al. 2021. SomaScan SMP normalization changes the sign of
effects at pleiotropic loci (Eldjarn et al. 2023, Supplementary Note SN4), and the evidence
state depends on sign. Effect-allele frequency is taken from `assocvariants.annotated.txt.gz`
(not `ImpMAF`), multiallelic rows use its corrected alleles, and variants in
`assocvariants.excluded.txt.gz` are dropped. The SMP-normalized release is a sensitivity
(S16) where its per-protein files can be retrieved.

**Outcome-GWAS candidates and tiers.** A candidate is a case–control, European-ancestry GWAS
with full summary statistics and effective N = 4 / (1/n_case + 1/n_control) ≥ 2,000 from the
GWAS Catalog, FinnGen R12 or IEU OpenGWAS, mapped to ontology IDs exactly as in
`count_all_indications.py`: GWAS Catalog `MAPPED_TRAIT_URI`; FinnGen through the Open Targets
study index; OpenGWAS `ebi-a-*` through the GWAS Catalog accession, `finn-b-*` through the
FinnGen endpoint, otherwise the OpenGWAS ontology field, otherwise an exact match of the trait
label to an Open Targets disease name or exact synonym. Relative to an indication, a candidate
is **tier 1** if a mapped ID equals the indication ID or an obsolete term Open Targets lists
for it, **tier 2** if it is an ontology descendant (subtype), **tier 3** if an ancestor
(broader phenotype). Case and control counts are those of the source record, stated per
candidate in the stage A table.

**Sample overlap.** For every candidate, stage A records whether it includes UK Biobank,
Icelandic (deCODE) or INTERVAL participants: *yes* if the source metadata names the cohort,
*no* if the source lists its cohorts and none matches (FinnGen by construction), otherwise
*unknown*. An OpenGWAS `consortium` field is not a cohort list, so OpenGWAS datasets other
than `finn-b-*` are never *no*. A candidate counts as a deCODE outcome GWAS when it includes
Icelandic participants. Only *no* counts as independent.

**Instrument and outcome GWAS, selected together.** For each hypothesis, instrument sources
are tried in order of discovery N (UKB-PPP, then deCODE, then INTERVAL), restricted to
sources with a genome-wide significant cis-pQTL for the gene. For each source, the outcome
GWAS is the first candidate by: (1) tier (1 before 2 before 3); (2) independence from that
source (*no* before *unknown* before *yes*); (3) larger effective N; (4) source order GWAS
Catalog, FinnGen, OpenGWAS; (5) lower accession string. **A deCODE instrument is never paired
with a deCODE outcome GWAS.** The first source whose selected candidate is tier 1 and *no*
is used. If none is, the hypothesis takes the first source's selection, even where a later
source has a tier-1 *unknown* candidate, and is flagged
`overlap_unknown` or `sample_overlap` accordingly, and tier 2 or 3 selections are flagged
`subtype_restricted` or `phenotype_broader`. The disease association is never consulted. The
earlier hand-built map in `code/disease_gwas.py` is not used.

**Genome build.** UKB-PPP and deCODE are GRCh38; OpenGWAS is GRCh37. Variants are matched by
rsID, using the UKB-PPP rsID maps (`syn51396727`) and deCODE's own rsID column; a variant
without an rsID in either file is dropped from colocalization. No liftover is applied.

**Colocalization (primary).** `coloc.abf` `@giambartolomei2014coloc` on the ±500 kb window
around the instrument's sentinel variant, priors p₁ = p₂ = 10⁻⁴, p₁₂ = 5 × 10⁻⁶
`@wallace2020coloc`, on variants present in both after harmonization. *Not run* if fewer
than 50 shared variants or either regional file is unavailable. PP.H3 and PP.H4 are recorded
for every run.

**Coverage.** Recorded for every run: shared variants; fraction of pQTL-window variants
retained; fraction of outcome-window variants retained; whether the sentinel variant, or a
variant at r² ≥ 0.8 with it in 1000 Genomes EUR, is retained. *Low coverage* means under 50%
of pQTL-window variants retained or the sentinel and its proxies absent.

**Harmonization.** Alleles aligned to the pQTL effect allele. Palindromic variants with
minor allele frequency 0.42–0.58 are dropped; others are aligned by frequency.

**Genetic effect direction.** Sign of the harmonized ratio β_outcome / β_protein at the
variant with the largest SNP.PP.H4 that survives harmonization. Outcome traits are
case–control disease, coded so positive β means higher risk; stage A records each selected
accession's trait coding, and an accession not coded that way sets its hypotheses' direction
to ambiguous.

**Evidence state (primary, binary).**

| state | condition |
|---|---|
| supportive | PP.H4 ≥ 0.80 and the genetic effect direction matches the intervention direction (decrease requires a risk-increasing protein; increase requires a protective one) |
| contradictory | PP.H4 ≥ 0.80 and the direction is opposite |
| inconclusive | PP.H4 < 0.80, colocalization not run, or intervention direction ambiguous |

S = 1 for supportive, 0 otherwise.

**Evidence score (secondary, continuous).** E = d × PP.H4, where d = +1 if the genetic effect
direction matches the intervention direction, −1 if opposite, and 0 if the direction is
ambiguous or colocalization was not run. E runs from −1 to +1.

**Covariates.** Instrument platform (Olink or SomaScan); log₁₀ effective N of the selected
outcome GWAS, standardized to mean 0 and SD 1 over the held-out set; oncology indication
(Open Targets `therapeuticAreas` contains MONDO_0045024, cancer or benign tumor). Trial phase
is not a covariate: every hypothesis in each outcome model starts at the same phase.

**Flags.** `heldout`; `sample_overlap`; `overlap_unknown`; `subtype_restricted`;
`phenotype_broader`; `low_coverage`; `protein_altering` (lead variant, or a
variant at r² ≥ 0.8 with it in 1000 Genomes EUR, has a missense, in-frame, stop, frameshift or
splice consequence in the target gene; or Eldjarn ST29 `PAV olink` / `PAV soma` = Y for the
instrument's platform); `platform_concordant` (Eldjarn ST29
`cis pQTL on both and high correlation (> 0.5)` = Y; proteins absent from ST29 or tested on
one platform are *untested*, not discordant); `splicing_candidate` (below); `karim_launched`
(gene–indication matches a Karim et al. 2026 launched pQTL-supported pair);
`pilot_indication` (one of the 24 V3.4 indications); `pre_pqtl_publication` (earliest Phase
II start of the hypothesis precedes the online publication of every source reporting a
cis-pQTL for the gene); `neuro_psych_indication` (neurologic or psychiatric: Open Targets
`therapeuticAreas` contains MONDO_0005071, nervous system disorder, or MONDO_0002025,
psychiatric disorder; see §Research questions for what the grouping does and does not
assert).

**Splicing flag.** For each instrument, in three GTEx v8 tissues — liver, whole blood, and the
gene's highest-median-TPM tissue in GTEx v8 — taken from the eQTL Catalogue: the lead sQTL
event is the leafcutter intron with the smallest sQTL p-value in the ±500 kb window, chosen
by the sQTL statistics alone, so one event is tested per gene per tissue. `splicing_candidate`
is set when coloc.abf (primary priors) gives PP.H4 ≥ 0.80 between the pQTL and the lead sQTL
in any of the three tissues and PP.H4 < 0.50 between the pQTL and the gene's eQTL in all
three. It is an instrument-validity flag; no hypothesis about splicing is registered here.

## Measured variables - File upload

N/A — the rule tables are in §Measured variables.

## Indices

Evidence state and evidence score are the indices (§Measured variables). Relative success
(RS), the advancement rate with supportive evidence divided by the rate without, is reported
beside each odds ratio.

## Indices - File upload

N/A — no index file.

## Statistical models

**Primary (H1).** Bayesian hierarchical logistic regression on the held-out aligned and
function-blocking hypotheses:

advance ~ S + A + S×A + platform + z(log₁₀N_eff) + oncology + (1 + S ‖ gene) + (1 | indication)
+ MM(program)

A = 1 for abundance-aligned, 0 for function-blocking. `(1 + S ‖ gene)` is a gene random
intercept and an uncorrelated gene random slope on S; the slope lets the supportive-evidence
association vary by gene, which the cross-level interaction S×A otherwise treats as fixed.
`MM(program)` is a multiple-membership drug-program intercept: each hypothesis carries the
mean of the program effects of every Phase II+ program it contains (weight 1/k over its k
programs), because the outcome is 1 when any of them advanced. Gene, indication and program
effects are crossed. Priors: Normal(0, 1.5) on fixed coefficients, Normal(0, 2.5) on the
intercept, Half-Normal(1) on the four random-effect SDs; all random effects non-centered. The
model is written directly in PyMC; `power/power_v9/power_v9_calibrate.py` (`fit_exact`) is the
reference implementation, and the stage D script adds the covariates. Four chains × 2,000
draws after 2,000 tuning steps, target acceptance 0.95, versions pinned in the stage D
lockfile. R̂ < 1.01, bulk ESS > 400 and zero divergent transitions are required; a fit that
misses them is rerun at target acceptance 0.99 with 4,000 draws. **If the rerun still fails,
the gene random slope is dropped, the model is refitted, and both fits are reported with their
diagnostics.**

**Identification.** Sampler diagnostics are necessary but not sufficient. If the posterior SD
of the focal coefficient (β_S×A for H1, β_S×C for H4) exceeds 0.8 times its prior SD of 1.5,
the estimate is reported as prior-dominated and no confirmatory decision is made. A gene-slope
SD near zero is not a reason to drop the slope; the slope then shrinks toward zero and the
model is kept. The slope is dropped only under the sampler-failure rule above, and both fits
are then reported. **The reduced model is descriptive only:** its decision rule was not
covered by the exact-model validation, so when the full model fails its diagnostics after the
rerun, no confirmatory decision is made for H1 (and therefore none for H4).

**Checks, every Bayesian model.** A prior-predictive check (distribution of stratum
advancement rates under the priors) and a posterior-predictive check (observed versus
replicated advancement rate per mechanism stratum × evidence state) are reported. Prior
sensitivity: the model is refitted with Normal(0, 1) and Student-t(3, 0, 1.5) on fixed
coefficients.

**H2 (key secondary).** The same model without A and S×A, on all held-out hypotheses
including *other*.

**H4.** The same model with A replaced by C (1 for a neurologic or psychiatric indication,
0 otherwise),
on all held-out hypotheses including *other*; the focal coefficient is β_S×C.

**Within-gene diagnostic (H1).** The H1 model restricted to the 11 genes with hypotheses in
both classes, with a gene fixed effect in place of the gene random intercept. Reported as an
estimate only, as a check on whether the H1 contrast survives within the same targets; it
cannot satisfy or fail H1. If its sign is opposite to the primary estimate, the paper states
that the H1 result may reflect the composition of the two gene sets.

**H3 (key secondary).** The same model as H2 with S replaced by the three-level state
(reference: inconclusive), restricted to held-out hypotheses; the contrast is
supportive − contradictory. Reported as an estimate only if stage B gives at least 20
supportive and 20 contradictory held-out hypotheses, a pragmatic reporting threshold for a
secondary estimate rather than a validated cutoff; otherwise the counts are reported and the
model is not fitted.

**Continuous evidence (secondary).** The H1 and H2 models with S replaced by E. The binary
state stays primary because PP.H4 as a continuous quantity depends on the coloc priors and
regional coverage, which the binary threshold is less sensitive to.

**Frequentist check, every model.** Logistic regression with the same fixed effects and
two-way cluster-robust standard errors, clusters being connected components of the gene–drug-
program graph and indications; its one-sided p-value is reported with the size-corrected
critical value from the matched no-effect simulation of `power/power_v9/` for the realized
structure, because the unadjusted test rejects too often (§Sample size). Beside it, a bootstrap
resampling whole gene–drug-program components (10,000 resamples) gives an interval for the
focal coefficient. Reported beside the Bayesian result, never instead of it.

**Descriptive, per mechanism stratum.** The 2×2 of supportive vs not by advanced vs not; RS
and OR; sensitivity, specificity, PPV, NPV and balanced accuracy with 90% intervals from the
bootstrap resampling whole gene–drug-program components.

## Statistical models - File upload

N/A — the model is specified in full above; the stage D script is committed before the join.

## Transformations

log₁₀ effective N standardized to mean 0 and SD 1 over the held-out set. Allele harmonization
as in §Measured variables. No outcome transformation.

## Inference criteria

| hypothesis | status | holds when |
|---|---|---|
| H1 | primary confirmatory | posterior Pr(β_S×A > 0) ≥ 0.95 in the primary model on the held-out set, and a confirmatory decision is permitted (below) |
| H4 | confirmatory, fixed sequence | H1 holds, posterior Pr(β_S×C < 0) ≥ 0.95 on the held-out set, and a confirmatory decision is permitted |
| H2 | key secondary | no confirmatory criterion; posterior Pr(β_S > 0), median and 90% credible interval reported |
| against H1 | reported | posterior Pr(β_S×A < 0) ≥ 0.95 |

**Fixed sequence.** H1 is the sole primary confirmatory test. H4 is tested confirmatorily only
if H1 holds; if H1 does not hold, or no confirmatory H1 decision is made, the H4 estimate is
reported as a secondary estimate and makes no confirmatory claim. H2 is a key secondary: its
model and report are unchanged, and it enters no confirmatory family. No other threshold
adjustment is applied. The 0.95 posterior threshold is a Bayesian decision rule; its
frequentist error rate is checked in the limited exact-model validation (§Sample size), not
guaranteed.

**What "confirmatory" means here.** A confirmatory test in this document is a preregistered
Bayesian decision rule applied once to the held-out set; no 5% familywise error control is
claimed across the design space. Under the fixed sequence a false confirmatory claim requires
either a false H1 decision, or a true H1 decision followed by a false H4 decision, so the
probability of any false confirmatory claim is at most the larger of the two rules'
false-positive rates. In the configurations checked (§Sample size) those rates were 2.0–4.0%
for H1 and 4.5% for H4; outside the checked configurations they are not established.

**Estimability and reliability, from stage B before outcomes are joined.** H1 is structurally
non-estimable only if the interaction lacks design variation: a mechanism stratum has no
supportive or no non-supportive hypothesis. It is then reported as not estimable, with the
counts. If either stratum has fewer than 10 genes with a supportive hypothesis or fewer than
20 supportive hypotheses, the registered model is fitted and reported as a low-information
estimate with its prior sensitivity, and no confirmatory H1 decision is made. The same two
levels apply to H4, with the neurologic-or-psychiatric and other strata in place of the
mechanism strata. The counts are fixed before any outcome is joined.

**Prior dominance.** No confirmatory decision is made for a focal coefficient whose posterior
SD exceeds 0.8 times its prior SD (§Statistical models, Identification); the estimate is
reported as prior-dominated.

**Magnitude.** An odds-ratio ratio of 1.5 (|β| = log 1.5) is the smallest difference of
interest, a convention fixed before any evidence state or outcome exists. For H1 and H4 the
report gives the posterior median, the 90% credible interval, Pr(β > 0), Pr(β > log 1.5) in
the predicted direction as graded evidence of magnitude, and the posterior mass within
±log 1.5, which is not evidence of equivalence. Neither of the last two is a confirmatory
criterion; power for a 1.5 odds-ratio ratio is at most 0.41 (§Sample size).

**Marginal quantities.** For H1 and H4 the report also gives, from the posterior, the
standardized advancement probability with and without supportive evidence in each stratum
and the risk-difference interaction: the difference between strata in the supportive-minus-not
risk difference. The standardized probabilities integrate over the random-effect
distribution: for each posterior draw, gene, indication and program effects are drawn as new
levels from normal distributions with that draw's SDs, and predictions are averaged over the
observed covariate distribution. A hypothesis with k drug programs receives k newly drawn
program effects, averaged with weight 1/k, as in the model.

**H1 is reported whether or not H2's posterior clears 0.95.** If it does not, a posterior below
0.95 for H1 is described as uninformative about moderation, and a posterior at or above 0.95 as
a difference between strata without an overall association, not as moderation of a positive
association.

**No analysis in §Other planned analysis, and no H2 or H3 result, can satisfy or fail H1 or
H4.**

## Data inclusion and exclusion

**Accession gate.** Before stage B, every selected outcome accession is checked against its
live source metadata: the mapped ID must still equal the indication ID (tier 1), a
descendant (tier 2) or an ancestor (tier 3) as selected. An accession that fails is replaced
by the rule's next candidate, and the replacement is logged, before any colocalization runs.

Excluded, each counted in the funnel:

- indications failing the indication rule;
- hypotheses with no Phase II or later drug for the indication;
- targets on no source's cis-pQTL list (§Study design, Instrument lists);
- hypotheses without a tier-1 outcome GWAS, and hypotheses flagged `sample_overlap` or
  `overlap_unknown` (each retained in a sensitivity set);
- indications dropped by the related-indication rule;
- hypotheses whose drug links to a different gene than the pQTL target (target mismatch, as in
  `protocol/DEVIATION_LOG.md`);
- for the primary outcome: active hypotheses, and hypotheses whose every stop is business;
- for H1 only: hypotheses classed *other*.

## Missing data

Colocalization not run, a regional file unavailable, or an ambiguous direction gives the
inconclusive state (S = 0, E = 0); no hypothesis is dropped for missing genetic data. A
missing covariate is imputed as the held-out median (effective N) or the modal category
(platform), with the count reported. A hypothesis with no retrievable phase is excluded and
counted. Missingness of each flag and of colocalization is tabulated by mechanism class and
evidence state.

## Other planned analysis

**Analysis sets.** Every row is run with the primary model and reported in one table. No set
is promoted on the basis of its result.

| set | definition |
|---|---|
| S1 (primary) | held-out; advancement outcome; 24-month maturation |
| S2 | pooled: held-out plus pilot keys re-derived under V8 rules |
| S3 | S1 restricted to `platform_concordant` proteins |
| S4 | S1 plus `overlap_unknown` hypotheses |
| S5 | S1 without `protein_altering` |
| S6 | S1 restricted to `pre_pqtl_publication` |
| S7 | S1 with no-observed-advancement restricted to efficacy-coded stops |
| S8 | S1 with antibody antagonists of targets that are not blood-secreted moved from *other* to function-blocking |
| S9 | S1 with UniProt Secreted in place of HPA blood-secreted |
| S10 | S1 plus tier-2 and tier-3 outcome GWAS (`subtype_restricted`, `phenotype_broader`) |
| S11 | S1 without `splicing_candidate` |
| S12 | one hypothesis per drug–gene program: within each gene, S1 hypotheses ordered by the indication's `specificity_rank` (most specific first), then indication ID, then lowest program ID; a hypothesis kept only if none of its programs is in a hypothesis already kept |
| S13 | S1 without `karim_launched` and without `pilot_indication` (foreknowledge) |
| S14a–c | S1 with 12-, 36- and 48-month maturation |
| S15a–g | colocalization: p₁₂ = 10⁻⁶; p₁₂ = 10⁻⁵; ±1 Mb window; PP.H4 threshold 0.75; 0.90; `low_coverage` excluded; `coloc.susie` with 1000 Genomes EUR LD where either trait has two or more credible sets (PP.H4 taken as the maximum over credible-set pairs) |
| S16 | S1 with deCODE SMP-normalized statistics where available |
| S17 | S1 with evidence redefined as outcome-association p < 0.05 at the sentinel, direction ignored (the pilot classifier) |
| S18 | S1 restricted to hypotheses with a single-protein mechanism row |
| S19 | H1 restricted to blood-secreted targets on both sides (neutralizing biologics vs small-molecule blockers of blood-secreted proteins) |
| S20 | S1 with source-native inclusive instrument lists (multi-gene assays split; V3.4 EpiGraphDB list included) |
| S21 | S1 with conflicted program–target mechanism rows restored, each row forming its own hypothesis |
| S22 | H4 with psychiatric-only indications (MONDO_0002025 without MONDO_0005071) removed |

**Further secondary analyses.** Time to Phase III: Cox proportional hazards with the H1 and
H2 fixed effects and two-way cluster-robust standard errors (components × indication). Phase III success: the H1 and H2 models
on hypotheses that reached Phase III, reported as estimates, not tests.

**Descriptive tables and figures, all reported.**

1. The funnel, including the indication-rule and outcome-GWAS-selection steps, with each
   exclusion and its count.
2. The outcome-GWAS selection: per indication, every candidate, tier, effective N, UK Biobank,
   Icelandic and INTERVAL content, and the accession selected.
3. Hypotheses by mechanism class: n, genes, indications, platform, effective N, oncology share;
   hypotheses per gene and per indication; supportive hypotheses, supportive genes,
   indications and drug programs per stratum (the reliability-gate counts).
4. Mechanism class × molecule type × action type × localization, with every *other*
   hypothesis's reason.
5. Evidence state by mechanism class and by instrument source; colocalization not-run counts
   by reason; distributions of PP.H4, PP.H3 and coverage.
6. Missingness of flags and colocalization by mechanism class and evidence state.
7. Outcome decomposition (advanced, efficacy, safety, business, unclassified, completed no
   successor, active) by mechanism class; follow-up time from last Phase II to the
   outcome-freeze date.
8. The per-stratum 2×2, RS, OR, sensitivity, specificity, PPV, NPV and balanced accuracy.
9. Model coefficients for H1, H2, H4, H3, the within-gene diagnostic and the
   continuous-evidence models, Bayesian and frequentist; standardized advancement
   probabilities and risk-difference interactions; prior sensitivity; prior- and
   posterior-predictive checks.
10. The analysis-set table above.
11. Per-indication counts and stratum ORs, indications with fewer than 5 hypotheses pooled into
    one row.
12. Platform concordance (Eldjarn ST29) by mechanism class: concordant, discordant, untested;
    and evidence-state agreement where a gene has cis-pQTLs in both UKB-PPP and deCODE.
13. Neurologic or psychiatric indications: the per-stratum 2×2 and RS split three ways —
    neurological (MONDO_0005071), psychiatric only (MONDO_0002025 without MONDO_0005071), and
    neither.
14. Held-out vs pilot: the per-stratum 2×2 on the pilot keys under V8 rules.
15. The hypotheses at the 12 Karim launched targets present in the universe, with class,
    state and outcome.

A forest plot shows the stratum ORs for S1, S2, S3, S4, S6, S7, S8 and S13; the remaining
sets appear in the supplementary forest plot.

**eQTL arm.** A cis-eQTL arm (eQTLGen `@vosa2021eqtlgen`) is registered separately, with its
own document and freeze, before any eQTL colocalization runs. It never pools with this study.

Any analysis not listed in this document is exploratory and labelled so wherever it appears.

## Context and additional information

Prior registrations: `protocol/v1`–`v7`, frozen by SHA-256 in this repository; V7
(`protocol/v7/PRESPEC_V7.md`) was registered and not run (`protocol/v7/NOT_RUN.md`). V8
supersedes V7's design for a new dataset and does not execute V7. This document is registered
on OSF at freeze.

**Study type.** Observational: a retrospective analysis of public genetic summary statistics and drug
development records.

**Causal interpretation.** No causal claim is made about drugs or mechanisms. Mechanism class is confounded with
modality, target biology and therapeutic area, and the design cannot separate them; because
class is nearly fixed per gene, an H1 difference may reflect the composition of the two gene
sets. The claims are associations between a genetic evidence state and a development outcome.

**Blinding of experimental treatments.** No treatment is assigned; the analysis-stage blinding is described under Additional blinding during research or analysis.

**Pre-stage amendments.** A post-freeze amendment, logged below the line before the stage it
concerns, may record availability, access terms, release numbers, hashes, software pins, or
the failure of a source named here. It may not change eligibility, source priority,
thresholds, model form or variable definitions, and it may not choose a replacement source or
dataset. Where a source fails, the affected hypotheses take the consequence registered below;
any unregistered replacement is exploratory or requires a new registration. Mechanism tokens,
data-source semantics, the stage D calibration protocol and outcome definitions are fixed in
this document. The permitted items and their frozen consequences:

| item | before stage | if the source fails |
|---|---|---|
| deCODE terms of use | B | deCODE is removed as an instrument source; the registered source order is applied without it (UKB-PPP, then INTERVAL); a hypothesis with no remaining source keeps its set membership and takes the inconclusive evidence state (§Missing data) |
| ChEMBL release and dump sha256 | A | ChEMBL is a cross-check only; if unavailable, the cross-check is reported as not run and no label changes |
| eQTL Catalogue regional queries (splicing flag) | B | `splicing_candidate` is missing for the affected instruments; S11 is run on the instruments where the flag was computed, with the missing count reported, or reported as not run if the flag is missing for all |
| retrieval of full summary statistics for a selected outcome GWAS | B | colocalization is not run for the affected hypotheses, which take the inconclusive state; no other candidate is selected (the accession gate in §Data inclusion and exclusion concerns metadata mismatch before stage B and is not a retrieval rule) |
| ClinicalTrials.gov dates | C | a non-advanced hypothesis with no retrievable Phase II completion or termination date is excluded from the primary outcome as *undated* and counted; one with no dated Phase II start is excluded from time to Phase III and counted |
| stage D environment | D | the Linux image digest is pinned before stage D with the package versions of the calibration image (PyMC 5.28.5, PyTensor 2.38.2); a version change after a failed fit is a deviation, logged as such |

The 1000 Genomes phase 3 EUR reference panel for LD is fixed in this document, not an
amendable item; it approximates the LD of the UK Biobank, Icelandic and other European samples
behind the instruments and outcome GWAS, and that approximation is a limitation of S15g and
of the LD-proxy steps.

**Maximum claim under this registration.** If H1 holds, the paper can state that in the held-out hypotheses from the indications
retained by the indication and outcome-GWAS rules,
colocalized cis-pQTL evidence pointing in the direction of the drug's action was associated
with Phase II → Phase III advancement more strongly for drugs that neutralize or lower a
blood-secreted target than for small-molecule function blockers, with the posterior
odds-ratio ratio, its interval, and the posterior probability that it exceeds 1.5; if H1 and
H4 hold, that the association was weaker for neurologic or psychiatric indications than for
other indications. The overall association (H2) is reported as a key secondary estimate, not
a confirmatory claim. It cannot claim that mechanism causes the difference, that the result holds for
approval (Phase III success is an estimate, not a test), for eQTL instruments, for antibodies
against targets that are not blood-secreted (S8 is a sensitivity analysis), or outside those
indications. The supportive-versus-contradictory comparison (H3) and everything in §Other
planned analysis is secondary or descriptive, and the pilot contributes to no confirmatory
claim.

---

## Log

Append only. Never edit above the line.

The last column is what distinguishes an amendment from a deviation, so you do not have to
decide which word to use: `nothing run`, `no results seen`, `results not opened`, `results seen`.

```
2026-09-29  created                              nothing run
2026-09-29  power_v8.py started, stopped on instruction before completion; 29 of 54 design scenarios in power/power_v8_checkpoint.jsonl, not opened; no study data involved   nothing run
2026-09-29  draft revised: outcome-GWAS selection rule replaces code/disease_gwas.py; deCODE non-normalized release; Eldjarn ST29 platform flag, S3 and tables 10, 13; never-read sheets; genome-build rule; data manifest   nothing run
2026-09-30  power_v8.py executed by Perplexity in its sandbox during round-2 review, before the script was reviewed; simulation only, no study data; summary in power/power_v8_perplexity_run_2026-09-30.json   nothing run
2026-09-30  draft revised for round 2: every indication with indication, related-indication and overlap rules; ontology-tier outcome-GWAS selection; outcome 0 renamed no observed advancement; maturation, coloc, prior and deCODE sensitivities; HPA blood-secreted primary; neutralizing BINDING AGENT rule; H3 demoted to key secondary; continuous evidence secondary; splicing flag; H1 threshold OR ratio 1.5; plasma-accessibility table; foreknowledge items 3, 4, 7; power_v8b design described, result not opened   nothing run
2026-09-30  draft revised: strict S1 counts and power_v8b / power_v8c results (run and read; simulation only, no study data) in Sample size; drug random intercept; instrument lists as source-reported cis lists; recount readings written as rules; S12 redefined, S18 added; threshold justification rewritten; foreknowledge items 3, 7, 8   nothing run
2026-09-30  draft revised after round 3: v4 outcome-blind recount (single-gene/single-protein instruments, conflicting mechanism rows to other, canonical parent-molecule programs, S12 order by ontology depth); gene random slope and multiple-membership program term in PyMC; H4 (CNS) in fixed sequence after H1; feasibility gates on supportive genes; no-relevant-difference row dropped; marginal risk-difference interaction; S19-S21; power_v9 grid and exact-model calibration run and read (simulation only, no study data)   nothing run
2026-09-30  draft revised for round 4: H1 sole primary confirmatory test, H4 in fixed sequence after H1, H2 key secondary; power reworded as size-corrected surrogate-test planning power and limited exact-model validation, representative exact-model cells specified (results pending, placeholder in Sample size); H4 stratum renamed neurologic or psychiatric with three-way decomposition (table 13) and psychiatric-excluded sensitivity S22; one specificity_rank for the related-indication rule and S12; two-level estimability/reliability gate; marginal quantities integrated over the random-effect distribution; prior-dominance identification rule; pre-stage amendment list   nothing run
2026-09-30  round-4 representative exact-model cells (1,040 NUTS fits on simulated data) run and read; table inserted in Sample size; simulation only, no study data   nothing run
2026-09-30  correction to the previous line: the representative exact-model table covers 1,300 simulated datasets (600 no-effect, 700 with an effect), one fit each, no reruns; 260 came from the earlier validation and 1,040 were new. Draft revised for round 5: reduced model after persistent sampler failure is descriptive only; post-freeze amendment rule with frozen failure consequences; meaning of confirmatory (no familywise error claim; fixed-sequence bound in checked configurations); multiple-membership term in marginal predictions; H3 20/20 labelled pragmatic   nothing run
2026-09-30  frozen at b946087070e4                nothing run  ·072508c4
2026-09-30  osf draft 6abd60246cf9e1ded2dbef1b of plan 6cbba0840b213f58  nothing run  ·d1aa4a36
```
