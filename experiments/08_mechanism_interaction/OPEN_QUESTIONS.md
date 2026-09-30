# Open questions before the V8 freeze

Started 2026-09-29; revised 2026-09-30 after rounds 2, 3 and 4. Items the repository does not settle,
or choices made in the draft `PREREG.md` that were not decided in conversation. Each needs
Elliot's call or a check before `prereg freeze`. The pre-revision copy is in
`planning/perplexity_v8_design/superseded_drafts/`.

## Needs a decision or a check

26. **H3 feasibility gate.** The draft reports H3 only with at least 20 supportive and 20
    contradictory held-out hypotheses. The number 20 is the drafter's.
27. **Target-form exceptions.** Round 2 accepted "frozen target-form exceptions". The draft
    registers none: every exception would be a per-target judgment, and S8 and S9 cover the
    boundary. Confirm or name the exceptions now.
28. **Neutralizing BINDING AGENT tokens.** The draft reads ChEMBL `mechanismOfAction` text for
    "inhibitor", "antagonist", "blocker", "neutraliz", "neutralis", "sequestr". Check the realized
    BINDING AGENT texts in the universe (mechanism text only, no outcomes) before the freeze.
29. **Splicing-flag tissues.** Liver, whole blood, and the gene's highest-median-TPM GTEx v8
    tissue; lead leafcutter intron chosen by its own sQTL p-value. Confirm the eQTL Catalogue
    serves GTEx v8 leafcutter sQTLs for these tissues by region query.
30. **Outcome source.** The draft takes eligibility and advancement from Open Targets 26.09
    `clinical_indication.maxClinicalStage` (drug × disease), with ChEMBL `max_phase_for_ind`
    as a reported cross-check; this settles item 13's per-disease question. Confirm.
31. **Time origin for time to Phase III.** Earliest dated Phase II start on ClinicalTrials.gov;
    hypotheses without one are excluded from that outcome. Many Open Targets phases come from
    sources other than ClinicalTrials.gov, so this secondary may cover a minority of
    hypotheses.
32. **Outcome-GWAS summary statistics at scale.** Selected GWAS Catalog studies need their full
    summary-statistics files (harmonised FTP) for regional extraction; with hundreds of
    indications this is a Modal job of its own, streamed per study and cut to the windows.
33. **LD reference.** 1000 Genomes phase 3 EUR (GRCh38) for the protein-altering r² proxy, the
    coverage proxy and `coloc.susie`. Reference LD differs from in-sample LD for UKB-PPP and
    deCODE; state this as a limit of S15g.
6. **Outcome rules (partly resolved).** Outcome 0 is *no observed advancement*; completed
   trials are *completed, no successor* and not read as efficacy failures; 12/36/48-month
   maturation is S14. Still the drafter's: the keyword table and its precedence
   (descriptive only now); engagement is dropped as having no systematic source.
21. **deCODE terms of use.** The form's terms have not been read. Confirm they permit this use
    and say how the data may be cited and stored, before stage B.
23. **Eldjarn ST29 platform scope.** ST29 compares Olink in UKB with SomaScan in Iceland;
    INTERVAL (SomaScan 1.1k) proteins fall in *untested*.
17. **Publication dates** for `pre_pqtl_publication` (INTERVAL 2018, deCODE 2021, UKB-PPP 2023)
    are taken from publisher records and frozen in stage A.
20. **Freeze mechanics.** `prereg freeze` needs `PREREG.md` committed alone. The plan freezes
    first; each stage commit is logged below the line. A rule a script cannot implement as
    written is an amendment, logged before its stage runs.
34. **Environment.** The all-indication count recreated the repository `.venv` under Python
    3.12; `uv sync` restores it.

42. **Trimming (not decided).** Round 3 recommended moving S6, S7, S10, S12, S16, H3, time to
    Phase III and Phase III success to the supplement, and cutting or making exploratory the
    splicing flag (S11), per-indication stratum ORs (table 11) and the Karim target table
    (table 15). Nothing has been cut or relabelled; Elliot has not decided. The plasma-
    accessibility question is now H4, and table 13 stays as its descriptive companion.
43. **Size correction at stage D.** The frequentist check's critical value comes from a
    no-effect simulation on the realized structure, which needs the realized supportive state
    per hypothesis (stage B) but no outcome. Confirm this is run after stage B and before the
    join, and logged.
44. **Gene–drug-program components.** One connected component holds 39% of H1 hypotheses
    (multi-target drugs link genes); cluster-robust inference is weak with one dominant
    cluster, which is why the Bayesian model is primary and the check is size-corrected.
45. **Within-gene diagnostic.** With a gene fixed effect on 11 crossover genes, the S×A term is
    identified only from genes with supportive hypotheses in both classes; it may be
    unestimable. It is reported as an estimate or as not estimable.
46. **Calibration omitted covariates.** The exact-model calibration fitted the registered
    random-effect structure without platform, N_eff and oncology, which the simulation does not
    generate. Adding them changes nothing in expectation but was not checked.
47. **Stage D environment.** PyMC 5.28.5 / PyTensor 2.38.2 do not install on the local Intel Mac
    (no llvmlite build); stage D runs on Linux (Modal), as the calibration did. Pin the stage D
    lockfile to the calibration image.
49. **Representative exact-model cells. RESOLVED 2026-09-30:** table inserted from `power/power_v9/results/calibrate_v9_round4_cells.json`. Previously: The round-4 NUTS cells are running; their
    table replaces `[[EXACT_MODEL_CELLS_TABLE]]` in §Sample size before the freeze.
50. **S12 count under `specificity_rank`.** The v4 recount ordered S12 by ontology depth (903
    hypotheses). S12 now uses `specificity_rank` (fewest descendants); its count is not
    recomputed and is not stated in the registration.

## Resolved

25. **H1 threshold** — kept as a graded magnitude (Pr > log 1.5) and a descriptive
    Pr(|β| < log 1.5); the "no relevant difference" confirmatory row is dropped (round 3).
35. **Drug dependence** — lead-drug rule removed; multiple-membership drug-program intercept in
    PyMC (round 3).
36. **Power with the final dependence** — power_v9 grid (324 scenarios × 1,000, size-corrected
    surrogate-test planning power) and limited exact-model validation (false-positive rates
    4.0% / 4.5% / 4.5% for H1 / H2 / H4); representative exact-model cells are item 49.
37. **Instrument lists** — single-gene / single-protein assays in every source for S1; inclusive
    lists are S20; UKB-PPP list pinned in `inputs/ukbppp/`.
39. **Conflicting mechanism rows** — program–target rows that conflict are *other*; 26
    program–targets in the universe, 1 in S1; S21 restores them.
40. **Related-indication rule** — by canonical parent-molecule program (S1 5,064; by raw
    `drugId` 5,126), applied before the held-out filter.
41. **S12 ordering** — `specificity_rank` (fewest Open Targets descendants, the same rank as
    the related-indication rule), then indication ID, then program ID; no GWAS effective N
    (round 4 unified the two specificity definitions).
48. **H4 grouping** — renamed neurologic or psychiatric (`neuro_psych_indication`), stated as
    an indication-level proxy, not a Mondo subclass relation; three-way decomposition in
    table 13; psychiatric-only indications removed in S22 (round 4).
51. **Multiplicity** — H1 sole primary confirmatory test, H4 in fixed sequence after H1, H2 key
    secondary (Elliot, round 4, option a).
52. **Feasibility gate** — split into structural non-estimability and a reliability gate that
    withholds a confirmatory decision but reports the fit (round 4).
53. **Identification** — prior-dominance rule (posterior SD > 0.8 × prior SD) added; near-zero
    slope SD does not trigger dropping the slope (round 4).
38. **OpenGWAS overlap** — kept strict (round 3): non-FinnGen OpenGWAS candidates are never *no*;
    195 S1 hypotheses use an OpenGWAS (FinnGen-mapped) outcome.
1. **Minimum detectable interaction** — power_v8c on the realized strict H1 structure is in
   §Sample size; power_v8b is the design-stage grid. The unopened
   `power/power_v8_checkpoint.jsonl` is superseded.
24. **Counts under the registered rules** — round-3 count in `feasibility/v4_round3/`: S1 5,064,
    H1 3,253 (582 aligned / 69 genes; 2,671 blocking / 170 genes), S4 5,355 (v3 counts
    superseded).
2. **H3** — demoted to key secondary with a pre-outcome feasibility gate (round 2).
3. **Primary evidence variable** — binary supportive state primary; continuous
   E = direction × PP.H4 registered as secondary (Option C), because continuous PP.H4 depends
   on the coloc priors and coverage.
4. **Model** — Bayesian hierarchical primary kept, with prior-predictive, posterior-predictive
   and prior-sensitivity checks; Firth + gene bootstrap + gene permutation beside it.
5. **Mechanism rule edges** — (a) antibodies against targets that are not blood-secreted are
   *other*, S8 moves them to function-blocking; (b) BINDING AGENT is aligned and *decrease* only
   when its mechanism text shows neutralization (item 28); (c) HPA "Secreted to blood" primary,
   UniProt Secreted is S9.
7. **Pilot definition** — the held-out set excludes all 161 keys of `classification_v5.csv`;
   the pooled set adds those keys that survive the V8 rules.
8. **Universe** — every indication in Open Targets 26.09 passing the indication rule
   (MONDO/EFO/Orphanet disease terms; HP, MP, GO, OBA, OTAR and non-disease therapeutic areas
   excluded).
9. **Sample overlap** — decided per candidate from source metadata (*yes* / *no* / *unknown*);
   only *no* is independent; unknown is excluded from S1 and returns in S4.
10. **Effective N without a source** — superseded: candidates are taken with the case and
    control counts of their source records; the old `EFFECTIVE_N` map is not used.
11. **Pilot allele harmonization** — V8 registers harmonization; bears on S2 and S17 only.
12. **Accession gate** — rewritten as an ID check against the selected tier.
13. **Open Targets stage semantics** — `clinical_indication` is drug × disease; see item 30.
14. **Genome build** — rsID matching, no liftover.
15. **Access to regional pQTL statistics** — UKB-PPP via `SYNAPSE_PAT`; deCODE download links
    received 2026-09-30 (valid 5 days, re-requested when stage B runs); annotation and
    excluded-variant files downloaded and hashed in `inputs/decode/`.
16. **Karim launched pairs** — workbook pinned by sha256 in `inputs/karim2026/`.
18. **Citation keys** — 12 verified records in `.citations/`, bib at `references.bib`.
19. **Uncommitted changes to frozen material** — the V6 hash-file edit was reverted to the
    committed text; the August `protocol/DEVIATION_LOG.md` addition is still uncommitted and
    is not part of the freeze commit.
22. **Outcome-GWAS rule inputs** — subtypes and broader phenotypes are the ontology descendants
    and ancestors of the indication ID (tiers 2 and 3), so no hand-written list is needed.
