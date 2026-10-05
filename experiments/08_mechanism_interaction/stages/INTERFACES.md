# V8 stage interfaces

The frozen plan (`../PREREG.md`, freeze commit b946087, OSF 9tzfk) is the specification. This
file fixes only what the plan leaves to implementation: the files each stage writes, their
columns, and who may read them. Where this file and the plan disagree, the plan wins and this
file is wrong.

## Rules for every stage

- Code lives in `stages/<stage>/`; pure logic in modules, Modal code only in a thin
  `modal_*.py` wrapper; tests in `stages/<stage>/tests/` on synthetic fixtures only.
- Imports at the top of the file; no `import sys`; Pydantic models for row schemas; raise
  exceptions, never print-and-continue; `uv run`; every dependency pinned with `==`.
- Every output is written to a file before anything is printed. Each output directory then gets
  `INPUTS.tsv` and `MANIFEST.tsv` in the one shared format of `run_guard/v8_manifest.py`, written
  by A, B, C and D and read by D. `INPUTS.tsv`: name, path, sha256, source, release,
  url_or_accession, download_date, one row per input, paths repo-relative
  (`experiments/08_mechanism_interaction/...`) or volume-relative (`pqtl-v8-inputs:...`), never
  absolute. `MANIFEST.tsv`: path (relative to the output directory), rows (data rows of a
  .csv/.tsv, else empty), sha256, script_sha256, utc_time, one row per output file and one for
  `INPUTS.tsv`. A stage's seal, `SEAL stage=X manifest_sha256=<h>`, is the sha256 of its
  MANIFEST.tsv; through the INPUTS.tsv row it covers the input provenance. D checks every listed
  file but never treats the inputs INPUTS.tsv names as outputs; B's and C's INPUTS.tsv must name
  A's `hypotheses.csv` sha256. Reading is strict (`v8_manifest.read_manifest`, `read_inputs`):
  every `sha256` and `script_sha256` (and every INPUTS.tsv `sha256`) is 64 lowercase hex
  characters; a listed path is relative, no component of it below the output directory is a
  symlink, and it resolves inside the resolved output directory, and the same holds for
  MANIFEST.tsv and INPUTS.tsv themselves; MANIFEST.tsv lists each path once and never itself,
  carries one `script_sha256` for the whole stage, a `rows` value that is empty or a non-negative
  integer, and a `utc_time` that is an ISO 8601 time at UTC offset zero; every row of either file
  has exactly the header's fields; INPUTS.tsv input names are non-empty and unique. The output
  directory itself must not be a symlink: `write_manifest` refuses a symlinked root before writing
  anything and `verify_output_dir` refuses one before reading (`v8_manifest.plain_root`).
  `verify_output_dir` checks every listed file, parses INPUTS.tsv, and refuses any file (or
  symlink) under the output directory that MANIFEST.tsv does not list, apart from MANIFEST.tsv and
  the side files the caller names; `verify_listed` checks named files only, for a stage that holds
  part of a predecessor's output. Stage D names the side files (`stage_d.guard.SIDE_FILES`):
  `unit_plan.json` in `B/output/` (written by `launch_stage_b.py plan`, and named with its sha256
  in B's INPUTS.tsv) and the repository's `.gitignore` in `C/output/`; nothing else.
- `script_sha256` is `v8_manifest.code_sha256` over the stage's source files keyed by a canonical
  relative POSIX path: for each file in the order of that path, the byte length of the path (8
  bytes, big-endian), the path, the byte length of the content, the content. Paths are relative to
  the stage directory (`run_stage_a.py`, `stage_a/build.py`, `stage_b/coloc_run.R`, ...), with the
  shared guard modules as `run_guard/v8_manifest.py` and `run_guard/v8_run_guard.py`, and stage D's
  fingerprint digest adds `power_v9/power_v9_fast.py`; the digest is the same in the repository
  and in a Modal image.
- Commit identity. A real run is started from a commit: every launcher (`modal_stage_a.py`,
  `launch_stage_b.py plan`, `collect`, `spawn` and `assemble`, `modal_stage_c.py`) first calls
  `v8_run_guard.clean_commit`, which refuses unless `git status --porcelain` is empty for
  `experiments/08_mechanism_interaction/stages` and `PREREG.md`. Each Modal image carries that
  commit as `V8_REPO_COMMIT` (`baked_commit`: empty when the tree is not clean), and every real-run
  function calls `run_commit`, which refuses in an image built from an unclean tree (stage D has no
  local launcher; its four real-run functions make this check). Each stage writes the commit to its
  run record (`run_info.json` in A, B and C: stage, run token, repository commit, code digest;
  `results.json` `meta.repo_commit` in D) and to INPUTS.tsv as the row `stage_code`: path
  `experiments/08_mechanism_interaction/stages/<X>`, sha256 the stage's `script_sha256`, source
  `git commit`, release the commit. The commit is not part of the B or D checkpoint fingerprints;
  the code digest is.
- Long loops checkpoint inside the unit of work and resume from disk. A checkpoint is resumed only
  under the fingerprint of the run that wrote it. Stage B binds each unit directory
  (`FINGERPRINT.json`, and the `fingerprint` field of `result.json`) to the sha256 of the full
  `InstrumentUnit` (which names each source file by identity: source, assay or SeqId and, for a
  deCODE file, its name, size and ETag in the pinned folder listing; never a URL or token, so a
  re-issued link or renewed token changes no fingerprint), the source pins (the deCODE annotation and excluded-variant files of
  `modal_inputs_manifest.json`, hashed in the container), the `stage_b/` code digest, the frozen
  plan hash, the tool versions read in the container when the unit starts
  (`stage_b.checkpoint.tool_versions`: Python and the Python packages, R, coloc, susieR, jsonlite,
  bcftools with its htslib, tabix) and the collect digest (`stage_b.checkpoint.collect_digest`, in
  the stage B section below), which carries the sha256 of every whole file the unit reads;
  `FINGERPRINT.json` stores those versions and the collect digest (`collect_sha256`) beside the
  fingerprint, and `result.json` carries the fingerprint and the collect digest. Stage D binds its work directory (`FINGERPRINT.json`), `plan.json`, every fit attempt,
  `final.json` and every frequentist result to the sha256 of the frozen plan hash, the sealed A, B
  and C manifest hashes, the script digest (`run_stage_d.py`, the `stage_d/` modules,
  `power_v9_fast.py`), the installed Python and package versions, the run token and the
  environment digest; fit and frequentist results also carry the sha256 of the design file they
  used. The environment digest is computed in the container (`stage_d.environment.collect`,
  hashed by `fingerprint.canonical_sha256`) over the Modal image id (`MODAL_IMAGE_ID`, where Modal
  sets it), the OS release (`/etc/os-release` PRETTY_NAME and VERSION_ID, `/etc/debian_version`),
  the machine architecture, the `libopenblas*` packages `dpkg-query` reports with version and
  status, the BLAS PyTensor links (`pytensor.config.blas__ldflags` and, for each `-l` library, its
  soname, loader-cache path, resolved file and sha256), the BLAS and LAPACK NumPy was built against
  (`numpy.show_config` Build Dependencies), and every distribution in the interpreter's
  site-packages as `name==version`. Values that differ between containers of one image or within
  one process (kernel release, CPU count, loaded libraries, distributions reached only through
  `sys.path`) are left out, so every phase in one image computes the same digest and a run resumed
  in a rebuilt image whose digest differs is refused. `FINGERPRINT.json` stores the fingerprint,
  its hashed components (`environment_sha256` among them) and the environment components. A
  mismatch raises
  (`StaleCheckpointError` in B, `StaleWork` in D); nothing is recomputed over it or resumed from it.
- Nothing runs on real study data until the OSF registration is approved and the run is logged
  below the line of `PREREG.md` (`prereg log`). Tests use synthetic data only. Every runner and
  every Modal function that reads study data calls the shared guard
  `run_guard/v8_run_guard.require_run(prereg, stage, token)` with its `--run-token`; it refuses
  unless `PREREG.md` passes `prereg check` (plan sha256
  `6cbba0840b213f586ef455b33e569a84a94fc676103b74adef6c46e2e0e70040`, log chain intact) and the log
  holds exactly one chained entry whose text is `RUN_START stage=<X> token=<token>`; stage A also
  needs OSF registration 9tzfk public with `pending_registration_approval` false.
- The stages run in the plan's order A → B → C → D, each sealed before the next starts. Before its
  RUN_START the log must hold a chained `SEAL stage=<P> manifest_sha256=<sha256>` for every
  predecessor: B needs A; C needs A and B; D needs A, B and C. The last SEAL of a stage before the
  RUN_START counts and a SEAL after it counts for nothing. `require_run(..., manifests=...)` takes
  the predecessor MANIFEST.tsv files the stage reads (B and C: A's; D: A's, B's and C's) and
  refuses unless each has the sealed sha256. The stage then checks the predecessor files it reads
  against that manifest, by sha256 and row count: B `hypotheses.csv` and `outcome_trait_coding.tsv`
  (`stage_b.launch.authorize`, at `plan`, `collect`, `spawn` and `assemble` and in every Modal call), C
  `hypotheses.csv` (`run_stage_c.run`), D every listed file of A, B and C. C needs B's seal in the
  log and opens no stage B file.
- Each Modal image bakes `PREREG.md` and the predecessor files its guard checks (B: A's
  `MANIFEST.tsv`, `hypotheses.csv`, `outcome_trait_coding.tsv`; C: A's `MANIFEST.tsv`,
  `hypotheses.csv`; D: the A, B and C output directories), so images are built after the SEAL and
  RUN_START entries are logged.
- Test suites run on Modal, on synthetic fixtures, in each stage's run image plus pytest:
  `modal run stages/<X>/modal_stage_<x>.py::tests` for A, B, C (report at `/tests/modal_tests.json`
  on the stage's volume) and D (`--what tests`, report in `D/work/`; `--what env` writes the
  environment components and digest of the image to `D/work/modal_env.json`), each under
  `uv run --project stages/<X> --with modal==1.4.3` because every wrapper imports the shared guard;
  `modal run stages/modal_shared_tests.py --what run_guard` for the shared guard and
  `--what cross_stage` for `D/tests/test_cross_stage_seal.py`, which needs all four stage projects
  in one image. A report is never overwritten: the earlier one is moved to `superseded/` first.
- Output directories `stages/<stage>/output/` hold small tables committed at the stage commit;
  raw downloads go under `experiments/08_mechanism_interaction/inputs/` (gitignored). The whole
  files stage B collects (deCODE per-SeqId files, UKB-PPP tars and rsID maps, GWAS Catalog files
  without an index) are held on the `pqtl-v8-stage-b` volume under `/stage_b/raw/` and are never
  copied into the repository; `collected_files.tsv` records each one's size and sha256. Stage B's
  regional extracts (each unit directory's `pqtl.tsv.gz`, `pqtl_smp.tsv.gz`, `outcome__*.tsv.gz`
  and `ld.npz`) hold rows of the downloaded files and are not committed either: they stay on the
  volume under `/stage_b/units/`, a local copy for `assemble` goes under the gitignored
  `inputs/stage_b/`, and `B/output/` carries only their hashes (`regional_manifest.tsv`).
- Pinned raw inputs live on the Modal volume `pqtl-v8-inputs` under `/08_mechanism_interaction/`,
  mirroring the experiment directory (`inputs/{decode,karim2026,ot_26.09,ukbppp}/`,
  `feasibility/v2_all_indications/inputs/`); `modal_inputs_manifest.json` pins every file by
  sha256 (`scripts/upload_inputs_to_modal.py` for the directories, `scripts/add_inputs_to_modal.py`
  for files added later, `scripts/modal_verify_inputs.py` to hash them on Modal). Among them are
  the five tables stage B's `plan` reads: the Sun 2023 (`inputs/ukbppp/sun2023_MOESM3_ESM.xlsx`),
  Ferkingstad 2021 (`inputs/decode/ferkingstad2021_MOESM4_ESM.xlsx`) and Sun 2018 supplementary
  workbooks, the OpenGWAS listing, and the listing of the deCODE proteomics folder
  (`inputs/decode/decode_proteomics_folder_listing_2026-09-30.json`: file name, size and ETag). Real runs mount it read-only at `/inputs`. ChEMBL 37 is read from
  `proteome-mr-claim-audit-inputs` at `/chembl_37/`, read-only, and never written.

## Who may read what

| stage | reads | never reads |
|---|---|---|
| A | Open Targets 26.09 tables, ChEMBL, HPA, UniProt, instrument lists, outcome-GWAS catalogues, Eldjarn ST20/ST29, pilot keys, Karim table | any phase value into an output (phase is read only to compute `eligible`); any disease-association statistic; Eldjarn ST16/17/18/38/39 |
| B | `A/output/hypotheses.csv`, `A/output/outcome_trait_coding.tsv`, `A/output/MANIFEST.tsv` (to check the two against A's seal), instrument and outcome regional statistics, VEP, 1000G EUR, eQTL Catalogue, Eldjarn ST20/ST29 | trial records, ClinicalTrials.gov, any phase or outcome field, `C/output/*` |
| C | `A/output/hypotheses.csv` (ids, drug programs, indication), `A/output/MANIFEST.tsv` (to check it against A's seal), Open Targets phase, ChEMBL `max_phase_for_ind`, ClinicalTrials.gov through the AACT daily export dated 2026-09-30 (`studies.txt`) | `B/output/*` (B's seal is read from the PREREG.md log, not from its manifest), any evidence or coloc value |
| D | `A/`, `B/`, `C/` outputs, only after A-C are committed and logged | nothing is off limits; runs last |

## Stage A → `A/output/`

`hypotheses.csv`, one row per therapeutic hypothesis in the full eligible universe (S1, S4,
S10, S20, S21 membership as flags so every set can be formed without re-running A):

| column | type | meaning |
|---|---|---|
| hypothesis_id | str | stable id: sha256 of `gene_ensembl|indication_id|direction|mechanism_class`, first 16 hex; rows with `variant` s8, s9 or s21 append `|s8`, `|s9` or `|s21` to the key |
| gene_symbol, gene_ensembl | str | target gene |
| indication_id, indication_name | str | Open Targets disease id and name |
| direction | decrease / increase / ambiguous | intervention direction (plan §Measured variables) |
| mechanism_class | aligned / blocking / other | plan mechanism rule |
| class_reason | str | which rule row fired, or why *other* |
| drug_program_ids | str | `;`-joined canonical parent-molecule ids, all Phase II+ programs of the hypothesis |
| eligible | bool | Phase II+ for the indication (the only phase-derived field; no phase value is written) |
| heldout | bool | gene+indication key absent from the pilot keys |
| specificity_rank | int | number of Open Targets descendants of the indication |
| instrument_source | ukbppp / decode / interval | selected by the plan's source order |
| instrument_assay_id | str | Olink assay id, SOMAmer SeqId, or INTERVAL SOMAmer |
| platform | Olink / SomaScan | |
| outcome_accession, outcome_source | str | selected outcome GWAS and catalogue |
| outcome_tier | 1 / 2 / 3 | |
| overlap | no / unknown / yes | vs the selected instrument source |
| outcome_n_case, outcome_n_control, outcome_neff | int, int, float | from the source record |
| blood_secreted_hpa, secreted_uniprot | bool | |
| neuro_psych, neuro_only, psych_only | bool | MONDO_0005071 / MONDO_0002025 membership |
| oncology | bool | MONDO_0045024 |
| pilot_indication, karim_launched | bool | foreknowledge flags |
| single_protein_row | bool | for S18 |
| conflicted_row_restored | bool | row exists only under S21's rule (restored conflicting rows); true exactly when `variant` is s21 |
| strict_instrument | bool | false only for rows present under S20's inclusive lists |
| in_s1, in_s4, in_s10, in_s12 | bool | set membership computable at stage A; `in_s1`, `in_s4`, `in_s10` mark passing the set's rules *before* the held-out filter (S1 = `in_s1 & heldout`, pooled S2 = `in_s1`), `in_s12` is formed on the held-out S1 |
| subtype_restricted, phenotype_broader, sample_overlap, overlap_unknown | bool | plan flags |
| pre_pqtl_publication_date | date | journal online publication date (Crossref published-online) of the earliest source reporting a cis-pQTL for the gene: UKB-PPP 2023-10-04, deCODE 2021-12-02, INTERVAL 2018-06-06 (`stage_a/flags.py` PUBLICATION_DATES). D sets `pre_pqtl_publication` when C's `earliest_phase2_start` is strictly before this date, which is to precede the online publication of every source reporting a cis-pQTL for the gene; missing when either date is missing |
| variant | primary / s8 / s9 / s21 | primary for every row formed under the registered rules; s8, s9 or s21 for a hypothesis that exists only when hypotheses are re-formed under that set's class table (S8, S9) or its restored-conflicts rule (S21) |
| in_s8, in_s9, in_s21 | bool | passes every S1 rule, before the held-out filter, with hypotheses re-formed under S8's table (antibody antagonists of non-blood-secreted targets are function-blocking), S9's (UniProt Secreted in place of HPA blood-secreted) or S21's rule (each conflicting program-target row forms its own hypothesis); S8 = `in_s8 & heldout`, likewise S9 and S21 |
| ukbppp_assay_ids, decode_assay_ids | str | the gene's assays on each source's list (the list the row was written under), `;`-joined, empty when none; `decode_assay_ids` also empty when the selected outcome GWAS includes Icelandic participants. Stage B uses them only for descriptive table 12 |
| s19_arm | neutralizing_biologic / small_molecule_blocker / "" | S19 arm, from the same molecule type, action type and localization fields as the class: *neutralizing_biologic* when a program's class was decided only by class-table row 1 (Antibody, Protein or Enzyme; INHIBITOR, ANTAGONIST, NEGATIVE ALLOSTERIC MODULATOR or a neutralizing BINDING AGENT; blood-secreted target), *small_molecule_blocker* when only by row 4 (Small molecule blocking row) and the target is blood-secreted; set for the hypothesis only when every one of its programs is in that arm, else empty. Aligned oligonucleotides (row 2) and degraders (row 3) are in neither arm. Variant rows use their own table's labels (S9: UniProt Secreted as the localization) |

S8, S9 and S21 re-form hypotheses. For S8 and S9 the drug-target rows are relabelled with the
set's class table; for S21 each conflicting program-target row keeps its own class and direction.
Either way the rows run through the whole stage A pipeline (conflicting-rows rule within the
table, hypothesis formation, indication rule, instruments, outcome selection, related-indication
rule), exactly as the primary. A re-formed hypothesis identical to a written row (same key,
programs, instrument and outcome selection, strict lists) is that row with `in_s8` / `in_s9` /
`in_s21` set; any other (one merged from two primary hypotheses, split from one, or sharing a
primary key with other programs) is written as its own row with `variant` s8, s9 or s21, `in_s1`,
`in_s4`, `in_s10`, `in_s12` false, so stages B and C code it like any row and it never enters
another set. No re-formed hypothesis is dropped because its key already exists.

Also: `funnel.csv` (step, remaining, excluded, reason; includes each re-formed set's funnel and its
shared and variant-only row counts, and two rows whose step starts `diagnostic:`, described below),
`run_info.json` (stage, run token, repository commit, code digest), `registered_count_check.json`
(below), `outcome_gwas_selection.csv` (per indication: every candidate, tier, n, overlap
content, selected flag — descriptive table 2), `mechanism_crosstab.csv` (descriptive table 4),
`outcome_trait_coding.tsv`, `INPUTS.tsv`, `MANIFEST.tsv`.

`outcome_trait_coding.tsv`, one row per selected outcome accession (PREREG §Genetic effect
direction: "stage A records each selected accession's trait coding"); stage B reads the first
two columns:

| column | meaning |
|---|---|
| outcome_accession | as in `hypotheses.csv` |
| risk_coded | true when the accession is a case–control GWAS (both case and control counts), whose effects are on case status, so a positive β means higher risk whatever the effect unit (log OR, or a linear model on the case indicator); false otherwise, which sets the genetic direction to 0 in stage B. The qualifying rule selects only case–control candidates, so every selected accession is true |
| outcome_source, trait, n_case, n_control | from the candidate record |
| coding_basis | the reason for `risk_coded` |

The instrument lists are those of `feasibility/v4_round3/count_v4.py`, as the plan defines them,
ported in `stage_a/instruments.py` (`ukbppp_lists`, `decode_lists`, `interval_lists`) to give the
same lists entry for entry: the workbooks are read as that script reads them (the recorded sheet
dimensions, no `reset_dimensions`), and its handling of edge values is kept and cited by line
(`UniProt2` not split on `_`; an empty ST02 gene cell read as the symbol `None`; an empty
EpiGraphDB token kept). Class and direction of a program-target come from the mechanism rows of its
Phase II+ drug × indication rows, and the class table reads the action type exactly as Open Targets
gives it, both as in that script. For deCODE a
SeqId is strict when its ST01 row (included in the analysis) has one Gene, one UniProt accession
and at most one Ensembl gene ID, and a strict SeqId lists its ST01 gene and every gene symbol its
ST02 cis rows name. The two `diagnostic:` funnel rows are outcome-blind counts on that list and
change no list (`stage_a.instruments.decode_strict_diagnostic`):
`diagnostic:decode_strict_seqids_st02_symbols_differ_from_st01_gene`, the strict SeqIds whose ST02
cis rows name gene symbols other than exactly the ST01 Gene, and
`diagnostic:decode_strict_symbols_only_through_st02_disagreement`, the gene symbols on the strict
list that are the ST01 Gene of no strict SeqId; each row's `reason` gives the denominator.

Before sealing, a run compares the outcome-blind counts PREREG §Sample size registers with the same
counts taken from `hypotheses.csv` (`stage_a/registered.py`): S1 hypotheses, gene-indication pairs,
genes and indications; H1 hypotheses, by class with genes, the genes in both classes and the *other*
hypotheses; the neurologic-or-psychiatric and other-indication hypotheses of S1; S4 hypotheses; and
the H1 hypotheses on HPA blood-secreted targets by class with genes (the count the plan states under
S19; the analysis set S19 itself is formed in stage D from `s19_arm`). `registered_count_check.json`
holds `source`, per quantity `registered`, `obtained` and `match`, and `all_match`. When any quantity
differs the run raises `RegisteredCountMismatch` after writing that file and before INPUTS.tsv and
MANIFEST.tsv, so the output is not sealed; stage A reads no outcome, so the difference is
investigated and the stage rerun under a new run token.

A run writes `run_info.json` first and refuses an output directory (or volume path) that is not
empty unless it holds the `run_info.json` of the same run token and no `MANIFEST.tsv`
(`run_stage_a.claim_output_dir`): files of another run, of no identifiable run, or of a finished
run are never overwritten.

Real runs go through `modal_stage_a.py`: inputs under `inputs/` and
`feasibility/v2_all_indications/inputs/` are read from the `pqtl-v8-inputs` volume (read-only,
mounted at `/inputs`, laid out as `/inputs/08_mechanism_interaction/<path in the experiment
directory>`); the repository files stage A reads are baked into the image. Output goes to the
`pqtl-v8-stage-a` volume.

## Stage B → `B/output/`

`evidence.csv`, one row per hypothesis in `hypotheses.csv`:

| column | meaning |
|---|---|
| hypothesis_id | from A |
| coloc_run | bool |
| not_run_reason | regional_file_unavailable / outcome_file_unavailable / fewer_than_50_shared / direction_ambiguous / "" |
| pp_h0 … pp_h4 | coloc.abf posterior probabilities (primary priors, ±500 kb) |
| n_shared, frac_pqtl_retained, frac_outcome_retained, sentinel_or_proxy_retained, low_coverage | coverage (plan §Coverage) |
| lead_variant, genetic_direction | variant with max SNP.PP.H4; +1 risk-increasing protein, −1 protective, 0 undetermined |
| evidence_state | supportive / contradictory / inconclusive |
| S, E | primary binary and secondary continuous evidence |
| protein_altering, platform_concordant, splicing_candidate | flags (platform_concordant: concordant / discordant / untested) |
| s15a_pp_h4 … s15g_pp_h4, s15f_low_coverage_excluded | sensitivity colocalizations (plan S15a–g) |
| s16_evidence_state | evidence state from the deCODE instrument's SMP-normalized statistics (S16); empty where no such colocalization exists: every instrument of another source, a deCODE instrument whose SMP-normalized file is not retrievable, and every row when no SMP-normalized listing is pinned (stage D then reports S16 as not run) |
| s17_sentinel_p | outcome-association p-value at the instrument's sentinel variant (rsID match in the outcome region, direction ignored), S17; empty where the outcome file is unavailable, lacks the sentinel, or lists it with two p-values |
| evidence_state_ukbppp, evidence_state_decode | the primary evidence rule applied with that source's instrument against the row's selected outcome GWAS (descriptive table 12, agreement where a gene has cis-pQTLs in both); equals `evidence_state` for the selected source; inconclusive where the source's regional file does not resolve; empty where the row names no assay in that source |

Also `regional_manifest.tsv` (one row per regional extract of a unit: source, protein/study, window,
variants, sha256; units that extract the same region, such as the units of an assay that serves
several genes, each give a row), `collected_files.tsv` (one row per collect record: source, key, status, name,
bytes, sha256, md5, etag, last_modified, source_url, detail, utc, layout, build, position_offset,
rsid_rule, meta_sha256, uncertainty_mode, header_sha256, validation_sha256; layout to meta_sha256 are
filled only for a GWAS Catalog file read without a harmonised copy, the last three only for a GWAS
Catalog whole file the pre-analysis validation passed, and `detail` holds the reason of an `unreadable` record; `source_url` holds only the host
of the record's address, for every source, such as `download.decode.is` or `ftp.ebi.ac.uk`, and
`synapse` for a Synapse entity, whose address has no host, while the record on the private volume
keeps the address it was written with: `stage_b.assemble.source_host`), `run_info.json` (stage, run
token, repository commit, code digest, and the tool versions the units ran under), `INPUTS.tsv`,
`MANIFEST.tsv`. For table 12, each hypothesis's outcome GWAS is also colocalized with its UKB-PPP
and deCODE instruments where it names assays in them; `unit_plan.json` maps each hypothesis to
those units (`hypothesis_source_units`). The deCODE annotation and excluded-variant files are
read from the `pqtl-v8-inputs` volume, mounted read-only.

`B/output/` holds exactly these files: `evidence.csv`, `regional_manifest.tsv`,
`collected_files.tsv`, `run_info.json`, `INPUTS.tsv`, `MANIFEST.tsv`, and the side file
`unit_plan.json`. They are tables of hashes, record metadata and derived summary results. No
regional extract and no whole source file is written there: `assemble` writes only the files
above, `verify_output_dir` (at `purge` and in stage D) refuses any other file in the directory,
and `assemble` refuses a `--units-dir` or `--collect-dir` inside the repository unless under the
gitignored `inputs/` directory (`stage_b.launch.private_copy`). The plan's data manifest rule is
the same: the manifest is committed, raw and extracted data files are not.

Stage B runs on Modal in three phases, each behind the run guard (`launch_stage_b.py`):

1. `plan` (`modal_stage_b.py::plan_remote`, `stage_b/plan.py`) reads the sealed stage A files and
   the five pinned instrument tables, each hashed against `modal_inputs_manifest.json` before it is
   read, and writes `units.jsonl` and `unit_plan.json`. A deCODE instrument is located by its SeqId
   in the pinned folder listing; a SeqId the listing does not name exactly once has no regional
   file. A unit holds identities only: `pqtl_locator` (UKB-PPP OID, deCODE file name, INTERVAL
   OpenGWAS id), `pqtl_listing` and `smp_listing` (name, size, ETag). No link or token is in
   `units.jsonl`, `unit_plan.json`, a fingerprint or a log.
   A unit is one instrument for one gene: a source, an assay and a gene, with every outcome GWAS
   that gene's hypotheses pair with the instrument. Its key is
   `<source>__<assay>__<gene_ensembl>`, each part with every character outside `[A-Za-z0-9_.-]`
   replaced by `_` (`stage_b.units.unit_key`); the key names the unit's directory
   `/stage_b/units/<unit_key>/` and is never parsed, and `unit_plan.json` gives each key's source,
   assay and gene (`unit_ids`). An assay that serves several genes (a SOMAmer for a protein
   complex, a SeqId that two genes' lists name) has one unit per gene. The sentinel is selected
   by source and assay alone, so those units carry the same sentinel, window and regional file,
   which `collect` downloads once, and differ in gene and outcomes; each runs the VEP lookup and
   the splicing flag for its own gene, and `assemble` refuses a hypothesis mapped to a unit of
   another gene.
2. `collect` (`collect_file`, `stage_b/collect.py`), one Modal call per whole file whatever the
   number of units that read it: deCODE per-SeqId files (and the SMP-normalized ones for S16 where
   a listing of that folder is pinned), UKB-PPP per-protein tars, UKB-PPP per-chromosome rsID maps,
   GWAS Catalog harmonised files without a tabix index and, for an accession without a harmonised
   file, the summary-statistics file of its study directory (below, "GWAS Catalog files without a
   harmonised copy"). A unit with a GWAS Catalog outcome also reads the rsID map of its sentinel's
   chromosome. A file is written to
   `/stage_b/raw/<source>/<name>` on the `pqtl-v8-stage-b` volume, resumed by HTTP Range, the
   volume committed every 256 MB, and gets a record `/stage_b/collect/<source>/<key>.json`:
   status, name, bytes, sha256, ETag, Last-Modified, the source address without query string or
   token, UTC time. It is accepted only with the length the source declared, the size and ETag of
   the pinned listing (deCODE), the MD5 the source declares (Synapse), and, for `.gz`, a stream
   that decompresses to its end. The record has status `absent` for a file the source does not
   hold, `unreadable` (with the reason in `detail`) for a GWAS Catalog file the source holds and
   no fixed rule reads, and `remote_indexed` for a GWAS Catalog file that has a tabix index
   (queried by region).
   A deCODE file is requested from `https://download.decode.is/s3/download?token=<token>&file=<Key>`
   (the form deCODE's download page builds), whose reply is followed whether it redirects to a
   signed address or carries one in a JSON body; the token and the signed address are never
   written. The files under `/stage_b/raw/` are working copies, held between collection and the
   stage B seal (a logged procedural deviation: the plan states that the full files are not
   stored). `launch_stage_b.py purge` (`purge_raw_files`) deletes them once stage B's SEAL is
   logged, its output verifies against the seal, every unit has its result and every planned file
   its record. The records and the regional extracts stay on the volume; the extracts stay private
   and only their hashes are published.
3. `spawn` (`run_unit`, `stage_b/pipeline.py`), one call per instrument unit, refused while a
   planned file has no record or while the pre-analysis validation report (below) is missing or
   does not cover the current GWAS Catalog collect records. Whole files are read only from the volume, after their size and
   sha256 match their record; a missing record raises `CollectError`. Queries by region stay
   remote: OpenGWAS associations, tabix on FinnGen, indexed GWAS Catalog files and the eQTL
   Catalogue (its EBI FTP paths read over HTTPS on the same host), bcftools on 1000 Genomes,
   Ensembl, GTEx. Each step of a unit is checkpointed.

Unavailable (`stage_b/remote.py`). A regional or outcome file is recorded as unavailable, and its
hypotheses take the registered consequence, only on a definitive absence (`SourceAbsent`): HTTP 404
or 410, a file its source's listing does not name, an accession with more than one harmonised file,
or a GWAS Catalog study directory with neither a harmonised file nor a summary-statistics file; or
when the file is `unreadable` (`SourceUnreadable`, a subclass): the source holds it and a fixed rule
finds it cannot be read without a guess. An accession without a harmonised file is not absent. A refused credential (401/403), a rate limit (429), a server error (5xx), any
other status, a timeout, a connection error, a failed tabix or bcftools call and a truncated or
corrupt gzip or tar stream raise `RetryableSourceError` with that class: nothing is written for the
step, the call fails and leaves a note under `/stage_b/errors/`, and the file or unit stays
unfinished until a later call finishes it. Messages carry a source label and the exception type,
never a URL.

GWAS Catalog files without a harmonised copy (`stage_b/outcome_files.py`, `stage_b/author_formats.py`).
The plan asks the GWAS Catalog for full summary statistics. An accession whose `harmonised/`
directory holds no `.h.tsv.gz` file (or does not exist) is resolved in its study directory:

- a reviewed author format, where `stage_b/author_formats.py` holds one for the accession: a fixed
  column map written from the deposited read-me and the data file's header line only, holding the
  sha256 of each; the header is read (HTTP Range, first 64 KB) and must be the reviewed one, else
  `unreadable`; an entry whose documents leave the effect allele, the effect scale, the rsID column
  or the build unstated is `unreadable` with that reason (GCST010514);
- else the study directory's one GWAS-SSF data file, `<acc>.tsv[.gz]` or
  `<acc>_build<GRCh37|GRCh38>.tsv[.gz]`, with its `-meta.yaml` and the directory's `md5sum.txt`.
  `genome_assembly` must be GRCh37 or GRCh38; `coordinate_system` 1-based, 0-based (+1 to each
  position) or absent (read as 1-based, a convention for legacy files: GWAS-SSF makes the key
  optional, and the metadata does not state it); `data_file_name` the data file; `data_file_md5sum`, where
  given, md5sum.txt's MD5. The header must hold `chromosome`, `base_pair_location`,
  `effect_allele`, `other_allele`, `standard_error` and `beta` or `odds_ratio` (both: column 4,
  the standard's effect column, decides). Any failure is `unreadable` with the reason, decided
  before a byte of data is downloaded. A downloaded data file whose MD5 is not md5sum.txt's raises
  `CollectError`: refused, not recorded. The MD5 md5sum.txt lists for the `-meta.yaml` is not a
  condition (the Catalog rewrites metadata in place without updating md5sum.txt); a difference is
  written to the record's `detail`;
- a study directory with neither is `absent`.

The file is downloaded whole and its record names how it is read: `layout` (`gwas_ssf`,
`author`), `build`, `position_offset`, `rsid_rule` and `meta_sha256`. The analyze phase reads it into
the canonical outcome table: beta, or ln(odds_ratio) where only an odds ratio is given; the
standard error as given; `effect_allele_frequency`, `p_value` (else 10^-`neg_log_10_p_value`) and
`n` where present; `#NA` missing; every field read with its surrounding whitespace removed. A row
is not read, for the first reason that applies (`outcome_files.REJECT_REASONS`), when it has
another number of fields than the header, no chromosome or position, a field that is not a number,
an effect or other allele that is not `[ACGT]+`, equal alleles, an effect allele outside the pair
(author formats), an odds ratio or standard error that is missing, infinite, zero or negative, an
ordinary p-value outside [0, 1], or a negative -log10 p-value; a missing p-value is kept as
missing. The `-meta.yaml` is read with PyYAML's SafeLoader (UTF-8 only, a byte-order mark allowed,
a repeated key refused) and its four keys checked as strings. The outcome window
is centred on the sentinel's position on the record's build (`VolumeFetcher.outcome_build`;
GRCh38 for every harmonised file, as before), which Ensembl's variation endpoint gives
(`positions`); no liftover. rsIDs (`rsid_rule`): `column` takes the `rsid` (else `rs_id`) column,
a value that is not `rs<digits>` giving none; `ukbppp_map`, for a file with neither column, takes
the rsID of the UKB-PPP rsID map row of the chromosome whose `POS19` (GRCh37) or `POS38` (GRCh38)
equals the variant's position and whose {REF, ALT} equals its {effect_allele, other_allele}, and
none where no row or more than one rsID matches. `variant_id` is never read as an rsID. An
`unreadable` record reads as `SourceUnreadable`, which the unit records as
`outcome_file_unavailable` with the reason in its detail; the record itself keeps status
`unreadable` (not `absent`) in the collect digest and in `collected_files.tsv`.

Pre-analysis validation (`stage_b/validate.py`, `modal_stage_b.py::validate_outcomes`,
`launch_stage_b.py validate`), after collect and before `spawn`: every collected GWAS Catalog
whole file (harmonised, GWAS-SSF or reviewed author format) is read whole on Modal by the reader
the analysis uses (`parsers.HarmonisedReader`, `outcome_files.SsfReader`, `outcome_files.AuthorReader`),
checkpointed every 1M rows. A `remote_indexed` file is queried by region through its tabix index
and is not validated; the report lists it with that reason. The harmonised reader applies no
allele, standard-error or p-value rule: a harmonised row is not read only for another number of
fields than the header, no chromosome or position, or a field that is not a number. `/stage_b/validation/outcome_validation.json` holds aggregate parser diagnostics only: per
file the header, rows, rows read, rows not read by reason (row-width and parse errors among them),
the fraction of rows read with an rsID, the number of rsIDs on more than one row read, and for
GCST008226 the OR/CI/SE integrity counts. Fixed rules: a file fails when its header disagrees with
the column map of its record, when fewer than 50% of its data rows are read, or (GCST008226) when
fewer than 99% of the rows with finite OR, confidence limits and SE give limits equal to
exp(ln(OR) ∓ 1.96·SE) within 1% or the printed rounding, whichever is looser. A failing file's
`collected` record is moved to `/stage_b/superseded/collect/` and replaced by an `unreadable`
record (reason in `detail`, name, size and sha256 kept), which takes the unavailable consequence.
`spawn` refuses unless the report exists, was made under the rules in force
(`stage_b.validate.RULES`), and covers every GWAS Catalog whole file of the current collect records:
each with a result for the sha256 its record names, no error, and status `unreadable` exactly where
the file failed (`stage_b.validate.require_validation`, called by `stage_b.launch.spawn_units` with
the report and records `stage_state` reads from the volume).

Uncertainty mode (`stage_b.outcome_files.choose_uncertainty_mode`). Before a GWAS-SSF or
author-format file is read, a first pass counts its `standard_error` field as missing or present
over the rows with a position and a valid pair of alleles; one mode is chosen for the whole file from
the header and those two counts, never row by row: `native_se` where it is present on every such
row; where it is missing on every such row, `or_ci_derived_se` for an odds-ratio file with
`ci_lower` and `ci_upper` (SE of ln(OR) = (ln U − ln L) / (2 × 1.96)) and `pvalue_coloc` for a beta
file with a p-value and an effect-allele frequency and no confidence limits (coloc.abf's p-value
form: p, MAF of the file's own frequency, N, `type = "cc"`, s; no beta or varbeta); a file with the
field present on some rows and missing on others, or with no route, fails. A harmonised file is
`native_se`. Mode-specific row rules (`outcome_files.REJECT_REASONS`): in `or_ci_derived_se` a
limit that is missing, infinite or not positive, a lower limit not below the upper, or an OR outside
the interval by more than half a unit in the last printed decimal of the OR and of the limit; in
`pvalue_coloc` a p that is missing or 0, or a frequency missing, 0 or 1. In every mode of the GWAS-SSF
and author readers a beta that is missing or infinite, and a frequency outside [0, 1], is not read.
An `or_ci_derived_se` file also fails unless at least 95% of its rows read with 0 < p < 1 give a
two-sided Wald p from ln(OR) / SE within 0.1 on the log10 scale of the published p, or within the
published p's printed rounding where that is looser (`outcome_files.CiPCheck`). A file that passes
has its `collected` record superseded by one naming `uncertainty_mode`, `header_sha256` and
`validation_sha256` (the sha256 of its result), which enter the collect digest; `spawn` also refuses
unless every `collected` record names the mode, header hash and result hash of its result.
`check-validation --expected <json>` compares the classification of a report (files by mode, failed,
error) with an expected one and writes the comparison to `inputs/stage_b/validation/`.
`restore-validated` moves the `collected` record of a file back only where its validation failed for
its standard error alone (`stage_b.validate.failed_only_for_standard_error`), after its bytes match
the recorded sha256; nothing is deleted, and the next `validate` reads it again.

Re-collection. `collect --run-token <token> --accessions <acc> ...` (or `--accessions-file`) spawns
`collect_file` for those planned GWAS Catalog tasks only, with `supersede_absent_record`: an
`absent` record is first moved to `/stage_b/superseded/collect/gwas_catalog/<acc>__<utc>.json`
(`stage_b.collect.supersede_absent`), and the task runs again; a record of any other status, or of
another source, is left in place and returned as it is.

Ensembl (`stage_b.fetch.ensembl_json`, `ensembl_unknown_id`). Ensembl's one definitive absence is
an HTTP 400 for which all of these hold: the body parses as its error JSON, an object whose only
key `error` holds a string; that string is, in full, its unknown-identifier message, `<rsID> not
found for human` on `GET /variation/human/<rsID>` or `No variant found with ID '<rsID>'` on `POST
/vep/human/id`; and the rsID it names is one that was requested. A 400 that says anything else or
anything more (request size, syntax, a backend fault, a page not found), a body that is not that
JSON, a message naming another rsID, and every other status raise; an Ensembl 404 or 410 raises as
`protocol`, because Ensembl answers 404 for a badly formed URL and 400 for an id it does not know.
An rsID Ensembl reports as unknown on the other build has no position there, like an rsID without
a single mapping on that chromosome. VEP is asked in batches of 200 rsIDs; Ensembl leaves an
unknown rsID out of a reply that holds records for others and answers with the 400 only when no
rsID of the batch gave a record, so such a batch adds no record, and the lookup is an absence
(`protein_altering` missing unless ST29 flags the protein) only when every batch was answered so.

Persistent outages. A source that keeps failing with a retryable error is never converted to
unavailable. No code path counts failures or writes an unavailable step, a collect record or a
result after any number of them, and Modal's `retries=3` only runs the same call again. The file
or unit stays unfinished, `status` shows it, and `assemble` refuses while a planned file has no
record or a unit no result.

`launch_stage_b.py status` writes `inputs/stage_b/status/status_<utc>.json`: per source, the files
collected, absent, queried by region, pending, failed (with the error class) and stranded, and the
units done, pending, failed and stranded. It reads the volume, the state Modal reports for the
latest call spawned for each file and unit, and whether the app is deployed, so unfinished work in
an app that is gone is `stranded`, not `pending`.

`unit_plan.json` also records what `plan` read and wrote: `a_outputs` (the sha256 of
`hypotheses.csv` and `outcome_trait_coding.tsv`), `tables` (the sha256 of each instrument table
read: `ukbppp_st9`, `decode_st02`, `interval_st4`, `opengwas_gwasinfo`, `decode_listing`, and
`decode_smp_listing` where used) and `units_sha256` (of `units.jsonl`). `collect`, `spawn`
and `assemble` take `--run-token` and refuse unless `a_outputs` equals the sealed stage A hashes
and `units.jsonl` is unchanged; `assemble` also refuses while a planned file has no collect
record. Each unit directory on the `pqtl-v8-stage-b` volume holds
`FINGERPRINT.json` (the fingerprint, and the tool versions and collect digest it covers); `run_unit`
refuses a directory written under another fingerprint or under none. `assemble` recomputes each
unit's fingerprint from the unit record, the pins, the code and the plan it holds, the tool versions
the unit recorded and the collect digest recomputed from the collect records in `--collect-dir`,
refuses a `FINGERPRINT.json` or `result.json` carrying any other fingerprint or collect digest, and
refuses units run under different tool versions.

Collect digest (`stage_b.checkpoint.collect_digest`). A unit's name, size and ETag do not fix the
bytes of its files, so the fingerprint also covers what the collect phase recorded. The digest is
the sha256 of the canonical JSON (sorted keys, no whitespace) of a list with one entry per whole
file the unit reads (`collect.collect_tasks` of the unit), ordered by source and key:

| record `status` | entry |
|---|---|
| `collected` | `source`, `key`, `status`, `name`, `bytes`, `sha256` of the record; for a GWAS Catalog file read without a harmonised copy also `layout`, `build`, `position_offset`, `rsid_rule`, `meta_sha256`; for a GWAS Catalog whole file the pre-analysis validation passed also `uncertainty_mode`, `header_sha256`, `validation_sha256` |
| `absent` | `source`, `key`, `status`, `name` as the pinned listing gives it (else empty), `bytes` null, `sha256` empty |
| `unreadable` | `source`, `key`, `status`, `name` empty, `bytes` null, `sha256` empty, `reason` (the record's `detail`) |
| `remote_indexed` | `source`, `key`, `status`, `name` as the source's listing gives it, `bytes` null, `sha256` empty, `release` (the pinned release of the source: GWAS Catalog `r2026-09-13`) |

No address, token, ETag, time or detail is in an entry. `run_unit` computes the digest from the
records under `/stage_b/collect/` before the first step and a missing record raises; the analyze
phase opens a whole file only after its size and sha256 match its record. A file collected again
with other bytes under the same name, size and ETag therefore has another sha256 in its record,
another collect digest and another fingerprint, and every checkpoint and result written from the
earlier bytes is refused by `run_unit` and by `assemble`. A file queried by region has no
whole-file hash; the sha256 of each extract taken from it is in `regional_manifest.tsv`.

## Stage C → `C/output/`

`outcomes.csv`, one row per hypothesis in `hypotheses.csv`:

| column | meaning |
|---|---|
| hypothesis_id | from A |
| status_24 | advanced / no_observed_advancement / active / business_only / undated / no_phase |
| advanced_24 | 1 / 0 / empty (empty when status is active, business_only, undated or no_phase; no_phase is PREREG §Missing data "no retrievable phase", excluded and counted) |
| status_12, status_36, status_48, advanced_12, advanced_36, advanced_48 | maturation sensitivities (S14) |
| stop_code | efficacy / safety / unclassified / business / completed_no_successor / "" |
| why_stopped_texts | original `why_stopped` strings, `||`-joined |
| efficacy_coded | bool, for S7 |
| phase2_start_date, phase3_start_date, time_to_phase3_days, phase3_event | time to Phase III (secondary) |
| phase3_success | 1 / 0 / empty (secondary) |
| approved | bool (descriptive) |
| chembl_max_phase_for_ind | cross-check |
| earliest_phase2_start | for the `pre_pqtl_publication` comparison in D |
| last_phase2_end_date | latest completion or termination date of a linked Phase II trial ended at the freeze (the date the maturation rule compares with the window); empty when none; follow-up column of descriptive table 7 |

Also `trials_snapshot/` (`studies_linked.tsv`: the archived `studies` rows of every linked NCT id;
`not_found.jsonl`; `manifest.json` with the source sha256), `trial_links.csv`,
`post_freeze_updates.csv` (hypothesis_id, n_linked_trials, n_not_in_archive,
n_first_submitted_after_freeze, n_last_update_after_freeze: a reported diagnostic that no label
reads), three more diagnostics described below (`linkage_audit.csv`, `partial_date_boundary.csv`,
`phase23_diagnostic.csv`), `run_info.json` (stage, run token, repository commit, code digest,
ChEMBL cross-check run or not, freeze date, AACT source), `INPUTS.tsv`, `MANIFEST.tsv`.

The three diagnostics are descriptive, lie outside the registered analyses, and change no label;
stage D does not read them. `linkage_audit.csv` (step, count, description;
`ot_phase.linkage_audit`) counts each step of the trial-linking rule over all hypotheses: the
hypotheses with and without a `clinical_indication` row, the distinct rows, their distinct report
ids, the report ids that are not ClinicalTrials.gov registrations, the distinct NCT ids, those with
and without a row in the AACT `studies` table, and the hypotheses with no linked trial and with no
trial in the snapshot. `partial_date_boundary.csv` and `phase23_diagnostic.csv` (hypothesis_id,
quantity, registered, alternative; `diagnostics.py`) hold one row per quantity of a hypothesis
whose value under the registered rules differs from its value under an alternative reading of the
same trials, with the rules themselves (`rules.code_hypothesis`) unchanged. In
`partial_date_boundary.csv` a date given as year-month or year takes the last day of that month or
year, and the quantities are `status_12`, `status_24`, `status_36`, `status_48`. In
`phase23_diagnostic.csv` a registration listing both PHASE2 and PHASE3 counts as Phase III and not
as Phase II, the Open Targets stages stay as registered, and the quantities are the four statuses,
`phase2_start_date`, `phase3_start_date`, `time_to_phase3_days` and `phase3_event`. A file with a
header and no row means no hypothesis differs.

ClinicalTrials.gov is read through CTTI's AACT pipe-delimited daily export dated 2026-09-30,
table `studies`, not from the live API. Its members were created at about 05:12 UTC on
2026-09-30, so it holds ClinicalTrials.gov as available at that extraction time, approximately
through the end of 2026-09-29, and not changes submitted during 2026-09-30. The author downloaded
it through a signed-in AACT session and uploaded it to the `pqtl-v8-inputs` volume at
`/aact_20260930/`; that sealed copy is the analysis source. `modal_stage_c.py::verify_uploaded_aact`
checks the zip's sha256 against the one recorded at download, extracts `studies.txt` and writes
`VERIFIED.json` (snapshot date, source, archive and `studies.txt` sha256); stage C reads
`studies.txt` only after its sha256 matches. `nct_id` must be the first column of `studies.txt`
(rows are selected on the text before the first `|`); any other position raises `AactFormatError`.

Real runs go through `modal_stage_c.py`: Open Targets tables and the AACT archive from the
`pqtl-v8-inputs` volume (read-only), ChEMBL 37 from `/chembl_37/chembl_37.db` on
`proteome-mr-claim-audit-inputs` (read-only), opened only after its sha256 matches
`/chembl_37/VERIFIED.json`; a missing or mismatched pin raises. Output goes to the
`pqtl-v8-stage-c` volume, committed after the linked trial rows are written and at the end.

## Stage D → `D/output/`

Joins A, B, C on `hypothesis_id` (after the run guard confirms the logged A, B and C seals);
forms sets S1–S22; applies gates (structural, reliability,
prior dominance, sampler fallback); fits H1, H4 (fixed sequence), H2, H3, continuous-evidence
and within-gene models in PyMC (reference implementation `power/power_v9/power_v9_calibrate.py`
`fit_exact`, with covariates added); frequentist check with the size correction; marginal
quantities; descriptive tables 1–15 and the forest plots. S8, S9 and S21 are `in_s8 & heldout`,
`in_s9 & heldout` and `in_s21 & heldout` with the ordinary `mechanism_class` (variant-only rows
enter no other set), S19 the held-out S1 rows with `blood_secreted_hpa` and an `s19_arm`
consistent with the class (aligned with neutralizing_biologic, blocking with
small_molecule_blocker), S17 the state from `s17_sentinel_p` < 0.05, table 7
the follow-up from `last_phase2_end_date`, a `no_phase` status a missing outcome counted in each
set's exclusion report, S16 the S1 rows with `s16_evidence_state` in place of the evidence state
where stage B wrote one (not formed, and reported as not run, when it is empty for every S1
hypothesis), and table 12 the cross-source agreement from
`evidence_state_ukbppp` / `evidence_state_decode`; the join raises if the selected source's state
differs from `evidence_state`. When S16 is formed it holds every held-out S1 row and not only the
rows with an `s16_evidence_state`: a row with one takes it as its state and every other row keeps
its primary `evidence_state` (`stage_d/sets.py`, `s16_evidence_state` filled from `evidence_state`
over the S1 mask; the set's `replaced` note counts the rows that took the SMP-normalized state). S5 keeps only `protein_altering` false, S11 only
`splicing_candidate` false and S15f only `low_coverage` false; a hypothesis whose flag is missing
is excluded from that set and counted in the set's `flag_missing_excluded` note. Work files live
in one directory per run (`--work`, default `D/work/run`; `/vol/work` on Modal) bound to the run
fingerprint: `prepare` refuses a non-empty directory under another fingerprint or under none and
does not rewrite a plan it already wrote under this one; the other phases refuse unless the
directory and `plan.json` carry this run's; `plan.json` records each design file's sha256;
`assemble` re-hashes every design and refuses a fit attempt, final fit or frequentist result with
another fingerprint or design sha256. `results.json` `meta.run_fingerprint` holds the fingerprint,
its hashed components and the environment components (as `FINGERPRINT.json` does), and `meta.repo_commit` the repository commit of the image; INPUTS.tsv names the
sealed A, B and C manifests and the `stage_code` row. Writes `results.json` (whole
structure: every coefficient, posterior probability, gate, diagnostic, set) before any table
or figure, and `tables/`, `figures/`, `INPUTS.tsv` (the sealed A, B, C manifests), `MANIFEST.tsv`.
