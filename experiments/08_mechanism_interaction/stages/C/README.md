# Stage C: outcome coding

Implements PREREG §Measured variables (primary outcome, decomposition, time to Phase III, Phase
III success, approval, ChEMBL cross-check) and writes `output/outcomes.csv` with the columns of
`../INTERFACES.md` §Stage C. Reads only `A/output/hypotheses.csv` (three columns), Open Targets
26.09 `clinical_indication`, `drug_molecule` and `clinical_report`, the ChEMBL 37 SQLite database, and
ClinicalTrials.gov as of the outcome-freeze date (2026-09-30): `studies.txt` of the AACT daily
flat-file snapshot dated 2026-09-30, whose members were created at about 05:12 UTC on 2026-09-30
(ClinicalTrials.gov as available at that extraction time, approximately through the end of
2026-09-29, without changes submitted during 2026-09-30), downloaded by the author through a signed-in AACT session and uploaded to
`pqtl-v8-inputs:/aact_20260930/`. AACT does not retain daily snapshots, so that copy is the record.
The live ClinicalTrials.gov API is not used. Never reads `stages/B`.

```
uv run pytest                                   # synthetic tests only
# once: check the uploaded snapshot's sha256, extract studies.txt and write VERIFIED.json in pqtl-v8-inputs:/aact_20260930/
uv run --with modal==1.4.3 modal run modal_stage_c.py::verify_uploaded_aact --name 20260930_export_ctgov.zip --expected-sha256 <sha256>
# real run, on Modal (inputs volume pqtl-v8-inputs, ChEMBL on proteome-mr-claim-audit-inputs):
uv run --with modal==1.4.3 modal run modal_stage_c.py --run-token <token>
```

The run refuses unless PREREG.md passes `prereg check` and its log holds the entry
`RUN_START stage=C token=<token>` (shared guard, `../run_guard/v8_run_guard.py`), checked locally and
again in the container against the PREREG.md baked into the image. The launcher also refuses unless
`git status --porcelain` is empty for `stages/` and `PREREG.md`; the commit goes into
`run_info.json` and the `stage_code` row of `INPUTS.tsv`.

ChEMBL is opened only after the sha256 of `/chembl_37/chembl_37.db` matches the pin in
`/chembl_37/VERIFIED.json` (`chembl.verified_chembl_db`); a missing or mismatched pin raises. The
AACT `studies.txt` is read only after its sha256 matches `/aact_20260930/VERIFIED.json`, which
`verify_uploaded_aact` writes (`ctgov.verified_studies`, `aact_download.verify_archive`).

## Rule to function

| plan rule | function |
|---|---|
| advanced: any drug at PHASE_3 / PREAPPROVAL / APPROVAL for the indication | `rules.code_hypothesis`, `models.ADVANCED_STAGES` |
| active: any trial recruiting, active, not yet recruiting, enrolling by invitation | `rules.primary_status`, `models.ACTIVE_STATUSES` |
| no observed advancement: last Phase II completion/termination ≥ 24 months before 2026-09-30 | `rules.primary_status`, `rules.last_end`, `rules.is_mature` |
| 12 / 36 / 48-month maturation (S14) | `models.MATURATION_WINDOWS_MONTHS`, `rules.code_hypothesis` |
| keyword table | `rules.classify_why_stopped`, `models.STOP_KEYWORDS` |
| precedence efficacy > safety > unclassified > business | `rules.combine_stop_codes`, `models.STOP_PRECEDENCE` |
| every stop business → excluded | `rules.stop_decomposition`, `rules.primary_status` |
| completed, no successor; never read as efficacy | `rules.stop_decomposition` |
| undated (pre-stage amendments table) | `rules.primary_status` |
| no retrievable phase → excluded and counted | `rules.code_hypothesis` (`no_phase`) |
| records dated after the freeze change no label | the AACT snapshot of 2026-09-30 (`ctgov`), `rules.registered_by_freeze`, `rules.status_at_freeze` |
| linked trials last updated after the freeze (reported diagnostic, no label) | `rules.post_freeze_counts` -> `post_freeze_updates.csv` |
| counts at each trial-linking step (diagnostic, no label) | `ot_phase.linkage_audit` -> `linkage_audit.csv` |
| statuses that differ with partial dates on the last day of the month or year (diagnostic, no label) | `diagnostics.partial_date_boundary`, `ctgov.parse_study(last_day=True)` -> `partial_date_boundary.csv` |
| statuses and time-to-Phase-III fields that differ with Phase II/III read as Phase III (diagnostic, no label) | `diagnostics.phase23_diagnostic`, `diagnostics.phase23_as_phase3` -> `phase23_diagnostic.csv` |
| time to Phase III, censored at the freeze | `rules.time_to_phase3` |
| Phase III success | `rules.phase3_success` |
| approval flag | `rules.code_hypothesis` |
| ChEMBL `max_phase_for_ind` | `chembl.verified_chembl_db`, `chembl.load_chembl_index`, `chembl.chembl_max_phase` |
| follow-up from the last Phase II (descriptive table 7) | `rules.code_hypothesis` (`last_phase2_end`, from `rules.last_end`) |
| drug program = parent molecule | `ot_phase.program_of` |
| trial linking (not in the plan) | `ot_phase.link_hypotheses`, `ot_phase.nct_of_report` |
| archived rows of every linked trial | `ctgov.AactTrials` -> `trials_snapshot/studies_linked.tsv` |
| archive check: zip sha256 against the recorded one, CRC, studies.txt header | `aact_download.verify_archive` (URL route: `aact_download.download_archive`) |
| C never reads B | `hypotheses.read_hypotheses`, `hypotheses.assert_not_stage_b`, `tests/test_blinding.py` |

## Choices needing a logged amendment before stage C runs

The plan does not fix these. Each is implemented as below and marked `IMPLEMENTATION CHOICE` in
the code.

1. **Trial linking.** A hypothesis's trials are the ClinicalTrials.gov registrations among the
   `clinicalReportIds` of the Open Targets 26.09 `clinical_indication` rows that give its phase:
   `diseaseId` equal to the indication id exactly (no descendants), drug program in
   `drug_program_ids` (salts through `parentId`). A report is a trial when
   `clinical_report.source` is "ClinicalTrials.gov" (its `id` is `nct` + 8 digits). This uses a
   table the plan does not list, `clinical_report` (26.09, sha256
   `ba0ea48b129182de119f5d676bb1483e6225ba65857a31b7319b65c0e6b8fa25`).
2. **Source of trial records.** ClinicalTrials.gov at the outcome-freeze date is read from the
   AACT daily flat-file snapshot of 2026-09-30, table `studies` (`studies.txt`), columns
   `nct_id`, `overall_status`, `why_stopped`, `phase`, `start_month_year`,
   `completion_month_year`, `study_first_submitted_date`, `last_update_submitted_date` (and the
   `*_date` / `*_date_type` columns, carried into `trials_snapshot/` but not read by a rule). A
   linked NCT id absent from the archive is a trial not found.
3. **Vocabulary.** `overall_status` and `phase` are mapped one-to-one to the ClinicalTrials.gov
   API v2 enums (`ACTIVE_NOT_RECRUITING`, `PHASE1/PHASE2` -> PHASE1, PHASE2) from either the enum
   or the legacy AACT display value ("Active, not recruiting", "Phase 1/Phase 2"); any other value
   raises.
4. **Phase of a registration.** Phase II = the phase list contains PHASE2 (Phase I/II and II/III
   included); Phase III = contains PHASE3 and not PHASE2, matching Open Targets, where
   PHASE_2_3 is not advanced.
5. **Completion date.** `completion_month_year` of COMPLETED or TERMINATED trials; a withdrawn
   or suspended trial has no completion or termination date.
6. **Partial dates.** The raw `start_month_year` / `completion_month_year` strings are parsed
   (`YYYY`, `YYYY-MM`, `YYYY-MM-DD`, or legacy `Month YYYY`, `Month D, YYYY`); a partial date takes
   the first day of the month or year. AACT's converted `start_date` / `completion_date` are not
   used, because the AACT dictionary describes their month-year conversion both as the first and
   as the last day of the month; a row with the raw value empty and the converted one set raises.
7. **Status at the freeze.** The snapshot is ClinicalTrials.gov as available at about 05:12 UTC on
   2026-09-30. Trials first submitted after 2026-09-30 are ignored; a trial whose completion
   date is after the freeze counts as ongoing, and its `why_stopped` is not read. Per hypothesis,
   `post_freeze_updates.csv` reports the linked trials first submitted after the freeze and those
   whose `last_update_submitted_date` is after it (none can be in a snapshot exported on the
   freeze date; the counts are kept as a check); no label reads these counts.
8. **One `why_stopped` matching several codes** takes the precedence order.
9. **Stops** are terminated, suspended or withdrawn trials of any phase. With no stop and at
   least one completed trial the code is `completed_no_successor`. A hypothesis with one
   business stop and otherwise only completed trials is `business_only`.
10. **Label order among exclusions:** ongoing trial -> `active`, then `business_only`, then
   `undated`.
11. **Time to Phase III.** Withdrawn trials and starts after the freeze are not starts; the event
   is the first Phase III start on or after the earliest Phase II start.
12. **Status value `no_phase`** implements PREREG §Missing data ("no retrievable phase", excluded
   and counted); it is the sixth status in INTERFACES.md, and stage D maps it to a missing outcome
   and counts it.
