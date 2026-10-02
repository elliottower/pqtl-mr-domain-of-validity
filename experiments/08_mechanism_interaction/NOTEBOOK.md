# V8 notebook

Dated provenance entries for `PREREG.md` and the stage runs. Append only.

## 2026-10-02 — results ledger started; stage A and the pre-freeze runs recorded

**Ledger.** `.results/` did not exist before today. `results init` was run at the repository root
after stage A had finished, so gate 1 (`results seal` before a run, `results run` after) was not
applied to stage A or to the three pre-freeze runs. They are recorded after the fact:

| run id | what | files |
|---|---|---|
| `v8-stage-a-f7014df548d04781` | stage A on Modal, commit dd7751f | `stages/A/output/*` |
| `v8-feasibility-count-v4` | outcome-blind count behind the registered sample sizes | `feasibility/v4_round3/coverage_v4.json`, `h1_membership.csv`, `s1_membership.csv` |
| `v8-power-v9-grid` | surrogate-test planning power, 324 scenarios | `power/power_v9/results/power_v9_summary.json` |
| `v8-calibrate-v9` | exact-model validation, 1,300 simulated datasets | `power/power_v9/results/calibrate_v9*.json`, `calibrate_v9/*.json` |

Stage A inputs: each of the 23 local input files was hashed and compared with the sha256 stage A
wrote to `stages/A/output/INPUTS.tsv` during the run; 23 of 23 match, and those files were then
sealed with role `input`. Stage A's own seal is `MANIFEST.tsv` sha256
`651ac6439cc724e2e4f34c1979c17c97a2e6486c6b4ba37e9cc2a17ba65a39ff`, logged as `SEAL stage=A` in
`PREREG.md`. The pre-freeze result files are unchanged since commit 0c168c6 (`git diff` empty).
`modal_calibrate_v9.py` was changed after the calibration runs (OpenBLAS added to the image), so it
was not sealed as the script of that run; `power_v9_calibrate.py`, which holds the model, is
unchanged and was sealed.

**Claims.** 25 claims bind the result-stating text of `PREREG.md` §Sample size to those runs: the
registered counts (to `v8-feasibility-count-v4`), the planning-power sentences and table rows (to
`v8-power-v9-grid`), and the exact-model sentences and table rows (to `v8-calibrate-v9`). All are
planning numbers, recorded as exploratory. `results coverage PREREG.md` still lists numbers that
are not results (dates, releases, design parameters such as the 324 scenarios).

**Values rechecked against the files.**

- Registered counts: stage A reproduced all 19 (`stages/A/output/registered_count_check.json`).
- Planning-power table: recomputed from `power_v9_summary.json`; all 12 cells and their ranges
  match. Two reference cells are exactly half-way values, 0.295 (H4, ratio 1/2, share 6%) and 0.735
  (H1, OR ratio 3, share 10%); the plan prints 0.30 and 0.74 (half up).
- Largest no-effect rejection rate of the surrogate test: 0.089 (H1), 0.194 (H2); H1 at an
  odds-ratio ratio of 1.5: 0.13–0.23 in the reference scenario, at most 0.41. All as printed.
- Exact-model table: all 10 rows match `calibrate_v9_round4_cells.json`. The "same decision" column
  prints 95%, 97% and 94% for 0.955, 0.975 and 0.935.

**Quotations.** The repository has no `claims/` directory, so `citations verify --claims claims/`
has nothing to check; `.citations/records/` holds 12 source records and no pinned quotation.
`PREREG.md` cites by key and quotes nothing. Any quotation in the paper is pinned with
`citations pin` before the sentence is written.

**`PREREG.md` in the ledger.** The file was sealed with role `prereg`. Its log grows with every
entry, so it is re-sealed after each `prereg log`; the frozen part is checked by `prereg check`.

**Inputs that cannot be sealed locally.** The AACT snapshot, ChEMBL 37 and the regional genetic
files exist only on Modal volumes. They are pinned by the sha256 values in each stage's
`INPUTS.tsv` and `VERIFIED.json`, and the ledger note of each run says so.
