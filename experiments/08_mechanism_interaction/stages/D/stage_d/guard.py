"""Input guard: stage D refuses to run unless stages A, B and C are sealed and the seals are logged.

`logged_seals` applies the shared run guard (stages/run_guard/v8_run_guard.py): PREREG.md passes
`prereg check`, its log holds `RUN_START stage=D token=<token>` and, before it, `SEAL stage=X
manifest_sha256=<sha256>` for A, B and C equal to the MANIFEST.tsv files under the stages root. The
logged values are the expected hashes `verify_inputs` receives.

A stage is sealed when `stages/<X>/output/MANIFEST.tsv` (the shared format of
stages/run_guard/v8_manifest.py) has the sha256 logged below the line of PREREG.md, lists every file
D reads and INPUTS.tsv, and every listed file's sha256 and row count match the file on disk; a
mismatch refuses the run. The inputs INPUTS.tsv names are provenance, not outputs, and are not
looked for; INPUTS.tsv itself is covered by the seal through its MANIFEST.tsv row. B's and C's
INPUTS.tsv must name A's `hypotheses.csv` sha256, so the three outputs are one chain.

A sealed output directory holds nothing its MANIFEST.tsv does not list, apart from MANIFEST.tsv
itself and the files of SIDE_FILES: `unit_plan.json` in B/output (written by `launch_stage_b.py
plan` before the run and named, with its sha256, in B's INPUTS.tsv) and the repository's
`.gitignore` in C/output. Any other file there refuses the run.
"""
from pathlib import Path

from pydantic import BaseModel, ConfigDict
from v8_manifest import (INPUTS_NAME, MANIFEST_NAME, ManifestEntry, ManifestError, read_inputs, read_manifest,
                         sha256_file, verify_output_dir)
from v8_run_guard import require_run

__all__ = ["ManifestEntry", "StageInputError", "StageSeal", "logged_seals", "read_manifest", "sha256_file",
           "verify_inputs", "verify_stage"]


class StageInputError(RuntimeError):
    pass


REQUIRED_FILES = {
    "A": ("hypotheses.csv", "funnel.csv", "outcome_gwas_selection.csv", "mechanism_crosstab.csv"),
    "B": ("evidence.csv",),
    "C": ("outcomes.csv",),
}
SIDE_FILES = {"A": (), "B": ("unit_plan.json",), "C": (".gitignore",)}
CHAINED = ("B", "C")


class StageSeal(BaseModel):
    model_config = ConfigDict(frozen=True)
    stage: str
    manifest_path: str
    manifest_sha256: str
    files: dict[str, str]


def logged_seals(prereg: Path, run_token: str, stages_root: Path) -> dict[str, str]:
    """The MANIFEST.tsv sha256 of A, B and C as sealed in the PREREG.md log (shared run guard)."""
    manifests = {s: stages_root / s / "output" / MANIFEST_NAME for s in ("A", "B", "C")}
    return require_run(prereg, "D", run_token, manifests=manifests).seals


def verify_stage(stages_root: Path, stage: str, expected_manifest_sha256: str) -> StageSeal:
    out = stages_root / stage / "output"
    manifest = out / MANIFEST_NAME
    if not manifest.is_file():
        raise StageInputError(f"{manifest} does not exist; the stage is not sealed")
    actual = sha256_file(manifest)
    if actual != expected_manifest_sha256:
        raise StageInputError(f"stage {stage}: MANIFEST.tsv sha256 {actual} != logged {expected_manifest_sha256}")
    try:
        files = verify_output_dir(out, REQUIRED_FILES[stage], SIDE_FILES[stage])
    except ManifestError as err:
        raise StageInputError(f"stage {stage}: {err}") from err
    return StageSeal(stage=stage, manifest_path=str(manifest), manifest_sha256=actual, files=files)


def verify_inputs(stages_root: Path, expected_manifest_sha256: dict[str, str]) -> dict[str, StageSeal]:
    if set(expected_manifest_sha256) != {"A", "B", "C"}:
        raise StageInputError(f"logged MANIFEST sha256 required for A, B and C; got {sorted(expected_manifest_sha256)}")
    seals = {s: verify_stage(stages_root, s, expected_manifest_sha256[s]) for s in ("A", "B", "C")}
    a_hyp = seals["A"].files["hypotheses.csv"]
    for stage in CHAINED:
        named = {r.sha256 for r in read_inputs(stages_root / stage / "output" / INPUTS_NAME)}
        if a_hyp not in named:
            raise StageInputError(f"stage {stage}: INPUTS.tsv does not name A hypotheses.csv ({a_hyp})")
    return seals
