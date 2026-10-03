"""Stage B behind the shared run guard (stages/run_guard/v8_run_guard.py), free of `modal` so the
guard seam is testable: `launch_stage_b.py` passes `modal.Function.spawn`, tests pass a fake.

`authorize` is the one check every entry point makes before reading anything: PREREG.md passes
`prereg check`; its log holds `SEAL stage=A manifest_sha256=<sha256>` and, after it, `RUN_START
stage=B token=<token>`; the seal equals the sha256 of stage A's MANIFEST.tsv; and the two stage A
files stage B reads (hypotheses.csv, outcome_trait_coding.tsv) have the sha256 and row count that
manifest lists. `launch_stage_b.py` calls it on the working tree, each Modal function on the copies
baked into its image.

`spawn_collect` and `spawn_units` then require that unit_plan.json was planned from those sealed
files and that units.jsonl is the file the plan wrote, so no file is collected and no unit analyzed
for another hypotheses.csv. `spawn_collect` starts one call per whole file (collect.collect_tasks);
`spawn_units` starts one call per unit and refuses while any of those files has no collect record.
Each call carries the token, so the Modal function applies the same guard.

`sealed_stage_b` is the check before the collected whole files are deleted (collect.purge_raw): the
log holds a chained `SEAL stage=B manifest_sha256=<sha256>` and, where the caller has stage B's
MANIFEST.tsv, that file has the sealed sha256.

`private_copy` is the check on where `assemble` reads the unit directories and collect records
from. A copy of the unit directories holds the regional extracts, which are rows of the downloaded
source files and are never committed, so the copy lies outside the repository or under the
gitignored inputs directory.
"""
import json
from collections.abc import Callable, Collection, Sequence
from pathlib import Path

from v8_manifest import MANIFEST_NAME, sha256_file, verify_listed
from v8_run_guard import RunNotAuthorized, check_plan, logged_seals, parse_log, require_run

from stage_b.collect import collect_tasks
from stage_b.schemas import CollectRecord, CollectTask, InputContractError, InstrumentUnit
from stage_b.validate import require_validation

A_FILES = ("hypotheses.csv", "outcome_trait_coding.tsv")


def authorize(prereg: Path, run_token: str, a_output: Path) -> dict[str, str]:
    """Raises unless stage B may run; returns the sealed sha256 of each file in A_FILES."""
    require_run(prereg, "B", run_token, manifests={"A": a_output / MANIFEST_NAME})
    return verify_listed(a_output, A_FILES)


def sealed_stage_b(prereg: Path, manifest: Path | None = None) -> str:
    """The sha256 stage B is sealed under; raises unless PREREG.md passes `prereg check` and its log
    holds a chained SEAL for stage B (and, when `manifest` is given, unless that MANIFEST.tsv is the
    sealed one)."""
    text = prereg.read_text(encoding="utf-8")
    check_plan(text)
    entries = parse_log(text)
    seal = logged_seals(entries, len(entries) + 1).get("B")
    if seal is None:
        raise RunNotAuthorized("the PREREG.md log holds no 'SEAL stage=B manifest_sha256=<sha256>'; the collected files "
                               "are kept until stage B is sealed")
    if manifest is not None and sha256_file(manifest) != seal:
        raise RunNotAuthorized(f"{manifest} has sha256 {sha256_file(manifest)}, the log seals stage B as {seal}")
    return seal


def private_copy(path: Path, repo: Path, private_root: Path) -> Path:
    """`path`, refused when it lies inside `repo` but not under `private_root` (the gitignored
    inputs directory of the experiment)."""
    resolved = path.resolve()
    if resolved.is_relative_to(repo.resolve()) and not resolved.is_relative_to(private_root.resolve()):
        raise InputContractError(f"{path} lies in the repository outside {private_root}; unit directories hold regional "
                                 "extracts, which are not committed: copy them outside the repository or under "
                                 f"{private_root}")
    return path


def check_plan_inputs(unit_plan: Path, units: Path, sealed: dict[str, str]) -> None:
    """unit_plan.json names the sealed stage A files as its inputs and units.jsonl as its output."""
    plan = json.loads(unit_plan.read_text())
    if plan.get("a_outputs") != sealed:
        raise InputContractError(f"{unit_plan} was planned from stage A files {plan.get('a_outputs')}, "
                                 f"the sealed ones are {sealed}; run `plan` again on the sealed stage A output")
    got = sha256_file(units)
    if plan.get("units_sha256") != got:
        raise InputContractError(f"{units} has sha256 {got}, {unit_plan} recorded {plan.get('units_sha256')}")


def planned_units(prereg: Path, run_token: str, units: Path, unit_plan: Path, a_output: Path) -> list[str]:
    """The lines of units.jsonl, after the guard and the plan check."""
    check_plan_inputs(unit_plan, units, authorize(prereg, run_token, a_output))
    return units.read_text().splitlines()


def planned_tasks(lines: list[str]) -> list[CollectTask]:
    return collect_tasks(InstrumentUnit.model_validate_json(line) for line in lines)


def spawn_collect(prereg: Path, run_token: str, units: Path, unit_plan: Path, a_output: Path,
                  spawn: Callable[[str, str], str], sources: Collection[str] = (),
                  accessions: Collection[str] = ()) -> list[dict]:
    """`spawn(task_json, run_token)` starts one collect call and returns its call id. `sources`
    restricts the calls to those sources (all when empty). `accessions` restricts them to those
    GWAS Catalog tasks, each of which must be planned (the re-collection of accessions recorded
    absent for want of a harmonised file; the caller's `spawn` asks the call to supersede an
    `absent` record, stage_b.collect.supersede_absent)."""
    tasks = planned_tasks(planned_units(prereg, run_token, units, unit_plan, a_output))
    if accessions:
        planned = {t.key for t in tasks if t.source == "gwas_catalog"}
        unknown = sorted(set(accessions) - planned)
        if unknown:
            raise InputContractError(f"{len(unknown)} accessions are not planned GWAS Catalog files (first: {unknown[0]})")
        tasks = [t for t in tasks if t.source == "gwas_catalog" and t.key in set(accessions)]
    return [{"source": t.source, "key": t.key, "call_id": spawn(t.model_dump_json(), run_token)}
            for t in tasks if not sources or t.source in sources]


def spawn_units(prereg: Path, run_token: str, units: Path, unit_plan: Path, a_output: Path,
                spawn: Callable[[str, str], str], recorded: Collection[tuple[str, str]] | None = None,
                validation: tuple[dict | None, Sequence[CollectRecord]] | None = None) -> list[dict]:
    """`spawn(unit_json, run_token)` starts one call and returns its call id. `recorded` is the set
    of (source, key) that have a collect record on the volume; when given, nothing is spawned
    unless it covers every file the units read. `validation` is the pre-analysis validation report
    on the volume (None when there is none) and the current GWAS Catalog collect records; when
    given, nothing is spawned unless the report covers them (validate.require_validation)."""
    lines = planned_units(prereg, run_token, units, unit_plan, a_output)
    if recorded is not None:
        missing = [(t.source, t.key) for t in planned_tasks(lines) if (t.source, t.key) not in recorded]
        if missing:
            raise InputContractError(f"{len(missing)} files have no collect record (first: {missing[0]}); "
                                     "the collect phase is not finished")
    if validation is not None:
        require_validation(*validation)
    calls = []
    for line in lines:
        u = InstrumentUnit.model_validate_json(line)
        calls.append({"unit_key": u.unit_key, "call_id": spawn(line, run_token)})
    return calls
