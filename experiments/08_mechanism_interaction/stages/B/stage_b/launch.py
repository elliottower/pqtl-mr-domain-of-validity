"""Stage B behind the shared run guard (stages/run_guard/v8_run_guard.py), free of `modal` so the
guard seam is testable: `launch_stage_b.py spawn` passes `modal.Function.spawn`, tests pass a fake.

`authorize` is the one check every entry point makes before reading anything: PREREG.md passes
`prereg check`; its log holds `SEAL stage=A manifest_sha256=<sha256>` and, after it, `RUN_START
stage=B token=<token>`; the seal equals the sha256 of stage A's MANIFEST.tsv; and the two stage A
files stage B reads (hypotheses.csv, outcome_trait_coding.tsv) have the sha256 and row count that
manifest lists. `launch_stage_b.py spawn` and `assemble` call it on the working tree, the Modal
function on the copies baked into its image.

`spawn_units` then requires that unit_plan.json was planned from those sealed files and that
units.jsonl is the file the plan wrote, so no unit derived from another hypotheses.csv is spawned.
Each unit carries the token so the Modal function applies the same guard.
"""
import json
from collections.abc import Callable
from pathlib import Path

from v8_manifest import MANIFEST_NAME, sha256_file, verify_listed
from v8_run_guard import require_run

from stage_b.schemas import InputContractError, InstrumentUnit

A_FILES = ("hypotheses.csv", "outcome_trait_coding.tsv")


def authorize(prereg: Path, run_token: str, a_output: Path) -> dict[str, str]:
    """Raises unless stage B may run; returns the sealed sha256 of each file in A_FILES."""
    require_run(prereg, "B", run_token, manifests={"A": a_output / MANIFEST_NAME})
    return verify_listed(a_output, A_FILES)


def check_plan_inputs(unit_plan: Path, units: Path, sealed: dict[str, str]) -> None:
    """unit_plan.json names the sealed stage A files as its inputs and units.jsonl as its output."""
    plan = json.loads(unit_plan.read_text())
    if plan.get("a_outputs") != sealed:
        raise InputContractError(f"{unit_plan} was planned from stage A files {plan.get('a_outputs')}, "
                                 f"the sealed ones are {sealed}; run `plan` again on the sealed stage A output")
    got = sha256_file(units)
    if plan.get("units_sha256") != got:
        raise InputContractError(f"{units} has sha256 {got}, {unit_plan} recorded {plan.get('units_sha256')}")


def spawn_units(prereg: Path, run_token: str, units: Path, unit_plan: Path, a_output: Path,
                spawn: Callable[[str, str], str]) -> list[dict]:
    """`spawn(unit_json, run_token)` starts one call and returns its call id."""
    check_plan_inputs(unit_plan, units, authorize(prereg, run_token, a_output))
    calls = []
    for line in units.read_text().splitlines():
        u = InstrumentUnit.model_validate_json(line)
        calls.append({"unit_key": u.unit_key, "call_id": spawn(line, run_token)})
    return calls
