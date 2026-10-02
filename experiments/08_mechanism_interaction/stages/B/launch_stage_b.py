"""Stage B launcher: plan the units, spawn them on the deployed Modal app, assemble the outputs.

    # 1. units from hypotheses.csv + the sources' sentinel tables (local, no network)
    uv run --project stages/B --with modal==1.4.3 python stages/B/launch_stage_b.py plan \
        --hypotheses stages/A/output/hypotheses.csv --trait-coding stages/A/output/outcome_trait_coding.tsv \
        --ukbppp-st9 <Sun 2023 supplementary.xlsx> --decode-st02 <Ferkingstad 2021 MOESM4.xlsx> \
        --interval-st4 <Sun 2018 MOESM4.xlsx> --opengwas-gwasinfo <gwasinfo_all.json> \
        --decode-urls <urls.txt> [--decode-smp-urls <urls.txt>]
    # 2. after `prereg log "SEAL stage=A manifest_sha256=<sha256>"` and then
    #    `prereg log "RUN_START stage=B token=<token>"`, then `modal deploy stages/B/modal_stage_b.py`
    #    (the image bakes PREREG.md and the sealed stage A files, so deploy after logging) and the
    #    smoke function passes
    uv run --project stages/B --with modal==1.4.3 python stages/B/launch_stage_b.py spawn --run-token <token>
    # 3. after every call finishes: modal volume get pqtl-v8-stage-b /stage_b/units <UNITS_DIR>
    uv run --project stages/B --with modal==1.4.3 python stages/B/launch_stage_b.py assemble \
        --run-token <token> --units-dir <UNITS_DIR> --hypotheses ... --st29 <Eldjarn 2023 MOESM3.xlsx>

The OpenGWAS listing and the INTERVAL supplement live on the `pqtl-v8-inputs` volume
(08_mechanism_interaction/feasibility/v2_all_indications/inputs/); `plan` reads local copies,
fetched with `modal volume get pqtl-v8-inputs <remote path> <local path>` and checked against
modal_inputs_manifest.json.

`spawn` and `assemble` refuse to start unless `git status --porcelain` is empty for stages/ and
PREREG.md (v8_run_guard.clean_commit; commit unit_plan.json after `plan`, and launch_log.jsonl
after `spawn`). `assemble` writes that commit to run_info.json and INPUTS.tsv. Both also refuse
unless PREREG.md passes `prereg check`, its log holds
`SEAL stage=A manifest_sha256=<sha256>` and then `RUN_START stage=B token=<token>` for the
`--run-token` given, the seal is the sha256 of stages/A/output/MANIFEST.tsv, and hypotheses.csv and
outcome_trait_coding.tsv match that manifest (stage_b.launch.authorize, the shared guard in
stages/run_guard); each Modal call makes the same check on its baked copies. `plan` records the
sha256 of the stage A files it read and of units.jsonl in unit_plan.json; `spawn` refuses units
planned from files other than the sealed ones. `assemble` accepts a unit directory only under the
fingerprint of this run (stage_b/checkpoint.py), recomputed from the unit record, the pins, the
code and the plan held here and the tool versions the unit recorded, and only when every unit ran
under the same tool versions.
Units carry deCODE links, which are credentials of a sort and expire, so units.jsonl lives in
the gitignored inputs directory, never in stages/B/output/.
"""
import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import modal

from v8_manifest import sha256_file
from v8_run_guard import PLAN_SHA256, clean_commit

from stage_b.assemble import build_evidence, collect_unit_dir, load_st29, regional_rows, write_outputs
from stage_b.checkpoint import common_tools, package_sha256, source_pins, unit_fingerprint, unit_tools
from stage_b.launch import A_FILES, authorize, check_plan_inputs, spawn_units
from stage_b.schemas import InputContractError, InstrumentUnit, Sentinel
from stage_b.sentinels import (decode_sentinels, interval_opengwas_id, interval_sentinels, load_decode_st02,
                               load_interval_st4, load_ukbppp_st9, ukbppp_sentinels)
from stage_b.units import build_units, load_hypotheses, load_trait_coding

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
REPO = EXP.parents[1]
PREREG = EXP / "PREREG.md"
A_OUTPUT = HERE.parent / "A" / "output"
INPUTS_MANIFEST = EXP / "modal_inputs_manifest.json"
UNITS = EXP / "inputs" / "stage_b" / "units.jsonl"
PLAN_JSON = HERE / "output" / "unit_plan.json"
LOG = HERE / "launch_log.jsonl"
APP = "pqtl-v8-stage-b"
DECODE_FILE = re.compile(r"/(\d+_\d+)_[^/?]+\.txt\.gz")


def decode_url_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    out = {}
    for line in path.read_text().split():
        m = DECODE_FILE.search(line)
        if not m:
            raise InputContractError(f"cannot read a SeqId from deCODE link {line.split('?')[0]}")
        out[m.group(1)] = line.strip()
    return out


def cmd_plan(a: argparse.Namespace) -> None:
    hyps = load_hypotheses(a.hypotheses)
    coding = load_trait_coding(a.trait_coding)
    sents: list[Sentinel] = (ukbppp_sentinels(load_ukbppp_st9(a.ukbppp_st9)) + decode_sentinels(load_decode_st02(a.decode_st02))
                             + interval_sentinels(load_interval_st4(a.interval_st4)))
    by_key = {(s.source, s.assay_id): s for s in sents}
    gwasinfo = json.loads(a.opengwas_gwasinfo.read_text())
    gwasinfo = list(gwasinfo.values()) if isinstance(gwasinfo, dict) else gwasinfo
    urls, smp = decode_url_map(a.decode_urls), decode_url_map(a.decode_smp_urls)

    def locate(s: Sentinel) -> str:
        if s.source == "ukbppp":
            return s.assay_id
        if s.source == "decode":
            if s.assay_id not in urls:
                raise InputContractError(f"no deCODE link supplied for SeqId {s.assay_id}")
            return urls[s.assay_id]
        return interval_opengwas_id(s.locator, gwasinfo)

    units, hyp_unit, unresolved, source_units = build_units(hyps, by_key, locate, coding, smp)
    UNITS.parent.mkdir(parents=True, exist_ok=True)
    UNITS.write_text("".join(u.model_dump_json() + "\n" for u in units))
    plan = {"generated_utc": datetime.now(timezone.utc).isoformat(), "hypotheses": len(hyps), "units": len(units),
            "a_outputs": dict(zip(A_FILES, (sha256_file(a.hypotheses), sha256_file(a.trait_coding)))),
            "units_sha256": sha256_file(UNITS),
            "hypothesis_unit": hyp_unit, "unresolved": unresolved, "hypothesis_source_units": source_units,
            "unit_ids": {u.unit_key: [u.source, u.assay_id] for u in units}}
    PLAN_JSON.parent.mkdir(parents=True, exist_ok=True)
    PLAN_JSON.write_text(json.dumps(plan, indent=1, sort_keys=True))
    print(f"wrote {UNITS} ({len(units)} units) and {PLAN_JSON}; {len(unresolved)} hypotheses without a regional file")


def cmd_spawn(a: argparse.Namespace) -> None:
    commit = clean_commit(REPO)
    fn = modal.Function.from_name(APP, "run_unit")
    calls = spawn_units(PREREG, a.run_token, UNITS, PLAN_JSON, A_OUTPUT,
                        lambda line, token: fn.spawn(line, token).object_id)
    entry = {"what": "spawn", "repo_commit": commit, "calls": calls, "utc": datetime.now(timezone.utc).isoformat()}
    with LOG.open("a") as fh:
        fh.write(json.dumps(entry) + "\n")
    print(f"{len(calls)} units spawned at {entry['utc']}")


def cmd_assemble(a: argparse.Namespace) -> None:
    commit = clean_commit(REPO)
    sealed = authorize(PREREG, a.run_token, A_OUTPUT)
    check_plan_inputs(PLAN_JSON, UNITS, sealed)
    if sha256_file(a.hypotheses) != sealed["hypotheses.csv"]:
        raise InputContractError(f"{a.hypotheses} is not the hypotheses.csv stage A sealed")
    hyps = load_hypotheses(a.hypotheses)
    plan = json.loads(PLAN_JSON.read_text())
    units = [InstrumentUnit.model_validate_json(line) for line in UNITS.read_text().splitlines()]
    pins, code = source_pins(INPUTS_MANIFEST), package_sha256(HERE / "stage_b")
    results, metas, tools = {}, {}, {}
    for u in units:
        tools[u.unit_key] = unit_tools(a.units_dir / u.unit_key)
        results[u.unit_key], metas[u.unit_key] = collect_unit_dir(
            u, a.units_dir / u.unit_key, unit_fingerprint(u, pins, code, PLAN_SHA256, tools[u.unit_key]))
    unit_ids = {k: tuple(v) for k, v in plan["unit_ids"].items()}
    rows = build_evidence(hyps, plan["hypothesis_unit"], unit_ids, results, load_st29(a.st29),
                          plan["hypothesis_source_units"])
    scripts = sorted((HERE / "stage_b").glob("*.py")) + [HERE / "stage_b" / "coloc_run.R", Path(__file__).resolve()]
    inputs = {"hypotheses": a.hypotheses, "st29_workbook": a.st29, "unit_plan": PLAN_JSON, "units": UNITS}
    manifest = write_outputs(HERE / "output", rows, regional_rows(metas), scripts, inputs, [("", REPO)],
                             script_root=HERE, run_token=a.run_token, repo_commit=commit, tools=common_tools(tools))
    print(f"wrote {HERE / 'output' / 'evidence.csv'} ({len(rows)} rows) and {manifest}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("plan")
    for name in ("hypotheses", "trait-coding", "ukbppp-st9", "decode-st02", "interval-st4", "opengwas-gwasinfo"):
        p.add_argument(f"--{name}", type=Path, required=True)
    p.add_argument("--decode-urls", type=Path, default=None)
    p.add_argument("--decode-smp-urls", type=Path, default=None)
    sp = sub.add_parser("spawn")
    sp.add_argument("--run-token", required=True, help="token of the PREREG.md entry 'RUN_START stage=B token=<token>'")
    s = sub.add_parser("assemble")
    s.add_argument("--run-token", required=True, help="token of the PREREG.md entry 'RUN_START stage=B token=<token>'")
    for name in ("hypotheses", "units-dir", "st29"):
        s.add_argument(f"--{name}", type=Path, required=True)
    a = ap.parse_args()
    {"plan": cmd_plan, "spawn": cmd_spawn, "assemble": cmd_assemble}[a.cmd](a)


if __name__ == "__main__":
    main()
