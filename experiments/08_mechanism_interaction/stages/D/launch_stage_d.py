"""Stage D launcher: start, status, assemble. Every step but `status` is part of the stage D run and
needs the seals of A, B and C, then `prereg log "RUN_START stage=D token=<token>"`, a commit, and
`modal deploy stages/D/modal_stage_d.py` after it (the image bakes PREREG.md and the sealed stage
outputs). Orchestration only: every remote function applies the shared guard and the run fingerprint
itself (modal_stage_d.py); nothing here reads a fit, a frequentist result or a decision.

    D="uv run --project stages/D --with modal==1.4.3 python stages/D/launch_stage_d.py"
    $D start --run-token <token>      # prepare on Modal, then one call per design (frequentist) and per fit
    $D status --run-token <token>     # until every fit and frequentist check is done
    $D assemble --run-token <token>   # decisions, results.json, tables, figures, MANIFEST.tsv -> D/output/

`start` is resumable: `prepare_remote` returns the existing plan when the work volume already holds
this run's, and a fit or frequentist check already done returns its checkpoint, so running `start`
again respawns only what is unfinished. `start` and `assemble` refuse unless `git status
--porcelain` is empty for stages/ and PREREG.md (commit launch_log.jsonl after each), and unless
PREREG.md authorizes the run for the token (stage_d.guard.logged_seals, the shared guard).

`status` writes work/status/status_<utc>.json: per fit and per design whether its result file
exists on the volume (`plan_status` returns booleans only) and the state of the last spawned call.

`assemble` writes nothing to the terminal from the result: the decisions the remote call returns are
discarded, and the files of /output on the volume are copied to D/output/, which must be empty.
"""
import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import modal
from stage_d.guard import logged_seals
from v8_run_guard import clean_commit

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
REPO = EXP.parents[1]
PREREG = EXP / "PREREG.md"
STAGES = HERE.parent
OUT = HERE / "output"
LOG = HERE / "launch_log.jsonl"
APP = "pqtl-v8-stage-d"
VOLUME = "pqtl-v8-stage-d"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log_entry(entry: dict) -> None:
    with LOG.open("a") as fh:
        fh.write(json.dumps(entry) + "\n")


def fn(name: str) -> modal.Function:
    return modal.Function.from_name(APP, name)


def cmd_start(a: argparse.Namespace) -> None:
    commit = clean_commit(REPO)
    logged_seals(PREREG, a.run_token, STAGES)
    plan = fn("prepare_remote").remote(a.run_token)
    freq = [{"design_id": d, "call_id": fn("frequentist").spawn(d, a.run_token).object_id} for d in plan["designs"]]
    fits = [{"fit_id": f, "call_id": fn("fit_one").spawn(f, a.run_token).object_id} for f in plan["fits"]]
    log_entry({"what": "start", "repo_commit": commit, "designs": len(freq), "fits": len(fits),
               "frequentist_calls": freq, "fit_calls": fits, "utc": utc_now()})
    print(f"{utc_now()} {len(freq)} frequentist calls and {len(fits)} fit calls spawned; commit {LOG.name} next")


def call_state(call_id: str) -> str:
    try:
        modal.FunctionCall.from_id(call_id).get(timeout=0)
    except TimeoutError:
        return "pending"
    except Exception as err:  # noqa: BLE001 -- whatever the call raised, remotely or when its result was fetched
        return f"failed: {type(err).__name__}: {str(err)[:300]}"
    return "done"


def latest_calls() -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {"frequentist": {}, "fits": {}}
    for line in LOG.read_text().splitlines() if LOG.exists() else []:
        entry = json.loads(line)
        out["frequentist"].update({c["design_id"]: c["call_id"] for c in entry.get("frequentist_calls", [])})
        out["fits"].update({c["fit_id"]: c["call_id"] for c in entry.get("fit_calls", [])})
    return out


def cmd_status(a: argparse.Namespace) -> None:
    done = fn("plan_status").remote(a.run_token)
    vol = modal.Volume.from_name(VOLUME)
    try:
        freq_done = {Path(e.path).stem for e in vol.listdir("/work/freq") if e.path.endswith(".json")}
    except modal.exception.NotFoundError:      # no frequentist check has finished yet
        freq_done = set()
    calls = latest_calls()
    fits = {f: "done" if ok else call_state(calls["fits"][f]) if f in calls["fits"] else "never spawned"
            for f, ok in done.items()}
    freq = {d: "done" if d in freq_done else call_state(c) for d, c in calls["frequentist"].items()}
    report = {"utc": utc_now(), "fits": fits, "frequentist": freq}
    out = HERE / "work" / "status" / f"status_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1))
    by_state = {k: dict(Counter(s.split(":")[0] for s in v.values())) for k, v in (("fits", fits), ("frequentist", freq))}
    print(f"fits {by_state['fits']}; frequentist {by_state['frequentist']}; full report: {out}")


def cmd_assemble(a: argparse.Namespace) -> None:
    commit = clean_commit(REPO)
    logged_seals(PREREG, a.run_token, STAGES)
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError(f"{OUT} is not empty; nothing is overwritten")
    call = fn("assemble_remote").spawn(a.run_token)
    log_entry({"what": "assemble", "repo_commit": commit, "call_id": call.object_id, "utc": utc_now()})
    call.get()      # the decisions it returns are not read here
    vol = modal.Volume.from_name(VOLUME)
    files = [e.path for e in vol.listdir("/output", recursive=True) if e.type == modal.volume.FileEntryType.FILE]
    for path in files:
        dest = OUT / Path(path).relative_to("output")
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"".join(vol.read_file(path)))
    print(f"{utc_now()} copied {len(files)} files to {OUT}; not opened")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("start", "status", "assemble"):
        sub.add_parser(name).add_argument("--run-token", required=True,
                                          help="token of the PREREG.md entry 'RUN_START stage=D token=<token>'")
    a = ap.parse_args()
    {"start": cmd_start, "status": cmd_status, "assemble": cmd_assemble}[a.cmd](a)


if __name__ == "__main__":
    main()
