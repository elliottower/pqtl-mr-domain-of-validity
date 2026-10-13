"""Launcher for the deployed `pqtl-v9-followup` app: start, status, collect. Nothing here reads a result.

    F="uv run --project ../08_mechanism_interaction/stages/D --with modal==1.4.3 python launch_followup.py"
    $F start --run-token <token>     # needs `prereg log "RUN_START followup token=<token>"`, a commit and a deploy
    $F status                        # the orchestrator call's state and the files on the volume
    $F collect                       # copies /results from the volume to results/ (refused if results/ holds files)
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import modal

HERE = Path(__file__).resolve().parent
LOG = HERE / "launch_log.jsonl"
OUT = HERE / "results"
APP = "pqtl-v9-followup"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def cmd_start(a: argparse.Namespace) -> None:
    call = modal.Function.from_name(APP, "run_all").spawn(a.run_token)
    with LOG.open("a") as fh:
        fh.write(json.dumps({"what": "start", "call_id": call.object_id, "utc": utc_now()}) + "\n")
    print(f"{utc_now()} run_all spawned: {call.object_id}; commit {LOG.name}")


def cmd_status(a: argparse.Namespace) -> None:
    last = [json.loads(line) for line in LOG.read_text().splitlines()][-1]
    try:
        modal.FunctionCall.from_id(last["call_id"]).get(timeout=0)
        state = "done"
    except TimeoutError:
        state = "running"
    except Exception as err:  # noqa: BLE001 -- whatever the call raised
        state = f"failed: {type(err).__name__}: {str(err)[:400]}"
    vol = modal.Volume.from_name(APP)
    counts = {}
    for top in ("work/q1", "work/fits", "work/q6/pass1", "work/q6/pass2", "results"):
        try:
            counts[top] = len([e for e in vol.listdir(top, recursive=True) if e.type == modal.volume.FileEntryType.FILE])
        except modal.exception.NotFoundError:
            counts[top] = 0
    print(f"{utc_now()} run_all {state}; files {counts}")


def cmd_collect(a: argparse.Namespace) -> None:
    if OUT.exists() and any(p.is_file() for p in OUT.rglob("*")):
        raise FileExistsError(f"{OUT} already holds files; nothing is overwritten")
    vol = modal.Volume.from_name(APP)
    files = [e.path for e in vol.listdir("results", recursive=True) if e.type == modal.volume.FileEntryType.FILE]
    for path in files:
        dest = OUT / Path(path).relative_to("results")
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"".join(vol.read_file(path)))
    print(f"{utc_now()} copied {len(files)} files to {OUT}; not opened")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("start").add_argument("--run-token", required=True)
    sub.add_parser("status")
    sub.add_parser("collect")
    a = ap.parse_args()
    {"start": cmd_start, "status": cmd_status, "collect": cmd_collect}[a.cmd](a)


if __name__ == "__main__":
    main()
