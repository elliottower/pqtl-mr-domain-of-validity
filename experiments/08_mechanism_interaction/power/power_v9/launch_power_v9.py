"""Spawn the deployed V8 power grid in shards, then the merge (run after `modal deploy`).

    uv run --with modal==1.4.3 python launch_power_v9.py shards   # start all shards
    uv run --with modal==1.4.3 python launch_power_v9.py merge    # after shards finish
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import modal

APP = "pqtl-v8-power-v9"
N_CELLS = {"h1": 432, "h2": 360, "h4": 378}
CELLS_PER_SHARD = 30
LOG = Path(__file__).resolve().parent / "launch_log.jsonl"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["shards", "merge"])
    args = ap.parse_args()
    if args.what == "merge":
        call = modal.Function.from_name(APP, "merge").spawn()
        entry = {"what": "merge", "call_id": call.object_id}
    else:
        run_shard = modal.Function.from_name(APP, "run_shard")
        ids = []
        for mode, n in N_CELLS.items():
            for start in range(0, n, CELLS_PER_SHARD):
                call = run_shard.spawn(mode, list(range(start, min(start + CELLS_PER_SHARD, n))))
                ids.append({"mode": mode, "start": start, "call_id": call.object_id})
        entry = {"what": "shards", "calls": ids}
    entry["utc"] = datetime.now(timezone.utc).isoformat()
    with LOG.open("a") as fh:
        fh.write(json.dumps(entry) + "\n")
    print(json.dumps(entry)[:400])


if __name__ == "__main__":
    main()
