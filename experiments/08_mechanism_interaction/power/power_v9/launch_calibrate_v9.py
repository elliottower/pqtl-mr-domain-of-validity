"""Spawn the exact-model check: 30 repetitions for each calibration scenario (after deploy)."""
import json
from datetime import datetime, timezone
from pathlib import Path

import modal

REPS = 30
LOG = Path(__file__).resolve().parent / "launch_log.jsonl"

fn = modal.Function.from_name("pqtl-v8-calibrate-v9", "run_rep")
calls = [{"label": lab, "rep": r, "call_id": fn.spawn(lab, r).object_id}
         for lab in ("null", "or2", "or3") for r in range(REPS)]
entry = {"what": "calibrate", "calls": calls, "utc": datetime.now(timezone.utc).isoformat()}
with LOG.open("a") as fh:
    fh.write(json.dumps(entry) + "\n")
print(len(calls), "calibration calls spawned at", entry["utc"])
