"""Thin Modal wrapper for the V8 power grid (logic in power_v9_grid.py).

Deploy once, then start work with spawn so it survives the client disconnecting:
    modal deploy experiments/08_mechanism_interaction/power/power_v9/modal_power_v9.py
    uv run --with modal==1.4.3 python experiments/08_mechanism_interaction/power/power_v9/launch_power_v9.py
Results land on the volume `pqtl-v8-power` under /vol/power_v9/; each cell's JSONL is
appended and committed every 100 repetitions, and a restarted shard resumes from disk.
"""
from pathlib import Path

import modal

image = modal.Image.debian_slim(python_version="3.12").pip_install(
    "numpy==2.1.3", "pandas==2.2.3", "scipy==1.14.1", "statsmodels==0.14.4", "tqdm==4.67.1"
).env({"PYTHONPATH": "/root/power"})
if modal.is_local():
    HERE = Path(__file__).resolve().parent
    FEAS = HERE.parents[1] / "feasibility" / "v4_round3"
    image = (
        image.add_local_file(str(HERE / "power_v9_logic.py"), "/root/power/power_v9_logic.py", copy=True)
        .add_local_file(str(HERE / "power_v9_fast.py"), "/root/power/power_v9_fast.py", copy=True)
        .add_local_file(str(HERE / "power_v9_grid.py"), "/root/power/power_v9_grid.py", copy=True)
        .add_local_file(str(FEAS / "h1_membership.csv"), "/root/data/h1_membership.csv", copy=True)
        .add_local_file(str(FEAS / "s1_membership.csv"), "/root/data/s1_membership.csv", copy=True)
    )

app = modal.App("pqtl-v8-power-v9", image=image)
vol = modal.Volume.from_name("pqtl-v8-power", create_if_missing=True)
OUT = "/vol/power_v9"
REPS = 1000


@app.function(cpu=2, memory=4096, timeout=6 * 3600, retries=3, volumes={"/vol": vol})
def run_shard(mode: str, cell_indices: list[int]):
    from datetime import datetime, timezone

    from power_v9_grid import cells, run_cell, structures

    out_dir = Path(OUT)
    out_dir.mkdir(parents=True, exist_ok=True)
    vol.reload()
    st, comp = structures(Path("/root/data"))[mode]
    all_cells = cells(mode)
    for i in cell_indices:
        kind, sc = all_cells[i]
        run_cell(out_dir, mode, kind, sc, st, comp, REPS, on_commit=vol.commit)
        print(f"{datetime.now(timezone.utc).isoformat()} {mode} cell {i} ({kind} {sc.scenario_id}) done")
    return {"mode": mode, "cells": len(cell_indices)}


@app.function(cpu=2, memory=4096, timeout=3600, volumes={"/vol": vol})
def merge():
    import json

    from power_v9_grid import summarize

    vol.reload()
    summary = summarize(Path(OUT), ["h1", "h2", "h4"], REPS)
    (Path(OUT) / "summary.json").write_text(json.dumps(summary, indent=2))
    vol.commit()
    return {m: len(rows) for m, rows in summary["modes"].items()}
