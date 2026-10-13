"""Thin Modal wrapper for the coloc.susie time-limit figure (logic and the registered rules in
susie_time_limit.py). Reads, on the stage B volume, the unit file-time snapshot taken before assemble
(/vol/diagnostics/<snapshot>), and of each unit the s15g_basis and susie_note of every coloc__<tag>.json and
the retrieved_utc of its outcome__<tag>.meta.json; no other field of a colocalization is read or returned.
The launch log is passed from the repository. Writes, under results/susie_time_limit/:
pairs.csv (one row per attempted pair: unit, tag, basis, tier, minute-resolution start and end, minutes, kind),
summary.json (per-tier counts and the curve, with the inputs' sha256) and susie_time_limit.png / .pdf.

    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \\
        modal run experiments/08_mechanism_interaction/stages/B/modal_susie_time_limit.py \\
        --snapshot unit_file_times_20261011T234239Z.json
"""
import hashlib
import io
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import modal

HERE = Path(__file__).resolve().parent
OUT_REL = Path("results") / "susie_time_limit"
EXTRACT = Path("/vol/diagnostics/susie_time_limit_extract")   # per-unit checkpoint of the fields read
SNAPSHOT_SHA256 = "4396657b443dfbde22aaf4e248b30d19f2123f52cc68d15b378c090de71d6e84"

image = modal.Image.debian_slim(python_version="3.12").pip_install("numpy==2.1.3", "pandas==2.2.3", "matplotlib==3.11.2")
if modal.is_local():
    image = image.add_local_file(str(HERE / "susie_time_limit.py"), "/root/susie_time_limit.py")
app = modal.App("pqtl-v8-susie-time-limit", image=image)
vol = modal.Volume.from_name("pqtl-v8-stage-b")


def extract_unit(unit_key: str, paths: list[str]) -> list[dict]:
    """s15g_basis and susie_note of each coloc__<tag>.json of one unit, and, for the pairs whose basis
    is attempted, the retrieved_utc of outcome__<tag>.meta.json; nothing else is read."""
    from susie_time_limit import ATTEMPTED
    out = []
    for path in paths:
        tag = Path(path).name.removeprefix("coloc__").removesuffix(".json")
        rec = json.loads((Path("/vol/stage_b/units") / path).read_text())
        basis = rec.get("s15g_basis") or ""
        if basis not in ATTEMPTED:
            continue
        meta = json.loads((Path("/vol/stage_b/units") / unit_key / f"outcome__{tag}.meta.json").read_text())
        out.append({"path": path, "tag": tag, "basis": basis, "note": rec.get("susie_note") or "",
                    "retrieved_utc": meta["retrieved_utc"]})
    return out


@app.function(volumes={"/vol": vol}, timeout=6 * 3600)
def pairs_and_figure(snapshot: str, launch_log: str) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    from susie_time_limit import LIMIT_MIN, pair_row, spawns_by_unit, summarize

    vol.reload()
    snap_path = Path("/vol/diagnostics") / snapshot
    if hashlib.sha256(snap_path.read_bytes()).hexdigest() != SNAPSHOT_SHA256:
        raise ValueError(f"{snap_path} is not the snapshot taken before assemble")
    mtimes = {f["path"]: f["mtime_utc"] for f in json.loads(snap_path.read_text())["files"]}
    spawns = spawns_by_unit(json.loads(line) for line in launch_log.splitlines() if line.strip())
    by_unit = defaultdict(list)
    for path in sorted(p for p in mtimes if Path(p).name.startswith("coloc__")):
        by_unit[path.split("/")[0]].append(path)
    # the three fields read per pair, checkpointed per unit so an interrupted call resumes
    EXTRACT.mkdir(parents=True, exist_ok=True)
    todo = [u for u in sorted(by_unit) if not (EXTRACT / f"{u}.json").is_file()]
    with ThreadPoolExecutor(max_workers=32) as pool:
        for i, (unit_key, fields) in enumerate(zip(todo, pool.map(lambda u: extract_unit(u, by_unit[u]), todo)), 1):
            (EXTRACT / f"{unit_key}.json").write_text(json.dumps(fields))
            if i % 25 == 0 or i == len(todo):
                vol.commit()
                print(f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} extracted {i}/{len(todo)} units")
    rows = []
    for unit_key in sorted(by_unit):
        for f in json.loads((EXTRACT / f"{unit_key}.json").read_text()):
            rows.append(pair_row(unit_key, f["tag"], f["basis"], f["note"], f["retrieved_utc"], mtimes[f["path"]], spawns))
    pairs = pd.DataFrame(rows)
    summary = summarize(pairs)
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    for tier, label in (("2cpu", "2 CPUs, 8 GiB"), ("8cpu", "8 CPUs, 32 GiB")):
        if tier in summary:
            c = summary[tier]["curve"]
            ax.step(c["minute"], c["share_completed"], where="post",
                    label=f"{label} (curve n = {summary[tier]['completed'] + summary[tier]['censored']})")
    ax.axvline(LIMIT_MIN, color="0.4", linestyle="--", linewidth=1)
    ax.set_xlim(0, LIMIT_MIN + 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("minutes from outcome region retrieved to colocalization written")
    ax.set_ylabel("share of coloc.susie runs completed")
    # curve n counts completed and censored runs; interrupted and stopped-in-R runs are counted in summary.json
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    figures = {}
    for ext in ("png", "pdf"):
        buf = io.BytesIO()
        fig.savefig(buf, format=ext, dpi=200)
        figures[ext] = buf.getvalue()
    return {"pairs_csv": pairs.to_csv(index=False), "summary": summary, "figures": figures}


@app.local_entrypoint()
def main(snapshot: str):
    out = HERE.parents[1] / OUT_REL
    if out.exists():
        raise FileExistsError(f"{out} exists; move it to results/superseded/ first, nothing is overwritten")
    log = (HERE / "launch_log.jsonl").read_text()
    got = pairs_and_figure.remote(snapshot, log)
    out.mkdir(parents=True)
    (out / "pairs.csv").write_text(got["pairs_csv"])
    summary = {"snapshot": snapshot, "snapshot_sha256": SNAPSHOT_SHA256,
               "launch_log_sha256": hashlib.sha256(log.encode()).hexdigest(),
               "logic_sha256": hashlib.sha256((HERE / "susie_time_limit.py").read_bytes()).hexdigest(),
               "tiers": got["summary"]}
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    for ext, data in got["figures"].items():
        (out / f"susie_time_limit.{ext}").write_bytes(data)
    print(f"wrote {out}: pairs.csv, summary.json, susie_time_limit.png/.pdf")
