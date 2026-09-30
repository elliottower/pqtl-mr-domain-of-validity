"""Grid runner for the V8 power simulation (pure logic; no Modal dependency).

Every alternative scenario has a null counterpart with the same variance structure and no
focal effect (H1: OR ratio 1; H4: CNS ratio 1; H2: supportive OR 1 in both strata). The null
runs give the size-corrected critical value (95th percentile of the one-sided z), because the
two-way cluster-robust test is anti-conservative when one gene-drug component dominates. Both
nominal and size-corrected power are reported.

Checkpointing: results for a (mode, scenario) cell are appended to a JSONL file after every
CHUNK repetitions; a restarted run reads the file and continues from the repetition count
already on disk. Each repetition's random stream is derived from (seed, mode, scenario, rep),
so a resumed run reproduces exactly what an uninterrupted run would have produced.

Usage (local, small):
    uv run python power_v9_grid.py --reps 50 --modes h1 --max-scenarios 3 --out /tmp/x
"""
import argparse
import json
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from tqdm import tqdm

from power_v9_fast import component_clusters, run_rep_fast
from power_v9_logic import load_structure, scenario_grid

SEED = 20260930
CHUNK = 100
MODE_CODE = {"h1": 1, "h2": 2, "h4": 4}


def null_of(sc, mode):
    if mode == "h1":
        return replace(sc, or_ratio=1.0)
    if mode == "h4":
        return replace(sc, cns_ratio=1.0)
    return replace(sc, or_blocking=1.0, or_ratio=1.0)


IRRELEVANT = {"h1": {"scenario_id", "cns_ratio"}, "h2": {"scenario_id", "cns_ratio"},
              "h4": {"scenario_id", "or_ratio", "aligned_support_mult"}}


def null_key(sc, mode):
    return tuple(v for k, v in asdict(null_of(sc, mode)).items() if k not in IRRELEVANT[mode])


def cells(mode, max_scenarios=None):
    """Alternative cells plus one null cell per distinct variance/support/blocking setting."""
    grid = scenario_grid()[:max_scenarios] if max_scenarios else scenario_grid()
    out, seen = [], set()
    for sc in grid:
        out.append(("alt", sc))
        nsc = null_of(sc, mode)
        key = null_key(sc, mode)
        if key not in seen:
            seen.add(key)
            out.append(("null", nsc))
    return out


def cell_path(out_dir: Path, mode, kind, sc):
    return out_dir / f"{mode}_{kind}_{sc.scenario_id:04d}.jsonl"


def run_cell(out_dir: Path, mode, kind, sc, st, comp, reps, on_commit=None):
    path = cell_path(out_dir, mode, kind, sc)
    done = sum(1 for _ in path.open()) if path.exists() else 0
    buf = []
    for rep in range(done, reps):
        rng = np.random.default_rng([SEED, MODE_CODE[mode], sc.scenario_id, int(kind == "null"), rep])
        r = run_rep_fast(rng, st, comp, sc, mode)
        buf.append(json.dumps({"rep": rep, "reject": r["reject"], "failed": r["failed"],
                               "z": (r["estimate"] / r["se"]) * (-1 if mode == "h4" else 1)
                               if not r["failed"] else None,
                               "reject_sesoi": r["reject_sesoi"],
                               "n_supportive": r["n_supportive"],
                               "n_supportive_genes_aligned": r["n_supportive_genes_aligned"]}))
        if len(buf) == CHUNK or rep == reps - 1:
            with path.open("a") as fh:
                fh.write("\n".join(buf) + "\n")
            buf = []
            if on_commit:
                on_commit()
    return path


def structures(membership_dir: Path):
    h1 = membership_dir / "h1_membership.csv"
    s1 = membership_dir / "s1_membership.csv"
    return {"h1": (load_structure(h1), component_clusters(h1)),
            "h2": (load_structure(s1), component_clusters(s1)),
            "h4": (load_structure(s1), component_clusters(s1))}


def summarize(out_dir: Path, modes, reps, max_scenarios=None) -> dict:
    result = {"generated_utc": datetime.now(timezone.utc).isoformat(), "seed": SEED,
              "reps_per_cell": reps, "modes": {}}
    for mode in modes:
        nulls, rows = {}, []
        for kind, sc in cells(mode, max_scenarios):
            recs = [json.loads(line) for line in cell_path(out_dir, mode, kind, sc).open()]
            zs = np.array([r["z"] for r in recs if r["z"] is not None])
            key = null_key(sc, mode)
            if kind == "null":
                nulls[key] = {"size_nominal": float(np.mean([r["reject"] for r in recs])),
                              "crit_z": float(np.quantile(zs, 0.95)) if len(zs) else float("nan")}
            else:
                rows.append((key, sc, recs, zs))
        out = []
        for key, sc, recs, zs in rows:
            n = len(recs)
            crit = nulls[key]["crit_z"]
            p_nom = float(np.mean([r["reject"] for r in recs]))
            p_adj = float(np.mean(zs >= crit)) if len(zs) else float("nan")
            mc = 1.96 * np.sqrt(p_adj * (1 - p_adj) / n)
            out.append({**asdict(sc), "reps": n, "power_nominal": p_nom,
                        "power_size_corrected": p_adj,
                        "power_size_corrected_mc95": [max(0.0, p_adj - mc), min(1.0, p_adj + mc)],
                        "null_size_nominal": nulls[key]["size_nominal"], "crit_z": crit,
                        "failed": sum(r["failed"] for r in recs),
                        "median_supportive_genes_aligned": float(np.median([r["n_supportive_genes_aligned"] for r in recs]))})
        result["modes"][mode] = out
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--membership-dir", default=str(Path(__file__).resolve().parents[2] / "feasibility" / "v4_round3"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--modes", nargs="+", default=["h1", "h2", "h4"])
    ap.add_argument("--max-scenarios", type=int, default=None)
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    st = structures(Path(args.membership_dir))
    for mode in args.modes:
        s, comp = st[mode]
        for kind, sc in tqdm(cells(mode, args.max_scenarios), desc=f"{mode}"):
            run_cell(out_dir, mode, kind, sc, s, comp, args.reps)
    summary = summarize(out_dir, args.modes, args.reps, args.max_scenarios)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"{datetime.now(timezone.utc).isoformat()} wrote {out_dir / 'summary.json'}")


if __name__ == "__main__":
    main()
