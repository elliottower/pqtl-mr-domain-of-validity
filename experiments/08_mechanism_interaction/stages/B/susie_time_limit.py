"""The coloc.susie time-limit figure registered on 2026-10-10 (PREREG log, "Amendment (B, two additions on the
coloc.susie time limit of S15g)", part 1). Pure logic; modal_susie_time_limit.py reads the volume and writes
the files.

Pairs: every instrument x outcome colocalization where coloc.susie was attempted, s15g_basis coloc.susie,
abf_susie_criterion_not_met or abf_susie_did_not_complete. Duration: from the outcome region's retrieved_utc
(outcome__<tag>.meta.json) to the modification time of coloc__<tag>.json, both read at the minute (seconds
dropped), the resolution the Modal volume listing reports. A pair is, in this order:
- interrupted, counted and left out of the curve, when a spawn of its unit in launch_log.jsonl falls inside
  (start, end], or the duration is over 70 minutes;
- stopped_in_r, counted and left out of the curve, when susie_note says R stopped the susie task;
- censored, at its duration, when susie_note says the call ran past the time limit;
- completed otherwise (coloc.susie or abf_susie_criterion_not_met: susie returned).
Tier: the unit's latest spawn at or before the start; before 2026-10-08T17:30:45Z the run_unit containers had
2 CPUs and 8 GiB, from then 8 CPUs and 32 GiB.
The curve is the Kaplan-Meier share completed by minute t, t = 0..60, per tier.
"""
from collections.abc import Iterable, Mapping
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd

ATTEMPTED = ("coloc.susie", "abf_susie_criterion_not_met", "abf_susie_did_not_complete")
EIGHT_CPU_FROM = datetime(2026, 10, 8, 17, 30, 45, tzinfo=timezone.utc)
INTERRUPTED_OVER_MIN = 70
LIMIT_MIN = 60
TIMEOUT_NOTE = "ran past"            # ColocTimeout: "Rscript ran past 3600 s on tasks [...]"
STOPPED_IN_R_NOTE = "stopped in R"   # SusieFailed: "susie did not complete: stopped in R: ..."
LATE_FROM_MIN = 45


class TimingError(ValueError):
    pass


def utc(text: str) -> datetime:
    t = datetime.fromisoformat(text)
    if t.tzinfo is None:
        raise TimingError(f"time without a zone: {text!r}")
    return t.astimezone(timezone.utc)


def to_minute(t: datetime) -> datetime:
    return t.replace(second=0, microsecond=0)


def spawns_by_unit(launch_log: Iterable[Mapping]) -> dict[str, list[datetime]]:
    """unit_key -> sorted spawn times, from the `spawn` entries of launch_log.jsonl."""
    out: dict[str, list[datetime]] = {}
    for entry in launch_log:
        if entry.get("what") == "spawn":
            when = utc(entry["utc"])
            for call in entry["calls"]:
                out.setdefault(call["unit_key"], []).append(when)
    return {k: sorted(v) for k, v in out.items()}


def classify(basis: str, note: str) -> str:
    if basis in ("coloc.susie", "abf_susie_criterion_not_met"):
        return "completed"
    if basis != "abf_susie_did_not_complete":
        raise TimingError(f"s15g_basis {basis!r} is not an attempted susie run")
    if STOPPED_IN_R_NOTE in note:
        return "stopped_in_r"
    if TIMEOUT_NOTE in note:
        return "censored"
    raise TimingError(f"abf_susie_did_not_complete with a note that names neither cause: {note!r}")


def pair_row(unit_key: str, tag: str, basis: str, note: str, retrieved_utc: str, coloc_mtime_utc: str,
             spawns: Mapping[str, list[datetime]]) -> dict:
    start, end = to_minute(utc(retrieved_utc)), to_minute(utc(coloc_mtime_utc))
    if end < start:
        raise TimingError(f"{unit_key}/{tag}: coloc file written before its outcome region")
    minutes = (end - start) / timedelta(minutes=1)
    unit_spawns = spawns.get(unit_key, [])
    before = [s for s in unit_spawns if s <= utc(retrieved_utc)]
    if not before:
        raise TimingError(f"{unit_key}/{tag}: no spawn of the unit before its outcome region was retrieved")
    restarted = any(utc(retrieved_utc) < s <= utc(coloc_mtime_utc) for s in unit_spawns)
    kind = classify(basis, note)
    if restarted:
        kind, why = "interrupted", "spawn inside the interval"
    elif minutes > INTERRUPTED_OVER_MIN:
        kind, why = "interrupted", f"over {INTERRUPTED_OVER_MIN} minutes"
    else:
        why = ""
    return {"unit_key": unit_key, "tag": tag, "s15g_basis": basis, "tier": "8cpu" if before[-1] >= EIGHT_CPU_FROM else "2cpu",
            "start_minute_utc": start.isoformat(), "end_minute_utc": end.isoformat(), "minutes": minutes,
            "kind": kind, "interrupted_because": why}


def kaplan_meier_completed(minutes: np.ndarray, completed: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Share completed by each time in `grid`: 1 - S(t), S the Kaplan-Meier survival of 'not yet
    completed', completions as events and the censored runs as censored at their duration."""
    order = np.argsort(minutes, kind="mergesort")
    t, e = minutes[order], completed[order].astype(bool)
    surv, s = [], 1.0
    times = np.unique(t[e])
    steps = {}
    for x in times:
        at_risk = int((t >= x).sum())
        events = int(((t == x) & e).sum())
        s *= 1 - events / at_risk
        steps[x] = s
    for g in grid:
        past = [steps[x] for x in times if x <= g]
        surv.append(past[-1] if past else 1.0)
    return 1 - np.array(surv)


def summarize(pairs: pd.DataFrame) -> dict:
    """Per tier: counts by kind, completions in minutes LATE_FROM_MIN..LIMIT_MIN, and the curve."""
    grid = np.arange(0, LIMIT_MIN + 1)
    out = {}
    for tier, g in pairs.groupby("tier", sort=True):
        kept = g[g["kind"].isin(["completed", "censored"])]
        curve = kaplan_meier_completed(kept["minutes"].to_numpy(float), (kept["kind"] == "completed").to_numpy(), grid)
        out[tier] = {"attempted": int(len(g)), **{k: int((g["kind"] == k).sum())
                                                  for k in ("completed", "censored", "stopped_in_r", "interrupted")},
                     "completed_in_minutes_45_to_60": int(((kept["kind"] == "completed")
                                                           & kept["minutes"].between(LATE_FROM_MIN, LIMIT_MIN)).sum()),
                     "curve": {"minute": grid.tolist(), "share_completed": curve.round(6).tolist()}}
    return out
