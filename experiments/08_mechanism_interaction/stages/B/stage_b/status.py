"""Where a stage B run stands: per source, the files collected, pending and failed; the units done,
pending and failed. Pure: `launch_stage_b.py status` gathers the three inputs and writes the report.

    volume   what is on the stage B volume (modal_stage_b.py::volume_state): the collect records,
             the partial downloads, the unit directories with a result.json, and the error markers
             the Modal functions leave when a call fails
    calls    the latest call the launcher spawned for each file and unit (launch_log.jsonl), with the
             state Modal reports for it: `done`, `failed` (with the exception class) or `pending`
    app      whether the deployed app answers, and each function's backlog and running containers

A file with a record is finished whatever its calls did. Without a record it is `failed` when its
latest call failed, `pending` when that call is still alive in a deployed app, `stranded` when the
call is pending but the app is not deployed (nothing will ever finish it: files on the volume alone
would look like slow progress), and `not_started` when no call was spawned. Units read the same way,
with result.json in place of the record.
"""
import json
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path

from stage_b.schemas import CollectTask, InstrumentUnit, RetryableSourceError, StageBError
from stage_b.validate import REPORT_NAME, VALIDATION_DIR


def error_of(err: BaseException) -> dict:
    """What the status report keeps of an exception: its class, the retryable kind, and the message
    only for stage B's own errors (their messages carry no URL)."""
    return {"error_class": type(err).__name__, "kind": err.kind if isinstance(err, RetryableSourceError) else "",
            "message": str(err)[:500] if isinstance(err, StageBError) else ""}


def error_marker(root: Path, kind: str, name: str) -> Path:
    """Where the failure of a collect task (`kind` "collect", name `<source>__<key>`) or of a unit
    (`kind` "units", name the unit key) is noted on the volume."""
    return root / "errors" / kind / f"{name}.json"


def marked(root: Path, kind: str, name: str, label: str, call: Callable[[], object], commit: Callable[[], None]):
    """`call()`. If it raises, the error is noted at `error_marker` (and re-raised); if it returns,
    an earlier note is removed. The volume then shows the last failure of everything unfinished."""
    path = error_marker(root, kind, name)
    try:
        out = call()
    except Exception as err:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"name": label, **error_of(err),
                                    "utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}, indent=1))
        commit()
        raise
    if path.exists():
        path.unlink()
        commit()
    return out


def volume_state(root: Path) -> dict:
    """What the stage B root on the volume holds: collect records, partial downloads, finished
    units, error markers, and, for the guard of `spawn` (validate.require_validation), the full
    GWAS Catalog records and the pre-analysis validation report (None when there is none). File
    names and small JSON records only; no regional file is opened."""
    records = [{k: json.loads(p.read_text())[k] for k in ("source", "key", "status")}
               for p in sorted((root / "collect").glob("*/*.json"))]
    partials = [{"source": p.parent.name, "name": p.name, "bytes": p.stat().st_size}
                for p in sorted((root / "raw").glob("*/*.part"))]
    done = sorted(d.name for d in (root / "units").glob("*") if (d / "result.json").is_file())
    started = sorted(d.name for d in (root / "units").glob("*") if d.is_dir())
    errors = {kind: {m["name"]: m for m in (json.loads(p.read_text()) for p in sorted((root / "errors" / kind).glob("*.json")))}
              for kind in ("collect", "units")}
    catalog = [json.loads(p.read_text()) for p in sorted((root / "collect" / "gwas_catalog").glob("*.json"))]
    report = root / VALIDATION_DIR / REPORT_NAME
    return {"records": records, "partials": partials, "units_done": done, "units_started": started, "errors": errors,
            "gwas_catalog_records": catalog, "outcome_validation": json.loads(report.read_text()) if report.is_file() else None}


def _state(finished: bool, call: Mapping | None, marker: Mapping | None, deployed: bool) -> dict:
    if finished:
        return {"state": "done"}
    if call is None:
        return {"state": "not_started", **({"last_error": marker} if marker else {})}
    if call["state"] == "failed":
        return {"state": "failed", "call_id": call["call_id"], **{k: call.get(k, "") for k in ("error_class", "kind", "message")}}
    if call["state"] == "done":       # the call returned but left nothing: the volume does not show its work
        return {"state": "failed", "call_id": call["call_id"], "error_class": "NoRecordAfterCall", "kind": "", "message": ""}
    if not deployed:
        return {"state": "stranded", "call_id": call["call_id"], **({"last_error": marker} if marker else {})}
    return {"state": "pending", "call_id": call["call_id"], **({"last_error": marker} if marker else {})}


def status_report(tasks: Sequence[CollectTask], units: Sequence[InstrumentUnit], volume: Mapping, calls: Mapping,
                  app: Mapping) -> dict:
    deployed = bool(app.get("deployed"))
    records = {(r["source"], r["key"]): r["status"] for r in volume["records"]}
    files, by_source = {}, {}
    for t in tasks:
        name = f"{t.source}/{t.key}"
        st = _state((t.source, t.key) in records, calls.get("collect", {}).get(name),
                    volume["errors"]["collect"].get(name), deployed)
        if st["state"] == "done":
            st["record"] = records[(t.source, t.key)]
        files[name] = st
        by_source.setdefault(t.source, Counter())[st.get("record") or st["state"]] += 1
    unit_states, units_by_source = {}, {}
    done = set(volume["units_done"])
    for u in units:
        st = _state(u.unit_key in done, calls.get("units", {}).get(u.unit_key), volume["errors"]["units"].get(u.unit_key),
                    deployed)
        unit_states[u.unit_key] = st
        units_by_source.setdefault(u.source, Counter())[st["state"]] += 1
    failed = {k: v for k, v in {**files, **unit_states}.items() if v["state"] in ("failed", "stranded")}
    return {"app": dict(app),
            "files_by_source": {s: dict(sorted(c.items())) for s, c in sorted(by_source.items())},
            "units_by_source": {s: dict(sorted(c.items())) for s, c in sorted(units_by_source.items())},
            "failures_by_error_class": dict(sorted(Counter(
                (v.get("error_class") or v["state"]) + (f":{v['kind']}" if v.get("kind") else "") for v in failed.values()).items())),
            "partial_downloads": volume.get("partials", []),
            "files": files, "units": unit_states}
