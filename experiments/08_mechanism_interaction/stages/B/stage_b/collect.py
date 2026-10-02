"""The collect phase: every whole file stage B needs, downloaded once to the stage B volume.

Whole files are the deCODE per-SeqId files (and their SMP-normalized release for S16), the UKB-PPP
per-protein tars, the UKB-PPP per-chromosome rsID maps and the GWAS Catalog harmonised files that
have no tabix index. `collect_tasks` lists them from the units, one task per file whatever the
number of units that read it. `collect_one` runs one task (one Modal call per file) and leaves:

    <root>/raw/<source>/<name>             the file
    <root>/collect/<source>/<key>.json     its record (schemas.CollectRecord): size, sha256, ETag,
                                           Last-Modified, the source URL without query or token,
                                           UTC time

A record is also written, without a file, when the source definitively does not hold the file
(`absent`) and for a GWAS Catalog file that has a tabix index (`remote_indexed`: queried by region
in the analyze phase, not downloaded). Any other fault raises (stage_b/remote.py) and writes no
record, so the task stays pending and a later call resumes it.

A download writes `<name>.part` and `<name>.part.json` (the validators of the response the bytes
came from). It commits the volume every CHECKPOINT_BYTES, and a later call continues with an HTTP
Range request from the bytes on disk; if the source then answers with other validators, or ignores
the range, the file starts over. A finished file must have the length the source declared, the
size and ETag of the pinned listing where the task carries them, the MD5 the source declares where
it declares one, and, if its name ends in .gz, must decompress to the end.

The analyze phase reads these files only through `read_record` and `collected_file`: a missing
record raises CollectError (it is not 'unavailable'), and a file is opened only after its size and
sha256 match its record.

The files under raw/ are working copies, held on the private stage B volume between collection
and the stage B seal. `purge_raw` deletes them once every unit has its result (the wrapper also
requires stage B's SEAL in the PREREG.md log). The records stay, with each file's size, sha256,
ETag and time, and so do the regional extracts of the unit directories; PURGED.json lists what was
deleted. The extracts hold rows of the downloaded files: they stay on the volume, as checkpoint
and provenance material, and are never written to the committed B/output/.

A record is what binds a unit's result to the bytes it was computed from: the sha256 here enters
the unit fingerprint through the collect digest (stage_b/checkpoint.py), so a file collected again
with other bytes under the same name, size and ETag invalidates every checkpoint that read it.
"""
import gzip
import hashlib
import json
import os
import re
from collections.abc import Callable, Iterable
from datetime import datetime, timezone
from pathlib import Path
from typing import Protocol

import requests
from v8_manifest import sha256_file

from stage_b.remote import CORRUPT, attempt
from stage_b.schemas import (CollectError, CollectRecord, CollectSource, CollectTask, InputContractError, InstrumentUnit,
                             RetryableSourceError, SourceAbsent)

CHECKPOINT_BYTES = 256 * 1024 * 1024
CHUNK_BYTES = 1024 * 1024
SAFE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.:-]*")
CONTENT_RANGE = re.compile(r"bytes (\d+)-(\d+)/(\d+)")


class Resolved(Protocol):
    """What a source says about one task's file, before any byte is downloaded."""

    name: str
    source_url: str     # without query string and without token
    indexed: bool       # a tabix index exists: the file is queried by region, not downloaded
    md5: str            # hex MD5 the source declares for the file, or ""

    def open_at(self, offset: int) -> requests.Response:
        """A 2xx streaming response for the file, from byte `offset` (an HTTP Range request when
        offset > 0); classified errors otherwise (stage_b.remote.http)."""


class Sources(Protocol):
    def resolve(self, task: CollectTask) -> Resolved:
        """The file of `task` at its source; SourceAbsent when the source does not hold it."""


def collect_tasks(units: Iterable[InstrumentUnit]) -> list[CollectTask]:
    """One task per whole file the units read, sorted by source and key."""
    tasks: dict[tuple[str, str], CollectTask] = {}

    def add(task: CollectTask) -> None:
        prev = tasks.setdefault((task.source, task.key), task)
        if prev != task:
            raise InputContractError(f"{task.source} {task.key} is listed with two records: {prev} and {task}")

    for u in units:
        if u.source == "decode":
            if u.pqtl_listing is None:
                raise InputContractError(f"deCODE unit {u.unit_key} carries no listing record for its file")
            add(CollectTask(source="decode", key=u.assay_id, **u.pqtl_listing.model_dump()))
            if u.smp_listing is not None:
                add(CollectTask(source="decode_smp", key=u.assay_id, **u.smp_listing.model_dump()))
        elif u.source == "ukbppp":
            add(CollectTask(source="ukbppp", key=u.assay_id))
            add(CollectTask(source="ukbppp_rsid_map", key=u.sentinel.chrom))
        for o in u.outcomes:
            if o.source == "gwas_catalog":
                add(CollectTask(source="gwas_catalog", key=o.accession))
    return [tasks[k] for k in sorted(tasks)]


def _safe(name: str, what: str) -> str:
    if SAFE_NAME.fullmatch(name) is None or ".." in name:
        raise CollectError(f"{what} {name!r} is not a plain file name")
    return name


def record_path(root: Path, source: CollectSource, key: str) -> Path:
    return root / "collect" / source / f"{_safe(key, 'collect key')}.json"


def read_record(root: Path, source: CollectSource, key: str) -> CollectRecord:
    """The record the collect phase left for a file; CollectError when there is none."""
    path = record_path(root, source, key)
    if not path.is_file():
        raise CollectError(f"no collect record for {source} {key} ({path}); the collect phase has not finished this "
                           "file, which is not the same as the file being unavailable")
    record = CollectRecord.model_validate_json(path.read_text())
    if (record.source, record.key) != (source, key):
        raise CollectError(f"{path} is the record of {record.source} {record.key}, not of {source} {key}")
    return record


def collected_file(root: Path, record: CollectRecord) -> Path:
    """The file of a `collected` record, after its size and sha256 are checked against the record."""
    if record.status != "collected":
        raise CollectError(f"{record.source} {record.key} is recorded as {record.status}; it has no file on the volume")
    path = root / record.path
    if not path.is_file():
        raise CollectError(f"{path} is missing although {record.source} {record.key} is recorded as collected")
    if path.stat().st_size != record.bytes:
        raise CollectError(f"{path} has {path.stat().st_size} bytes, its record {record.bytes}")
    got = sha256_file(path)
    if got != record.sha256:
        raise CollectError(f"{path} has sha256 {got}, its record {record.sha256}")
    return path


def _write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=True))
    os.replace(tmp, path)


def _validators(r: requests.Response, offset: int, what: str) -> dict:
    """ETag, Last-Modified and total length of the file a response is serving."""
    if r.status_code == 206:
        m = CONTENT_RANGE.fullmatch(r.headers.get("Content-Range", "").strip())
        if m is None or int(m.group(1)) != offset:
            raise RetryableSourceError("protocol", f"{what}: a range response without a usable Content-Range")
        size: int | None = int(m.group(3))
    else:
        length = r.headers.get("Content-Length", "")
        size = int(length) if length.isdigit() else None
    return {"etag": r.headers.get("ETag", "").strip().removeprefix("W/").strip('"'),
            "last_modified": r.headers.get("Last-Modified", ""), "size": size}


def download(open_at: Callable[[int], requests.Response], part: Path, what: str, commit: Callable[[], None],
             checkpoint_bytes: int = CHECKPOINT_BYTES) -> dict:
    """Bring `part` to the full length of the remote file and return the validators of the file it
    holds. One connection: a fault propagates with the bytes received so far flushed and committed,
    and the next call continues from them."""
    state_path = part.with_name(part.name + ".json")
    state = json.loads(state_path.read_text()) if state_path.is_file() and part.is_file() else None
    have = part.stat().st_size if state is not None else 0
    if have and state["size"] == have:
        return state                       # every byte is already on disk; a range request past the end would fail
    try:
        r = open_at(have)
    except RetryableSourceError as err:
        if not have or not str(err).endswith("HTTP 416"):
            raise
        have = 0                           # the bytes on disk reach past the end of the file now served: start over
        r = open_at(0)
    try:
        found = _validators(r, have, what)
        if have and (r.status_code != 206 or found != state):
            r.close()                      # the range was ignored, or the file changed: start over
            have = 0
            r = open_at(0)
            found = _validators(r, 0, what)
        if not have:
            part.parent.mkdir(parents=True, exist_ok=True)
            part.write_bytes(b"")
            _write_json(state_path, found)
        since = 0
        try:
            with part.open("ab") as f:
                for chunk in r.iter_content(CHUNK_BYTES):
                    f.write(chunk)
                    since += len(chunk)
                    if since >= checkpoint_bytes:
                        f.flush()
                        commit()
                        since = 0
        finally:
            commit()
    finally:
        r.close()
    got = part.stat().st_size
    if found["size"] is not None and got != found["size"]:
        raise RetryableSourceError("connection", f"{what}: the stream ended at {got} of {found['size']} bytes")
    return found


def gzip_intact(path: Path) -> bool:
    """True when the whole gzip file decompresses (what `gzip -t` checks)."""
    try:
        with gzip.open(path, "rb") as f:
            while f.read(CHUNK_BYTES):
                pass
    except CORRUPT:
        return False
    return True


def md5_file(path: Path) -> str:
    h = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(CHUNK_BYTES), b""):
            h.update(chunk)
    return h.hexdigest()


def _put_record(root: Path, record: CollectRecord, commit: Callable[[], None]) -> CollectRecord:
    _write_json(record_path(root, record.source, record.key), record.model_dump())
    commit()
    return record


def collect_one(task: CollectTask, sources: Sources, root: Path, commit: Callable[[], None] = lambda: None,
                checkpoint_bytes: int = CHECKPOINT_BYTES) -> CollectRecord:
    """Run one collect task to its record. A task that already has a record returns it and touches
    nothing, so a re-issued link or renewed token cannot change a collected file."""
    path = record_path(root, task.source, task.key)
    if path.is_file():
        return read_record(root, task.source, task.key)
    what = f"{task.source} {task.key}"
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    base = {"source": task.source, "key": task.key, "utc": now}
    try:
        resolved = attempt(lambda: sources.resolve(task), what)
    except SourceAbsent as err:
        return _put_record(root, CollectRecord(status="absent", detail=str(err), **base), commit)
    name = _safe(resolved.name, "source file name")
    if task.name and name != task.name:
        raise CollectError(f"{what}: the source names the file {name}, the pinned listing {task.name}")
    if resolved.indexed:
        return _put_record(root, CollectRecord(status="remote_indexed", name=name, source_url=resolved.source_url,
                                               **base), commit)
    final = root / "raw" / task.source / name
    part = final.with_name(name + ".part")
    state_path = part.with_name(part.name + ".json")
    if not final.is_file():
        try:
            found = attempt(lambda: download(resolved.open_at, part, what, commit, checkpoint_bytes), what)
        except SourceAbsent as err:
            return _put_record(root, CollectRecord(status="absent", name=name, detail=str(err), **base), commit)
        problems = []
        if task.size is not None and part.stat().st_size != task.size:
            problems.append(f"{part.stat().st_size} bytes, the pinned listing gives {task.size}")
        if task.etag and found["etag"] and found["etag"] != task.etag:
            problems.append("an ETag other than the pinned listing's")
        if resolved.md5 and md5_file(part) != resolved.md5:
            problems.append("an MD5 other than the one the source declares")
        if problems:
            raise CollectError(f"{what}: the downloaded file has " + "; ".join(problems))
        if name.endswith(".gz") and not gzip_intact(part):
            part.unlink()
            state_path.unlink()
            commit()
            raise RetryableSourceError("corrupt", f"{what}: the downloaded file is not a complete gzip stream")
        os.replace(part, final)
        commit()
    found = json.loads(state_path.read_text()) if state_path.is_file() else {"etag": "", "last_modified": ""}
    record = CollectRecord(status="collected", name=name, path=final.relative_to(root).as_posix(),
                           bytes=final.stat().st_size, sha256=sha256_file(final), md5=resolved.md5, etag=found["etag"],
                           last_modified=found["last_modified"], source_url=resolved.source_url, **base)
    _put_record(root, record, commit)
    if state_path.is_file():
        state_path.unlink()
        commit()
    return record


PURGED_NAME = "PURGED.json"


def purge_raw(root: Path, units: Iterable[InstrumentUnit], commit: Callable[[], None] = lambda: None) -> dict:
    """Delete every whole file under `root`/raw, after checking that the run no longer needs one:
    each unit has its result.json, each planned file its collect record, and each file to delete is
    the one a `collected` record names, with the recorded size. Nothing is deleted unless all of
    that holds. Records, unit directories and extracts are not touched. Returns (and writes to
    PURGED.json, extending an earlier one) the files deleted."""
    units = list(units)
    unfinished = sorted(u.unit_key for u in units if not (root / "units" / u.unit_key / "result.json").is_file())
    if unfinished:
        raise CollectError(f"{len(unfinished)} units have no result.json (first: {unfinished[0]}); raw files are kept")
    records = [read_record(root, t.source, t.key) for t in collect_tasks(units)]
    recorded = {r.path: r for r in records if r.status == "collected"}
    raw = root / "raw"
    found = sorted(p for p in raw.rglob("*") if p.is_file()) if raw.is_dir() else []
    doomed = []
    for p in found:
        rel = p.relative_to(root).as_posix()
        leftover = p.name.endswith((".part", ".part.json"))
        if not leftover and rel not in recorded:
            raise CollectError(f"{p} is the file of no collect record; nothing is deleted")
        if not leftover and p.stat().st_size != recorded[rel].bytes:
            raise CollectError(f"{p} has {p.stat().st_size} bytes, its record {recorded[rel].bytes}; nothing is deleted")
        doomed.append({"path": rel, "bytes": p.stat().st_size, "sha256": "" if leftover else recorded[rel].sha256})
    log_path = root / PURGED_NAME
    earlier = json.loads(log_path.read_text())["deleted"] if log_path.is_file() else []
    for d in doomed:
        (root / d["path"]).unlink()
    out = {"utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "units": len(units), "records_kept": len(records),
           "deleted": earlier + doomed, "bytes_deleted": sum(d["bytes"] for d in earlier + doomed)}
    _write_json(log_path, out)
    commit()
    return out
