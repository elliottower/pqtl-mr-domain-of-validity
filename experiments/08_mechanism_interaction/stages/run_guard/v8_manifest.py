"""The one MANIFEST.tsv / INPUTS.tsv format, written by stages A, B, C (and D) and read by D.

An output directory holds its tables, then:

INPUTS.tsv    one row per input: name, path, sha256, source, release, url_or_accession,
              download_date (the last four empty where not recorded). Paths are repo-relative
              (`experiments/08_mechanism_interaction/...`) or volume-relative
              (`pqtl-v8-inputs:aact_20260930/studies.txt`), never absolute.
MANIFEST.tsv  one row per file of the directory, INPUTS.tsv included: path (relative to the
              directory), rows (data rows of a .csv/.tsv, empty otherwise), sha256, script_sha256,
              utc_time.

The SEAL entry of a stage (`SEAL stage=X manifest_sha256=<h>`) is the sha256 of its MANIFEST.tsv.
Because MANIFEST.tsv lists INPUTS.tsv with its sha256, the seal covers the input provenance too;
stage D checks INPUTS.tsv's hash like any listed file but never treats the inputs it names as
outputs. INPUTS.tsv is written first, MANIFEST.tsv last.

The output directory itself is a directory, not a symlink to one: `write_manifest` and
`verify_output_dir` refuse a symlinked root, so a seal names the directory it was written in.

Reading is strict (`read_manifest`, `read_inputs`, `verify_output_dir`):
- every sha256 and script_sha256 is 64 lowercase hex characters;
- a listed path is relative, has no symlink anywhere below the output directory and resolves inside
  it; the same holds for MANIFEST.tsv and INPUTS.tsv themselves;
- MANIFEST.tsv lists each path once, never itself, carries one script_sha256 for the whole stage,
  and every utc_time is an ISO 8601 time at UTC offset zero;
- INPUTS.tsv parses with the columns above, and every input name is non-empty and unique;
- `verify_output_dir` refuses a file under the output directory that MANIFEST.tsv does not list,
  unless the caller names it in `side_files` (a file a stage is known to leave beside its outputs).

`code_sha256` is the digest written as script_sha256: over source files keyed by a canonical
relative POSIX path, sorted by that path, each as (path length, path, file length, bytes).
`code_record` is the INPUTS.tsv row that carries the digest and the repository commit of the run.
"""
import csv
import hashlib
import os
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, fields
from datetime import datetime, timedelta, timezone
from pathlib import Path, PurePosixPath

MANIFEST_NAME = "MANIFEST.tsv"
INPUTS_NAME = "INPUTS.tsv"
RUN_INFO_NAME = "run_info.json"
STAGES_REPO_PATH = "experiments/08_mechanism_interaction/stages"
CODE_RECORD_NAME = "stage_code"
MANIFEST_COLUMNS = ("path", "rows", "sha256", "script_sha256", "utc_time")
TABULAR = {".csv": ",", ".tsv": "\t"}
SHA256_RE = re.compile(r"[0-9a-f]{64}")
COMMIT_RE = re.compile(r"[0-9a-f]{40}")


class ManifestError(RuntimeError):
    """A MANIFEST.tsv or INPUTS.tsv that does not follow the format, or a directory it does not describe."""


@dataclass(frozen=True)
class InputRecord:
    name: str
    path: str
    sha256: str
    source: str = ""
    release: str = ""
    url_or_accession: str = ""
    download_date: str = ""


INPUTS_COLUMNS = tuple(f.name for f in fields(InputRecord))


@dataclass(frozen=True)
class ManifestEntry:
    path: str
    rows: int | None
    sha256: str
    script_sha256: str
    utc_time: str


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def relative_files(root: Path, files: Iterable[Path]) -> dict[str, Path]:
    """relative POSIX path -> file, for files under `root`. A file outside `root`, or two files with
    one relative path, raises."""
    out: dict[str, Path] = {}
    base = root.resolve()
    for f in files:
        try:
            rel = f.resolve().relative_to(base).as_posix()
        except ValueError as err:
            raise ManifestError(f"{f} is not under {root}") from err
        if rel in out:
            raise ManifestError(f"{rel} is listed twice")
        out[rel] = f
    return out


def guard_files() -> dict[str, Path]:
    """The shared guard modules a stage runs with, keyed as they sit in the repository
    (`run_guard/<name>`), wherever this copy was loaded from (stages/run_guard, or beside the
    stage's modules in a Modal image)."""
    here = Path(__file__).resolve().parent
    return {f"run_guard/{name}": here / name for name in ("v8_manifest.py", "v8_run_guard.py")}


def code_sha256(files: Mapping[str, Path]) -> str:
    """One hash over source files keyed by canonical relative POSIX path, in the order of that
    path. Each file contributes the byte length of its path (8 bytes, big-endian), the path, the
    byte length of its content (8 bytes, big-endian) and the content, so neither a directory layout
    nor a boundary between path and bytes is ambiguous."""
    h = hashlib.sha256()
    for rel in sorted(files):
        _relative(rel)
        if PurePosixPath(rel).as_posix() != rel:
            raise ManifestError(f"{rel!r} is not a canonical relative POSIX path")
        name, data = rel.encode("utf-8"), files[rel].read_bytes()
        h.update(len(name).to_bytes(8, "big"))
        h.update(name)
        h.update(len(data).to_bytes(8, "big"))
        h.update(data)
    return h.hexdigest()


def code_record(stage: str, digest: str, repo_commit: str) -> InputRecord:
    """The INPUTS.tsv row naming the code of a run: the stage directory, its `code_sha256`, and the
    repository commit the run was started from (in `release`)."""
    if COMMIT_RE.fullmatch(repo_commit) is None:
        raise ManifestError(f"repository commit {repo_commit!r} is not 40 lowercase hex characters")
    return InputRecord(name=CODE_RECORD_NAME, path=f"{STAGES_REPO_PATH}/{stage}", sha256=_sha256(digest, "code digest"),
                       source="git commit", release=repo_commit)


def count_rows(path: Path) -> int | None:
    """Data rows of a .csv / .tsv (header excluded, quoted newlines respected); None otherwise."""
    sep = TABULAR.get(path.suffix)
    if sep is None:
        return None
    with path.open(newline="", encoding="utf-8") as f:
        return max(sum(1 for _ in csv.reader(f, delimiter=sep)) - 1, 0)


def _relative(path: str) -> str:
    p = PurePosixPath(path)
    if not path or p.is_absolute() or ".." in p.parts or path.startswith("\\") or ":\\" in path:
        raise ManifestError(f"path {path!r} must be relative, without '..'")
    return path


def _sha256(value: str, what: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise ManifestError(f"{what} {value!r} is not 64 lowercase hex characters")
    return value


def listed_file(out_dir: Path, rel: str) -> Path:
    """`out_dir / rel` for a listed path: no component below `out_dir` is a symlink, and the
    resolved file lies inside the resolved `out_dir`."""
    _relative(rel)
    p = out_dir
    for part in PurePosixPath(rel).parts:
        p = p / part
        if p.is_symlink():
            raise ManifestError(f"{out_dir}: listed path {rel} goes through the symlink {p}")
    if not p.resolve().is_relative_to(out_dir.resolve()):
        raise ManifestError(f"{out_dir}: listed path {rel} resolves outside the output directory")
    return p


def plain_root(out_dir: Path) -> Path:
    """`out_dir`, refused when it is itself a symlink."""
    if out_dir.is_symlink():
        raise ManifestError(f"output directory {out_dir} is a symlink to {out_dir.resolve()}; "
                            "a sealed output directory is not a symlink")
    return out_dir


def _portable(path: str) -> str:
    """An INPUTS.tsv path: repo-relative, or `<volume>:<relative path>`."""
    label, sep, rest = path.partition(":")
    if sep and label and "/" not in label:
        _relative(rest)
    else:
        _relative(path)
    return path


def portable_path(path: Path, roots: Sequence[tuple[str, Path]]) -> str:
    """`path` relative to the first root containing it, prefixed by that root's label: a label
    ending in ':' is a volume (`pqtl-v8-inputs:` + rel), any other a repository directory
    (`label/rel`), the empty label the repository root. Raises when no root contains it."""
    resolved = path.resolve()
    for label, root in roots:
        try:
            rel = resolved.relative_to(root.resolve()).as_posix()
        except ValueError:
            continue
        return f"{label}{rel}" if label.endswith(":") else (f"{label}/{rel}" if label else rel)
    raise ManifestError(f"{path} is under none of the roots {[label or '<repo>' for label, _ in roots]}")


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def write_manifest(out_dir: Path, outputs: Sequence[Path], script_sha256: str, inputs: Sequence[InputRecord]) -> Path:
    """INPUTS.tsv, then MANIFEST.tsv listing every file in `outputs` and INPUTS.tsv. An `out_dir`
    that is a symlink is refused before anything is written."""
    plain_root(out_dir)
    _sha256(script_sha256, "script_sha256")
    for r in inputs:
        _portable(r.path)
        _sha256(r.sha256, f"input {r.name} sha256")
    names = [r.name for r in inputs]
    if len(set(names)) != len(names):
        raise ManifestError(f"duplicate input names in {names}")
    _write_tsv(out_dir / INPUTS_NAME, INPUTS_COLUMNS, [asdict(r) for r in inputs])
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    files = [*outputs, out_dir / INPUTS_NAME]
    rows = []
    for p in files:
        rel = _relative(p.resolve().relative_to(out_dir.resolve()).as_posix())
        if rel == MANIFEST_NAME:
            raise ManifestError("MANIFEST.tsv cannot list itself")
        n = count_rows(p)
        rows.append({"path": rel, "rows": "" if n is None else n, "sha256": sha256_file(p),
                     "script_sha256": script_sha256, "utc_time": now})
    if len({r["path"] for r in rows}) != len(rows):
        raise ManifestError("a file is listed twice")
    path = out_dir / MANIFEST_NAME
    _write_tsv(path, MANIFEST_COLUMNS, rows)
    return path


def _rows(path: Path, columns: Sequence[str]) -> list[dict[str, str]]:
    """The rows of one of the two files: not a symlink, inside its directory, the exact header,
    and every row with exactly the header's fields."""
    listed_file(path.parent, path.name)
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(columns):
            raise ManifestError(f"{path}: header {reader.fieldnames} != {list(columns)}")
        rows = list(reader)
    for i, r in enumerate(rows, start=2):
        if None in r or any(v is None for v in r.values()):
            raise ManifestError(f"{path}: line {i} does not have {len(columns)} tab-separated fields")
    return rows


def _row_count(value: str, what: str) -> int | None:
    if value == "":
        return None
    if not value.isascii() or not value.isdigit():
        raise ManifestError(f"{what} rows {value!r} is not a non-negative integer")
    return int(value)


def _utc_time(value: str, what: str) -> str:
    try:
        t = datetime.fromisoformat(value)
    except ValueError as err:
        raise ManifestError(f"{what} utc_time {value!r} is not an ISO 8601 time") from err
    if t.utcoffset() != timedelta(0):
        raise ManifestError(f"{what} utc_time {value!r} is not at UTC offset zero")
    return value


def read_manifest(path: Path) -> list[ManifestEntry]:
    if not path.is_file():
        raise ManifestError(f"{path} does not exist; the stage is not sealed")
    entries = [ManifestEntry(path=_relative(r["path"]), rows=_row_count(r["rows"], f"{path}: {r['path']}"),
                             sha256=_sha256(r["sha256"], f"{path}: {r['path']} sha256"),
                             script_sha256=_sha256(r["script_sha256"], f"{path}: {r['path']} script_sha256"),
                             utc_time=_utc_time(r["utc_time"], f"{path}: {r['path']}"))
               for r in _rows(path, MANIFEST_COLUMNS)]
    if not entries:
        raise ManifestError(f"{path} lists no files")
    paths = [e.path for e in entries]
    if MANIFEST_NAME in paths:
        raise ManifestError(f"{path} lists itself")
    twice = sorted({p for p in paths if paths.count(p) > 1})
    if twice:
        raise ManifestError(f"{path} lists {twice} more than once")
    scripts = sorted({e.script_sha256 for e in entries})
    if len(scripts) != 1:
        raise ManifestError(f"{path} carries {len(scripts)} different script_sha256 values; a stage has one")
    return entries


def read_inputs(path: Path) -> list[InputRecord]:
    if not path.is_file():
        raise ManifestError(f"{path} does not exist")
    records = [InputRecord(**r) for r in _rows(path, INPUTS_COLUMNS)]
    for r in records:
        if not r.name:
            raise ManifestError(f"{path}: an input has no name")
        _portable(r.path)
        _sha256(r.sha256, f"{path}: input {r.name} sha256")
    names = [r.name for r in records]
    twice = sorted({n for n in names if names.count(n) > 1})
    if twice:
        raise ManifestError(f"{path}: input names {twice} appear more than once")
    return records


def files_under(out_dir: Path) -> list[str]:
    """Relative POSIX path of every file and every symlink under `out_dir`; symlinks are not followed."""
    found = []
    for base, dirs, files in os.walk(out_dir, followlinks=False):
        here = Path(base)
        links = [d for d in dirs if (here / d).is_symlink()]
        found += [(here / name).relative_to(out_dir).as_posix() for name in (*files, *links)]
    return sorted(found)


def verify_output_dir(out_dir: Path, required: Sequence[str], side_files: Sequence[str] = ()) -> dict[str, str]:
    """Check every file MANIFEST.tsv lists (INPUTS.tsv among them) against its sha256 and row
    count; `required` files must be listed; INPUTS.tsv must parse. A file under `out_dir` that is
    neither listed nor MANIFEST.tsv nor named in `side_files` is refused, and so is an `out_dir`
    that is itself a symlink. Returns path -> sha256."""
    entries = read_manifest(plain_root(out_dir) / MANIFEST_NAME)
    listed = {e.path: e for e in entries}
    absent = [f for f in (*required, INPUTS_NAME) if f not in listed]
    if absent:
        raise ManifestError(f"{out_dir / MANIFEST_NAME} does not list {absent}")
    for name in side_files:
        _relative(name)
    hashes = {e.path: _check_entry(out_dir, e) for e in entries}
    read_inputs(out_dir / INPUTS_NAME)
    unlisted = sorted(set(files_under(out_dir)) - set(listed) - {MANIFEST_NAME} - set(side_files))
    if unlisted:
        raise ManifestError(f"{out_dir} holds files MANIFEST.tsv does not list: {unlisted}")
    return hashes


def verify_listed(out_dir: Path, names: Sequence[str]) -> dict[str, str]:
    """Check only the files `names` against MANIFEST.tsv: each must be listed, and its sha256 and
    row count must match. For a stage that consumes some of a predecessor's outputs (B and C read
    A's) and has only those and the MANIFEST.tsv at hand. Returns name -> sha256."""
    listed = {e.path: e for e in read_manifest(out_dir / MANIFEST_NAME)}
    absent = [n for n in names if n not in listed]
    if absent:
        raise ManifestError(f"{out_dir / MANIFEST_NAME} does not list {absent}")
    return {n: _check_entry(out_dir, listed[n]) for n in names}


def _check_entry(out_dir: Path, e: ManifestEntry) -> str:
    p = listed_file(out_dir, e.path)
    if not p.is_file():
        raise ManifestError(f"{out_dir}: listed file {e.path} is missing")
    digest = sha256_file(p)
    if digest != e.sha256:
        raise ManifestError(f"{out_dir}: {e.path} sha256 {digest} != manifest {e.sha256}")
    if count_rows(p) != e.rows:
        raise ManifestError(f"{out_dir}: {e.path} has {count_rows(p)} rows, manifest says {e.rows}")
    return digest
