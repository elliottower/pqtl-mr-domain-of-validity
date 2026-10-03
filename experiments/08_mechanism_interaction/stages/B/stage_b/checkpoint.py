"""What a unit checkpoint is bound to. A unit directory on the stage B volume outlives the call
that wrote it, so each one carries the fingerprint of the run that may resume it:

    sha256 of {the full InstrumentUnit, the source pins, the stage B code digest, the frozen plan
               hash, the tool versions, the collect digest}

`DirStore.bind` (pipeline.py) writes it to FINGERPRINT.json before the first step and refuses a
directory holding another one, so a checkpoint written for a different unit record (another
outcome list, another source file), other pinned inputs, other code, another plan, another image
or other collected bytes is never resumed. `result.json` carries the same value and assembly
checks it.

The unit record names each source file by identity: the source, the assay or SeqId, and for a
deCODE file its name, size and ETag in the pinned folder listing (schemas.InstrumentUnit, which
refuses a URL). No link and no token is in it, so a re-issued deCODE folder link or a renewed
Synapse or OpenGWAS token leaves every fingerprint, and so every finished checkpoint, as it was.

Name, size and ETag do not fix a file's bytes, so the fingerprint also covers what the collect
phase recorded. The collect digest (`collect_digest`) is the sha256 of the canonical JSON of one
entry per whole file the unit reads (`collect.collect_tasks` of the unit), sorted by source and key:

    collected       source, key, status, name, bytes, sha256 of the file on the volume; for a GWAS
                    Catalog file without a harmonised copy (a record with a `layout`), also layout,
                    build, position_offset, rsid_rule and meta_sha256, which say how it is read
    absent          source, key, status, and the name the pinned listing gives the file (or "")
    unreadable      source, key, status, name, bytes, sha256 and the reason (`detail`), which a
                    fixed rule of stage_b/outcome_files.py or stage_b/author_formats.py gave; name
                    "", bytes null and sha256 "" for a file found unreadable at collect, the file's
                    for one the pre-analysis validation found unreadable (stage_b/validate.py)
    remote_indexed  source, key, status, the file's name in the source's listing, and the pinned
                    release of the source (REMOTE_RELEASE); the file is queried by region and has
                    no whole-file hash

No address and no token is in an entry. The analyze phase opens a whole file only after its size
and sha256 match its record (collect.collected_file), so a result is bound to the bytes it was
computed from: a file collected again with other bytes under the same name, size and ETag has
another sha256 in its record, another collect digest, and another fingerprint. FINGERPRINT.json and
result.json carry the digest (`collect_sha256`); assembly recomputes it from the collect records it
is given.

Tool versions (`tool_versions`) are read where the unit runs, not copied from the image definition:
Python and the Python packages, R, coloc, susieR and jsonlite as R reports them, bcftools with the
htslib it was built against, and tabix. FINGERPRINT.json stores them beside the fingerprint, so
assembly, which runs where those tools are not installed, recomputes each unit's fingerprint from
the unit record, the pins, the code and the plan it holds and the tool versions the unit recorded;
an edited record no longer matches the fingerprint. Assembly refuses units run under different
tool versions (`common_tools`).

Source pins are the sha256 values modal_inputs_manifest.json records for the pinned files stage B
reads from the `pqtl-v8-inputs` volume (the deCODE annotation and excluded-variant files);
`verified_source_pins` hashes the mounted files and raises unless they match, so the pins in the
fingerprint are those of the bytes read. The remote endpoints are constants of fetch.py and are
covered by the code digest.
"""
import hashlib
import json
import platform
import subprocess
from collections.abc import Callable, Mapping
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from v8_manifest import SHA256_RE, code_sha256, guard_files, relative_files, sha256_file

from stage_b.collect import collect_tasks, read_record
from stage_b.schemas import (CollectError, CollectRecord, CollectSource, CollectTask, InputContractError, InstrumentUnit,
                             StaleCheckpointError)

FINGERPRINT_NAME = "FINGERPRINT.json"
# The pinned release of a source whose files are queried by region and not downloaded (PREREG §Data
# collection: "GWAS Catalog r2026-09-13").
REMOTE_RELEASE: dict[str, str] = {"gwas_catalog": "r2026-09-13"}
PINNED_SOURCES = ("inputs/decode/assocvariants.annotated.txt.gz", "inputs/decode/assocvariants.excluded.txt.gz")
PY_PACKAGES = ("numpy", "pandas", "pydantic", "requests", "openpyxl", "synapseclient", "prereg", "provenance-core", "PyYAML")
R_PACKAGES = ("coloc", "susieR", "jsonlite")
R_VERSIONS = ("cat(as.character(getRversion()), "
              + ", ".join(f'as.character(packageVersion("{p}"))' for p in R_PACKAGES) + ', sep = "\\n")')


def package_files(package_dir: Path) -> dict[str, Path]:
    """The stage B sources by relative path (`stage_b/<module>.py`, `stage_b/coloc_run.R`) and the
    shared guard modules (`run_guard/...`)."""
    return {**relative_files(package_dir.parent, [*sorted(package_dir.glob("*.py")), package_dir / "coloc_run.R"]),
            **guard_files()}


def package_sha256(package_dir: Path) -> str:
    """One hash over `package_files`, each by relative path, length and bytes (v8_manifest.code_sha256)."""
    return code_sha256(package_files(package_dir))


def run_tool(cmd: list[str]) -> str:
    """stdout of a version command; a tool that is absent or fails raises."""
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as err:
        raise InputContractError(f"{cmd[0]} could not be run for its version: {err}") from err
    if res.returncode != 0:
        raise InputContractError(f"{' '.join(cmd[:2])} exited {res.returncode}: {res.stderr[-500:]}")
    return res.stdout


def tool_versions(run: Callable[[list[str]], str] = run_tool,
                  package_version: Callable[[str], str] = version) -> dict[str, str]:
    """The versions installed where this is called: Python, PY_PACKAGES, R and R_PACKAGES as R
    reports them, bcftools and the htslib it uses, and tabix. Anything unreadable raises."""
    out = {"python": platform.python_version()}
    for p in PY_PACKAGES:
        try:
            out[f"python:{p}"] = package_version(p)
        except PackageNotFoundError as err:
            raise InputContractError(f"Python package {p} is not installed; its version cannot be recorded") from err
    r = run(["Rscript", "-e", R_VERSIONS]).split()
    if len(r) != 1 + len(R_PACKAGES):
        raise InputContractError(f"Rscript reported {r}, not the versions of R and {list(R_PACKAGES)}")
    out.update({"R": r[0], **{f"R:{p}": v for p, v in zip(R_PACKAGES, r[1:])}})
    bcftools = run(["bcftools", "--version"]).splitlines()
    tabix = run(["tabix", "--version"]).splitlines()
    if len(bcftools) < 2 or not bcftools[0].startswith("bcftools ") or not bcftools[1].startswith("Using htslib "):
        raise InputContractError(f"unexpected `bcftools --version` output: {bcftools[:2]}")
    if not tabix or not tabix[0].startswith("tabix (htslib) "):
        raise InputContractError(f"unexpected `tabix --version` output: {tabix[:1]}")
    out.update({"bcftools": bcftools[0].split(" ", 1)[1], "bcftools:htslib": bcftools[1].removeprefix("Using htslib "),
                "tabix:htslib": tabix[0].removeprefix("tabix (htslib) ")})
    return out


def source_pins(inputs_manifest: Path) -> dict[str, str]:
    """path -> pinned sha256 of PINNED_SOURCES, from modal_inputs_manifest.json."""
    listed = {f["path"]: f["sha256"] for f in json.loads(inputs_manifest.read_text())["files"]}
    missing = [p for p in PINNED_SOURCES if p not in listed]
    if missing:
        raise InputContractError(f"{inputs_manifest} pins no sha256 for {missing}")
    return {p: listed[p] for p in PINNED_SOURCES}


def verified_source_pins(inputs_manifest: Path, experiment_root: Path) -> dict[str, str]:
    """`source_pins`, after each pinned file under `experiment_root` is hashed and matches."""
    pins = source_pins(inputs_manifest)
    for rel, pinned in pins.items():
        got = sha256_file(experiment_root / rel)
        if got != pinned:
            raise InputContractError(f"{experiment_root / rel}: sha256 {got} differs from the pin {pinned}")
    return pins


def collect_entry(task: CollectTask, record: CollectRecord) -> dict:
    """What the collect digest keeps of one collect record: identity, terminal status and, for a
    collected file, its bytes and sha256. Never an address."""
    what = f"{task.source} {task.key}"
    if (record.source, record.key) != (task.source, task.key):
        raise CollectError(f"the record of {record.source} {record.key} was given for {what}")
    entry = {"source": task.source, "key": task.key, "status": record.status}
    if record.status == "collected":
        if not record.name or record.bytes is None or SHA256_RE.fullmatch(record.sha256) is None:
            raise CollectError(f"{what}: the record of a collected file lacks its name, size or sha256")
        if task.name and record.name != task.name:
            raise CollectError(f"{what}: the record names the file {record.name}, the pinned listing {task.name}")
        entry = {**entry, "name": record.name, "bytes": record.bytes, "sha256": record.sha256}
        if record.layout:
            entry.update({k: getattr(record, k) for k in ("layout", "build", "position_offset", "rsid_rule", "meta_sha256")})
        return entry
    if record.status == "absent":
        return {**entry, "name": task.name, "bytes": None, "sha256": ""}
    if record.status == "unreadable":
        return {**entry, "name": record.name, "bytes": record.bytes, "sha256": record.sha256, "reason": record.detail}
    if record.source not in REMOTE_RELEASE or not record.name:
        raise CollectError(f"{what}: a file queried by region needs its name and the pinned release of its source")
    return {**entry, "name": record.name, "bytes": None, "sha256": "", "release": REMOTE_RELEASE[record.source]}


def collect_digest(unit: InstrumentUnit, record_of: Callable[[CollectSource, str], CollectRecord]) -> str:
    """The sha256 over the collect records of every whole file `unit` reads (module docstring).
    `record_of(source, key)` returns a file's record and raises when there is none."""
    entries = [collect_entry(t, record_of(t.source, t.key)) for t in collect_tasks([unit])]
    return hashlib.sha256(json.dumps(entries, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def volume_collect_digest(root: Path, unit: InstrumentUnit) -> str:
    """`collect_digest` from the records under `root`/collect; CollectError when one is missing."""
    return collect_digest(unit, lambda source, key: read_record(root, source, key))


def unit_fingerprint(unit: InstrumentUnit, pins: Mapping[str, str], code_digest: str, plan_sha256: str,
                     tools: Mapping[str, str], collect_sha256: str) -> str:
    """`tools` is `tool_versions()` of the container the unit runs in (or, at assembly, the record
    the unit left in FINGERPRINT.json); `collect_sha256` is `collect_digest` of the unit."""
    if not tools:
        raise InputContractError("a unit fingerprint needs the tool versions of the image")
    if SHA256_RE.fullmatch(collect_sha256) is None:
        raise InputContractError("a unit fingerprint needs the collect digest of the unit's files")
    record = {"unit": unit.model_dump(mode="json"), "source_pins": dict(sorted(pins.items())),
              "code_sha256": code_digest, "plan_sha256": plan_sha256, "tools": dict(sorted(tools.items())),
              "collect_sha256": collect_sha256}
    return hashlib.sha256(json.dumps(record, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def unit_tools(unit_dir: Path) -> dict[str, str]:
    """The tool versions a unit directory's FINGERPRINT.json records."""
    path = unit_dir / FINGERPRINT_NAME
    tools = json.loads(path.read_text()).get("tools") if path.is_file() else None
    if not isinstance(tools, dict) or not tools:
        raise StaleCheckpointError(f"{path} records no tool versions; the unit was not written by this pipeline")
    return tools


def common_tools(tools_by_unit: Mapping[str, Mapping[str, str]]) -> dict[str, str]:
    """The one tool-version record every unit shares; units run under different versions raise."""
    distinct = {json.dumps(dict(sorted(t.items())), sort_keys=True) for t in tools_by_unit.values()}
    if len(distinct) != 1:
        raise StaleCheckpointError(f"{len(tools_by_unit)} units were run under {len(distinct)} different sets of tool "
                                   "versions; stage B assembles units of one image only")
    return json.loads(distinct.pop())
