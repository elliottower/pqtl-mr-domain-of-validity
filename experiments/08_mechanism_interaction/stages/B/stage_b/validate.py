"""Pre-analysis validation of every GWAS Catalog whole file on the stage B volume.

After the collect phase and before any unit runs, every `collected` GWAS Catalog record is read
whole by the reader the analysis uses, on Modal (modal_stage_b.py::validate_outcomes): a harmonised
file (no `layout`) by parsers.HarmonisedReader, a GWAS-SSF file by outcome_files.SsfReader, a
reviewed author format by outcome_files.AuthorReader. A `remote_indexed` file is queried by region
through its tabix index and never held whole, so it is not validated; the report lists it.
The result per file holds aggregate parser diagnostics only: the header, the number of rows, the
rows read, the rows not read by reason (outcome_files.REJECT_REASONS, `row_width` and
`unparseable_number` among them), the fraction of rows read that carry an rsID, the number of rsIDs
on more than one row read, and, for an author format with an integrity rule, that rule's counts and
pass/fail (outcome_files.OrCiCheck). No per-variant value and no hypothesis is in it. The
harmonised reader applies no allele, standard-error or p-value rule, so only `row_width`,
`no_position` and `unparseable_number` occur for a harmonised file.

Pass/fail rules, fixed before any file is read:
- the header must agree with the column map the file's record was resolved under: for a harmonised
  file, HarmonisedReader finds every column it requires; for GWAS-SSF, `ssf_columns` reads it and
  gives the record's rsID rule; for an author format, it is the header the reviewed map was written
  for (author_formats.check_header);
- at least MIN_KEPT_FRACTION (50%) of the data rows must be read;
- an author format with an integrity rule must pass it (for GCST008226: outcome_files.OrCiCheck).
A file that fails is made `unreadable` (`apply_validation`): its `collected` record is moved, not
deleted, to <root>/superseded/collect/gwas_catalog/, and a record with status `unreadable` takes its
place, with the reasons in `detail` after VALIDATION_DETAIL and the file's name, path, size and
sha256 kept. It then takes the registered consequence of an unavailable outcome file, as an
`unreadable` record found at collect does.

A file is read in one call, with a checkpoint every CHECKPOINT_ROWS rows: the counts so far, and the
rsID, chromosome, position and alleles of each row read, appended to a shard. A restarted call
resumes at the row after the checkpoint. The result is written to <root>/validation/files/<key>.json
and reused while the record names the same sha256; the shard is then deleted. The report of a run is
<root>/validation/outcome_validation.json (REPORT_NAME).

`require_validation` is the guard of `spawn` (stage_b/launch.py): no unit is started unless the
report exists, was made under the rules in force (RULES), and covers every GWAS Catalog whole file
of the current collect records, each with a result for the sha256 its record names, no error, and
a record status that agrees with the result (`unreadable` exactly where the file failed).
"""
import gzip
import json
import os
from collections import Counter
from collections.abc import Callable, Iterable
from datetime import datetime, timezone
from itertools import islice
from pathlib import Path

import pandas as pd

from stage_b.author_formats import AUTHOR_FORMATS, check_header
from stage_b.collect import SUPERSEDED_DIR, collected_file, put_record, read_record, record_path, supersede_record
from stage_b.outcome_files import (REJECT_REASONS, AuthorReader, OrCiCheck, SsfReader, attach_map_rsids, map_rsids,
                                   ssf_columns)
from stage_b.parsers import HarmonisedReader, RowReader, split_lines
from stage_b.schemas import CollectError, CollectRecord, InputContractError, SourceUnreadable

MIN_KEPT_FRACTION_PERCENT = 50
CHECKPOINT_ROWS = 1_000_000
VALIDATION_DIR = "validation"
VALIDATION_DETAIL = "pre-analysis validation: "
SHARD_COLUMNS = ["rsid", "chrom", "pos", "ea", "oa"]
REPORT_NAME = "outcome_validation.json"
RULES = {"min_kept_fraction_percent": MIN_KEPT_FRACTION_PERCENT,
         "header": "agrees with the column map of the file's record",
         "integrity": "author formats with an integrity rule must pass it (outcome_files.OrCiCheck)",
         "reject_reasons": list(REJECT_REASONS)}


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=True))
    os.replace(tmp, path)


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.name.endswith(".gz") else path.open("rt", encoding="utf-8")


def gwas_catalog_records(root: Path) -> list[CollectRecord]:
    """Every GWAS Catalog collect record under `root`, sorted by key."""
    paths = sorted((root / "collect" / "gwas_catalog").glob("*.json"))
    return [CollectRecord.model_validate_json(p.read_text()) for p in paths]


def in_scope(record: CollectRecord) -> bool:
    """A GWAS Catalog whole file on the volume: collected, or made unreadable here before."""
    return record.source == "gwas_catalog" and (
        record.status == "collected" or (record.status == "unreadable" and record.detail.startswith(VALIDATION_DETAIL)))


def result_path(root: Path, key: str) -> Path:
    return root / VALIDATION_DIR / "files" / f"{key}.json"


def _header_check(record: CollectRecord, header: list[str]) -> tuple[RowReader | None, OrCiCheck | None, list[str]]:
    """The reader of the file under its record's column map, its integrity check, and the reasons the
    header disagrees with that map (empty when it agrees)."""
    what = f"GWAS Catalog {record.key}"
    try:
        if record.layout == "":
            return HarmonisedReader(header), None, []
        if record.layout == "author":
            fmt = AUTHOR_FORMATS.get(record.key)
            if fmt is None:
                return None, None, [f"no reviewed author-format map for {record.key}"]
            check_header(fmt, header)
            reader = AuthorReader(fmt, header)
            return reader, OrCiCheck(fmt, header) if fmt.integrity == "or_ci_se" else None, []
        rule = "column" if ssf_columns(header, what).rsid else "ukbppp_map"
        if rule != record.rsid_rule:
            return None, None, [f"the header gives the rsID rule {rule}, the record {record.rsid_rule}"]
        return SsfReader(header, record.position_offset, what), None, []
    except (SourceUnreadable, InputContractError) as err:
        return None, None, [str(err)]


def _scan(record: CollectRecord, path: Path, work: Path, commit: Callable[[], None], checkpoint_rows: int) -> dict:
    """Read every data row, resuming from the checkpoint in `work`; the counts and the header check."""
    state_path, shard = work / "state.json", work / "variants.tsv"
    with _open_text(path) as fh:
        rows = split_lines(fh, "\t")
        header = next(rows, None)
        if header is None:
            return {"header": [], "header_problems": ["the file has no header line"]}
        reader, integrity, problems = _header_check(record, header)
        if reader is None:
            return {"header": header, "header_problems": problems}
        state = json.loads(state_path.read_text()) if state_path.is_file() else {"rows": 0, "kept": 0, "rejected": {},
                                                                                     "shard_bytes": 0, "integrity": None}
        done, kept, rejected = state["rows"], state["kept"], Counter(state["rejected"])
        if integrity is not None and state["integrity"]:
            integrity.checked, integrity.agree, integrity.not_checkable = (state["integrity"][k] for k in
                                                                            ("checked", "agree", "not_checkable"))
        work.mkdir(parents=True, exist_ok=True)
        with shard.open("ab") as out:
            out.truncate(state["shard_bytes"])          # rows written after the checkpoint are written again
        with shard.open("ab") as out:
            def checkpoint(n: int) -> None:
                out.flush()
                _write_json(state_path, {"rows": n, "kept": kept, "rejected": dict(rejected), "shard_bytes": out.tell(),
                                         "integrity": None if integrity is None else {
                                             "checked": integrity.checked, "agree": integrity.agree,
                                             "not_checkable": integrity.not_checkable}})
                commit()

            n = done
            for r in islice(rows, done, None):
                n += 1
                got = reader.row(r)
                if integrity is not None:
                    integrity.add(r)
                if isinstance(got, str):
                    rejected[got] += 1
                else:
                    kept += 1
                    out.write(f"{got['rsid']}\t{got['chrom']}\t{got['pos']}\t{got['ea']}\t{got['oa']}\n".encode())
                if n % checkpoint_rows == 0:
                    checkpoint(n)
            checkpoint(n)
    return {"header": header, "header_problems": [], "rows": n, "kept": kept, "rejected": rejected,
            "integrity": None if integrity is None else integrity.result()}


def _map_file(root: Path, chrom: str) -> Path | None:
    """The collected UKB-PPP rsID map of `chrom`, or None where none was collected."""
    if not record_path(root, "ukbppp_rsid_map", chrom).is_file():
        return None
    record = read_record(root, "ukbppp_rsid_map", chrom)
    return collected_file(root, record) if record.status == "collected" else None


def _rsid_stats(root: Path, record: CollectRecord, shard: Path) -> dict:
    """rsID recovery and duplicate-rsID counts over the rows read (the shard)."""
    df = pd.read_csv(shard, sep="\t", header=None, names=SHARD_COLUMNS, keep_default_na=False,
                     dtype={"rsid": str, "chrom": str, "pos": "int64", "ea": str, "oa": str})
    no_map = 0
    if record.rsid_rule == "ukbppp_map":
        parts = []
        for chrom, part in df.groupby("chrom", sort=True):
            path = _map_file(root, str(chrom))
            if path is None:
                no_map += len(part)
                parts.append(part)
                continue
            with gzip.open(path, "rt") as fh:
                rows = split_lines(fh, "\t")
                header = next(rows)
                mapping = map_rsids(rows, header, record.build, {int(p) for p in part["pos"]})
            parts.append(attach_map_rsids(part, mapping))
        df = pd.concat(parts) if parts else df
    with_rsid = df["rsid"][df["rsid"] != ""]
    counts = with_rsid.value_counts()
    return {"rows_with_rsid": len(with_rsid),
            "rsid_recovery_rate": len(with_rsid) / len(df) if len(df) else None,
            "duplicate_rsids": int((counts > 1).sum()), "rows_with_duplicate_rsid": int(counts[counts > 1].sum()),
            "rows_on_a_chromosome_without_a_collected_rsid_map": no_map}


def validate_file(root: Path, record: CollectRecord, commit: Callable[[], None] = lambda: None,
                  checkpoint_rows: int = CHECKPOINT_ROWS) -> dict:
    """The validation result of one in-scope record (module docstring), written to `result_path`."""
    if not in_scope(record):
        raise CollectError(f"{record.source} {record.key} is not a GWAS Catalog whole file on the volume")
    out = result_path(root, record.key)
    if out.is_file():
        done = json.loads(out.read_text())
        if done["sha256"] == record.sha256:
            return done
        stale = root / VALIDATION_DIR / SUPERSEDED_DIR / f"{record.key}__{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}.json"
        stale.parent.mkdir(parents=True, exist_ok=True)
        os.replace(out, stale)
    if record.status != "collected":
        raise CollectError(f"{record.key}: recorded {record.status} with no validation result for its sha256")
    path = collected_file(root, record)
    work = root / VALIDATION_DIR / "work" / f"{record.key}__{record.sha256[:16]}"
    scan = _scan(record, path, work, commit, checkpoint_rows)
    result = {"key": record.key, "layout": record.layout, "build": record.build, "rsid_rule": record.rsid_rule,
              "position_offset": record.position_offset, "name": record.name, "sha256": record.sha256,
              "header": scan["header"], "header_agrees": not scan["header_problems"], "utc": _utc()}
    reasons = list(scan["header_problems"])
    if not reasons:
        rows, kept, rejected = scan["rows"], scan["kept"], scan["rejected"]
        result.update({"rows": rows, "rows_kept": kept, "row_width_errors": rejected["row_width"],
                       "parse_errors": rejected["unparseable_number"],
                       "rejected": {reason: rejected[reason] for reason in REJECT_REASONS},
                       "integrity": scan["integrity"], **_rsid_stats(root, record, work / "variants.tsv")})
        if rows == 0 or 100 * kept < MIN_KEPT_FRACTION_PERCENT * rows:
            reasons.append(f"{kept} of {rows} data rows read, fewer than {MIN_KEPT_FRACTION_PERCENT}%")
        integrity = scan["integrity"]
        if integrity is not None and not integrity["passed"]:
            reasons.append(f"OR/CI/SE integrity: {integrity['rows_agree']} of {integrity['rows_checked']} checked rows "
                           f"agree, fewer than {integrity['min_agreement_percent']}%")
    result.update({"passed": not reasons, "reasons": reasons})
    _write_json(out, result)
    commit()
    for p in sorted(work.glob("*")) if work.is_dir() else []:
        p.unlink()
    if work.is_dir():
        work.rmdir()
    commit()
    return result


def apply_validation(root: Path, record: CollectRecord, result: dict, commit: Callable[[], None] = lambda: None) -> CollectRecord:
    """A collected record whose file failed: superseded by an `unreadable` record (module docstring).
    Any other record is returned unchanged."""
    if result["passed"] or record.status != "collected":
        return record
    supersede_record(root, record.source, record.key, commit)
    detail = VALIDATION_DETAIL + "; ".join(result["reasons"]) + (f" ({record.detail})" if record.detail else "")
    return put_record(root, record.model_copy(update={"status": "unreadable", "detail": detail, "utc": _utc()}), commit)


def validation_report(records: Iterable[CollectRecord], results: dict[str, dict | str]) -> dict:
    """The one report of a validation run: the rules, and per file its result (or the error that
    stopped it); the GWAS Catalog records not validated, with why. Aggregate diagnostics only."""
    records = list(records)
    files = {k: results[k] for k in sorted(results)}
    failed = sorted(k for k, r in files.items() if isinstance(r, dict) and not r["passed"])
    errors = sorted(k for k, r in files.items() if not isinstance(r, dict))
    return {"utc": _utc(), "rules": RULES,
            "files_validated": len(files) - len(errors), "passed": len(files) - len(failed) - len(errors),
            "failed": failed, "errors": errors,
            "unreadable_at_collect": sorted(r.key for r in records if r.status == "unreadable" and not in_scope(r)),
            "absent": sorted(r.key for r in records if r.status == "absent"),
            "not_validated_remote_indexed": sorted(r.key for r in records if r.status == "remote_indexed"),
            "not_validated_remote_indexed_reason": "queried by region through its tabix index, never held whole",
            "files": files}


def require_validation(report: dict | None, records: Iterable[CollectRecord]) -> None:
    """Raises InputContractError unless `report` (the validation report on the volume, None when
    there is none) covers the current GWAS Catalog collect `records` (module docstring)."""
    if report is None:
        raise InputContractError(f"no pre-analysis validation report ({VALIDATION_DIR}/{REPORT_NAME}); run "
                                 "`launch_stage_b.py validate` before `spawn`")
    if report.get("rules") != RULES:
        raise InputContractError("the pre-analysis validation report was made under other rules than the ones in force; "
                                 "run `launch_stage_b.py validate` again")
    files = report.get("files", {})
    problems = []
    for r in sorted((r for r in records if in_scope(r)), key=lambda r: r.key):
        got = files.get(r.key)
        if not isinstance(got, dict):
            problems.append(f"{r.key}: " + ("not in the report" if got is None else "its validation ended in an error"))
        elif got.get("sha256") != r.sha256:
            problems.append(f"{r.key}: the report validated other bytes than the record names")
        elif (r.status == "unreadable") == bool(got.get("passed")):
            problems.append(f"{r.key}: the record is {r.status} and the report says the file "
                            f"{'passed' if got.get('passed') else 'failed'}")
    if problems:
        raise InputContractError(f"the pre-analysis validation report does not cover the current collect records "
                                 f"({len(problems)} files; first: {problems[0]}); run `launch_stage_b.py validate` again")
