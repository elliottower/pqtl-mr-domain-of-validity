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

Uncertainty mode (outcome_files module docstring). For a GWAS-SSF or author-format file, a first pass
counts the `standard_error` field over the rows with a position and a valid pair of alleles
(`se_census`: missing, present, not counted), and the mode is chosen once for the whole file from
the header and those two counts (outcome_files.choose_uncertainty_mode); a harmonised file is read in
`native_se`. The second pass reads every row in that mode. No value of a row enters the choice.

Pass/fail rules, fixed before any file is read:
- the header must agree with the column map the file's record was resolved under: for a harmonised
  file, HarmonisedReader finds every column it requires; for GWAS-SSF, `ssf_columns` reads it and
  gives the record's rsID rule; for an author format, it is the header the reviewed map was written
  for (author_formats.check_header);
- an uncertainty mode must be chosen: a file whose standard_error is present on some rows and missing
  on others fails closed, as does one with no route;
- at least MIN_KEPT_FRACTION (50%) of the data rows must be read;
- an author format with an integrity rule must pass it (for GCST008226: outcome_files.OrCiCheck);
- an `or_ci_derived_se` file must pass the CI-versus-p check (outcome_files.CiPCheck) under the rule
  version outcome_files.CI_P_RULE selects (frozen or rounding-aware; under the rounding-aware rule
  the frozen rule's counts are reported beside its own).
A file that fails is made `unreadable` (`apply_validation`): its `collected` record is moved, not
deleted, to <root>/superseded/collect/gwas_catalog/, and a record with status `unreadable` takes its
place, with the reasons in `detail` after VALIDATION_DETAIL and the file's name, path, size and
sha256 kept. It then takes the registered consequence of an unavailable outcome file, as an
`unreadable` record found at collect does. A file that passes has its `collected` record superseded
the same way by one that names its `uncertainty_mode`, bound to the sha256 of the header line
(`header_sha256`, author_formats.header_sha256) and of the validation result (`validation_sha256`,
`result_sha256`); the record's sha256 and meta_sha256 are kept, and all of it enters the collect
digest (stage_b/checkpoint.py).

A file is read in one call, with a checkpoint every CHECKPOINT_ROWS rows in each pass: the counts so
far, and the rsID, chromosome, position and alleles of each row read, appended to a shard. A
restarted call resumes at the row after the checkpoint. The work directory is named for the file's
sha256 and the rules (RULES_SHA256), so work begun under other rules is never resumed. The result is
written to <root>/validation/files/<key>.json and reused while the record names the same sha256 and
the result was made under RULES; the shard is then deleted. A file recorded `unreadable` by an
earlier validation is read again, its bytes checked against the record, when its result is stale.
The report of a run is <root>/validation/outcome_validation.json (REPORT_NAME).

`restorable` and `restore_validated` undo `apply_validation` for a file whose validation failed for
its standard error alone (`failed_only_for_standard_error`: one reason, too few rows read, which
would not have held had the rows rejected for their standard error been read; launch_stage_b.py
`restore-validated`): refused unless the file's bytes still have the recorded sha256; the result and
the `unreadable` record are moved to superseded/, the `collected` record moved back, nothing
deleted. The next `validate` reads the file again under the rules in force.

`require_validation` is the guard of `spawn` (stage_b/launch.py): no unit is started unless the
report exists, was made under the rules in force (RULES), and covers every GWAS Catalog whole file
of the current collect records, each with a result for the sha256 its record names, made under those
rules, no error, a record status that agrees with the result (`unreadable` exactly where the file
failed), and, for every `collected` record, an uncertainty mode, header hash and result hash that are
those of its result.
"""
import gzip
import hashlib
import json
import os
from collections import Counter
from collections.abc import Callable, Iterable, Iterator
from contextlib import contextmanager
from datetime import datetime, timezone
from itertools import islice
from pathlib import Path

import pandas as pd

from stage_b.author_formats import AUTHOR_FORMATS, header_sha256
from stage_b.collect import SUPERSEDED_DIR, collected_file, put_record, read_record, record_path, supersede_record
from stage_b.outcome_files import (CI_LEVEL, CI_P_LOG10_TOLERANCE, CI_P_MIN_AGREEMENT_PERCENT, CI_P_RULE, CI_P_RULE_FROZEN,
                                   CI_Z, REJECT_REASONS, AuthorReader, CiPCheck, OrCiCheck, SsfReader, attach_map_rsids,
                                   choose_uncertainty_mode, map_rsids, ssf_columns)
from stage_b.parsers import HarmonisedReader, RowReader, split_lines
from stage_b.schemas import (UNCERTAINTY_MODES, CollectError, CollectRecord, InputContractError, SourceUnreadable,
                             UncertaintyMode)

MIN_KEPT_FRACTION_PERCENT = 50
CHECKPOINT_ROWS = 1_000_000
VALIDATION_DIR = "validation"
VALIDATION_DETAIL = "pre-analysis validation: "
SHARD_COLUMNS = ["rsid", "chrom", "pos", "ea", "oa"]
REPORT_NAME = "outcome_validation.json"
RULES = {"min_kept_fraction_percent": MIN_KEPT_FRACTION_PERCENT,
         "header": "agrees with the column map of the file's record",
         "integrity": "author formats with an integrity rule must pass it (outcome_files.OrCiCheck)",
         "reject_reasons": list(REJECT_REASONS),
         "uncertainty_mode": ("chosen once per file from the header and the standard_error census of the rows with a "
                              "position and a valid pair of alleles (outcome_files.choose_uncertainty_mode); a file with "
                              "standard_error present on some rows and missing on others fails"),
         "or_ci_derived_se": {"ci_level": CI_LEVEL, "z": CI_Z, "se": "(ln U - ln L) / (2 z)",
                              "or_inside_ci": "within half a unit in the last printed decimal of OR plus of the limit"},
         "ci_p_check": None}


def ci_p_rules(rule: str) -> dict:
    """The RULES entry of the CI-versus-p check under rule version `rule` (outcome_files.CI_P_RULES);
    the frozen entry is the one the frozen rule was validated under, unchanged."""
    if rule == CI_P_RULE_FROZEN:
        return {"rows": "read, with 0 < p < 1", "log10_tolerance": CI_P_LOG10_TOLERANCE,
                "or_within_the_printed_rounding_of_p": True, "min_agreement_percent": CI_P_MIN_AGREEMENT_PERCENT}
    return {"rule_version": rule,
            "rows": ("read, with OR, L and U finite and positive, L < U and 0 < p < 1 as printed; not checkable where the "
                     "rounding interval of L reaches 0, counted apart"),
            "intervals": ("each printed value +/- half a unit in its last printed digit (mantissa digit in scientific "
                          "notation); the implied two-sided p over the exact extremes of |ln OR| 2 z / (ln U - ln L) on "
                          "the box of the three intervals with ln L <= ln OR <= ln U (outcome_files.wald_ratio_extremes); "
                          "no feasible point disagrees; published p over its own interval"),
            "rounding_uninformative": ("a zero interval width reachable at a feasible OR other than 1, and a feasible OR "
                                       "range holding 1: not checkable, counted apart, outside the agreement denominator"),
            "agrees": "the published-p interval meets the implied-p interval, each widened by log10_tolerance in log10",
            "log10_tolerance": CI_P_LOG10_TOLERANCE, "min_agreement_percent": CI_P_MIN_AGREEMENT_PERCENT,
            "no_informative_row": "the file fails",
            "row_crosstab": "every row read by its frozen and rounding-aware verdicts, counts only",
            "frozen_rule_reported_beside": ci_p_rules(CI_P_RULE_FROZEN)}


RULES["ci_p_check"] = ci_p_rules(CI_P_RULE)
RULES_SHA256 = hashlib.sha256(json.dumps(RULES, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


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


def result_sha256(result: dict) -> str:
    """sha256 of a validation result as canonical JSON (sorted keys, no whitespace): the hash a
    `collected` record's `validation_sha256` binds its uncertainty mode to."""
    return hashlib.sha256(json.dumps(result, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _reader(record: CollectRecord, header: list[str], mode: UncertaintyMode) -> RowReader:
    """The analysis reader of the file under its record's column map, in `mode`; SourceUnreadable or
    InputContractError, with the reason, where the header disagrees with that map."""
    what = f"GWAS Catalog {record.key}"
    if record.layout == "":
        return HarmonisedReader(header)
    if record.layout == "author":
        fmt = AUTHOR_FORMATS.get(record.key)
        if fmt is None:
            raise SourceUnreadable(f"no reviewed author-format map for {record.key}")
        return AuthorReader(fmt, header, mode)          # check_header: the header the map was reviewed against
    rule = "column" if ssf_columns(header, what).rsid else "ukbppp_map"
    if rule != record.rsid_rule:
        raise SourceUnreadable(f"the header gives the rsID rule {rule}, the record {record.rsid_rule}")
    return SsfReader(header, record.position_offset, what, mode)


# The checks a file passes or fails beside its rows: name -> the check of a reader, or None where none applies.
CHECKS: dict[str, Callable[[CollectRecord, list[str], RowReader], OrCiCheck | CiPCheck | None]] = {
    "integrity": lambda record, header, reader: (
        OrCiCheck(AUTHOR_FORMATS[record.key], header)
        if record.layout == "author" and AUTHOR_FORMATS[record.key].integrity == "or_ci_se" else None),
    "ci_p_check": lambda record, header, reader: (
        CiPCheck(reader, CI_P_RULE) if isinstance(reader, SsfReader) and reader.mode == "or_ci_derived_se" else None),
}


@contextmanager
def _data(path: Path) -> Iterator[tuple[list[str] | None, Iterator[list[str]]]]:
    """The header line and an iterator of the data rows of a collected file."""
    with _open_text(path) as fh:
        rows = split_lines(fh, "\t")
        yield next(rows, None), rows


def _census(reader: RowReader, rows: Iterator[list[str]], work: Path, commit: Callable[[], None],
            checkpoint_rows: int) -> dict:
    """The standard_error census of every data row (RowReader.se_state): rows whose field is missing,
    present, or not counted (no position or no valid pair of alleles). Checkpointed in
    work/census.json every `checkpoint_rows` rows and resumed from there."""
    path = work / "census.json"
    state = json.loads(path.read_text()) if path.is_file() else {"rows": 0, "missing": 0, "present": 0,
                                                                     "not_counted": 0, "done": False}
    if state["done"]:
        return state
    work.mkdir(parents=True, exist_ok=True)
    counts = Counter({k: state[k] for k in ("missing", "present", "not_counted")})
    n = state["rows"]
    for r in islice(rows, n, None):
        n += 1
        counts[reader.se_state(r) or "not_counted"] += 1
        if n % checkpoint_rows == 0:
            _write_json(path, {"rows": n, **counts, "done": False})
            commit()
    _write_json(path, {"rows": n, "missing": counts["missing"], "present": counts["present"],
                       "not_counted": counts["not_counted"], "done": True})
    commit()
    return json.loads(path.read_text())


def _read_pass(reader: RowReader, checks: dict, rows: Iterator[list[str]], work: Path, commit: Callable[[], None],
               checkpoint_rows: int) -> dict:
    """Read every data row in the reader's mode, resuming from the checkpoint in `work`: the rows, the
    rows read and the rows not read by reason; each row read appended to the shard; each check fed."""
    state_path, shard = work / "state.json", work / "variants.tsv"
    state = json.loads(state_path.read_text()) if state_path.is_file() else {"rows": 0, "kept": 0, "rejected": {},
                                                                                 "shard_bytes": 0, "checks": {}}
    done, kept, rejected = state["rows"], state["kept"], Counter(state["rejected"])
    for name, check in checks.items():
        for k, v in state["checks"].get(name, {}).items():
            setattr(check, k, v)
    work.mkdir(parents=True, exist_ok=True)
    with shard.open("ab") as out:
        out.truncate(state["shard_bytes"])          # rows written after the checkpoint are written again
    with shard.open("ab") as out:
        def checkpoint(n: int) -> None:
            out.flush()
            _write_json(state_path, {"rows": n, "kept": kept, "rejected": dict(rejected), "shard_bytes": out.tell(),
                                     "checks": {name: {k: getattr(c, k) for k in c.count_names} for name, c in checks.items()}})
            commit()

        integrity, ci_p = checks.get("integrity"), checks.get("ci_p_check")
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
                if ci_p is not None:
                    ci_p.add(got, r)
                out.write(f"{got['rsid']}\t{got['chrom']}\t{got['pos']}\t{got['ea']}\t{got['oa']}\n".encode())
            if n % checkpoint_rows == 0:
                checkpoint(n)
        checkpoint(n)
    return {"rows": n, "kept": kept, "rejected": rejected}


def _scan(record: CollectRecord, path: Path, work: Path, commit: Callable[[], None], checkpoint_rows: int) -> dict:
    """The header check, the standard_error census and the uncertainty mode, then the read pass in
    that mode with its checks (module docstring)."""
    with _data(path) as (header, rows):
        if header is None:
            return {"header": [], "header_problems": ["the file has no header line"]}
        try:
            reader = _reader(record, header, "native_se")
        except (SourceUnreadable, InputContractError) as err:
            return {"header": header, "header_problems": [str(err)]}
        if record.layout == "":
            census, mode, problem = None, "native_se", ""
        else:
            census = _census(reader, rows, work, commit, checkpoint_rows)
            cols = reader.cols if isinstance(reader, SsfReader) else None
            mode, problem = choose_uncertainty_mode(cols, census["missing"], census["present"])
    out = {"header": header, "header_problems": [], "uncertainty_mode": mode or "",
           "se_census": None if census is None else {k: census[k] for k in ("rows", "missing", "present", "not_counted")}}
    if mode is None:
        return {**out, "mode_problem": problem}
    with _data(path) as (_, rows):
        try:
            reader = _reader(record, header, mode)
            checks = {name: c for name, make in CHECKS.items() if (c := make(record, header, reader)) is not None}
        except (SourceUnreadable, InputContractError) as err:
            return {**out, "header_problems": [str(err)]}
        read = _read_pass(reader, checks, rows, work, commit, checkpoint_rows)
    return {**out, **read, **{name: (checks[name].result() if name in checks else None) for name in CHECKS}}


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
        if done["sha256"] == record.sha256 and done.get("rules") == RULES:
            return done
        stale = root / VALIDATION_DIR / SUPERSEDED_DIR / f"{record.key}__{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}.json"
        stale.parent.mkdir(parents=True, exist_ok=True)
        os.replace(out, stale)
    # the bytes the record names, whatever its status: a file an earlier validation made unreadable is read again
    path = collected_file(root, record.model_copy(update={"status": "collected"}))
    work = root / VALIDATION_DIR / "work" / f"{record.key}__{record.sha256[:16]}__{RULES_SHA256[:12]}"
    scan = _scan(record, path, work, commit, checkpoint_rows)
    result = {"key": record.key, "layout": record.layout, "build": record.build, "rsid_rule": record.rsid_rule,
              "position_offset": record.position_offset, "name": record.name, "sha256": record.sha256,
              "header": scan["header"], "header_sha256": header_sha256(scan["header"]),
              "header_agrees": not scan["header_problems"], "rules": RULES, "utc": _utc()}
    reasons = list(scan["header_problems"])
    if not reasons:
        result.update({"uncertainty_mode": scan["uncertainty_mode"], "se_census": scan["se_census"]})
        if "mode_problem" in scan:
            reasons.append(f"uncertainty mode: {scan['mode_problem']}")
    if not reasons:
        rows, kept, rejected = scan["rows"], scan["kept"], scan["rejected"]
        result.update({"rows": rows, "rows_kept": kept, "row_width_errors": rejected["row_width"],
                       "parse_errors": rejected["unparseable_number"],
                       "rejected": {reason: rejected[reason] for reason in REJECT_REASONS},
                       "integrity": scan["integrity"], "ci_p_check": scan["ci_p_check"],
                       **_rsid_stats(root, record, work / "variants.tsv")})
        if rows == 0 or 100 * kept < MIN_KEPT_FRACTION_PERCENT * rows:
            reasons.append(f"{kept} of {rows} data rows read, fewer than {MIN_KEPT_FRACTION_PERCENT}%")
        integrity, ci_p = scan["integrity"], scan["ci_p_check"]
        if integrity is not None and not integrity["passed"]:
            reasons.append(f"OR/CI/SE integrity: {integrity['rows_agree']} of {integrity['rows_checked']} checked rows "
                           f"agree, fewer than {integrity['min_agreement_percent']}%")
        if ci_p is not None and not ci_p["passed"]:
            label = "CI-versus-p check" if ci_p["rule"] == "or_ci_vs_p" else "rounding-aware CI-versus-p check"
            if ci_p["rows_checked"] == 0 and ci_p["rule"] != "or_ci_vs_p":
                reasons.append(f"{label}: no informative row ({ci_p['rows_rounding_uninformative']} rounding-uninformative, "
                               f"{ci_p['rows_not_checkable']} not checkable in all)")
            else:
                reasons.append(f"{label}: {ci_p['rows_agree']} of {ci_p['rows_checked']} checkable rows agree, "
                               f"fewer than {ci_p['min_agreement_percent']}%")
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
    """A collected record whose file failed: superseded by an `unreadable` record. A collected record
    whose file passed: superseded by one naming the uncertainty mode of `result`, bound to its header
    hash and its result hash, unless it already names exactly those (module docstring). Any other
    record is returned unchanged."""
    if record.status != "collected":
        return record
    if result["passed"]:
        bound = record.model_copy(update={"uncertainty_mode": result["uncertainty_mode"],
                                          "header_sha256": result["header_sha256"], "validation_sha256": result_sha256(result)})
        if bound == record:
            return record
        supersede_record(root, record.source, record.key, commit)
        return put_record(root, bound, commit)
    supersede_record(root, record.source, record.key, commit)
    detail = VALIDATION_DETAIL + "; ".join(result["reasons"]) + (f" ({record.detail})" if record.detail else "")
    return put_record(root, record.model_copy(update={"status": "unreadable", "detail": detail, "utc": _utc(),
                                                      "uncertainty_mode": "", "header_sha256": "", "validation_sha256": ""}),
                      commit)


def failed_only_for_standard_error(result: dict) -> bool:
    """True when a validation result failed for one reason only, fewer than MIN_KEPT_FRACTION_PERCENT
    of its rows read, and at least that fraction would have been read had the rows rejected as
    `se_nonpositive_or_nonfinite` been read: the header agreed, no integrity rule failed, and every
    other rule held."""
    rows, kept = result.get("rows", 0), result.get("rows_kept", 0)
    se = (result.get("rejected") or {}).get("se_nonpositive_or_nonfinite", 0)
    return (result.get("passed") is False and result.get("header_agrees") is True and result.get("integrity") is None
            and result.get("reasons") == [f"{kept} of {rows} data rows read, fewer than {MIN_KEPT_FRACTION_PERCENT}%"]
            and se > 0 and 100 * (kept + se) >= MIN_KEPT_FRACTION_PERCENT * rows)


def _failed_result(root: Path, key: str, sha256: str) -> dict:
    """The validation result that made `key` unreadable: the current one, or, after a restore stopped
    once it had archived the result, the archived one; CollectError when neither names `sha256`."""
    found = [result_path(root, key)] if result_path(root, key).is_file() else sorted(
        (root / VALIDATION_DIR / SUPERSEDED_DIR).glob(f"{key}__*__restored.json"))
    results = [r for r in (json.loads(p.read_text()) for p in found) if r.get("sha256") == sha256]
    if not results:
        raise CollectError(f"{key}: no validation result for its sha256; nothing is restored")
    return results[-1]


def restorable(root: Path, key: str) -> dict:
    """What `restore_validated` would move for GWAS Catalog `key`, after every check, moving nothing.
    Refused (CollectError) unless: the current record is `unreadable` by this validation, or is gone
    after an interrupted restore; exactly one superseded `collected` record of `key` names the same
    name, path, size and sha256; the file on the volume still has that size and sha256; and the
    validation result that made it unreadable failed for its standard error alone
    (`failed_only_for_standard_error`)."""
    path = record_path(root, "gwas_catalog", key)
    current = read_record(root, "gwas_catalog", key) if path.is_file() else None
    if current is not None and not (current.status == "unreadable" and current.detail.startswith(VALIDATION_DETAIL)):
        raise CollectError(f"{key}: the record is {current.status}, not made unreadable by the pre-analysis validation")
    found = []
    for p in sorted((root / SUPERSEDED_DIR / "collect" / "gwas_catalog").glob(f"{key}__*.json")):
        r = CollectRecord.model_validate_json(p.read_text())
        if r.source == "gwas_catalog" and r.key == key and r.status == "collected" and (
                current is None or (r.name, r.path, r.bytes, r.sha256) == (current.name, current.path, current.bytes,
                                                                           current.sha256)):
            found.append((p, r))
    if len(found) != 1:
        raise CollectError(f"{key}: {len(found)} superseded collected records of the same file, not exactly one")
    superseded, collected = found[0]
    collected_file(root, collected)          # CollectError unless the bytes on the volume are the recorded ones
    if not failed_only_for_standard_error(_failed_result(root, key, collected.sha256)):
        raise CollectError(f"{key}: its validation did not fail for the standard error alone; it is not restored")
    result = result_path(root, key)
    return {"key": key, "sha256": collected.sha256, "current": current, "superseded": superseded,
            "result": result if result.is_file() else None}


def restore_validated(root: Path, plan: dict, commit: Callable[[], None] = lambda: None) -> dict:
    """Undo `apply_validation` for a file whose failure was a reader defect (plan from `restorable`):
    its validation result is moved to <root>/validation/superseded/ (so the next validation reads the
    file again), the `unreadable` record to <root>/superseded/collect/gwas_catalog/, and the superseded
    `collected` record back where read_record finds it. Nothing is deleted."""
    key, stamp = plan["key"], f"{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}"
    moved = {}

    def archive(src: Path, dest: Path, what: str) -> None:
        if dest.exists():
            raise CollectError(f"{dest} already exists; an archived file is never overwritten")
        dest.parent.mkdir(parents=True, exist_ok=True)
        os.replace(src, dest)
        commit()
        moved[what] = str(dest.relative_to(root))

    if plan["result"] is not None:
        archive(plan["result"], root / VALIDATION_DIR / SUPERSEDED_DIR / f"{key}__{stamp}__restored.json",
                "validation_result")
    if plan["current"] is not None:
        archive(record_path(root, "gwas_catalog", key),
                root / SUPERSEDED_DIR / "collect" / "gwas_catalog" / f"{key}__{stamp}__unreadable.json", "unreadable_record")
    os.replace(plan["superseded"], record_path(root, "gwas_catalog", key))
    commit()
    moved["collected_record"] = str(plan["superseded"].relative_to(root))
    return {"key": key, "sha256": plan["sha256"], "moved": moved}


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


def frozen_ci_p_verdicts(report: dict) -> dict[str, bool]:
    """GWAS Catalog accession -> whether its file passed the frozen CI-versus-p rule, for every file of
    `report` the check ran on (an or_ci_derived_se file): under the frozen rule its own verdict, under
    the rounding-aware rule the frozen verdict reported beside it, whether the file passed or failed
    the rule in force. Stage B writes it per hypothesis (evidence.csv `outcome_file_frozen_ci_p_pass`)
    so that stage D can form the frozen-rule sensitivity set without the validation report."""
    out = {}
    for key, got in report["files"].items():
        check = got.get("ci_p_check") if isinstance(got, dict) else None
        if check is not None:
            out[key] = bool((check["frozen_rule"] if check["rule"] == "or_ci_vs_p_rounding_aware" else check)["passed"])
    return out


CLASSES = (*UNCERTAINTY_MODES, "failed", "error")


def classification(report: dict) -> dict[str, list[str]]:
    """The validated files of a report by outcome: the uncertainty mode a file passed in, `failed`, or
    `error`. Accessions only."""
    out: dict[str, list[str]] = {c: [] for c in CLASSES}
    for key, got in sorted(report["files"].items()):
        out["error" if not isinstance(got, dict) else got["uncertainty_mode"] if got["passed"] else "failed"].append(key)
    return out


def check_classification(report: dict, expected: dict[str, list[str]]) -> dict:
    """The classification of `report` against `expected` (class -> accessions; a class not given is
    expected empty): the counts of both and, per class, the accessions in one and not the other."""
    unknown = sorted(set(expected) - set(CLASSES))
    if unknown:
        raise InputContractError(f"the expected classification names classes {unknown}; the classes are {list(CLASSES)}")
    got = classification(report)
    mismatches = {}
    for c in CLASSES:
        want, have = set(expected.get(c, [])), set(got[c])
        if want != have:
            mismatches[c] = {"expected_not_found": sorted(want - have), "found_not_expected": sorted(have - want)}
    return {"counts": {c: len(got[c]) for c in CLASSES}, "expected_counts": {c: len(expected.get(c, [])) for c in CLASSES},
            "mismatches": mismatches, "agrees": not mismatches}


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
        elif got.get("rules") != RULES:
            problems.append(f"{r.key}: its result was made under other rules than the ones in force")
        elif (r.status == "unreadable") == bool(got.get("passed")):
            problems.append(f"{r.key}: the record is {r.status} and the report says the file "
                            f"{'passed' if got.get('passed') else 'failed'}")
        elif r.status == "collected" and not (r.uncertainty_mode and r.uncertainty_mode == got.get("uncertainty_mode")
                                              and r.header_sha256 == got.get("header_sha256")
                                              and r.validation_sha256 == result_sha256(got)):
            problems.append(f"{r.key}: the record's uncertainty mode ({r.uncertainty_mode or 'none'}) is not bound to its "
                            f"validation result (mode {got.get('uncertainty_mode') or 'none'})")
    if problems:
        raise InputContractError(f"the pre-analysis validation report does not cover the current collect records "
                                 f"({len(problems)} files; first: {problems[0]}); run `launch_stage_b.py validate` again")
