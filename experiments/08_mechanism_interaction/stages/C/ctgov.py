"""ClinicalTrials.gov at the outcome-freeze date, from the AACT daily flat-file snapshot dated
2026-09-30, not the live API. Its members were created at about 05:12 UTC on 2026-09-30, so it
holds ClinicalTrials.gov as available at that extraction time, approximately through the end of
2026-09-29, and not changes submitted during 2026-09-30.

PREREG §Starting and stopping rules fixes the outcome-freeze date to the freeze commit date,
2026-09-30, and "trial records dated after it do not change any label". A live record can have
been updated after the freeze; the snapshot of 2026-09-30 is fixed at that date. AACT does not
retain daily snapshots, so the copy the author downloaded and uploaded to the `pqtl-v8-inputs`
volume at /aact_20260930/ is the record; `modal_stage_c.py::verify_uploaded_aact` (logic in
aact_download.py) checks its sha256 and writes VERIFIED.json, and this module reads only after the
sha256 of `studies.txt` matches that file.

Table and columns (AACT data dictionary, `ctgov.studies`): nct_id, overall_status, why_stopped,
phase, start_month_year, start_date, start_date_type, completion_month_year, completion_date,
completion_date_type, study_first_submitted_date, last_update_submitted_date. `calculated_values`
is not needed. Flat-file format (AACT "Pipe-Delimited Flat Files Instructions"): UTF-8, `|`
between fields, header row of column names, consecutive pipes for null, fields not quoted except
a string holding an embedded `|`, which is enclosed in double quotes; line feeds inside fields
removed.

Vocabulary. The AACT dictionary lists legacy display values ("Phase 1/Phase 2", "Active, not
recruiting"); AACT loaded from ClinicalTrials.gov API v2 writes the API's enum values
("PHASE1/PHASE2", "ACTIVE_NOT_RECRUITING"). Both are mapped one-to-one to the API v2 enums the
rules use (STATUS_VALUES, PHASE_TOKENS); any other value raises, so an unexpected vocabulary stops
the run instead of changing a label.

Dates. The rules read the raw `*_month_year` strings (the date exactly as ClinicalTrials.gov gives
it) and apply the registered partial-date choice (first day of the month or year) themselves. The
AACT `start_date` / `completion_date` columns are not used: the dictionary describes their
month-year conversion both as the first and as the last day of the month. A row with an empty
`*_month_year` and a non-empty `*_date` raises. `study_first_submitted_date` and
`last_update_submitted_date` are full dates. `last_day=True` reads a partial date as the last day
of its month or year instead; only the partial_date_boundary.csv diagnostic uses it, and no label does.
"""
import calendar
import csv
import hashlib
import io
import json
import re
from datetime import date, datetime
from pathlib import Path

from models import StageCError, Trial

AACT_SNAPSHOT_DATE = "2026-09-30"
AACT_DIR_NAME = "aact_20260930"
STUDIES_FILE = "studies.txt"
VERIFIED_NAME = "VERIFIED.json"
STUDIES_COLUMNS = ("nct_id", "overall_status", "why_stopped", "phase", "start_month_year", "start_date",
                   "start_date_type", "completion_month_year", "completion_date", "completion_date_type",
                   "study_first_submitted_date", "last_update_submitted_date")

# ClinicalTrials.gov API v2 OverallStatus enums, and the legacy AACT display values of each.
STATUS_VALUES = {
    "ACTIVE_NOT_RECRUITING": "Active, not recruiting", "COMPLETED": "Completed",
    "ENROLLING_BY_INVITATION": "Enrolling by invitation", "NOT_YET_RECRUITING": "Not yet recruiting",
    "RECRUITING": "Recruiting", "SUSPENDED": "Suspended", "TERMINATED": "Terminated", "WITHDRAWN": "Withdrawn",
    "AVAILABLE": "Available", "NO_LONGER_AVAILABLE": "No longer available",
    "TEMPORARILY_NOT_AVAILABLE": "Temporarily not available", "APPROVED_FOR_MARKETING": "Approved for marketing",
    "WITHHELD": "Withheld", "UNKNOWN": "Unknown status",
}
STATUS_OF = {**{k: k for k in STATUS_VALUES}, **{v: k for k, v in STATUS_VALUES.items()}}
# API v2 Phase enums (one per token of AACT's `/`-joined `phase`), and the legacy display values.
PHASE_TOKENS = {"EARLY_PHASE1": "Early Phase 1", "PHASE1": "Phase 1", "PHASE2": "Phase 2", "PHASE3": "Phase 3",
                "PHASE4": "Phase 4"}
PHASE_OF = {**{k: k for k in PHASE_TOKENS}, **{v: k for k, v in PHASE_TOKENS.items()}}
NO_PHASE = frozenset({"", "NA", "Not Applicable", "N/A"})
MONTHS = {m: i for i, m in enumerate(("january", "february", "march", "april", "may", "june", "july", "august",
                                      "september", "october", "november", "december"), start=1)}
ISO_PARTIAL = re.compile(r"(\d{4})(?:-(\d{2}))?(?:-(\d{2}))?")
LEGACY_DATE = re.compile(r"([A-Za-z]+)(?: (\d{1,2}))?,? (\d{4})")
NCT = re.compile(r"NCT\d{8}")


class AactFormatError(StageCError):
    """An AACT row or value outside the documented format or vocabulary."""


class AactNotVerified(StageCError):
    """The AACT archive is missing, or studies.txt is not the file VERIFIED.json pins."""


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verified_studies(aact_dir: Path) -> tuple[Path, str]:
    """(studies.txt, sha256) after checking it against VERIFIED.json written by the download."""
    verified, studies = aact_dir / VERIFIED_NAME, aact_dir / STUDIES_FILE
    if not verified.is_file():
        raise AactNotVerified(f"{verified} is missing; the AACT archive is not read")
    if not studies.is_file():
        raise AactNotVerified(f"{studies} is missing")
    v = json.loads(verified.read_text())
    if v.get("snapshot_date") != AACT_SNAPSHOT_DATE:
        raise AactNotVerified(f"{verified}: snapshot_date {v.get('snapshot_date')!r}, not {AACT_SNAPSHOT_DATE}")
    pinned = v.get("files", {}).get(STUDIES_FILE, {}).get("sha256")
    if not isinstance(pinned, str) or len(pinned) != 64:
        raise AactNotVerified(f"{verified} pins no sha256 for {STUDIES_FILE}")
    got = sha256_file(studies)
    if got != pinned:
        raise AactNotVerified(f"{studies}: sha256 {got} differs from {pinned} in {verified}")
    return studies, got


def split_record(line: str, n_fields: int) -> list[str]:
    """One flat-file line into fields. A plain split is used when it gives the header's field
    count (embedded quotes are then literal text); otherwise the line is read with `"` quoting,
    for a field enclosed in quotes because it holds `|`. Anything else raises."""
    fields = line.split("|")
    if len(fields) == n_fields:
        return fields
    quoted = next(csv.reader(io.StringIO(line), delimiter="|", quotechar='"', strict=True))
    if len(quoted) != n_fields:
        raise AactFormatError(f"row has {len(fields)} fields ({len(quoted)} with quoting), header has {n_fields}: "
                              f"{line[:80]!r}")
    return quoted


def read_studies(path: Path, nct_ids: set[str]) -> dict[str, dict[str, str]]:
    """The raw `studies` rows of `nct_ids` (every column), streamed; other rows are skipped.

    Rows are selected on the text before the first `|`, read without quote handling, so `nct_id`
    must be the first column: a quoted field holding `|` in an earlier column would shift a plain
    split. Any other position raises."""
    out: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as f:
        header = f.readline().rstrip("\r\n").split("|")
        missing = [c for c in STUDIES_COLUMNS if c not in header]
        if missing:
            raise AactFormatError(f"{path}: header lacks {missing}")
        if header[0] != "nct_id":
            raise AactFormatError(f"{path}: nct_id is column {header.index('nct_id')}, not column 0")
        for line in f:
            line = line.rstrip("\r\n")
            if not line:
                continue
            nct = line.split("|", 1)[0].strip('"')
            if nct not in nct_ids:
                continue
            row = dict(zip(header, split_record(line, len(header))))
            if row["nct_id"] in out:
                raise AactFormatError(f"{path}: {row['nct_id']} appears twice")
            out[row["nct_id"]] = row
    return out


def _day(year: int, month: int | None, day: int | None, last_day: bool) -> date:
    """The date of a (possibly partial) year / month / day: the first day of the stated month or
    year, or with `last_day` the last."""
    if month is None:
        return date(year, 12, 31) if last_day else date(year, 1, 1)
    if day is None:
        return date(year, month, calendar.monthrange(year, month)[1] if last_day else 1)
    return date(year, month, day)


def parse_partial_date(raw: str, last_day: bool = False) -> date | None:
    """ClinicalTrials.gov dates as given: 'YYYY', 'YYYY-MM', 'YYYY-MM-DD' (API v2) or 'Month YYYY',
    'Month, YYYY', 'Month D, YYYY' (legacy). IMPLEMENTATION CHOICE: a partial date takes the first
    day of its month (or year). `last_day=True` takes the last day instead (diagnostic only)."""
    raw = (raw or "").strip()
    if not raw:
        return None
    try:
        m = ISO_PARTIAL.fullmatch(raw)
        if m:
            y, mo, d = m.groups()
            return _day(int(y), int(mo) if mo else None, int(d) if d else None, last_day)
        m = LEGACY_DATE.fullmatch(raw)
        if m and m.group(1).lower() in MONTHS:
            return _day(int(m.group(3)), MONTHS[m.group(1).lower()], int(m.group(2)) if m.group(2) else None, last_day)
    except ValueError as err:
        raise AactFormatError(f"impossible ClinicalTrials.gov date {raw!r}") from err
    raise AactFormatError(f"unparseable ClinicalTrials.gov date {raw!r}")


def parse_full_date(raw: str) -> date | None:
    raw = (raw or "").strip()
    if not raw:
        return None
    try:
        return datetime.strptime(raw, "%Y-%m-%d").date()
    except ValueError as err:
        raise AactFormatError(f"not a YYYY-MM-DD date: {raw!r}") from err


def raw_date(row: dict[str, str], stem: str, last_day: bool = False) -> date | None:
    """`<stem>_month_year` as given; `<stem>_date` is never converted (module docstring)."""
    given, converted = row.get(f"{stem}_month_year", ""), row.get(f"{stem}_date", "")
    if not given.strip() and converted.strip():
        raise AactFormatError(f"{row.get('nct_id')}: {stem}_month_year empty but {stem}_date = {converted!r}")
    return parse_partial_date(given, last_day)


def parse_status(raw: str) -> str:
    value = (raw or "").strip()
    if value not in STATUS_OF:
        raise AactFormatError(f"overall_status {raw!r} is not a ClinicalTrials.gov status")
    return STATUS_OF[value]


def parse_phases(raw: str) -> tuple[str, ...]:
    value = (raw or "").strip()
    if value in NO_PHASE:
        return ()
    tokens = [t.strip() for t in value.split("/")]
    if not all(t in PHASE_OF for t in tokens):
        raise AactFormatError(f"phase {raw!r} is not a ClinicalTrials.gov phase list")
    return tuple(PHASE_OF[t] for t in tokens)


def parse_study(row: dict[str, str], last_day: bool = False) -> Trial:
    """One AACT `studies` row -> Trial. `last_day` reads partial start and completion dates as the
    last day of the stated month or year (diagnostic only; the registered reading is the first)."""
    if not NCT.fullmatch(row.get("nct_id", "")):
        raise AactFormatError(f"not an NCT id: {row.get('nct_id')!r}")
    return Trial(
        nct_id=row["nct_id"],
        phases=parse_phases(row["phase"]),
        overall_status=parse_status(row["overall_status"]),
        why_stopped=row.get("why_stopped", ""),
        start_date=raw_date(row, "start", last_day),
        completion_date=raw_date(row, "completion", last_day),
        first_submit_date=parse_full_date(row["study_first_submitted_date"]),
        last_update_submit_date=parse_full_date(row["last_update_submitted_date"]),
    )


class AactTrials:
    """The linked trials of stage C, read from the verified archive. `load` writes the raw rows of
    every linked NCT id to `<out>/studies_linked.tsv`, the ids absent from the archive to
    `<out>/not_found.jsonl`, and `<out>/manifest.json` (source sha256, counts) before returning.
    `raw` then holds those rows by NCT id, for the diagnostics that re-read their dates."""

    def __init__(self, aact_dir: Path, out: Path):
        self.aact_dir = aact_dir
        self.out = out
        self.rows = out / "studies_linked.tsv"
        self.missing = out / "not_found.jsonl"
        self.manifest = out / "manifest.json"
        self.raw: dict[str, dict[str, str]] = {}

    def load(self, nct_ids: list[str]) -> dict[str, Trial | None]:
        studies, sha = verified_studies(self.aact_dir)
        wanted = set(nct_ids)
        rows = read_studies(studies, wanted)
        self.raw = rows
        self.out.mkdir(parents=True, exist_ok=True)
        with self.rows.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(STUDIES_COLUMNS)
            for n in sorted(rows):
                w.writerow([rows[n][c] for c in STUDIES_COLUMNS])
        absent = sorted(wanted - set(rows))
        self.missing.write_text("".join(json.dumps({"nct_id": n, "source": STUDIES_FILE}) + "\n" for n in absent))
        self.manifest.write_text(json.dumps({"source": str(studies), "source_sha256": sha,
                                             "snapshot_date": AACT_SNAPSHOT_DATE, "requested": len(wanted),
                                             "found": len(rows), "not_found": len(absent)}, indent=1))
        trials: dict[str, Trial | None] = {n: None for n in absent}
        trials.update({n: parse_study(r) for n, r in rows.items()})
        return trials
