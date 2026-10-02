import hashlib
import io
import json
import zipfile
from datetime import date

import pytest

from aact_download import HEADERS, UPLOADED_SOURCE, ArchiveHashMismatch, DownloadError, download_archive, verify_archive
from ctgov import (AactFormatError, AactNotVerified, AactTrials, parse_partial_date, parse_phases, parse_status,
                   parse_study, read_studies, split_record, verified_studies)
from tests.conftest import AACT_HEADER, aact_row, write_aact


# ---- values -------------------------------------------------------------------------------------

def test_parse_study_reads_the_aact_columns():
    t = parse_study(aact_row("NCT00000001", "TERMINATED", "PHASE1/PHASE2", "2015-03", "2017-06-12", "futility",
                             submitted="2014-12-01", updated="2026-10-01"))
    assert (t.nct_id, t.phases, t.overall_status, t.why_stopped) == ("NCT00000001", ("PHASE1", "PHASE2"),
                                                                     "TERMINATED", "futility")
    assert (t.start_date, t.completion_date) == (date(2015, 3, 1), date(2017, 6, 12))
    assert (t.first_submit_date, t.last_update_submit_date) == (date(2014, 12, 1), date(2026, 10, 1))


def test_legacy_and_api_v2_vocabularies_give_the_same_trial():
    enum = parse_study(aact_row("NCT00000001", "ACTIVE_NOT_RECRUITING", "PHASE2/PHASE3", "2015-03", "2027-01"))
    legacy = parse_study(aact_row("NCT00000001", "Active, not recruiting", "Phase 2/Phase 3", "March 2015",
                                  "January, 2027"))
    assert enum == legacy


@pytest.mark.parametrize("raw,expected", [
    ("2024", date(2024, 1, 1)), ("2024-09", date(2024, 9, 1)), ("2024-09-17", date(2024, 9, 17)),
    ("September 2024", date(2024, 9, 1)), ("September, 2024", date(2024, 9, 1)),
    ("September 17, 2024", date(2024, 9, 17)), ("", None)])
def test_partial_dates_take_the_first_day(raw, expected):
    assert parse_partial_date(raw) == expected


@pytest.mark.parametrize("raw", ["Sept 2024", "2024/09", "Q3 2024", "2024-13"])
def test_unparseable_dates_raise(raw):
    with pytest.raises(AactFormatError):
        parse_partial_date(raw)


@pytest.mark.parametrize("raw,expected", [
    ("PHASE2", ("PHASE2",)), ("PHASE1/PHASE2", ("PHASE1", "PHASE2")), ("Phase 1/Phase 2", ("PHASE1", "PHASE2")),
    ("EARLY_PHASE1", ("EARLY_PHASE1",)), ("Early Phase 1", ("EARLY_PHASE1",)), ("NA", ()), ("Not Applicable", ()),
    ("", ())])
def test_phase_lists(raw, expected):
    assert parse_phases(raw) == expected


@pytest.mark.parametrize("raw", ["Phase 5", "PHASE2 PHASE3", "phase2"])
def test_unknown_phase_raises(raw):
    with pytest.raises(AactFormatError):
        parse_phases(raw)


@pytest.mark.parametrize("raw", ["Stopped", "", "completed"])
def test_unknown_status_raises(raw):
    with pytest.raises(AactFormatError):
        parse_status(raw)


def test_converted_date_without_the_raw_value_raises():
    row = aact_row("NCT00000001", "COMPLETED", "PHASE2", "2015-03", "")
    row["completion_date"] = "2017-06-30"
    with pytest.raises(AactFormatError, match="completion_month_year empty"):
        parse_study(row)


# ---- the flat file --------------------------------------------------------------------------------

def test_split_record_reads_a_quoted_field_holding_a_pipe_and_keeps_literal_quotes():
    assert split_record('NCT1|"a|b"|c', 3) == ["NCT1", "a|b", "c"]
    assert split_record('NCT1|the "x" trial|c', 3) == ["NCT1", 'the "x" trial', "c"]
    with pytest.raises(AactFormatError):
        split_record("NCT1|a|b|c", 3)


def test_read_studies_returns_only_the_requested_rows(tmp_path):
    d = write_aact(tmp_path / "a", [aact_row("NCT00000001", "COMPLETED", "PHASE2"),
                                    aact_row("NCT00000002", "COMPLETED", "PHASE3", title='x "y"'),
                                    aact_row("NCT00000003", "WITHDRAWN", "PHASE2", title="a|b")])
    text = (d / "studies.txt").read_text().replace("|a|b|", '|"a|b"|')
    (d / "studies.txt").write_text(text)
    rows = read_studies(d / "studies.txt", {"NCT00000002", "NCT00000003", "NCT00000009"})
    assert set(rows) == {"NCT00000002", "NCT00000003"}
    assert rows["NCT00000002"]["brief_title"] == 'x "y"' and rows["NCT00000003"]["brief_title"] == "a|b"
    assert list(rows["NCT00000003"]) == list(AACT_HEADER)


def test_read_studies_refuses_a_duplicate_or_a_missing_column(tmp_path):
    d = write_aact(tmp_path / "a", [aact_row("NCT00000001", "COMPLETED", "PHASE2")] * 2)
    with pytest.raises(AactFormatError, match="twice"):
        read_studies(d / "studies.txt", {"NCT00000001"})
    (d / "studies.txt").write_text("nct_id|phase\nNCT00000001|PHASE2\n")
    with pytest.raises(AactFormatError, match="header lacks"):
        read_studies(d / "studies.txt", {"NCT00000001"})


def test_read_studies_refuses_nct_id_in_any_column_but_the_first(tmp_path):
    row = aact_row("NCT00000001", "COMPLETED", "PHASE2")
    first = tmp_path / "first.txt"
    first.write_text("|".join(AACT_HEADER) + "\n" + "|".join(row[c] for c in AACT_HEADER) + "\n")
    assert set(read_studies(first, {"NCT00000001"})) == {"NCT00000001"}
    # the same row with brief_title moved in front of nct_id and holding a quoted pipe: a plain
    # split at column 1 would read `b"` as the id and skip the row
    header = ("brief_title", "nct_id", *AACT_HEADER[2:])
    second = tmp_path / "second.txt"
    second.write_text("|".join(header) + "\n" + "|".join('"a|b"' if c == "brief_title" else row[c] for c in header) + "\n")
    with pytest.raises(AactFormatError, match="nct_id is column 1, not column 0"):
        read_studies(second, {"NCT00000001"})


# ---- VERIFIED.json ---------------------------------------------------------------------------------

def test_verified_studies_accepts_the_pinned_file(tmp_path):
    d = write_aact(tmp_path / "a", [aact_row("NCT00000001", "COMPLETED", "PHASE2")])
    path, sha = verified_studies(d)
    assert sha == hashlib.sha256(path.read_bytes()).hexdigest()


def test_verified_studies_refuses_without_verified_json(tmp_path):
    with pytest.raises(AactNotVerified, match="missing"):
        verified_studies(write_aact(tmp_path / "a", [], pin=False))


def test_verified_studies_refuses_an_edited_file(tmp_path):
    d = write_aact(tmp_path / "a", [aact_row("NCT00000001", "TERMINATED", "PHASE2", why="business")])
    (d / "studies.txt").write_text((d / "studies.txt").read_text().replace("business", "futility"))
    with pytest.raises(AactNotVerified, match="differs"):
        verified_studies(d)


def test_verified_studies_refuses_another_snapshot_date(tmp_path):
    with pytest.raises(AactNotVerified, match="snapshot_date"):
        verified_studies(write_aact(tmp_path / "a", [], snapshot_date="2026-11-01"))


def test_aact_trials_writes_the_rows_it_used(tmp_path):
    d = write_aact(tmp_path / "a", [aact_row("NCT00000001", "COMPLETED", "PHASE2"),
                                    aact_row("NCT00000002", "RECRUITING", "PHASE3")])
    src = AactTrials(d, tmp_path / "snap")
    trials = src.load(["NCT00000002", "NCT00000005"])
    assert trials["NCT00000005"] is None and trials["NCT00000002"].overall_status == "RECRUITING"
    assert src.rows.read_text().splitlines()[1].startswith("NCT00000002\t")
    assert len(src.rows.read_text().splitlines()) == 2
    assert [json.loads(x)["nct_id"] for x in src.missing.read_text().splitlines()] == ["NCT00000005"]
    m = json.loads(src.manifest.read_text())
    assert (m["requested"], m["found"], m["not_found"], m["snapshot_date"]) == (2, 1, 1, "2026-09-30")


# ---- download (synthetic zip, fake session) ---------------------------------------------------------

def zip_bytes(studies: str) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_STORED) as z:
        z.writestr("studies.txt", studies)
        z.writestr("interventions.txt", "id|nct_id\n" + "".join(f"{i}|NCT{i:08d}\n" for i in range(2000)))
    return buf.getvalue()


class FakeResponse:
    def __init__(self, status: int, body: bytes, fail_after: int | None = None):
        self.status_code, self.body, self.fail_after = status, body, fail_after

    def iter_content(self, n):
        for i in range(0, len(self.body), 1000):
            if self.fail_after is not None and i >= self.fail_after:
                raise ConnectionError("dropped")
            yield self.body[i:i + 1000]

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class FakeSession:
    def __init__(self, body: bytes, honour_range: bool = True, fail_after: int | None = None):
        self.body, self.honour_range, self.fail_after, self.requests = body, honour_range, fail_after, []

    def get(self, url, headers, stream, timeout):
        self.requests.append(headers)
        start = int(headers["Range"].split("=")[1].rstrip("-")) if "Range" in headers else 0
        fail, self.fail_after = self.fail_after, None
        if start and self.honour_range:
            return FakeResponse(206, self.body[start:], fail)
        return FakeResponse(200, self.body, fail)


STUDIES = "|".join(AACT_HEADER) + "\n" + "|".join(aact_row("NCT00000001", "COMPLETED", "PHASE2")[c]
                                                    for c in AACT_HEADER) + "\n"
URL = "https://example.org/static/exported_files/daily/20260930_export_ctgov.zip?X-Signature=abc"


def test_download_writes_a_verified_json_stage_c_accepts(tmp_path):
    body = zip_bytes(STUDIES)
    commits = []
    rec = download_archive(URL, tmp_path / "aact_20260930", FakeSession(body), commit=lambda: commits.append(1))
    assert rec["archive"]["sha256"] == hashlib.sha256(body).hexdigest()
    assert rec["source_url"].endswith("20260930_export_ctgov.zip") and "Signature" not in rec["source_url"]
    assert rec["files"]["studies.txt"]["rows"] == 1 and rec["snapshot_date"] == "2026-09-30"
    path, sha = verified_studies(tmp_path / "aact_20260930")
    assert path.read_text() == STUDIES and commits


def test_download_resumes_from_the_bytes_on_disk_with_a_range_request(tmp_path):
    body = zip_bytes(STUDIES)
    dest = tmp_path / "aact_20260930"
    s1 = FakeSession(body, fail_after=5000)
    commits = []
    with pytest.raises(ConnectionError):
        download_archive(URL, dest, s1, commit=lambda: commits.append((dest / "20260930_export_ctgov.zip.part")
                                                                     .stat().st_size), checkpoint_bytes=2000)
    assert commits and all(c % 1000 == 0 for c in commits)   # committed inside the unit, before the failure
    partial = (dest / "20260930_export_ctgov.zip.part").stat().st_size
    assert 0 < partial < len(body)
    s2 = FakeSession(body)
    rec = download_archive(URL, dest, s2)
    assert s2.requests[0]["Range"] == f"bytes={partial}-"
    assert rec["archive"]["sha256"] == hashlib.sha256(body).hexdigest()


def test_download_restarts_when_the_server_ignores_range(tmp_path):
    body = zip_bytes(STUDIES)
    dest = tmp_path / "aact_20260930"
    dest.mkdir()
    (dest / "20260930_export_ctgov.zip.part").write_bytes(b"garbage")
    rec = download_archive(URL, dest, FakeSession(body, honour_range=False))
    assert rec["archive"]["sha256"] == hashlib.sha256(body).hexdigest()


def test_download_skips_everything_once_verified(tmp_path):
    body = zip_bytes(STUDIES)
    download_archive(URL, tmp_path / "a", FakeSession(body))
    again = FakeSession(body)
    download_archive(URL, tmp_path / "a", again)
    assert again.requests == []


def test_download_refuses_a_non_zip_and_a_zip_without_the_columns(tmp_path):
    with pytest.raises(DownloadError):
        download_archive(URL, tmp_path / "a", FakeSession(b"<html>sign in</html>"))
    assert not list((tmp_path / "a").glob("*.zip"))
    with pytest.raises(AactFormatError, match="header lacks"):
        download_archive(URL, tmp_path / "b", FakeSession(zip_bytes("nct_id|phase\nNCT00000001|PHASE2\n")))
    assert not (tmp_path / "a" / "VERIFIED.json").exists() and not (tmp_path / "b" / "VERIFIED.json").exists()


def test_download_request_headers_carry_no_email_or_personal_field(tmp_path):
    s = FakeSession(zip_bytes(STUDIES))
    download_archive(URL, tmp_path / "a", s)
    for h in s.requests + [HEADERS]:
        assert set(h) <= {"User-Agent", "Range"}
        assert "@" not in " ".join(h.values())


def test_download_with_a_wrong_expected_sha_raises_and_writes_no_verified_json(tmp_path):
    with pytest.raises(DownloadError, match="expected"):
        download_archive(URL, tmp_path / "a", FakeSession(zip_bytes(STUDIES)), expected_sha256="0" * 64)
    assert not (tmp_path / "a" / "VERIFIED.json").exists() and not (tmp_path / "a" / "studies.txt").exists()


# ---- uploaded archive -------------------------------------------------------------------------------

def test_verify_archive_pins_an_uploaded_zip_stage_c_accepts(tmp_path):
    body = zip_bytes(STUDIES)
    dest = tmp_path / "aact_20260930"
    dest.mkdir()
    (dest / "20260930_export_ctgov.zip").write_bytes(body)
    rec = verify_archive(dest / "20260930_export_ctgov.zip", dest, UPLOADED_SOURCE, hashlib.sha256(body).hexdigest())
    assert (rec["source"], rec["archive"]["sha256"]) == (UPLOADED_SOURCE, hashlib.sha256(body).hexdigest())
    assert "source_url" not in rec and rec["files"]["studies.txt"]["rows"] == 1
    path, sha = verified_studies(dest)
    assert path.read_text() == STUDIES and sha == rec["files"]["studies.txt"]["sha256"]


def test_verify_archive_refuses_a_zip_whose_sha_differs_and_keeps_it(tmp_path):
    body = zip_bytes(STUDIES)
    dest = tmp_path / "aact_20260930"
    dest.mkdir()
    (dest / "x.zip").write_bytes(body)
    wrong = hashlib.sha256(body + b"\0").hexdigest()
    with pytest.raises(ArchiveHashMismatch, match=wrong):
        verify_archive(dest / "x.zip", dest, UPLOADED_SOURCE, wrong)
    assert (dest / "x.zip").read_bytes() == body
    assert not (dest / "VERIFIED.json").exists() and not (dest / "studies.txt").exists()
