"""The AACT daily flat-file snapshot of 2026-09-30 on `pqtl-v8-inputs`, and its VERIFIED.json.

The snapshot (AACT "Latest" daily export, every member dated 2026-09-30 about 01:00) was
downloaded by the author through a signed-in AACT session and uploaded to the volume at
/aact_20260930/. AACT deletes daily snapshots at month end, so the copy on the volume is the
record. Two routes reach the same check, `verify_archive`:

- uploaded (`modal_stage_c.py::verify_uploaded_aact`): the zip is already on the volume; its
  sha256 must equal the one recorded when it was downloaded.
- URL (`modal_stage_c.py::download_aact`): the zip is streamed to `<dest>/<name>.part`. Every
  CHECKPOINT_BYTES the part file is flushed and `commit()` (the volume commit) runs; a restarted
  call resumes with an HTTP Range request from the bytes already on disk. A server that ignores
  Range restarts the file. Without an expected sha256 the sha256 of the fetched file is pinned.

`verify_archive`, resumable: if `<dest>/VERIFIED.json` exists nothing is done. The zip's sha256
must equal the expected one; its members are checked with `ZipFile.testzip`; `studies.txt` is
extracted beside it and its header must hold every column stage C reads (ctgov.STUDIES_COLUMNS).
VERIFIED.json records the snapshot date, the source as a sentence (and, for the URL route, the URL
without its query string; a signed link is not written down), the zip's and studies.txt's sha256
and sizes, the studies.txt member's zip timestamp, header and row count, and the UTC time. Stage C
reads studies.txt only when its sha256 matches. No row content is read or written anywhere else.

Requests carry a fixed User-Agent and no email or other personal header.
"""
import json
import shutil
import zipfile
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit

import requests

from ctgov import AACT_SNAPSHOT_DATE, STUDIES_COLUMNS, STUDIES_FILE, VERIFIED_NAME, AactFormatError, sha256_file
from models import StageCError

HEADERS = {"User-Agent": "pqtl-mr-v8-stage-c/1.0"}
CHECKPOINT_BYTES = 256 * 1024 * 1024
CHUNK_BYTES = 1024 * 1024
UPLOADED_SOURCE = ("AACT daily flat-file snapshot 2026-09-30, downloaded by the author from "
                   "aact.ctti-clinicaltrials.org (signed-in session) on 2026-09-30, uploaded to pqtl-v8-inputs")


class DownloadError(StageCError):
    """The archive could not be retrieved or is not a readable AACT flat-file zip."""


class ArchiveHashMismatch(DownloadError):
    """The zip on disk is not the file whose sha256 was recorded."""


def fetch_zip(url: str, part: Path, session: requests.Session, commit: Callable[[], None],
              checkpoint_bytes: int = CHECKPOINT_BYTES) -> None:
    have = part.stat().st_size if part.exists() else 0
    headers = dict(HEADERS, **({"Range": f"bytes={have}-"} if have else {}))
    with session.get(url, headers=headers, stream=True, timeout=120) as r:
        if r.status_code == 200 and have:
            have = 0                       # Range ignored: start over
        elif r.status_code not in (200, 206):
            raise DownloadError(f"{urlsplit(url).netloc}{urlsplit(url).path}: HTTP {r.status_code}")
        since = 0
        with part.open("ab" if have else "wb") as f:
            for chunk in r.iter_content(CHUNK_BYTES):
                f.write(chunk)
                since += len(chunk)
                if since >= checkpoint_bytes:
                    f.flush()
                    commit()
                    since = 0
    commit()


def studies_rows(path: Path) -> tuple[list[str], int]:
    with path.open(encoding="utf-8") as f:
        header = f.readline().rstrip("\r\n").split("|")
        n = sum(1 for line in f if line.strip())
    missing = [c for c in STUDIES_COLUMNS if c not in header]
    if missing:
        raise AactFormatError(f"{path}: header lacks {missing}")
    return header, n


def verify_archive(final: Path, dest: Path, source: str, expected_sha256: str,
                   commit: Callable[[], None] = lambda: None, extra: dict | None = None) -> dict:
    """Returns the VERIFIED.json content (module docstring). Never deletes `final`."""
    verified = dest / VERIFIED_NAME
    if verified.exists():
        return json.loads(verified.read_text())
    if not final.is_file():
        raise DownloadError(f"{final} is missing")
    archive_sha = sha256_file(final)
    if archive_sha != expected_sha256.lower():
        raise ArchiveHashMismatch(f"{final}: sha256 {archive_sha}, expected {expected_sha256}")
    studies, part = dest / STUDIES_FILE, dest / f"{STUDIES_FILE}.part"
    try:
        with zipfile.ZipFile(final) as z:
            bad = z.testzip()
            if bad is not None:
                raise DownloadError(f"{final}: member {bad} fails its CRC")
            member = next((m for m in z.namelist() if Path(m).name == STUDIES_FILE), None)
            if member is None:
                raise DownloadError(f"{final} holds no {STUDIES_FILE}")
            member_time = datetime(*z.getinfo(member).date_time).isoformat()
            with z.open(member) as src, part.open("wb") as out:
                shutil.copyfileobj(src, out)
    except zipfile.BadZipFile as err:
        raise DownloadError(f"{final} is not a readable zip: {err}") from err
    part.replace(studies)
    commit()
    header, n = studies_rows(studies)
    record = {
        "snapshot_date": AACT_SNAPSHOT_DATE,
        "source": source,
        **(extra or {}),
        "archive": {"name": final.name, "sha256": archive_sha, "bytes": final.stat().st_size},
        "files": {STUDIES_FILE: {"sha256": sha256_file(studies), "bytes": studies.stat().st_size,
                                 "zip_member_time": member_time, "rows": n, "header": header}},
        "verified_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    tmp = dest / f"{VERIFIED_NAME}.part"
    tmp.write_text(json.dumps(record, indent=1))
    tmp.replace(verified)
    commit()
    return record


def download_archive(url: str, dest: Path, session: requests.Session | None = None,
                     commit: Callable[[], None] = lambda: None, checkpoint_bytes: int = CHECKPOINT_BYTES,
                     expected_sha256: str | None = None) -> dict:
    """URL route: fetch (resumable), then `verify_archive`. Returns the VERIFIED.json content."""
    verified = dest / VERIFIED_NAME
    if verified.exists():
        return json.loads(verified.read_text())
    dest.mkdir(parents=True, exist_ok=True)
    name = Path(urlsplit(url).path).name or f"{AACT_SNAPSHOT_DATE}_pipe-delimited-export.zip"
    part, final = dest / f"{name}.part", dest / name
    if not final.exists():
        fetch_zip(url, part, session or requests.Session(), commit, checkpoint_bytes)
        part.replace(final)
        commit()
    s = urlsplit(url)
    source_url = urlunsplit((s.scheme, s.netloc, s.path, "", ""))
    try:
        return verify_archive(final, dest, f"AACT flat-file archive downloaded from {source_url}",
                              expected_sha256 or sha256_file(final), commit,
                              extra={"source_url": source_url, "sha256_pinned_before_download": bool(expected_sha256)})
    except DownloadError as err:
        final.unlink()                     # a rerun downloads again rather than reusing a bad file
        commit()
        raise DownloadError(f"{final} is not a readable AACT archive: {err}") from err
