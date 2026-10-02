"""Run guard shared by stages A, B, C and D (PREREG §Additional blinding; INTERFACES.md §Rules).

A real run of stage X is allowed only when all of these hold, checked in this order:

1. `prereg check` semantics on PREREG.md: the frozen plan hash is 6cbba084...0040, the plan above
   the log hashes to it (prereg's own `plan_of` / `sha256_of`), nothing sits outside the hash
   (`unhashed_content`), and the log chain is intact (`log_problems`).
2. The log below `## Log` holds exactly one entry whose event text is exactly
   `RUN_START stage=X token=<token>`, written by `prereg log` (it carries a chain value), and the
   token equals the one the runner was given. Log lines are parsed one by one into date, event
   and access level; the log is never substring-searched, so a mention such as
   "Pre-stage note (B, ...)" or "stage B committed" matches nothing.
3. Stage A only: OSF registration 9tzfk is approved (`pending_registration_approval` false and
   `public` true on https://api.osf.io/v2/registrations/9tzfk/).
4. The stages run in the plan's order A -> B -> C -> D, each sealed before the next starts. Before
   the stage's RUN_START entry the log holds a chained `SEAL stage=P manifest_sha256=<h>` for every
   predecessor P (PREDECESSORS: B needs A; C needs A and B; D needs A, B and C). The last such
   entry per stage before the RUN_START counts; a SEAL logged after it counts for nothing.
5. Each predecessor MANIFEST.tsv the stage reads (MANIFEST_READS: B and C read A's; D reads A's,
   B's and C's) is given to the guard, and its sha256 equals the logged seal. The stage then checks
   the predecessor files it consumes against that manifest (v8_manifest.verify_listed or
   verify_output_dir), so the seal covers the bytes read, not only the manifest.

Tokens are 8 or more characters of letters, digits, `.`, `_` or `-`. Log a start with, e.g.:
    prereg log "RUN_START stage=A token=$(openssl rand -hex 8)" --access "nothing run"

Commit identity. A real run is identified by one repository commit:
- `clean_commit(repo)` raises unless `git status --porcelain` is empty for the stage code
  (experiments/08_mechanism_interaction/stages) and PREREG.md, and returns HEAD. Each launcher calls
  it before it starts anything.
- `baked_commit(repo)` is the same value for an image build, or "" when the tree is not clean; each
  Modal wrapper writes it to the image environment as V8_REPO_COMMIT.
- `run_commit(repo)` is what a runner records in run_info.json and INPUTS.tsv: V8_REPO_COMMIT where
  the environment sets it (so a container built from an unclean tree refuses), else
  `clean_commit(repo)`.
"""
import hashlib
import json
import os
import re
import subprocess
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from urllib.request import Request, urlopen

from prereg.log import ACCESS, LOG_MARK, log_lines, log_problems
from prereg.plan import plan_of, sha256_of, unhashed_content

PLAN_SHA256 = "6cbba0840b213f586ef455b33e569a84a94fc676103b74adef6c46e2e0e70040"
OSF_REGISTRATION = "9tzfk"
OSF_API = "https://api.osf.io/v2/registrations/{rid}/"
STAGES = ("A", "B", "C", "D")
PREDECESSORS = {"A": (), "B": ("A",), "C": ("A", "B"), "D": ("A", "B", "C")}
MANIFEST_READS = {"A": (), "B": ("A",), "C": ("A",), "D": ("A", "B", "C")}
TOKEN_RE = r"[A-Za-z0-9._-]{8,}"
RUN_START_RE = re.compile(rf"RUN_START stage=([ABCD]) token=({TOKEN_RE})")
SEAL_RE = re.compile(r"SEAL stage=([ABC]) manifest_sha256=([0-9a-f]{64})")
FROZEN_HASH_RE = re.compile(r"^\*\*Plan sha256:\*\* `([0-9a-f]{64})`", re.M)
ENTRY_RE = re.compile(r"(\d{4}-\d{2}-\d{2})  (.*?)\s{2,}(" + "|".join(re.escape(a) for a in ACCESS) + r")")
HEADERS = {"Accept": "application/vnd.api+json", "User-Agent": "pqtl-mr-v8-run-guard/1.0"}
COMMIT_ENV = "V8_REPO_COMMIT"
COMMIT_RE = re.compile(r"[0-9a-f]{40}")
CLEAN_PATHS = ("experiments/08_mechanism_interaction/stages", "experiments/08_mechanism_interaction/PREREG.md")


class RunNotAuthorized(RuntimeError):
    """A stage was started without the logged, approved run the plan requires."""


@dataclass(frozen=True)
class LogEntry:
    index: int          # 1-based position in the log
    date: str
    event: str
    access: str
    chained: bool       # written by `prereg log` (carries a chain value)


@dataclass(frozen=True)
class Authorization:
    stage: str
    token: str
    run_start_index: int
    seals: dict[str, str]   # predecessor stage -> logged MANIFEST.tsv sha256 (empty for stage A)


def parse_log(text: str) -> list[LogEntry]:
    """Every log entry below the line, parsed into date, event and access level. A line that does
    not have that shape is kept with an empty event, so it can never match a RUN_START or SEAL."""
    out = []
    for i, line in enumerate(log_lines(text), start=1):
        entry, _, recorded = line.rpartition(LOG_MARK)
        chained = bool(entry) and bool(re.fullmatch(r"[0-9a-f]{8}", recorded.strip()))
        body = (entry if entry else line).rstrip()
        m = ENTRY_RE.fullmatch(body)
        if m is None:
            out.append(LogEntry(index=i, date="", event="", access="", chained=chained))
        else:
            out.append(LogEntry(index=i, date=m.group(1), event=m.group(2).strip(), access=m.group(3),
                                chained=chained))
    return out


def check_plan(text: str, plan_sha256: str = PLAN_SHA256) -> None:
    """`prereg check` semantics: frozen, the recorded hash is the registered one, the plan still
    hashes to it, nothing is outside the hash, and the log chain is intact."""
    m = FROZEN_HASH_RE.search(text)
    if m is None:
        raise RunNotAuthorized("PREREG.md is not frozen (no **Plan sha256:** line)")
    if m.group(1) != plan_sha256:
        raise RunNotAuthorized(f"PREREG.md records plan sha256 {m.group(1)}, not the registered {plan_sha256}")
    now = sha256_of(plan_of(text))
    if now != plan_sha256:
        raise RunNotAuthorized(f"the plan above the log hashes to {now}, not the frozen {plan_sha256}")
    problems = unhashed_content(text) + log_problems(text)
    if problems:
        raise RunNotAuthorized("PREREG.md fails prereg check: " + " | ".join(problems))


def find_run_start(entries: list[LogEntry], stage: str, token: str) -> LogEntry:
    """The single RUN_START entry carrying `token`; it must name `stage` and be chained."""
    if stage not in STAGES:
        raise ValueError(f"unknown stage {stage!r}")
    if not re.fullmatch(TOKEN_RE, token or ""):
        raise RunNotAuthorized(f"run token {token!r} is not 8+ characters of [A-Za-z0-9._-]")
    hits = []
    for e in entries:
        m = RUN_START_RE.fullmatch(e.event)
        if m is not None and m.group(2) == token:
            hits.append((e, m.group(1)))
    if not hits:
        raise RunNotAuthorized(f"the PREREG.md log has no entry 'RUN_START stage={stage} token={token}'; "
                               "log the run with `prereg log` first")
    if len(hits) > 1:
        raise RunNotAuthorized(f"token {token} appears in {len(hits)} RUN_START entries; tokens must be unique")
    entry, logged_stage = hits[0]
    if logged_stage != stage:
        raise RunNotAuthorized(f"token {token} was logged for stage {logged_stage}, not stage {stage}")
    if not entry.chained:
        raise RunNotAuthorized(f"log entry {entry.index} is not chained; RUN_START must be written by `prereg log`")
    return entry


def logged_seals(entries: list[LogEntry], before: int) -> dict[str, str]:
    """stage -> manifest sha256 of the last chained SEAL entry for that stage before entry `before`."""
    seals: dict[str, str] = {}
    for e in entries:
        if e.index >= before:
            break
        m = SEAL_RE.fullmatch(e.event)
        if m is not None and e.chained:
            seals[m.group(1)] = m.group(2)
    return seals


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_seals(entries: list[LogEntry], run_start: LogEntry, stage: str,
                manifests: Mapping[str, Path]) -> dict[str, str]:
    """Every predecessor of `stage` is sealed in the log before its RUN_START, and each MANIFEST.tsv
    the stage is given (exactly those of MANIFEST_READS[stage]) has the logged sha256."""
    if set(manifests) != set(MANIFEST_READS[stage]):
        raise ValueError(f"stage {stage} needs the MANIFEST.tsv of {list(MANIFEST_READS[stage])}; "
                         f"got {sorted(manifests)}")
    seals = logged_seals(entries, run_start.index)
    for s in PREDECESSORS[stage]:
        if s not in seals:
            raise RunNotAuthorized(f"no 'SEAL stage={s} manifest_sha256=<sha256>' entry before the stage {stage} "
                                   f"RUN_START (log entry {run_start.index})")
    for s, path in sorted(manifests.items()):
        got = sha256_file(path)
        if got != seals[s]:
            raise RunNotAuthorized(f"stage {s}: {path} has sha256 {got}, the log seals {seals[s]}")
    return {s: seals[s] for s in PREDECESSORS[stage]}


def fetch_osf_registration(registration: str = OSF_REGISTRATION, token: str | None = None) -> dict:
    """The registration's JSON:API document. The token (env OSF_PAT) is optional: a public
    registration answers unauthenticated. No email or personal header is sent."""
    headers = dict(HEADERS)
    if token:
        headers["Authorization"] = f"Bearer {token}"
    try:
        with urlopen(Request(OSF_API.format(rid=registration), headers=headers), timeout=30) as r:
            return json.loads(r.read())
    except OSError as err:  # HTTPError (401 while the registration is not public) and URLError
        raise RunNotAuthorized(f"OSF registration {registration} could not be read: {err}") from err


def check_osf_approved(doc: dict, registration: str = OSF_REGISTRATION) -> None:
    try:
        data = doc["data"]
        attrs = data["attributes"]
    except (KeyError, TypeError) as err:
        raise RunNotAuthorized(f"OSF response for {registration} has no data.attributes") from err
    if data.get("id") != registration:
        raise RunNotAuthorized(f"OSF returned registration {data.get('id')!r}, not {registration}")
    if attrs.get("pending_registration_approval") is not False:
        raise RunNotAuthorized(f"OSF {registration}: pending_registration_approval is "
                               f"{attrs.get('pending_registration_approval')!r}, not false")
    if attrs.get("public") is not True:
        raise RunNotAuthorized(f"OSF {registration}: public is {attrs.get('public')!r}, not true")


def require_run(prereg: Path, stage: str, token: str, *, manifests: Mapping[str, Path] | None = None,
                osf_fetch: Callable[[], dict] | None = None) -> Authorization:
    """Raise RunNotAuthorized unless stage `stage` may run with `token` (module docstring).
    `manifests` maps each stage of MANIFEST_READS[stage] to the MANIFEST.tsv the caller reads.
    `osf_fetch` returns the OSF registration document; the default queries the OSF API."""
    text = prereg.read_text(encoding="utf-8")
    check_plan(text)
    entries = parse_log(text)
    start = find_run_start(entries, stage, token)
    if stage == "A":
        fetch = osf_fetch or (lambda: fetch_osf_registration(OSF_REGISTRATION, os.environ.get("OSF_PAT")))
        check_osf_approved(fetch())
    seals = check_seals(entries, start, stage, manifests or {})
    return Authorization(stage=stage, token=token, run_start_index=start.index, seals=seals)


# ---- commit identity ---------------------------------------------------------------------------

def git_output(repo: Path, *args: str) -> str:
    """stdout of `git -C <repo> <args>`; a missing git, a directory that is no repository or a
    failing command refuses the run."""
    try:
        res = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as err:
        raise RunNotAuthorized(f"git {' '.join(args)} could not run in {repo}: {err}") from err
    if res.returncode != 0:
        raise RunNotAuthorized(f"git {' '.join(args)} failed in {repo}: {res.stderr.strip()[-500:]}")
    return res.stdout


def clean_commit(repo: Path, git: Callable[..., str] = git_output) -> str:
    """HEAD of `repo`, after `git status --porcelain` shows nothing modified, staged or untracked
    under CLEAN_PATHS. `git(repo, *args)` returns the command's stdout."""
    status = git(repo, "status", "--porcelain", "--", *CLEAN_PATHS)
    if status.strip():
        lines = status.strip().splitlines()
        raise RunNotAuthorized(f"{len(lines)} uncommitted path(s) under {list(CLEAN_PATHS)} (first: {lines[0].strip()}); "
                               "commit the stage code and PREREG.md before a real run")
    sha = git(repo, "rev-parse", "HEAD").strip()
    if COMMIT_RE.fullmatch(sha) is None:
        raise RunNotAuthorized(f"git rev-parse HEAD gave {sha!r}, not a 40-character commit sha")
    return sha


def baked_commit(repo: Path, git: Callable[..., str] = git_output) -> str:
    """`clean_commit(repo)` for an image build, or "" when the tree is not clean, so that building
    an image (which the synthetic test suites also do) never fails and `run_commit` refuses in it."""
    try:
        return clean_commit(repo, git)
    except RunNotAuthorized:
        return ""


def run_commit(repo: Path | None, env: Mapping[str, str] | None = None, git: Callable[..., str] = git_output) -> str:
    """The commit a real run records. Where `env` (default: the process environment) sets
    V8_REPO_COMMIT, that value, which must be a commit sha: an image built from an unclean tree
    carries "" and is refused. Otherwise `clean_commit(repo)`; with no `repo` either, refused."""
    env = os.environ if env is None else env
    if COMMIT_ENV in env:
        sha = env[COMMIT_ENV]
        if COMMIT_RE.fullmatch(sha) is None:
            raise RunNotAuthorized(f"{COMMIT_ENV} is {sha!r}: this image was not built from a clean commit of "
                                   f"{list(CLEAN_PATHS)}")
        return sha
    if repo is None:
        raise RunNotAuthorized(f"{COMMIT_ENV} is not set and no repository was given; the run has no commit identity")
    return clean_commit(repo, git)
