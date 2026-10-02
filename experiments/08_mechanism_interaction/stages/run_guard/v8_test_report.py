"""Run a stage's synthetic pytest suite inside its Modal image and record the summary.

Each stage's Modal wrapper has a `tests` function that calls `run_pytest` in the stage directory
baked into the image and `write_report` on the stage's volume. A report already at the target path
is moved to `superseded/` beside it, named by the time it was written, before the new one is
written; no report is overwritten. Synthetic fixtures only: the suites read no study data.
"""
import json
import os
import re
import subprocess
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path

COUNT_RE = re.compile(r"(\d+) (passed|failed|skipped|errors?|xfailed|xpassed|deselected|warnings?)")


def run_pytest(cwd: Path, pythonpath: Sequence[Path], args: Sequence[str] = ("tests",),
               command: Sequence[str] = ("python", "-m", "pytest"), timeout_s: int = 3 * 3600) -> dict:
    """pytest in `cwd` with PYTHONPATH set to `pythonpath`; the summary line, its counts and the
    output tails. A non-zero exit code is reported, not raised."""
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(str(p) for p in pythonpath)}
    res = subprocess.run([*command, "-q", "-p", "no:cacheprovider", *args], cwd=cwd, env=env, capture_output=True,
                         text=True, timeout=timeout_s)
    lines = [line for line in res.stdout.splitlines() if line.strip()]
    summary = lines[-1] if lines else ""
    counts = {kind.rstrip("s") if kind.startswith(("error", "warning")) else kind: int(n)
              for n, kind in COUNT_RE.findall(summary)}
    return {"synthetic": True, "cwd": str(cwd), "args": list(args), "returncode": res.returncode, "summary": summary,
            "counts": counts, "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "stdout_tail": res.stdout[-6000:], "stderr_tail": res.stderr[-2000:]}


def write_report(path: Path, report: dict) -> Path:
    """Write `report` to `path`, first moving any report already there to superseded/."""
    if path.exists():
        written = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        old = path.parent / "superseded" / f"{path.stem}_{written}{path.suffix}"
        old.parent.mkdir(parents=True, exist_ok=True)
        if old.exists():
            raise FileExistsError(f"{old} already exists; not overwriting an archived report")
        path.rename(old)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2))
    return path
