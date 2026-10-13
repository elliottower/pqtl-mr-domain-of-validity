"""Shared constants, the run guard and small helpers for the V9 follow-up (PREREG.md, frozen
2026-10-13, plan sha256 ed1285b8…).

The guard: the plan must be byte-identical to its freeze (`.prereg/PREREG.md.json`), and its
`PREREG.log` must hold `RUN_START followup token=<token>` for the token the caller passes. V8's
sealed stage outputs are checked separately by `stage_d.guard.logged_seals` against the V8 log.
"""
import hashlib
import json
import re
from math import sqrt
from pathlib import Path

SEED = 20261013
FROZEN_PLAN_SHA256 = "ed1285b842a9a8d0fe5e72d78dbed7a5ae420f756a37b4b9dd1042151513f1bd"
V8_STAGE_D_TOKEN = "a23dc31477e8f272"   # V8's RUN_START stage=D token; its log carries the A-C seals
RUN_START_RE = re.compile(r"RUN_START followup token=([0-9a-f]{16})\b")
DECLARATION = ("These analyses were selected after inspection of the frozen V8 stage D results. They are post hoc "
               "and exploratory, not independent tests of V8 H1 or H4. Their estimates and intervals do not alter "
               "the registered no-decision status of any V8 hypothesis.")
BOOT_REPS = 10_000
BOOT_FAIL_MAX = 0.05


class FollowupError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def require_run(plan: Path, log: Path, token: str) -> None:
    """Refuse unless `plan` is the frozen plan and `log` holds the RUN_START line for `token`."""
    got = sha256_file(plan)
    if got != FROZEN_PLAN_SHA256:
        raise FollowupError(f"{plan} has sha256 {got}, the frozen plan is {FROZEN_PLAN_SHA256}")
    tokens = {m.group(1) for m in RUN_START_RE.finditer(log.read_text(encoding="utf-8"))}
    if token not in tokens:
        raise FollowupError(f"{log} has no `RUN_START followup token={token}` entry")


def wilson(successes: int, n: int, z: float = 1.959963984540054) -> tuple[float, float] | None:
    """Wilson score 95% interval for a binomial proportion; None when n = 0."""
    if n == 0:
        return None
    p = successes / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def write_json(path: Path, obj: dict) -> None:
    """Atomic write; every output carries the declaration."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps({"declaration": DECLARATION, **obj}, indent=1, default=_default))
    tmp.replace(path)


def _default(x):
    if hasattr(x, "item"):
        return x.item()
    if hasattr(x, "isoformat"):
        return x.isoformat()
    raise TypeError(f"not JSON serializable: {type(x).__name__}")
