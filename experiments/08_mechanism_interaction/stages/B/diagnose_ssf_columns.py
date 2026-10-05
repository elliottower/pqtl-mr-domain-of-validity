"""Diagnostic, not part of the stage B run: why the GWAS-SSF files that failed the pre-analysis
validation had every row rejected as `se_nonpositive_or_nonfinite`.

For each file named, read on Modal from the stage B volume (mounted read-only) after its sha256 is
checked against its collect record, every data row is classified field by field into aggregate
categories only. No value of any row is returned or written. Per column of interest
(standard_error, the effect column, the p column, effect_allele_frequency, ci_lower, ci_upper, n):

    width_mismatch    the row has another number of fields than the header (not classified further)
    missing:<token>   a token in stage_b.parsers.MISSING, named (`NA`, `#NA`, empty, ...)
    zero, negative, positive_finite, infinite, nan_parsed (a token float() reads as NaN that is not in MISSING)
    non_numeric       a token float() refuses: without a digit the token itself is counted (it is not a
                      number); with a digit only its shape, digits replaced by 9 (e.g. `9,999`), is counted
    padded            the field has surrounding whitespace (counted in addition to its category)

plus the joint count of (standard_error, effect, p) coarse categories, and the allele tokens that are
not [ACGT]+ (allele labels, not association values). Checkpoint every CHECKPOINT_ROWS rows to the
volume `pqtl-v8-diagnostics` (the stage B volume is never written); a restarted call resumes there.

    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \
        modal run experiments/08_mechanism_interaction/stages/B/diagnose_ssf_columns.py \
        --report experiments/08_mechanism_interaction/inputs/stage_b/validation/outcome_validation.json

`--report` is the validation report; its `failed` list is diagnosed (or pass `--keys GCST1,GCST2`).
The result is written to inputs/stage_b/diagnostics/ssf_columns_<utc>.json (gitignored) and to
/ssf_columns/report_<utc>.json on the diagnostics volume.
"""
import gzip
import hashlib
import json
import math
import os
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import modal

MISSING = {"", "NA", "#NA", "nan", "NaN", "None", ".", "-"}   # stage_b.parsers.MISSING, copied: no stage_b import here
ALLELE = re.compile(r"[ACGT]+")
DIGIT = re.compile(r"[0-9]")
COLUMNS = ("standard_error", "beta", "odds_ratio", "p_value", "neg_log_10_p_value", "effect_allele_frequency",
           "ci_lower", "ci_upper", "n")
MAX_TOKENS = 50
CHECKPOINT_ROWS = 2_000_000
HERE = Path(__file__).resolve().parent

app = modal.App("pqtl-v8-stage-b-diagnose", image=modal.Image.debian_slim(python_version="3.12"))
stage_vol = modal.Volume.from_name("pqtl-v8-stage-b").with_mount_options(read_only=True)
diag_vol = modal.Volume.from_name("pqtl-v8-diagnostics", create_if_missing=True)
ROOT = Path("/stage/stage_b")
OUT = Path("/diag/ssf_columns")


def field_class(raw: str) -> str:
    """The aggregate category of one field (module docstring); never the value of a number."""
    x = raw.strip()
    if x in MISSING:
        return f"missing:{x}"
    try:
        v = float(x)
    except ValueError:
        return f"non_numeric_shape:{DIGIT.sub('9', x)[:20]}" if DIGIT.search(x) else f"non_numeric_token:{x[:40]}"
    if math.isnan(v):
        return f"nan_parsed:{x[:20]}"
    if math.isinf(v):
        return "infinite"
    return "zero" if v == 0 else "negative" if v < 0 else "positive_finite"


def coarse(category: str) -> str:
    return category.split(":", 1)[0] if category.startswith(("missing", "non_numeric", "nan_parsed")) else category


def _cap(counter: Counter) -> dict:
    """Categories with counts; distinct non-numeric tokens and shapes capped at MAX_TOKENS each."""
    kept, dropped = {}, Counter()
    by_kind = Counter()
    for k, n in counter.most_common():
        kind = k.split(":", 1)[0]
        if kind in ("non_numeric_token", "non_numeric_shape", "allele_token"):
            by_kind[kind] += 1
            if by_kind[kind] > MAX_TOKENS:
                dropped[kind] += n
                continue
        kept[k] = n
    return {**kept, **{f"{kind}:<{MAX_TOKENS}+ others>": n for kind, n in dropped.items()}}


def classify_rows(lines, header: list[str], state: dict, every: int, checkpoint) -> dict:
    """Classify every data row after the `state['rows']` already done; `checkpoint(state)` every `every` rows."""
    idx = {h: i for i, h in enumerate(header)}
    cols = [c for c in COLUMNS if c in idx]
    effect = "beta" if "beta" in idx and (header[4] == "beta" or "odds_ratio" not in idx) else "odds_ratio"
    p = "p_value" if "p_value" in idx else "neg_log_10_p_value" if "neg_log_10_p_value" in idx else None
    per = {c: Counter(state["columns"].get(c, {})) for c in cols}
    padded = Counter(state["padded"])
    joint = Counter({tuple(k.split("|")): n for k, n in state["joint"].items()})
    alleles = Counter(state["alleles"])
    width = len(header)
    rows, mismatch, skip = state["rows"], state["width_mismatch"], state["rows"]

    def snapshot() -> dict:
        return {"rows": rows, "width_mismatch": mismatch, "columns": {c: dict(v) for c, v in per.items()},
                "padded": dict(padded), "joint": {"|".join(k): n for k, n in joint.items()}, "alleles": dict(alleles)}

    for line in lines:
        line = line.rstrip("\r\n")
        if not line or line.startswith("##"):
            continue
        if skip:
            skip -= 1
            continue
        rows += 1
        r = line.split("\t")
        if len(r) != width:
            mismatch += 1
        else:
            got = {}
            for c in cols:
                raw = r[idx[c]]
                got[c] = field_class(raw)
                per[c][got[c]] += 1
                if raw != raw.strip():
                    padded[c] += 1
            joint[(coarse(got["standard_error"]), coarse(got[effect]) if effect in got else "absent",
                   coarse(got[p]) if p else "absent")] += 1
            for c in ("effect_allele", "other_allele"):
                a = r[idx[c]].strip().upper()
                if ALLELE.fullmatch(a) is None:
                    alleles[f"allele_token:{c}={a[:20]}"] += 1
        if rows % every == 0:
            checkpoint(snapshot())
    final = snapshot()
    checkpoint(final)
    return {"header": header, "effect_column": effect, "p_column": p, "rows": final["rows"],
            "width_mismatch": final["width_mismatch"],
            "columns": {c: _cap(Counter(v)) for c, v in final["columns"].items()}, "padded": final["padded"],
            "joint_se_effect_p": {"|".join(k): n for k, n in joint.most_common()}, "alleles": _cap(alleles)}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


@app.function(volumes={"/stage": stage_vol, "/diag": diag_vol}, cpu=1, memory=4096, timeout=6 * 3600, retries=2,
              max_containers=8)
def diagnose_ssf_columns(key: str) -> dict:
    diag_vol.reload()
    record = json.loads((ROOT / "collect" / "gwas_catalog" / f"{key}.json").read_text())
    path = ROOT / record["path"]
    tag = f"{key}__{record['sha256'][:16]}"
    done = OUT / "files" / f"{tag}.json"
    if done.is_file():
        return json.loads(done.read_text())
    got = _sha256(path)
    if got != record["sha256"]:
        raise RuntimeError(f"{key}: the file on the volume has sha256 {got}, its record {record['sha256']}")
    state_path = OUT / "work" / f"{tag}.json"
    state = json.loads(state_path.read_text()) if state_path.is_file() else {
        "rows": 0, "width_mismatch": 0, "columns": {}, "padded": {}, "joint": {}, "alleles": {}}

    def checkpoint(s: dict) -> None:
        state_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = state_path.with_name(f".{state_path.name}.tmp")
        tmp.write_text(json.dumps(s))
        os.replace(tmp, state_path)
        diag_vol.commit()
        print(f"{datetime.now(timezone.utc):%H:%M:%S} {key}: {s['rows']} rows", flush=True)

    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        header = next(l for l in (x.rstrip("\r\n") for x in fh) if l and not l.startswith("##")).split("\t")
        result = classify_rows(fh, header, state, CHECKPOINT_ROWS, checkpoint)
    result = {"key": key, "name": record["name"], "sha256": record["sha256"], "record_status": record["status"],
              "resumed_from_row": state["rows"], **result}
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(json.dumps(result, indent=1))
    diag_vol.commit()
    return result


@app.local_entrypoint()
def main(report: str = "", keys: str = ""):
    names = [k for k in keys.split(",") if k] or json.loads(Path(report).read_text())["failed"]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    results = {}
    for key, got in zip(names, diagnose_ssf_columns.map(names, return_exceptions=True), strict=True):
        results[key] = got if isinstance(got, dict) else f"{type(got).__name__}: {got}"[:500]
    out = {"utc": stamp, "what": "aggregate category counts of GWAS-SSF columns; no association value",
           "source_report": report, "keys": names, "files": results}
    local = HERE.parents[1] / "inputs" / "stage_b" / "diagnostics" / f"ssf_columns_{stamp}.json"
    local.parent.mkdir(parents=True, exist_ok=True)
    local.write_text(json.dumps(out, indent=1))
    print(f"wrote {local}")
    for key, r in results.items():
        if isinstance(r, str):
            print(key, "ERROR", r)
            continue
        se = r["columns"].get("standard_error", {})
        print(key, r["rows"], "se:", dict(list(se.items())[:4]), "joint:", dict(list(r["joint_se_effect_p"].items())[:3]))
