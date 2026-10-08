"""Diagnostic, not part of the stage B run: the CI-versus-p disagreement in the Verma et al. 2024
MVP files (PubMed 39024449), for the round-17 review. Association-blind: aggregate counts only.

For each file named, read on Modal from the stage B volume (mounted read-only) after its sha256 is
checked against its collect record, the data rows are sampled deterministically across the genome:
a row is kept when crc32("<chromosome>:<base_pair_location>") % SAMPLE_MOD == 0, so the same
positions are sampled in every file. On each sampled row with OR, ci_lower, ci_upper and p_value:

    verdict       the production rounding-aware row rule (stage_b.outcome_files.ci_p_rounding_aware_verdict,
                  z = 1.96, 0.1 log10), crossed by reported-p bin x MAF bin
    gap           log10 p_Wald - log10 p_reported, p_Wald from SE = (ln U - ln L) / (2 x 1.96),
                  histogrammed in 0.05 bins by reported-p bin (positive = reported p smaller)
    asymmetry     |ln OR - (ln L + ln U)/2| / (ln U - ln L), in bins: is ln OR the CI's midpoint?

and over every row (not only sampled ones): num_cases and num_controls (distinct values, min and
max), and the minor-allele frequency from effect_allele_frequency against the one from
(num_cases x case_af + num_controls x control_af) / n (author columns the harmoniser does not flip,
so only the orientation-free MAFs are compared). The file's case count is the median of
num_cases over sampled rows; files are binned by it locally. No row value, position or association
value of any hypothesis is returned. Checkpoint every CHECKPOINT_ROWS rows to the volume
`pqtl-v8-diagnostics` (the stage B volume is never written); a restarted call resumes there.

    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \
        modal run experiments/08_mechanism_interaction/stages/B/diagnose_mvp_ci_p.py \
        --studies experiments/08_mechanism_interaction/inputs/stage_b/diagnostics/20261007/gwas_catalog_study_pubmed.json

Every file of PubMed 39024449 in the studies map is diagnosed, passing and failing alike. The result
is written to inputs/stage_b/diagnostics/mvp_ci_p_<utc>.json (gitignored) and to
/mvp_ci_p/report_<utc>.json on the diagnostics volume.
"""
import gzip
import hashlib
import json
import math
import os
import statistics
import zlib
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import modal

PUBMED_ID = 39024449
SAMPLE_MOD = 20
VERSION = "v2"          # v1 compared EAF with the author columns' frequency without regard to orientation
CHECKPOINT_ROWS = 2_000_000
REQUIRED = ("chromosome", "base_pair_location", "odds_ratio", "ci_lower", "ci_upper", "p_value",
            "effect_allele_frequency", "num_cases", "num_controls", "case_af", "control_af")
P_BINS = (1e-8, 1e-5, 1e-3, 0.05, 1.0)
MAF_BINS = (0.001, 0.005, 0.01, 0.05, 0.5)
ASYM_BINS = (0.001, 0.01, 0.05, math.inf)
EAF_GAP_BINS = (0.001, 0.01, 0.05, math.inf)
HERE = Path(__file__).resolve().parent

image = (modal.Image.debian_slim(python_version="3.12")
         .pip_install("numpy==2.1.3", "pandas==2.2.3", "pydantic==2.13.5", "PyYAML==6.0.3")
         .env({"PYTHONPATH": "/root/stageb"})
         .add_local_dir(str(HERE / "stage_b"), "/root/stageb/stage_b", copy=True))
app = modal.App("pqtl-v8-stage-b-diagnose-mvp", image=image)
with image.imports():
    from stage_b.outcome_files import ci_p_rounding_aware_verdict, log10_two_sided_p
stage_vol = modal.Volume.from_name("pqtl-v8-stage-b").with_mount_options(read_only=True)
diag_vol = modal.Volume.from_name("pqtl-v8-diagnostics", create_if_missing=True)
ROOT = Path("/stage/stage_b")
OUT = Path("/diag/mvp_ci_p")


def bin_label(x: float, edges: tuple[float, ...]) -> str:
    """The first edge x lies below or at, as '<=edge'; '>last' beyond the last."""
    for e in edges:
        if x <= e:
            return f"<={e:g}"
    return f">{edges[-1]:g}"


def fnum(s: str) -> float | None:
    try:
        x = float(s)
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def sampled(chrom: str, pos: str) -> bool:
    return zlib.crc32(f"{chrom}:{pos}".encode()) % SAMPLE_MOD == 0


def new_state() -> dict:
    return {"rows": 0, "sampled": 0, "width_mismatch": 0, "verdict": {}, "gap": {}, "asymmetry": {},
            "num_cases": {}, "num_controls": {}, "eaf_gap": {}, "not_parsed": {}}


def scan(lines, header: list[str], state: dict, every: int, checkpoint) -> dict:
    """Every data row after the `state['rows']` already done; `checkpoint(state)` every `every` rows."""
    ix = {h: i for i, h in enumerate(header)}
    missing = [c for c in REQUIRED if c not in ix]
    if missing:
        raise RuntimeError(f"header lacks {missing}")
    c = {k: Counter(state[k]) for k in ("verdict", "gap", "asymmetry", "num_cases", "num_controls", "eaf_gap", "not_parsed")}
    rows, kept, mismatch, skip, width = state["rows"], state["sampled"], state["width_mismatch"], state["rows"], len(header)

    def snapshot() -> dict:
        return {"rows": rows, "sampled": kept, "width_mismatch": mismatch, **{k: dict(v) for k, v in c.items()}}

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
            nca, nco = r[ix["num_cases"]], r[ix["num_controls"]]
            c["num_cases"][nca] += 1
            c["num_controls"][nco] += 1
            eaf, caf, coaf, a, b = (fnum(r[ix[k]]) for k in ("effect_allele_frequency", "case_af", "control_af",
                                                                "num_cases", "num_controls"))
            if None in (eaf, caf, coaf, a, b) or a + b <= 0:
                c["eaf_gap"]["not_parsed"] += 1
            else:
                w = (a * caf + b * coaf) / (a + b)      # author columns, not flipped by the harmoniser: compare MAFs
                c["eaf_gap"][bin_label(abs(min(eaf, 1 - eaf) - min(w, 1 - w)), EAF_GAP_BINS)] += 1
            if sampled(r[ix["chromosome"]], r[ix["base_pair_location"]]):
                kept += 1
                texts = (r[ix["odds_ratio"]], r[ix["ci_lower"]], r[ix["ci_upper"]])
                p_text = r[ix["p_value"]]
                o, lo, up, p = (fnum(t) for t in (*texts, p_text))
                if None in (o, lo, up, p, eaf) or min(o, lo, up) <= 0 or up <= lo or not 0 < p < 1:
                    c["not_parsed"]["or_ci_p_eaf"] += 1
                else:
                    pb = bin_label(p, P_BINS)
                    mb = bin_label(min(eaf, 1 - eaf), MAF_BINS)
                    v = ci_p_rounding_aware_verdict(texts, p_text, False, 1.96)
                    c["verdict"][f"{pb}|{mb}|{v}"] += 1
                    ln_o, ln_l, ln_u = math.log(o), math.log(lo), math.log(up)
                    gap = log10_two_sided_p(ln_o / ((ln_u - ln_l) / (2 * 1.96))) - math.log10(p)
                    c["gap"][f"{pb}|{math.floor(gap / 0.05) * 0.05:.2f}"] += 1
                    c["asymmetry"][bin_label(abs(ln_o - (ln_l + ln_u) / 2) / (ln_u - ln_l), ASYM_BINS)] += 1
        if rows % every == 0:
            checkpoint(snapshot())
    final = snapshot()
    checkpoint(final)
    return final


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def summarize(s: dict) -> dict:
    """Per-file aggregates from the counts: agreement by reported-p bin, the median gap by reported-p
    bin (from the 0.05 histogram), the median num_cases over all rows, and the count of distinct
    num_cases and num_controls values."""
    by_p: dict[str, Counter] = {}
    for k, n in s["verdict"].items():
        pb, _, v = k.split("|")
        by_p.setdefault(pb, Counter())[v] += n
    gaps: dict[str, list[tuple[float, int]]] = {}
    for k, n in s["gap"].items():
        pb, g = k.split("|")
        gaps.setdefault(pb, []).append((float(g), n))
    med = {}
    for pb, xs in gaps.items():
        xs.sort()
        half, run = sum(n for _, n in xs) / 2, 0
        for g, n in xs:
            run += n
            if run >= half:
                med[pb] = g + 0.025
                break
    cases = [(fnum(k), n) for k, n in s["num_cases"].items() if fnum(k) is not None]
    cases.sort()
    half, run, median_cases = sum(n for _, n in cases) / 2, 0, None
    for x, n in cases:
        run += n
        if run >= half:
            median_cases = x
            break
    return {"agree_by_p_bin": {pb: {"agree": v.get("agree", 0), "disagree": v.get("disagree", 0)} for pb, v in by_p.items()},
            "median_gap_by_p_bin": med, "median_num_cases": median_cases,
            "distinct_num_cases": len(s["num_cases"]), "distinct_num_controls": len(s["num_controls"])}


@app.function(volumes={"/stage": stage_vol, "/diag": diag_vol}, cpu=1, memory=4096, timeout=6 * 3600, retries=2,
              max_containers=48)
def diagnose(key: str) -> dict:
    diag_vol.reload()
    record = json.loads((ROOT / "collect" / "gwas_catalog" / f"{key}.json").read_text())
    path = ROOT / record["path"]
    tag = f"{key}__{record['sha256'][:16]}__{VERSION}"
    done = OUT / "files" / f"{tag}.json"
    if done.is_file():
        return json.loads(done.read_text())
    got = _sha256(path)
    if got != record["sha256"]:
        raise RuntimeError(f"{key}: the file on the volume has sha256 {got}, its record {record['sha256']}")
    state_path = OUT / "work" / f"{tag}.json"
    state = json.loads(state_path.read_text()) if state_path.is_file() else new_state()
    resumed = state["rows"]

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
        final = scan(fh, header, state, CHECKPOINT_ROWS, checkpoint)
    validated = json.loads((ROOT / "validation" / "files" / f"{key}.json").read_text())["rows"]
    if final["rows"] != validated:
        raise RuntimeError(f"{key}: read {final['rows']} data rows, the validation read {validated}")
    result = {"key": key, "sha256": record["sha256"], "record_status": record["status"], "resumed_from_row": resumed,
              "header": header, "sample_mod": SAMPLE_MOD, "counts": final, "summary": summarize(final)}
    done.parent.mkdir(parents=True, exist_ok=True)
    done.write_text(json.dumps(result, indent=1))
    diag_vol.commit()
    return result


@app.local_entrypoint()
def main(studies: str):
    names = sorted(k for k, v in json.loads(Path(studies).read_text())["studies"].items() if v["pubmed_id"] == PUBMED_ID)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    results = {}
    for key, got in zip(names, diagnose.map(names, return_exceptions=True), strict=True):
        results[key] = got if isinstance(got, dict) else f"{type(got).__name__}: {got}"[:500]
    code = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out = {"utc": stamp, "what": "aggregate counts of the CI-versus-p disagreement in PubMed 39024449 files; no association value",
           "script_sha256": code, "studies_map": studies, "studies_map_sha256": hashlib.sha256(Path(studies).read_bytes()).hexdigest(),
           "sample_mod": SAMPLE_MOD, "keys": names, "files": results}
    local = HERE.parents[1] / "inputs" / "stage_b" / "diagnostics" / f"mvp_ci_p_{stamp}.json"
    local.parent.mkdir(parents=True, exist_ok=True)
    local.write_text(json.dumps(out, indent=1))
    print(f"wrote {local}")
    errors = {k: r for k, r in results.items() if isinstance(r, str)}
    print(f"{len(names)} files, {len(errors)} errors: {errors}")
