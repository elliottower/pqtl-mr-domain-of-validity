"""Thin Modal wrapper for the V9 follow-up (logic in followup/; plan PREREG.md, frozen 2026-10-13).

Images: `analysis` is the V8 stage D image definition (PyMC 5.28.5, OpenBLAS, the same pins) plus
openpyxl 3.1.5 and pyarrow 21.0.0; `power` is the V8 power image (numpy 2.1.3, statsmodels 0.14.4);
`rq` is Debian bookworm with r-base-core and r-cran-survival from apt. Every real-run function first
applies the guard: this plan unchanged since its freeze, `RUN_START followup token=<token>` in its log,
V8's A-C seals (stage_d.guard.logged_seals) and a clean committed tree baked into the image.

    uv run --project ../08_mechanism_interaction/stages/D --with modal==1.4.3 modal deploy modal_followup.py
    uv run --project ../08_mechanism_interaction/stages/D --with modal==1.4.3 python launch_followup.py start --run-token <t>
    ... status / collect
    modal run modal_followup.py --what tests      # synthetic test suite, before RUN_START
"""
import json
import os
import subprocess
import urllib.request
from pathlib import Path

import modal

HERE = Path(__file__).resolve().parent
V8 = HERE.parent / "08_mechanism_interaction"
OPENBLAS_APT = ("libopenblas-dev=0.3.21+ds-4", "libopenblas-pthread-dev=0.3.21+ds-4", "libopenblas0-pthread=0.3.21+ds-4")
COMMIT_ENV = "V9_REPO_COMMIT"
OT = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output"
OT_FILES = {
    **{f"evidence_gwas_credible_sets/part-0000{i}-c9f0fb93-8f5e-4da3-93c3-703223d7d0ae-c000.zstd.parquet":
       f"ot_gcs_part{i}.parquet" for i in range(3)},
    "disease/00000000.parquet": "ot_disease.parquet",
}
DOCS_URL = ("https://raw.githubusercontent.com/opentargets/platform-docs/"
            "02978a52a708d8105a04c63cd89f9abe10a4d6cb/evidence.md")
PINNED = {"sun2023_MOESM3_ESM.xlsx": ("/inputs/08_mechanism_interaction/inputs/ukbppp/sun2023_MOESM3_ESM.xlsx",
                                      "e174964b81e0d16f84e485593d4b349e9c24790f60983724d4f9f37681938cc7"),
          "ferkingstad2021_MOESM4_ESM.xlsx": ("/inputs/08_mechanism_interaction/inputs/decode/ferkingstad2021_MOESM4_ESM.xlsx",
                                              "e3f907a11b452f979db8241ed852c1542d00ed5373c39b53c6a79702f8b80ef6")}


def baked_commit() -> str:
    """HEAD when the follow-up and V8 stage code are committed and clean, else "" (the run refuses)."""
    repo = HERE.parents[1]
    paths = ["experiments/09_v8_followup", "experiments/08_mechanism_interaction/stages",
             "experiments/08_mechanism_interaction/power/power_v9", "experiments/08_mechanism_interaction/PREREG.md"]
    dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "--", *paths], capture_output=True, text=True)
    if dirty.returncode != 0 or dirty.stdout.strip():
        return ""
    return subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()


analysis = modal.Image.debian_slim(python_version="3.12").apt_install(*OPENBLAS_APT).pip_install(
    "pymc==5.28.5", "pytensor==2.38.2", "arviz==0.23.4", "numpy==2.5.3", "scipy==1.18.1", "pandas==3.0.6",
    "statsmodels==0.15.0", "tqdm==4.70.1", "matplotlib==3.11.2", "pydantic==2.13.5", "pytest==9.1.1",
    "prereg==0.4.2", "provenance-core==0.4.2", "openpyxl==3.1.5", "pyarrow==21.0.0",
).env({"PYTHONPATH": "/root/v9:/root/stage_d_root:/root/power", "PYTENSOR_FLAGS": "blas__ldflags=-lopenblas"})
power = modal.Image.debian_slim(python_version="3.12").pip_install(
    "numpy==2.1.3", "pandas==2.2.3", "scipy==1.14.1", "statsmodels==0.14.4", "tqdm==4.67.1", "pytest==9.1.1",
).env({"PYTHONPATH": "/root/v9:/root/power"})
rq = modal.Image.debian_slim(python_version="3.12").apt_install("r-base-core", "r-cran-survival")

if modal.is_local():
    commit = baked_commit()
    stage_d = V8 / "stages" / "D"
    common = lambda img: (img.add_local_dir(str(HERE / "followup"), "/root/v9/followup", copy=True, ignore=["__pycache__"])  # noqa: E731
                          .add_local_dir(str(HERE / "tests"), "/root/v9/tests", copy=True, ignore=["__pycache__"])
                          .add_local_file(str(HERE / "PREREG.md"), "/root/v9/PREREG.md", copy=True)
                          .add_local_file(str(HERE / "PREREG.log"), "/root/v9/PREREG.log", copy=True)
                          .env({COMMIT_ENV: commit}))
    analysis = common(analysis.add_local_dir(str(stage_d / "stage_d"), "/root/stage_d_root/stage_d", copy=True,
                                             ignore=["__pycache__"])
                      .add_local_file(str(V8 / "stages" / "run_guard" / "v8_run_guard.py"), "/root/stage_d_root/v8_run_guard.py", copy=True)
                      .add_local_file(str(V8 / "stages" / "run_guard" / "v8_manifest.py"), "/root/stage_d_root/v8_manifest.py", copy=True)
                      .add_local_file(str(V8 / "power" / "power_v9" / "power_v9_fast.py"), "/root/power/power_v9_fast.py", copy=True)
                      .add_local_file(str(V8 / "PREREG.md"), "/root/v8/PREREG.md", copy=True)
                      .add_local_dir(str(V8 / "stages" / "A" / "output"), "/root/v8/stages/A/output", copy=True)
                      .add_local_dir(str(V8 / "stages" / "B" / "output"), "/root/v8/stages/B/output", copy=True)
                      .add_local_dir(str(V8 / "stages" / "C" / "output"), "/root/v8/stages/C/output", copy=True))
    power = common(power.add_local_file(str(V8 / "power" / "power_v9" / "power_v9_logic.py"), "/root/power/power_v9_logic.py", copy=True))
    rq = common(rq)

app = modal.App("pqtl-v9-followup")
vol = modal.Volume.from_name("pqtl-v9-followup", create_if_missing=True)
inputs_vol = modal.Volume.from_name("pqtl-v8-inputs").with_mount_options(read_only=True)
ROOT = Path("/vol")
VOLS = {"/vol": vol}


def guard(token: str) -> str:
    from followup.common import FollowupError, require_run
    require_run(Path("/root/v9/PREREG.md"), Path("/root/v9/PREREG.log"), token)
    commit = os.environ.get(COMMIT_ENV, "")
    if len(commit) != 40:
        raise FollowupError("this image was not built from a clean committed tree")
    return commit


def v8_seals() -> dict[str, str]:
    from followup.common import V8_STAGE_D_TOKEN
    from stage_d.guard import logged_seals
    return logged_seals(Path("/root/v8/PREREG.md"), V8_STAGE_D_TOKEN, Path("/root/v8/stages"))


def fingerprint(token: str) -> str:
    from followup.run import run_fingerprint
    return run_fingerprint(token, Path("/root/v9/followup"))


@app.function(image=analysis, cpu=4, memory=16384, timeout=3 * 3600, volumes={**VOLS, "/inputs": inputs_vol})
def prepare_remote(token: str) -> dict:
    from followup.common import FollowupError, sha256_file, write_json
    from followup.run import prepare
    commit = guard(token)
    vol.reload()
    hashes = {}
    for name, (path, sha) in PINNED.items():
        got = sha256_file(Path(path))
        if got != sha:
            raise FollowupError(f"{path} has sha256 {got}, pinned {sha}")
        hashes[name] = got
    out = prepare(Path("/root/v8/stages"), v8_seals(), Path(PINNED["sun2023_MOESM3_ESM.xlsx"][0]),
                  Path(PINNED["ferkingstad2021_MOESM4_ESM.xlsx"][0]), ROOT)
    write_json(ROOT / "results" / "run_info.json", {"repo_commit": commit, "run_token": token,
                                                   "fingerprint": fingerprint(token), "pinned_spreadsheets": hashes})
    vol.commit()
    return out


def _download(url: str, dest: Path) -> None:
    if dest.exists():
        return
    tmp = dest.with_suffix(".part")
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fh:
                while block := r.read(1 << 20):
                    fh.write(block)
            tmp.replace(dest)
            return
        except OSError:
            if attempt == 4:
                raise


@app.function(image=analysis, cpu=4, memory=32768, timeout=6 * 3600, volumes=VOLS)
def prepare_q4_remote(token: str) -> dict:
    import pandas as pd
    import pyarrow.parquet as pq

    from followup.common import sha256_file, write_json
    from followup.q4_q5 import check_evidence
    from followup.run import load_frame, prepare_q4
    guard(token)
    vol.reload()
    d = ROOT / "inputs"
    d.mkdir(parents=True, exist_ok=True)
    for rel, name in OT_FILES.items():
        _download(f"{OT}/{rel}", d / name)
        vol.commit()
    _download(DOCS_URL, d / "platform_docs_evidence_02978a52.md")
    hashes = {p.name: sha256_file(p) for p in sorted(d.iterdir()) if p.suffix in (".parquet", ".md")}
    write_json(ROOT / "results" / "inputs_sha256.json", {
        "files": hashes, "sources": {**{v: f"{OT}/{k}" for k, v in OT_FILES.items()},
                                     "platform_docs_evidence_02978a52.md": DOCS_URL}})
    vol.commit()
    genes = set(load_frame(ROOT)["gene_ensembl"])
    cols = ["studyLocusId", "targetId", "diseaseId", "resourceScore", "score", "publicationDate"]
    parts, total = [], 0
    for i in range(3):
        t = pq.read_table(d / f"ot_gcs_part{i}.parquet", columns=cols).to_pandas()
        check_evidence(t)                                   # every row, before filtering
        total += len(t)
        parts.append(t.loc[t["targetId"].isin(genes)])
    ev = pd.concat(parts, ignore_index=True)
    disease = pq.read_table(d / "ot_disease.parquet", columns=["id", "descendants"]).to_pandas()
    designs = prepare_q4(ev, disease, ROOT)
    write_json(ROOT / "results" / "q4_evidence_rows.json", {"rows_all_parts": total, "rows_for_v8_genes": int(len(ev))})
    vol.commit()
    return designs


@app.function(image=analysis, cpu=2, memory=8192, timeout=24 * 3600, retries=2, volumes=VOLS)
def q1_boot(model: str, token: str) -> str:
    from followup.common import write_json
    from followup.q1 import run_bootstrap, summarize_bootstrap
    from followup.run import q1_point
    guard(token)
    vol.reload()
    pop, spec, point = q1_point(ROOT, model)
    if spec is None:
        write_json(ROOT / "results" / f"q1_model_{model}.json", {"fitted": False, "reason": f"{len(pop)} runs < 50"})
        vol.commit()
        return model
    path = ROOT / "work" / "q1" / f"{model}.jsonl"
    run_bootstrap(pop, spec, model, path, on_commit=vol.commit)
    write_json(ROOT / "results" / f"q1_model_{model}.json", {
        "fitted": True, "n": int(len(pop)), "genes": int(pop["gene_ensembl"].nunique()), "platform_dropped": spec.drop_platform,
        "levels": spec.levels, "knots": spec.knots, **summarize_bootstrap(point, path)})
    vol.commit()
    return model


@app.function(image=analysis, cpu=4, memory=16384, timeout=24 * 3600, retries=2, volumes=VOLS)
def fit_one(fit_id: str, token: str) -> str:
    from followup.common import sha256_file
    from followup.run import FITS
    from stage_d.bayes import fit
    from stage_d.designs import load_design
    from stage_d.fitting import run_fit
    guard(token)
    vol.reload()
    spec = next(s for s in FITS if s.fit_id == fit_id)
    path = ROOT / "work" / "designs" / f"{spec.design_id}.npz"
    final = run_fit(spec, load_design(path), ROOT / "work" / "fits", fit, vol.commit,
                    fingerprint=fingerprint(token), design_sha256=sha256_file(path))
    vol.commit()
    return final["outcome"]


@app.function(image=power, cpu=2, memory=8192, timeout=24 * 3600, retries=2, volumes=VOLS)
def q6_cell(k_idx: int, r_idx: int, pass_no: int, n: int, token: str) -> str:
    from followup.q6 import run_cell
    from power_v9_logic import load_structure
    guard(token)
    vol.reload()
    st = load_structure(ROOT / "work" / "q6_structure.csv")
    run_cell(st, k_idx, r_idx, pass_no, n, ROOT / "work" / "q6" / f"pass{pass_no}" / f"{k_idx}_{r_idx}.jsonl", vol.commit)
    vol.commit()
    return f"{k_idx}_{r_idx}_{pass_no}"


@app.function(image=rq, cpu=2, memory=8192, timeout=3 * 3600, volumes=VOLS)
def q7_remote(token: str) -> str:
    guard(token)
    vol.reload()
    w = ROOT / "work" / "q7"
    subprocess.run(["Rscript", "/root/v9/followup/q7.R", str(w / "survival.csv"), str(w / "terms_H1.txt"),
                    str(w / "terms_H2.txt"), str(ROOT / "results" / "q7")], check=True)
    vol.commit()
    return "q7"


@app.function(image=power, cpu=2, memory=4096, timeout=3600, volumes=VOLS)
def q6_summary(pass_no: int, token: str) -> dict:
    """After pass 1: the second-pass cells. After pass 2: the summary over both passes."""
    from followup.common import write_json
    from followup.q6 import K_GRID, RATIOS, second_pass_cells, summarize_cell
    guard(token)
    vol.reload()
    cells = {}
    for k in range(len(K_GRID)):
        for r in range(len(RATIOS)):
            recs = []
            for p in range(1, pass_no + 1):
                f = ROOT / "work" / "q6" / f"pass{p}" / f"{k}_{r}.jsonl"
                if f.exists():
                    recs += [json.loads(line) for line in f.read_text().splitlines()]
            cells[(k, r)] = summarize_cell(recs)
    if pass_no == 1:
        return {"second_pass": second_pass_cells({c: v["decision"]["p"] for c, v in cells.items()})}
    first_k = {}
    for r, ratio in enumerate(RATIOS):
        hit = [K_GRID[k] for k in range(len(K_GRID)) if cells[(k, r)]["decision"]["p"] >= 0.8]
        first_k[str(ratio)] = hit[0] if hit else None
    write_json(ROOT / "results" / "q6_calibration.json", {
        "label": "estimated decision probabilities under the V8 planning surrogate; scenarios, not forecasts",
        "cells": [{"k": K_GRID[k], "or_ratio": RATIOS[r], **v} for (k, r), v in sorted(cells.items())],
        "smallest_k_reaching_0.8": first_k})
    vol.commit()
    return first_k


@app.function(image=analysis, cpu=2, memory=8192, timeout=3600, volumes=VOLS)
def assemble_remote(token: str) -> list[str]:
    from followup.common import write_json
    from followup.run import FITS
    guard(token)
    vol.reload()
    fits = {}
    for spec in FITS:
        final = json.loads((ROOT / "work" / "fits" / spec.fit_id / "final.json").read_text())
        att = next((a for a in final["attempts"] if a["passed"]), final["attempts"][-1])
        fits[spec.fit_id] = {"spec": final["spec"], "outcome": final["outcome"], "reported_attempt": att["attempt"],
                             "diagnostics": att["diagnostics"], "summary": {k: att["summary"][k] for k in (
                                 "focal", "focal_posterior_sd", "focal_prior_sd", "prior_dominated", "coefficients",
                                 "stratum_odds_ratios", "marginal")}}
    write_json(ROOT / "results" / "fits.json", fits)
    vol.commit()
    return sorted(str(p.relative_to(ROOT / "results")) for p in (ROOT / "results").rglob("*") if p.is_file())


@app.function(image=analysis, cpu=1, memory=2048, timeout=24 * 3600, volumes=VOLS)
def run_all(token: str) -> list[str]:
    """The whole run in order; every step resumes from the volume, so a restart repeats nothing done."""
    from followup.q1 import MODELS
    from followup.q6 import FIRST_PASS, K_GRID, RATIOS, SECOND_PASS
    from followup.run import FITS
    guard(token)
    prepare_remote.remote(token)
    prepare_q4_remote.remote(token)
    calls = [q1_boot.spawn(m, token) for m in MODELS] + [fit_one.spawn(s.fit_id, token) for s in FITS]
    calls.append(q7_remote.spawn(token))
    first = [q6_cell.spawn(k, r, 1, FIRST_PASS, token) for k in range(len(K_GRID)) for r in range(len(RATIOS))]
    for c in first:
        c.get()
    second = q6_summary.remote(1, token)["second_pass"]
    for c in [q6_cell.spawn(k, r, 2, SECOND_PASS, token) for k, r in second]:
        c.get()
    q6_summary.remote(2, token)
    for c in calls:
        c.get()
    return assemble_remote.remote(token)


@app.function(image=analysis, cpu=4, memory=16384, timeout=3 * 3600)
def tests_analysis() -> str:
    r = subprocess.run(["python", "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests", "--ignore=tests/test_q6_power.py"],
                       cwd="/root/v9", capture_output=True, text=True)
    return r.stdout[-8000:] + r.stderr[-2000:]


@app.function(image=power, cpu=2, memory=4096, timeout=3600)
def tests_power() -> str:
    r = subprocess.run(["python", "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests/test_q6_power.py"],
                       cwd="/root/v9", capture_output=True, text=True)
    return r.stdout[-8000:] + r.stderr[-2000:]


@app.function(image=rq, cpu=2, memory=4096, timeout=1800)
def tests_r() -> str:
    """Q7's R script on a synthetic survival table (no study data)."""
    import csv
    import random
    rnd = random.Random(7)
    w = Path("/tmp/q7")
    w.mkdir(parents=True, exist_ok=True)
    with open(w / "survival.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["hypothesis_id", "gene_ensembl", "time", "status", "cls", "S", "A", "SxA", "platform",
                     "z_log10_neff", "oncology", "group"])
        for i in range(600):
            cls = rnd.choice(["aligned", "blocking", "other"])
            s = int(rnd.random() < 0.2)
            a = int(cls == "aligned")
            wr.writerow([f"h{i}", f"ENSG{rnd.randint(0, 80):011d}", rnd.randint(30, 6000), rnd.choice([0, 1, 2]), cls, s, a,
                         s * a, rnd.randint(0, 1), rnd.gauss(0, 1), rnd.randint(0, 1), f"{cls}|S={s}"])
    (w / "terms_H1.txt").write_text("S\nA\nSxA\nplatform\nz_log10_neff\noncology\n")
    (w / "terms_H2.txt").write_text("S\nplatform\nz_log10_neff\noncology\n")
    r = subprocess.run(["Rscript", "/root/v9/followup/q7.R", str(w / "survival.csv"), str(w / "terms_H1.txt"),
                        str(w / "terms_H2.txt"), str(w / "out")], capture_output=True, text=True)
    files = sorted(p.name for p in (w / "out").iterdir()) if (w / "out").exists() else []
    head = (w / "out" / "fine_gray_H1.csv").read_text()[:600] if (w / "out" / "fine_gray_H1.csv").exists() else ""
    return f"exit {r.returncode}\n{r.stderr[-3000:]}\nfiles {files}\n{head}"


@app.function(image=analysis, cpu=2, memory=8192, timeout=1800, volumes={"/inputs": inputs_vol})
def probe() -> dict:
    """WIRING ONLY, no outcome read: the pinned spreadsheets' hashes and parse, V8's seals, the imports."""
    from followup.common import sha256_file
    from followup.data import decode_strength, ukbppp_strength
    from stage_d.bayes import fit  # noqa: F401  -- the sampler module imports in this image
    hashes = {n: sha256_file(Path(p)) == sha for n, (p, sha) in PINNED.items()}
    sun = ukbppp_strength(Path(PINNED["sun2023_MOESM3_ESM.xlsx"][0]))
    ferk = decode_strength(Path(PINNED["ferkingstad2021_MOESM4_ESM.xlsx"][0]))
    return {"hashes_match": hashes, "ukbppp_assays": len(sun), "decode_assays": len(ferk),
            "ukbppp_truncated": sum(v["truncated"] for v in sun.values()),
            "decode_truncated": sum(v["truncated"] for v in ferk.values()),
            "v8_seals": v8_seals(), "commit_baked": len(os.environ.get(COMMIT_ENV, "")) == 40}


@app.local_entrypoint()
def main(what: str = "tests"):
    """Synthetic checks (`tests`) and the wiring probe (`probe`); the real run goes through launch_followup.py."""
    if what == "probe":
        report = probe.remote()
        (HERE / "work").mkdir(exist_ok=True)
        (HERE / "work" / "modal_probe.json").write_text(json.dumps(report, indent=1))
        print(json.dumps(report, indent=1))
        return
    if what != "tests":
        raise SystemExit("--what is tests or probe")
    out = {"analysis": tests_analysis.spawn(), "power": tests_power.spawn(), "r": tests_r.spawn()}
    report = {k: v.get() for k, v in out.items()}
    path = HERE / "work" / "modal_tests.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(report, indent=1))
    for k, v in report.items():
        print(f"== {k}\n{v[-1500:]}")
