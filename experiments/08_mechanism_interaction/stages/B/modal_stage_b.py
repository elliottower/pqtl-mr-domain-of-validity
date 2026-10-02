"""Thin Modal wrapper for stage B (logic in stage_b/). Nothing here decides anything.

Log `SEAL stage=A manifest_sha256=<sha256>` and then `RUN_START stage=B token=<token>` below the
line of PREREG.md first (`prereg log`): the image bakes PREREG.md and stage A's MANIFEST.tsv,
hypotheses.csv and outcome_trait_coding.tsv, and `run_unit` refuses unless both entries are in the
baked log, the seal is the baked manifest's sha256 and the two files match it (stage_b.launch.
authorize over the shared guard, stages/run_guard/v8_run_guard.py). Commit the stage code and
PREREG.md before deploying: the image carries the commit of a clean tree as V8_REPO_COMMIT and
`run_unit` refuses in an image built from an unclean one. Then deploy, and start work
with spawn so it survives the client disconnecting:
    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \
        modal deploy experiments/08_mechanism_interaction/stages/B/modal_stage_b.py
    uv run --with modal==1.4.3 python -c "import modal; print(modal.Function.from_name('pqtl-v8-stage-b','smoke').remote())"
    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \
        python experiments/08_mechanism_interaction/stages/B/launch_stage_b.py spawn --run-token <token>

Synthetic test suite (no study data), in the run image plus pytest, with this directory and
stages/run_guard baked at /root/exp/stages/; the R coloc tests skipped on the laptop run here. The
summary is printed and written to /vol/tests/modal_tests.json on pqtl-v8-stage-b, the previous
report moved to /vol/tests/superseded/ first. `tests` belongs to its own app, pqtl-v8-stage-b-tests,
because the run app references the secret `pqtl-stage-b`, which exists only once the real run is
being prepared, and the suite needs no credentials:
    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \
        modal run experiments/08_mechanism_interaction/stages/B/modal_stage_b.py::tests

Before the first real run (each step is logged below the line of PREREG.md):
- Modal secret `pqtl-stage-b` with SYNAPSE_PAT and OPENGWAS (read from the 1Password pipe).

Inputs: the deCODE annotation and excluded-variant files are read from the `pqtl-v8-inputs`
volume, mounted read-only at /inputs (08_mechanism_interaction/inputs/decode/, uploaded and
sha256-verified by scripts/upload_inputs_to_modal.py and scripts/modal_verify_inputs.py).

One call per instrument unit. Each step of a unit (positions, pQTL region, LD, each outcome
region, each colocalization, each VEP lookup, splicing, S16) is written to
/vol/stage_b/units/<unit_key>/ and the volume committed before the next step; a restarted
call resumes at the first missing step, and only in a directory bound to this run's fingerprint
(stage_b/checkpoint.py: the unit record, the pinned deCODE files as hashed in the container, the
stage_b/ code digest, the frozen plan hash, and the versions of Python, the Python packages, R,
coloc, susieR, jsonlite, bcftools/htslib and tabix read in the container at the start of the call);
any other directory raises. Pinned: R 4.4.3 (rocker/r-ver, dated CRAN snapshot for
transitive packages), coloc 5.2.3, susieR 0.12.35, jsonlite 2.0.0, htslib/bcftools 1.21, and the
Python wheels below.
"""
import os
import subprocess
from pathlib import Path

import modal

from v8_run_guard import COMMIT_ENV, baked_commit

R_VERSION = "4.4.3"
HTSLIB_VERSION = "1.21"
R_PACKAGES = {"jsonlite": "2.0.0", "susieR": "0.12.35", "coloc": "5.2.3"}
PY_PACKAGES = ["numpy==2.1.3", "pandas==2.2.3", "pydantic==2.13.5", "requests==2.32.3", "openpyxl==3.1.5",
               "synapseclient==4.12.0", "prereg==0.4.2", "provenance-core==0.4.2"]
HERE = Path(__file__).resolve().parent
PREREG_IMAGE = Path("/root/exp/PREREG.md")
INPUTS_MANIFEST_IMAGE = Path("/root/exp/modal_inputs_manifest.json")
A_IMAGE = Path("/root/stages/A/output")
PACKAGE_IMAGE = Path("/root/stageb/stage_b")
A_BAKED = ("MANIFEST.tsv", "hypotheses.csv", "outcome_trait_coding.tsv")
TEST_PACKAGES = ["pytest==9.1.1"]
TEST_IGNORE = [".venv", "**/__pycache__", "**/.pytest_cache", "output", "work", "launch_log.jsonl"]
STAGES_IMAGE = Path("/root/exp/stages")

_r_install = "; ".join(
    f"remotes::install_version('{p}', version='{v}', upgrade='never')" for p, v in R_PACKAGES.items())
image = (
    modal.Image.from_registry(f"rocker/r-ver:{R_VERSION}", add_python="3.12")
    .apt_install("wget", "bzip2", "make", "gcc", "libcurl4-openssl-dev", "libbz2-dev", "liblzma-dev", "zlib1g-dev",
                 "libssl-dev", "libncurses-dev")
    .run_commands(
        f"wget -q https://github.com/samtools/bcftools/releases/download/{HTSLIB_VERSION}/bcftools-{HTSLIB_VERSION}.tar.bz2"
        f" && tar xjf bcftools-{HTSLIB_VERSION}.tar.bz2"
        f" && cd bcftools-{HTSLIB_VERSION}/htslib-{HTSLIB_VERSION} && ./configure --enable-libcurl && make -j4 && make install"
        f" && cd .. && ./configure --with-htslib=system && make -j4 && make install && ldconfig",
        "Rscript -e \"install.packages('remotes')\"",
        f"Rscript -e \"{_r_install}\"",
        "Rscript -e \"stopifnot(packageVersion('coloc') == '5.2.3', packageVersion('susieR') == '0.12.35')\"",
    )
    .pip_install(*PY_PACKAGES)
    .env({"PYTHONPATH": "/root/stageb"})
)
test_image = image
if modal.is_local():
    EXP = HERE.parents[1]   # local only: in the container this file sits directly under /root
    image = (image.add_local_dir(str(HERE / "stage_b"), "/root/stageb/stage_b", copy=True)
             .add_local_file(str(HERE.parent / "run_guard" / "v8_run_guard.py"), "/root/stageb/v8_run_guard.py", copy=True)
             .add_local_file(str(HERE.parent / "run_guard" / "v8_manifest.py"), "/root/stageb/v8_manifest.py", copy=True)
             .add_local_file(str(EXP / "PREREG.md"), str(PREREG_IMAGE), copy=True)
             .add_local_file(str(EXP / "modal_inputs_manifest.json"), str(INPUTS_MANIFEST_IMAGE), copy=True))
    for _name in A_BAKED:
        if (HERE.parent / "A" / "output" / _name).is_file():
            image = image.add_local_file(str(HERE.parent / "A" / "output" / _name), str(A_IMAGE / _name), copy=True)
    image = image.env({COMMIT_ENV: baked_commit(EXP.parents[1])})
    test_image = (image.pip_install(*TEST_PACKAGES)
                  .env({"PYTHONPATH": f"/root/stageb:{STAGES_IMAGE / 'run_guard'}"})
                  .add_local_dir(str(HERE), str(STAGES_IMAGE / "B"), copy=True, ignore=TEST_IGNORE)
                  .add_local_dir(str(HERE.parent / "run_guard"), str(STAGES_IMAGE / "run_guard"), copy=True,
                                 ignore=TEST_IGNORE)
                  .add_local_file(str(HERE.parent / "INTERFACES.md"), str(STAGES_IMAGE / "INTERFACES.md"), copy=True))

app = modal.App("pqtl-v8-stage-b", image=image)
test_app = modal.App("pqtl-v8-stage-b-tests")
vol = modal.Volume.from_name("pqtl-v8-stage-b", create_if_missing=True)
inputs_vol = modal.Volume.from_name("pqtl-v8-inputs").with_mount_options(read_only=True)
secret = modal.Secret.from_name("pqtl-stage-b")
ROOT = Path("/vol/stage_b")
EXP_INPUTS = Path("/inputs/08_mechanism_interaction")
DECODE_INPUTS = EXP_INPUTS / "inputs" / "decode"

with image.imports():
    import numpy as np

    from stage_b.checkpoint import package_sha256, tool_versions, unit_fingerprint, verified_source_pins
    from stage_b.coloc_backend import ColocDataset, ColocTask, RscriptColoc
    from stage_b.fetch import RemoteFetcher
    from stage_b.launch import authorize
    from stage_b.pipeline import DirStore, process_unit, reset_unavailable
    from stage_b.schemas import InstrumentUnit, StageBError
    from v8_run_guard import PLAN_SHA256, run_commit

with test_image.imports():
    from v8_test_report import run_pytest, write_report


@app.function(cpu=2, memory=8192, timeout=24 * 3600, retries=3, volumes={"/vol": vol, "/inputs": inputs_vol},
              secrets=[secret], ephemeral_disk=20 * 1024)
def run_unit(unit_json: str, run_token: str) -> dict:
    authorize(PREREG_IMAGE, run_token, A_IMAGE)
    run_commit(None)
    vol.reload()
    unit = InstrumentUnit.model_validate_json(unit_json)
    tools = tool_versions()
    fingerprint = unit_fingerprint(unit, verified_source_pins(INPUTS_MANIFEST_IMAGE, EXP_INPUTS),
                                   package_sha256(PACKAGE_IMAGE), PLAN_SHA256, tools)
    store = DirStore(ROOT / "units" / unit.unit_key, commit=vol.commit)
    fetcher = RemoteFetcher(os.environ["SYNAPSE_PAT"], os.environ["OPENGWAS"],
                            DECODE_INPUTS / "assocvariants.annotated.txt.gz",
                            DECODE_INPUTS / "assocvariants.excluded.txt.gz", ROOT / "cache")
    result = process_unit(unit, fetcher, RscriptColoc(), store, fingerprint, tools)
    vol.commit()
    return {"unit_key": unit.unit_key, "pqtl_available": result["pqtl_available"], "outcomes": len(result["outcomes"])}


@app.function(timeout=600, volumes={"/vol": vol})
def reset_unit(unit_key: str) -> list[str]:
    vol.reload()
    return reset_unavailable(DirStore(ROOT / "units" / unit_key, commit=vol.commit))


@test_app.function(image=test_image, cpu=2, memory=8192, timeout=3600, volumes={"/vol": vol})
def tests() -> dict:
    """SYNTHETIC: the stage B pytest suite on the copy of this directory baked at /root/exp/stages/B."""
    report = run_pytest(STAGES_IMAGE / "B", [STAGES_IMAGE / "B", STAGES_IMAGE / "B" / "tests", STAGES_IMAGE / "run_guard"])
    vol.reload()
    write_report(Path("/vol/tests/modal_tests.json"), report)
    vol.commit()
    print(f"stage B tests: {report['summary']} (exit {report['returncode']})")
    if report["returncode"] != 0:
        print(report["stdout_tail"], report["stderr_tail"])
    return report


@app.function(cpu=2, memory=4096, timeout=1800)
def smoke() -> dict:
    """coloc.abf and coloc.susie on a tiny synthetic region with a known answer, plus tool
    versions. Region: 60 variants, AR(1) LD r = 0.6^|i-j|; both traits' z-scores are the LD
    column of variant 30 scaled to z = 10 (one shared causal variant), so PP.H4 must dominate
    and the lead variant must be v30. A second pair puts the outcome signal on v10 (distinct
    causal variants), so PP.H3 must dominate."""
    m, causal, other = 60, 30, 10
    idx = np.arange(m)
    ld = 0.6 ** np.abs(idx[:, None] - idx[None, :])
    snp = [f"v{i}" for i in range(m)]
    se = np.full(m, 0.02)

    def ds(zc: np.ndarray, kind: str) -> ColocDataset:
        extra = {"sdY": 1.0} if kind == "quant" else {"s": 0.3}
        return ColocDataset(snp=snp, beta=list(zc * se), varbeta=list(se ** 2), N=20000.0, type=kind,
                            MAF=[0.3] * m, **extra)

    z_same, z_other = 10 * ld[causal], 10 * ld[other]
    tasks = [ColocTask(id="shared", method="abf", p1=1e-4, p2=1e-4, p12=5e-6, d1=ds(z_same, "quant"), d2=ds(z_same, "cc")),
             ColocTask(id="distinct", method="abf", p1=1e-4, p2=1e-4, p12=5e-6, d1=ds(z_same, "quant"), d2=ds(z_other, "cc")),
             ColocTask(id="susie", method="susie", p1=1e-4, p2=1e-4, p12=5e-6, d1=ds(z_same, "quant"),
                       d2=ds(z_same, "cc"), LD=ld.tolist())]
    backend = RscriptColoc()
    res = backend.run(tasks)
    shared, distinct = res["shared"], res["distinct"]
    if not (shared.pp[4] > 0.9 and shared.lead_variant() == f"v{causal}" and distinct.pp[3] > 0.9):
        raise StageBError(f"coloc smoke failed: shared {shared.pp}, lead {shared.lead_variant()}, distinct {distinct.pp}")
    tools = {t: subprocess.run([t, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
             for t in ("tabix", "bcftools")}
    return {"shared_pp": shared.pp, "lead": shared.lead_variant(), "distinct_pp": distinct.pp,
            "susie": res["susie"].model_dump(), "session": backend.session, "tools": tools}
