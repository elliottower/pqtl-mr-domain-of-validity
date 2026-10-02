"""The real stage A, B and C writers, each run in its own uv project on that stage's synthetic test
fixtures, under one copy of the frozen PREREG.md whose log grows in the plan's order: RUN_START A,
SEAL A, RUN_START B, SEAL B, RUN_START C, SEAL C, RUN_START D. Each writer passes its own run guard
on that log (B and C against stage A's sealed MANIFEST.tsv and the files they read from it), and
stage D's guard accepts the three outputs; one changed output byte is refused.

Needs `uv` and the sibling stage directories, so it is skipped in the stage D Modal image, which
carries stage D only. It runs in the combined image of stages/modal_shared_tests.py
(`cross_stage`), which holds all four stage projects."""
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from prereg.log import append
from v8_manifest import read_inputs, read_manifest, sha256_file
from v8_run_guard import parse_log

from stage_d.guard import StageInputError, logged_seals, verify_inputs

STAGES = Path(__file__).resolve().parents[2]
PREREG = Path(os.environ.get("V8_PREREG", STAGES.parent / "PREREG.md"))
UV = shutil.which("uv")
COMMIT = "0123456789abcdef0123456789abcdef01234567"

WRITE_A = """
import os, shutil
from pathlib import Path
from prereg.log import append
import run_stage_a
from tests.synth import DATES, WORLD_COUNTS, World, write_world
root, token = Path(os.environ["V8_TMP"]), "crossstage-a-01"
run_stage_a.PUBLICATION_DATES = dict(DATES)
run_stage_a.REGISTERED_COUNTS = dict(WORLD_COUNTS)     # the synthetic world's counts, not the plan's
prereg = Path(os.environ["V8_CHAIN"])
append(prereg, "2026-10-02", f"RUN_START stage=A token={token}", "nothing run")
approved = {"data": {"id": "9tzfk", "attributes": {"pending_registration_approval": False, "public": True}}}
run_stage_a.main(["--paths-json", str(write_world(World(), root / "in_a")), "--out-dir",
                  str(root / "stages" / "A" / "output"), "--run-token", token, "--prereg", str(prereg),
                  "--repo-root", str(root)], osf_fetch=lambda: approved, repo_commit=os.environ["V8_TEST_COMMIT"])
"""

WRITE_B = """
import os
from pathlib import Path
from prereg.log import append
from stage_b.assemble import build_evidence, collect_unit_dir, write_outputs
from stage_b.launch import authorize
from stage_b.pipeline import DirStore
from tests.test_pipeline import H4, TOOLS, FakeFetcher, StubBackend, fp, hyp, run, unit
root, token = Path(os.environ["V8_TMP"]), "crossstage-b-01"
prereg = Path(os.environ["V8_CHAIN"])
append(prereg, "2026-10-02", f"RUN_START stage=B token={token}", "nothing run")
sealed = authorize(prereg, token, root / "stages" / "A" / "output")
assert set(sealed) == {"hypotheses.csv", "outcome_trait_coding.tsv"}
run(unit(("F_ok",)), FakeFetcher(), StubBackend(H4), DirStore(root / "b_units"))
res, _ = collect_unit_dir(unit(("F_ok",)), root / "b_units", fp(unit(("F_ok",))))
rows = build_evidence([hyp("h1", "decrease", "F_ok")], {"h1": "decode__1_1"}, {"decode__1_1": ("decode", "1_1")},
                      {"decode__1_1": res}, {})
(root / "stages" / "B" / "output").mkdir(parents=True)
(root / "stages" / "B" / "output" / "unit_plan.json").write_text("{}")     # the side file `plan` leaves there
write_outputs(root / "stages" / "B" / "output", rows, [], sorted(Path("stage_b").resolve().glob("*.py")),
              {"hypotheses": root / "stages" / "A" / "output" / "hypotheses.csv"}, [("", root)],
              script_root=Path.cwd(), run_token=token, repo_commit=os.environ["V8_TEST_COMMIT"], tools=TOOLS)
"""

WRITE_C = """
import hashlib, json, os, shutil
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
from prereg.log import append
import run_stage_c
from tests.conftest import AACT_ROWS, CI, write_aact
root, token = Path(os.environ["V8_TMP"]), "crossstage-c-01"
inputs = root / "c_inputs"
ot = inputs / run_stage_c.OT_DIR_REL
ot.mkdir(parents=True)
pq.write_table(pa.table({"drugId": [r[0] for r in CI], "diseaseId": [r[1] for r in CI],
                         "maxClinicalStage": [r[2] for r in CI], "clinicalReportIds": [r[3] for r in CI]}),
               ot / "clinical_indication.parquet")
pq.write_table(pa.table({"id": ["P1", "S1"], "parentId": [None, "P1"]}), ot / "drug_molecule.parquet")
cr = inputs / run_stage_c.CLINICAL_REPORT_REL
cr.parent.mkdir(parents=True)
reports = sorted({r for row in CI for r in row[3]})
pq.write_table(pa.table({"id": reports, "source": ["ClinicalTrials.gov"] * len(reports)}), cr)
pins = root / "pins.json"
pins.write_text(json.dumps({"inputs": {f"ot_{n}": {"sha256": hashlib.sha256((ot / f"{n}.parquet").read_bytes()).hexdigest()}
                                       for n in ("clinical_indication", "drug_molecule")}}))
prereg = Path(os.environ["V8_CHAIN"])
append(prereg, "2026-10-02", f"RUN_START stage=C token={token}", "nothing run")
run_stage_c.run(token, root / "stages" / "A" / "output" / "hypotheses.csv", inputs,
                write_aact(root / "aact_20260930", AACT_ROWS), None, root / "stages" / "C" / "output",
                prereg=prereg, pins=pins, roots=[("", root)], repo_commit=os.environ["V8_TEST_COMMIT"])
"""

pytestmark = pytest.mark.skipif(UV is None or not (STAGES / "A" / "pyproject.toml").is_file(),
                                reason="needs uv and the stage A, B, C projects")


def run_writer(stage: str, code: str, tmp: Path) -> None:
    env = {**os.environ, "V8_TMP": str(tmp), "V8_PREREG": str(PREREG), "V8_CHAIN": str(tmp / "PREREG_CHAIN.md"),
           "V8_TEST_COMMIT": COMMIT,
           "PYTHONPATH": os.pathsep.join([str(STAGES / stage), str(STAGES / stage / "tests")])}
    env.pop("VIRTUAL_ENV", None)
    res = subprocess.run([UV, "run", "--project", str(STAGES / stage), "python", "-c", code], cwd=STAGES / stage,
                         env=env, capture_output=True, text=True, timeout=600)
    if res.returncode != 0:
        raise AssertionError(f"stage {stage} writer failed:\n{res.stdout[-3000:]}\n{res.stderr[-3000:]}")


@pytest.fixture(scope="module")
def written(tmp_path_factory) -> Path:
    tmp = tmp_path_factory.mktemp("cross_stage")
    chain = tmp / "PREREG_CHAIN.md"
    shutil.copyfile(PREREG, chain)
    for stage, code in (("A", WRITE_A), ("B", WRITE_B), ("C", WRITE_C)):
        run_writer(stage, code, tmp)      # the writer logs its own RUN_START and passes its guard
        append(chain, "2026-10-02", f"SEAL stage={stage} manifest_sha256="
               f"{sha256_file(tmp / 'stages' / stage / 'output' / 'MANIFEST.tsv')}", "results not opened")
    return tmp


def sealed_copy(written: Path, tmp_path: Path) -> tuple[Path, Path]:
    root = tmp_path / "stages"
    shutil.copytree(written / "stages", root)
    prereg = tmp_path / "PREREG.md"
    shutil.copyfile(written / "PREREG_CHAIN.md", prereg)
    append(prereg, "2026-10-09", "RUN_START stage=D token=crossstage-d-01", "results not opened")
    return root, prereg


def test_the_log_holds_the_chain_in_the_plans_order(written):
    events = [e.event for e in parse_log((written / "PREREG_CHAIN.md").read_text()) if e.chained]
    chain = [e.split(" token=")[0].split(" manifest_sha256=")[0] for e in events if e.startswith(("RUN_START", "SEAL"))]
    assert chain == ["RUN_START stage=A", "SEAL stage=A", "RUN_START stage=B", "SEAL stage=B", "RUN_START stage=C",
                     "SEAL stage=C"]


def test_the_three_writers_produce_the_shared_format(written):
    for s, required in (("A", "hypotheses.csv"), ("B", "evidence.csv"), ("C", "outcomes.csv")):
        out = written / "stages" / s / "output"
        paths = {e.path for e in read_manifest(out / "MANIFEST.tsv")}
        assert {required, "INPUTS.tsv"} <= paths
        records = read_inputs(out / "INPUTS.tsv")
        for r in records:
            assert not Path(r.path).is_absolute() and ".." not in Path(r.path).parts, (s, r.path)
        code = records[-1]
        assert (code.name, code.path, code.release) == ("stage_code", f"experiments/08_mechanism_interaction/stages/{s}", COMMIT)
        assert {e.script_sha256 for e in read_manifest(out / "MANIFEST.tsv")} == {code.sha256}
        info = json.loads((out / "run_info.json").read_text())
        assert (info["stage"], info["repo_commit"], info["code_sha256"]) == (s, COMMIT, code.sha256)
        assert info["run_token"] == f"crossstage-{s.lower()}-01" and "run_info.json" in paths


def test_d_guard_accepts_the_sealed_outputs_of_the_real_writers(written, tmp_path):
    root, prereg = sealed_copy(written, tmp_path)
    seals = logged_seals(prereg, "crossstage-d-01", root)
    verified = verify_inputs(root, seals)
    assert {s: v.manifest_sha256 for s, v in verified.items()} == seals
    assert "INPUTS.tsv" in verified["C"].files and "trials_snapshot/studies_linked.tsv" in verified["C"].files
    assert {"linkage_audit.csv", "partial_date_boundary.csv", "phase23_diagnostic.csv"} <= set(verified["C"].files)
    assert (root / "B" / "output" / "unit_plan.json").is_file() and "unit_plan.json" not in verified["B"].files
    check = json.loads((root / "A" / "output" / "registered_count_check.json").read_text())
    assert check["all_match"] is True and "registered_count_check.json" in verified["A"].files
    (root / "A" / "output" / "left_over.csv").write_text("x")
    with pytest.raises(StageInputError, match="does not list"):
        verify_inputs(root, seals)


@pytest.mark.parametrize("stage,name", [("A", "hypotheses.csv"), ("B", "evidence.csv"), ("C", "outcomes.csv"),
                                        ("C", "trials_snapshot/studies_linked.tsv")])
def test_one_changed_output_byte_is_refused(written, tmp_path, stage, name):
    root, prereg = sealed_copy(written, tmp_path)
    seals = logged_seals(prereg, "crossstage-d-01", root)
    p = root / stage / "output" / name
    data = bytearray(p.read_bytes())
    i = len(data) // 2
    data[i] = data[i] ^ 0x01
    p.write_bytes(bytes(data))
    with pytest.raises(StageInputError, match=f"{name} sha256"):
        verify_inputs(root, seals)
