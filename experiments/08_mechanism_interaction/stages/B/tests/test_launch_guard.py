import ast
import hashlib
import json
import shutil
from pathlib import Path

import pytest
from prereg.log import append
from v8_manifest import ManifestError, write_manifest
from v8_run_guard import RunNotAuthorized

from stage_b.launch import A_FILES, authorize, spawn_units
from stage_b.schemas import InputContractError, InstrumentUnit, Sentinel

HERE = Path(__file__).resolve().parents[1]
PREREG = HERE.parents[1] / "PREREG.md"
TOKEN = "stageb-token-0001"
START = f"RUN_START stage=B token={TOKEN}"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prereg_with(tmp_path: Path, *events: str) -> Path:
    p = tmp_path / "PREREG.md"
    shutil.copyfile(PREREG, p)
    for e in events:
        append(p, "2026-10-02", e, "nothing run")
    return p


@pytest.fixture
def stage_a(tmp_path) -> Path:
    """A sealed stage A output directory holding the two files stage B reads."""
    out = tmp_path / "stages" / "A" / "output"
    out.mkdir(parents=True)
    (out / "hypotheses.csv").write_text("hypothesis_id,gene_symbol\nh1,G1\nh2,G1\n")
    (out / "outcome_trait_coding.tsv").write_text("outcome_accession\trisk_coded\nGCST1\ttrue\n")
    write_manifest(out, [out / n for n in A_FILES], "c" * 64, [])
    return out


def seal_a(stage_a: Path) -> str:
    return f"SEAL stage=A manifest_sha256={sha(stage_a / 'MANIFEST.tsv')}"


@pytest.fixture
def planned(tmp_path, stage_a) -> tuple[Path, Path]:
    """units.jsonl and the unit_plan.json `plan` writes for it from the sealed stage A files."""
    s = Sentinel(source="ukbppp", assay_id="OID1", rsid="rs1", chrom="1", pos=1000, build="GRCh38", neg_log10_p=20.0)
    lines = [InstrumentUnit(unit_key=f"u{i}", source="ukbppp", assay_id="OID1", gene_symbol="G1", gene_ensembl="ENSG1",
                            platform="Olink", sentinel=s, pqtl_locator="syn1", outcomes=()).model_dump_json()
             for i in range(3)]
    units = tmp_path / "units.jsonl"
    units.write_text("\n".join(lines) + "\n")
    plan = tmp_path / "unit_plan.json"
    plan.write_text(json.dumps({"a_outputs": {n: sha(stage_a / n) for n in A_FILES}, "units_sha256": sha(units)}))
    return units, plan


class FakeSpawn:
    def __init__(self):
        self.calls = []

    def __call__(self, line: str, token: str) -> str:
        self.calls.append(token)
        return f"fc-{len(self.calls)}"


@pytest.mark.parametrize("events", [
    (),
    ("Pre-stage note (B, eQTL Catalogue access): stage B route",),
    ("<SEAL_A>", f"stage B started token={TOKEN}"),
    ("<SEAL_A>", f"RUN_START stage=C token={TOKEN}"),
    (START,),                                                          # stage A not sealed
    ("stage A sealed and committed", START),                           # a mention is not a seal
    (START, "<SEAL_A>"),                                               # sealed after the start
    (f"SEAL stage=A manifest_sha256={'0' * 64}", START),               # another manifest sealed
    ("<SEAL_A>", f"SEAL stage=A manifest_sha256={'0' * 64}", START),   # the last seal before the start counts
])
def test_spawn_refuses_without_a_seal_of_a_before_an_exact_stage_b_run_start_and_spawns_nothing(tmp_path, stage_a,
                                                                                                planned, events):
    spawn = FakeSpawn()
    prereg = prereg_with(tmp_path, *(e.replace("<SEAL_A>", seal_a(stage_a)) for e in events))
    with pytest.raises(RunNotAuthorized):
        spawn_units(prereg, TOKEN, *planned, stage_a, spawn)
    assert spawn.calls == []


def test_the_frozen_prereg_with_its_stage_b_note_does_not_authorize_stage_b(stage_a, planned):
    assert "Pre-stage note (B" in PREREG.read_text()
    spawn = FakeSpawn()
    with pytest.raises(RunNotAuthorized):
        spawn_units(PREREG, TOKEN, *planned, stage_a, spawn)
    assert spawn.calls == []


@pytest.mark.parametrize("before", [(), (f"SEAL stage=A manifest_sha256={'0' * 64}",)])
def test_spawn_passes_with_the_seal_then_the_exact_line_and_hands_every_unit_the_token(tmp_path, stage_a, planned, before):
    spawn = FakeSpawn()
    calls = spawn_units(prereg_with(tmp_path, *before, seal_a(stage_a), START), TOKEN, *planned, stage_a, spawn)
    assert [c["unit_key"] for c in calls] == ["u0", "u1", "u2"]
    assert spawn.calls == [TOKEN] * 3


def test_authorize_returns_the_sealed_hashes_of_the_files_stage_b_reads(tmp_path, stage_a):
    assert authorize(prereg_with(tmp_path, seal_a(stage_a), START), TOKEN, stage_a) == {n: sha(stage_a / n) for n in A_FILES}


@pytest.mark.parametrize("name", A_FILES)
def test_a_stage_a_file_changed_after_the_seal_is_refused_and_nothing_is_spawned(tmp_path, stage_a, planned, name):
    prereg = prereg_with(tmp_path, seal_a(stage_a), START)
    (stage_a / name).write_text((stage_a / name).read_text().replace("1", "2", 1))
    spawn = FakeSpawn()
    with pytest.raises(ManifestError, match=f"{name} sha256"):
        spawn_units(prereg, TOKEN, *planned, stage_a, spawn)
    assert spawn.calls == []


def test_a_manifest_that_omits_a_file_stage_b_reads_is_refused(tmp_path, stage_a):
    write_manifest(stage_a, [stage_a / "hypotheses.csv"], "c" * 64, [])
    with pytest.raises(ManifestError, match="does not list"):
        authorize(prereg_with(tmp_path, seal_a(stage_a), START), TOKEN, stage_a)


def test_units_planned_from_other_stage_a_files_are_not_spawned(tmp_path, stage_a, planned):
    units, plan = planned
    record = json.loads(plan.read_text())
    record["a_outputs"]["hypotheses.csv"] = "d" * 64
    plan.write_text(json.dumps(record))
    spawn = FakeSpawn()
    with pytest.raises(InputContractError, match="planned from stage A files"):
        spawn_units(prereg_with(tmp_path, seal_a(stage_a), START), TOKEN, units, plan, stage_a, spawn)
    assert spawn.calls == []


def test_a_units_file_edited_after_planning_is_not_spawned(tmp_path, stage_a, planned):
    units, plan = planned
    units.write_text(units.read_text().replace("syn1", "syn2"))
    spawn = FakeSpawn()
    with pytest.raises(InputContractError, match="units.jsonl has sha256"):
        spawn_units(prereg_with(tmp_path, seal_a(stage_a), START), TOKEN, units, plan, stage_a, spawn)
    assert spawn.calls == []


def _function(path: Path, name: str) -> ast.FunctionDef:
    tree = ast.parse(path.read_text())
    return next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name)


def test_modal_run_unit_guards_on_the_baked_prereg_and_stage_a_files_before_anything_else():
    fn = _function(HERE / "modal_stage_b.py", "run_unit")
    first = fn.body[0]
    assert isinstance(first, ast.Expr) and isinstance(first.value, ast.Call)
    assert ast.unparse(first.value) == "authorize(PREREG_IMAGE, run_token, A_IMAGE)"
    src = (HERE / "modal_stage_b.py").read_text()
    assert 'PREREG_IMAGE = Path("/root/exp/PREREG.md")' in src and '"PREREG.md"), str(PREREG_IMAGE)' in src
    assert 'A_BAKED = ("MANIFEST.tsv", "hypotheses.csv", "outcome_trait_coding.tsv")' in src
    assert "process_unit(unit, fetcher, RscriptColoc(), store, fingerprint, tools)" in src
    assert ast.unparse(fn.body[1]) == "run_commit(None)"          # an image built from an unclean tree refuses
    assert "image = image.env({COMMIT_ENV: baked_commit(EXP.parents[1])})" in src


def test_launcher_spawns_and_assembles_only_through_the_guard():
    src = (HERE / "launch_stage_b.py").read_text()
    assert "spawn_units(PREREG, a.run_token, UNITS, PLAN_JSON, A_OUTPUT" in src
    for name in ("cmd_spawn", "cmd_assemble"):
        assert ast.unparse(_function(HERE / "launch_stage_b.py", name).body[0]) == "commit = clean_commit(REPO)"
    assemble = ast.unparse(_function(HERE / "launch_stage_b.py", "cmd_assemble").body[1])
    assert assemble == "sealed = authorize(PREREG, a.run_token, A_OUTPUT)"
    assert "repo_commit=commit, tools=common_tools(tools))" in src
    assert "re.search" not in src and "require_logged_start" not in src
