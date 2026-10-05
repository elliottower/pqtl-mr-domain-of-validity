import ast
import hashlib
import json
import shutil
from pathlib import Path

import pytest
from prereg.log import append
from v8_manifest import ManifestError, write_manifest
from v8_run_guard import RunNotAuthorized

from stage_b.launch import A_FILES, authorize, private_copy, sealed_stage_b, spawn_collect, spawn_units
from stage_b.schemas import CollectRecord, InputContractError, InstrumentUnit, OutcomeSpec, Sentinel
from stage_b.validate import RULES, result_sha256

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
    shared = OutcomeSpec(accession="GCST1", source="gwas_catalog", n_case=10, n_control=90, risk_coded=True)
    lines = [InstrumentUnit(unit_key=f"u{i}", source="ukbppp", assay_id=f"OID{i}", gene_symbol="G1", gene_ensembl="ENSG1",
                            platform="Olink", sentinel=s, pqtl_locator="syn1", outcomes=(shared,)).model_dump_json()
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


FILES = {("ukbppp", "OID0"), ("ukbppp", "OID1"), ("ukbppp", "OID2"), ("ukbppp_rsid_map", "1"), ("gwas_catalog", "GCST1")}


def test_collect_spawns_one_call_per_file_behind_the_same_guard(tmp_path, stage_a, planned):
    spawn = FakeSpawn()
    with pytest.raises(RunNotAuthorized):
        spawn_collect(prereg_with(tmp_path, seal_a(stage_a)), TOKEN, *planned, stage_a, spawn)
    assert spawn.calls == []
    prereg = prereg_with(tmp_path, seal_a(stage_a), START)
    calls = spawn_collect(prereg, TOKEN, *planned, stage_a, spawn)
    assert {(c["source"], c["key"]) for c in calls} == FILES and len(calls) == 5 == len(spawn.calls)   # three units, five files
    only = spawn_collect(prereg, TOKEN, *planned, stage_a, FakeSpawn(), sources=("gwas_catalog",))
    assert [(c["source"], c["key"]) for c in only] == [("gwas_catalog", "GCST1")]
    units, plan = planned
    units.write_text(units.read_text().replace("GCST1", "GCST2"))
    with pytest.raises(InputContractError, match="units.jsonl has sha256"):
        spawn_collect(prereg, TOKEN, units, plan, stage_a, spawn)
    assert len(spawn.calls) == 5


def test_units_are_not_spawned_while_a_file_has_no_collect_record(tmp_path, stage_a, planned):
    prereg, spawn = prereg_with(tmp_path, seal_a(stage_a), START), FakeSpawn()
    with pytest.raises(InputContractError, match="1 files have no collect record"):
        spawn_units(prereg, TOKEN, *planned, stage_a, spawn, recorded=FILES - {("gwas_catalog", "GCST1")})
    assert spawn.calls == []
    assert len(spawn_units(prereg, TOKEN, *planned, stage_a, spawn, recorded=FILES)) == 3


def catalog_record(key: str = "GCST1", sha: str = "a" * 64, status: str = "collected", detail: str = "") -> CollectRecord:
    return CollectRecord(status=status, source="gwas_catalog", key=key, name=f"{key}.h.tsv.gz", path=f"raw/gwas_catalog/{key}",
                         bytes=5, sha256=sha, detail=detail)


def report_for(*records: CollectRecord, mode: str = "native_se") -> dict:
    return {"rules": RULES, "files": {r.key: {"sha256": r.sha256, "passed": r.status == "collected", "rules": RULES,
                                              "uncertainty_mode": mode, "header_sha256": "c" * 64} for r in records}}


def bound(record: CollectRecord, report: dict) -> CollectRecord:
    """`record` as apply_validation leaves it when its file passed: bound to the mode, header and result."""
    got = report["files"][record.key]
    return record.model_copy(update={"uncertainty_mode": got["uncertainty_mode"], "header_sha256": got["header_sha256"],
                                     "validation_sha256": result_sha256(got)})


def test_units_are_not_spawned_without_a_validation_report_or_with_one_that_is_stale(tmp_path, stage_a, planned):
    prereg, spawn = prereg_with(tmp_path, seal_a(stage_a), START), FakeSpawn()
    record = catalog_record()
    with pytest.raises(InputContractError, match="no pre-analysis validation report"):
        spawn_units(prereg, TOKEN, *planned, stage_a, spawn, recorded=FILES, validation=(None, [record]))
    recollected = catalog_record(sha="b" * 64)                      # the file collected again after the report
    with pytest.raises(InputContractError, match="GCST1: the report validated other bytes than the record names"):
        spawn_units(prereg, TOKEN, *planned, stage_a, spawn, recorded=FILES, validation=(report_for(record), [recollected]))
    with pytest.raises(InputContractError, match="GCST2: not in the report"):
        spawn_units(prereg, TOKEN, *planned, stage_a, spawn, recorded=FILES,
                    validation=(report_for(record), [bound(record, report_for(record)), catalog_record("GCST2")]))
    with pytest.raises(InputContractError, match=r"GCST1: the record's uncertainty mode \(none\) is not bound"):
        spawn_units(prereg, TOKEN, *planned, stage_a, spawn, recorded=FILES, validation=(report_for(record), [record]))
    other_mode = bound(record, report_for(record, mode="or_ci_derived_se"))
    with pytest.raises(InputContractError, match=r"GCST1: the record's uncertainty mode \(or_ci_derived_se\) is not bound"):
        spawn_units(prereg, TOKEN, *planned, stage_a, spawn, recorded=FILES, validation=(report_for(record), [other_mode]))
    assert spawn.calls == []
    indexed = CollectRecord(status="remote_indexed", source="gwas_catalog", key="GCST3", name="x.h.tsv.gz")
    record = bound(record, report_for(record))
    calls = spawn_units(prereg, TOKEN, *planned, stage_a, spawn, recorded=FILES, validation=(report_for(record), [record, indexed]))
    assert len(calls) == 3
    src = (HERE / "launch_stage_b.py").read_text()
    assert 'validation=(state["outcome_validation"],' in src and 'state["gwas_catalog_records"]' in src


def test_the_whole_files_are_purged_only_after_stage_b_is_sealed(tmp_path, stage_a):
    seal_b = f"SEAL stage=B manifest_sha256={sha(stage_a / 'MANIFEST.tsv')}"      # any manifest stands in for stage B's
    for events in ((), (seal_a(stage_a), START), (seal_a(stage_a), START, "stage B sealed and committed")):
        with pytest.raises(RunNotAuthorized, match="no 'SEAL stage=B"):
            sealed_stage_b(prereg_with(tmp_path, *events))
    prereg = prereg_with(tmp_path, seal_a(stage_a), START, seal_b)
    assert sealed_stage_b(prereg) == sealed_stage_b(prereg, stage_a / "MANIFEST.tsv") == sha(stage_a / "MANIFEST.tsv")
    (stage_a / "MANIFEST.tsv").write_text((stage_a / "MANIFEST.tsv").read_text() + "\n")
    with pytest.raises(RunNotAuthorized, match="the log seals stage B as"):
        sealed_stage_b(prereg, stage_a / "MANIFEST.tsv")
    purge = _function(HERE / "modal_stage_b.py", "purge_raw_files")
    assert ast.unparse(purge.body[1]) == "sealed_stage_b(PREREG_IMAGE)"                # the docstring, then the guard
    launcher = _function(HERE / "launch_stage_b.py", "cmd_purge")
    assert [ast.unparse(n) for n in launcher.body[1:3]] == [
        "commit = clean_commit(REPO)", "seal = sealed_stage_b(PREREG, HERE / 'output' / MANIFEST_NAME)"]


@pytest.mark.parametrize("name", ["run_unit", "collect_file", "plan_remote", "validate_outcomes", "validate_outcome_file",
                                  "restore_validated_records"])
def test_every_real_run_modal_function_guards_before_anything_else(name):
    fn = _function(HERE / "modal_stage_b.py", name)
    assert ast.unparse(fn.body[0]) == "authorize(PREREG_IMAGE, run_token, A_IMAGE)"
    assert ast.unparse(fn.body[1]) == "run_commit(None)"          # an image built from an unclean tree refuses


def test_no_token_or_link_reaches_a_unit_a_record_or_a_log():
    wrapper = (HERE / "modal_stage_b.py").read_text()
    unit_call = ast.unparse(_function(HERE / "modal_stage_b.py", "run_unit"))
    assert "DECODE_FOLDER_TOKEN" not in unit_call and "SYNAPSE_PAT" not in unit_call      # analyze holds no download credential
    assert wrapper.count('os.environ.get("DECODE_FOLDER_TOKEN", "")') == 1                 # read once, by collect_file
    launcher = (HERE / "launch_stage_b.py").read_text()
    assert "os.environ" not in launcher and "decode-urls" not in launcher and "decode_url" not in launcher


def test_modal_run_unit_guards_on_the_baked_prereg_and_stage_a_files_before_anything_else():
    fn = _function(HERE / "modal_stage_b.py", "run_unit")
    first = fn.body[0]
    assert isinstance(first, ast.Expr) and isinstance(first.value, ast.Call)
    assert ast.unparse(first.value) == "authorize(PREREG_IMAGE, run_token, A_IMAGE)"
    src = (HERE / "modal_stage_b.py").read_text()
    assert 'PREREG_IMAGE = Path("/root/exp/PREREG.md")' in src and '"PREREG.md"), str(PREREG_IMAGE)' in src
    assert 'A_BAKED = ("MANIFEST.tsv", "hypotheses.csv", "outcome_trait_coding.tsv")' in src
    assert "process_unit(unit, fetcher, RscriptColoc(), store, fingerprint, tools, collected)" in src
    assert "collected = volume_collect_digest(fetcher.root, unit)" in src       # the records the fetcher reads from
    assert ast.unparse(fn.body[1]) == "run_commit(None)"          # an image built from an unclean tree refuses
    assert "image = image.env({COMMIT_ENV: baked_commit(EXP.parents[1])})" in src


def test_launcher_spawns_and_assembles_only_through_the_guard():
    src = (HERE / "launch_stage_b.py").read_text()
    assert "spawn_units(PREREG, a.run_token, UNITS, PLAN_JSON, A_OUTPUT" in src
    assert "spawn_collect(PREREG, a.run_token, UNITS, PLAN_JSON, A_OUTPUT" in src
    for name in ("cmd_plan", "cmd_collect", "cmd_validate", "cmd_restore_validated", "cmd_spawn", "cmd_assemble"):
        assert ast.unparse(_function(HERE / "launch_stage_b.py", name).body[0]) == "commit = clean_commit(REPO)"
    for name in ("cmd_validate", "cmd_restore_validated"):
        assert ast.unparse(_function(HERE / "launch_stage_b.py", name).body[1]) == "authorize(PREREG, a.run_token, A_OUTPUT)"
    assemble = ast.unparse(_function(HERE / "launch_stage_b.py", "cmd_assemble").body[1])
    assert assemble == "sealed = authorize(PREREG, a.run_token, A_OUTPUT)"
    assert "repo_commit=commit, tools=common_tools(tools)," in src and "collected=collected)" in src
    assert 'units_dir = private_copy(a.units_dir, REPO, EXP / "inputs")' in src
    assert 'load_collect_records(private_copy(a.collect_dir, REPO, EXP / "inputs"))' in src
    assert "digest = collect_digest(u, lambda source, key: records[(source, key)])" in src
    assert "unit_fingerprint(u, pins, code, PLAN_SHA256, tools[u.unit_key], digest), digest)" in src
    assert ast.unparse(_function(HERE / "launch_stage_b.py", "cmd_plan").body[1]) == "authorize(PREREG, a.run_token, A_OUTPUT)"
    assert "re.search" not in src and "require_logged_start" not in src


def test_unit_directories_are_copied_only_outside_the_repository_or_under_the_gitignored_inputs_directory(tmp_path):
    repo = tmp_path / "repo"
    private = repo / "experiments" / "08_mechanism_interaction" / "inputs"
    inside = repo / "experiments" / "08_mechanism_interaction" / "stages" / "B" / "output" / "units"
    for d in (private / "stage_b" / "units", inside, tmp_path / "elsewhere" / "units"):
        d.mkdir(parents=True)
    assert private_copy(private / "stage_b" / "units", repo, private) == private / "stage_b" / "units"
    assert private_copy(tmp_path / "elsewhere" / "units", repo, private) == tmp_path / "elsewhere" / "units"
    for path in (inside, repo / "units", repo / "experiments" / "08_mechanism_interaction" / "inputs_copy"):
        with pytest.raises(InputContractError, match="regional extracts, which are not committed"):
            private_copy(path, repo, private)
