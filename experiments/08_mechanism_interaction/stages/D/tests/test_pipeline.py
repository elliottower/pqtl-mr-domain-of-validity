import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from helpers import stub_fit_factory
from v8_manifest import code_sha256, read_inputs

from stage_d import frequentist
from stage_d.designs import save_design
from stage_d.fingerprint import FINGERPRINT_NAME, StaleWork, canonical_sha256
from stage_d.fitting import load_attempts, run_fit
from stage_d.guard import read_manifest, sha256_file
from stage_d.pipeline import (SCRIPT_FILES, IncompleteRun, assemble, load_plan, load_planned_design, prepare,
                              run_frequentist, script_files, script_sha256)
from stage_d.synthetic import make_tables, synthetic_fingerprint, write_stage_outputs

RUNNER = Path(__file__).resolve().parent.parent / "run_stage_d.py"
COMMIT = "0123456789abcdef0123456789abcdef01234567"


class Interrupted(RuntimeError):
    pass


def fit_one(work, plan, fp, spec, fit_fn, fits_dir=None):
    """run_fit on the planned design, as run_stage_d.run_fits calls it."""
    design, design_sha = load_planned_design(work, plan, spec.design_id)
    return run_fit(spec, design, fits_dir or work / "fits", fit_fn, fingerprint=fp.digest, design_sha256=design_sha)


def test_run_fit_sequence_reaches_reduced_model_and_reports_all_attempts(sealed_stages, tmp_path):
    root, expected = sealed_stages
    fp = synthetic_fingerprint(expected)
    plan = prepare(root, expected, tmp_path / "work", fp)
    spec = next(s for s in plan.fits if s.fit_id == "h1__S1__normal15")
    calls = []
    final = fit_one(tmp_path / "work", plan, fp, spec, stub_fit_factory(fail_full=True, calls=calls))
    assert final["outcome"] == "reduced"
    assert [a["attempt"] for a in final["attempts"]] == ["full", "full_rerun", "reduced"]
    assert [(c[1], c[2]) for c in calls] == [(0.95, True), (0.99, True), (0.95, False)]
    assert final["attempts"][1]["sampler"]["draws"] == 4000
    assert "slope" not in final["attempts"][2]["random_effects"]


def test_run_fit_resumes_from_attempts_on_disk(sealed_stages, tmp_path):
    root, expected = sealed_stages
    fp, work = synthetic_fingerprint(expected), tmp_path / "work"
    plan = prepare(root, expected, work, fp)
    spec = next(s for s in plan.fits if s.fit_id == "h4__S1__normal15")
    design_sha = plan.designs[spec.design_id]["sha256"]
    inner = stub_fit_factory(fail_full=True)

    def dies_on_second(design, prior, sampler, keep_slope, seed):
        if sampler.target_accept == 0.99:
            raise Interrupted
        return inner(design, prior, sampler, keep_slope, seed)

    with pytest.raises(Interrupted):
        fit_one(work, plan, fp, spec, dies_on_second)
    assert [a["attempt"] for a in load_attempts(work / "fits" / spec.fit_id, fp.digest, design_sha)] == ["full"]
    calls = []
    final = fit_one(work, plan, fp, spec, stub_fit_factory(fail_full=True, calls=calls))
    assert [c[1:] for c in calls] == [(0.99, True), (0.95, False)]
    assert final["outcome"] == "reduced"


def test_prepare_registers_every_model_and_set(sealed_stages, tmp_path):
    root, expected = sealed_stages
    plan = prepare(root, expected, tmp_path / "work", synthetic_fingerprint(expected))
    roles = {s.role for s in plan.fits}
    assert {"H1", "H4", "H2", "continuous_H1", "continuous_H2", "phase3_success_H1", "phase3_success_H2"} <= roles
    for prefix in ("h1", "h4", "h2", "h1_E", "h2_E"):
        assert {s.prior for s in plan.fits if s.design_id == f"{prefix}__S1"} == {"normal15", "normal1", "student_t3"}
    set_rows = {s.set_id for s in plan.fits if s.role == "set_row"}
    formed = {s["set_id"] for s in plan.sets if s["formed"]} - {"S1"}
    assert set_rows | {k.split("__")[1] for k in plan.skipped if "__S" in k} >= formed
    assert plan.designs["h4__S22"]["kind"] == "h4"
    assert "A" not in plan.designs["h2__S1"]["columns"]
    assert plan.designs["h1__S1"]["columns"][:3] == ["S", "A", "SxA"]
    assert plan.designs["h4__S1"]["columns"][:3] == ["S", "C", "SxC"]


def _stub(spec, focal_mean, fail_full):
    mean = {"h1__S1": focal_mean, "h4__S1": -focal_mean}.get(spec.design_id, 0.0)
    return stub_fit_factory(focal_mean=mean, focal_sd=0.3, fail_full=fail_full and spec.design_id == "h1__S1")


def _fit_all(work, plan, fp, focal_mean, fail_full=False):
    for spec in plan.fits:
        fit_one(work, plan, fp, spec, _stub(spec, focal_mean, fail_full))


def _compute(root, expected, work, focal_mean, fail_full=False, fp=None, reps=20, n_boot=30):
    """prepare, frequentist checks and stub fits: everything `assemble` reads."""
    fp = fp or synthetic_fingerprint(expected)
    plan = prepare(root, expected, work, fp)
    run_frequentist(work, fp, reps=reps, n_boot=n_boot)
    _fit_all(work, plan, fp, focal_mean, fail_full)
    return fp


def _assemble(root, expected, work, out, fp):
    return assemble(root, expected, work, out, fp=fp, runner_path=RUNNER, package_versions={"synthetic": "test"},
                    repo_commit=COMMIT, n_boot=200)


def _run_all(root, expected, work, out, focal_mean, fail_full=False):
    return _assemble(root, expected, work, out, _compute(root, expected, work, focal_mean, fail_full))


def test_end_to_end_on_synthetic_data(tmp_path):
    root = tmp_path / "stages"
    expected = write_stage_outputs(root, *make_tables(seed=5, n_genes=150))
    out = tmp_path / "out"
    res = _run_all(root, expected, tmp_path / "work", out, focal_mean=1.2)
    assert res["gates"]["h1"]["reliable"] and res["gates"]["h4"]["reliable"]
    assert res["decisions"]["H1"]["status"] == "holds"
    assert res["decisions"]["H4"]["status"] == "holds" and res["decisions"]["H4"]["confirmatory"]
    assert res["H2"]["confirmatory"] is False
    results_mtime = (out / "results.json").stat().st_mtime_ns
    produced = [p for p in out.rglob("*") if p.is_file() and p.name not in ("results.json", "MANIFEST.tsv")]
    assert any(p.parent.name == "figures" for p in produced)
    assert all(p.stat().st_mtime_ns >= results_mtime for p in produced)
    entries = read_manifest(out / "MANIFEST.tsv")
    assert {e.path for e in entries} == {str(p.relative_to(out)) for p in produced} | {"results.json"}
    for e in entries:
        assert sha256_file(out / e.path) == e.sha256
    meta = json.loads((out / "results.json").read_text())["meta"]
    assert meta["plan_sha256"].startswith("6cbba084")
    fp = synthetic_fingerprint(expected)
    assert meta["run_fingerprint"] == fp.record() == json.loads((tmp_path / "work" / FINGERPRINT_NAME).read_text())
    assert meta["run_fingerprint"]["environment"] == {"synthetic": "test"}
    assert meta["run_fingerprint"]["components"]["environment_sha256"] == fp.environment_sha256
    assert len(res["descriptive"]["table10_analysis_sets"]) == 30
    assert (meta["repo_commit"], meta["run_token"]) == (COMMIT, fp.run_token)
    assert meta["script_sha256"] == script_sha256(RUNNER.parent / "stage_d", RUNNER) == entries[0].script_sha256
    code = read_inputs(out / "INPUTS.tsv")[-1]
    assert (code.name, code.path, code.sha256, code.source, code.release) == (
        "stage_code", "experiments/08_mechanism_interaction/stages/D", meta["script_sha256"], "git commit", COMMIT)


def test_reduced_model_gives_no_confirmatory_decision_end_to_end(sealed_stages, tmp_path):
    root, expected = sealed_stages
    res = _run_all(root, expected, tmp_path / "work", tmp_path / "out", focal_mean=1.2, fail_full=True)
    assert res["fits"]["h1__S1__normal15"]["outcome"] == "reduced"
    assert res["decisions"]["H1"]["status"] in ("no_confirmatory_decision",)
    assert res["decisions"]["H4"]["status"] == "secondary_estimate"


def test_assemble_refuses_with_missing_fits(sealed_stages, tmp_path):
    root, expected = sealed_stages
    fp = synthetic_fingerprint(expected)
    prepare(root, expected, tmp_path / "work", fp)
    run_frequentist(tmp_path / "work", fp, reps=5, n_boot=5)
    with pytest.raises(IncompleteRun, match="final.json"):
        assemble(root, expected, tmp_path / "work", tmp_path / "out", fp=fp, runner_path=RUNNER, package_versions={},
                 repo_commit=COMMIT)
    assert not (tmp_path / "out").exists()


# ---- work bound to the run fingerprint ---------------------------------------------------------------

def tree(work) -> dict[str, bytes]:
    return {str(p.relative_to(work)): p.read_bytes() for p in sorted(work.rglob("*")) if p.is_file()}


def test_the_fingerprint_changes_with_every_component(sealed_stages):
    _, expected = sealed_stages
    fp = synthetic_fingerprint(expected)
    assert fp.digest == synthetic_fingerprint(dict(expected)).digest and len(fp.digest) == 64
    others = [synthetic_fingerprint(expected, plan_sha256="1" * 64), synthetic_fingerprint({**expected, "A": "2" * 64}),
              synthetic_fingerprint({**expected, "B": "2" * 64}), synthetic_fingerprint({**expected, "C": "2" * 64}),
              synthetic_fingerprint(expected, script_sha256="3" * 64),
              synthetic_fingerprint(expected, pins={"synthetic": "other"}),
              synthetic_fingerprint(expected, run_token="another-token-01"),
              synthetic_fingerprint(expected, environment={"synthetic": "rebuilt"})]
    assert len({fp.digest, *(o.digest for o in others)}) == len(others) + 1
    assert fp.environment_sha256 == canonical_sha256({"synthetic": "test"}) == fp.components()["environment_sha256"]
    assert "environment" not in fp.components() and fp.record()["environment"] == {"synthetic": "test"}


@pytest.mark.parametrize("change", [{"run_token": "another-token-01"}, {"script_sha256": "3" * 64},
                                    {"pins": {"synthetic": "other"}}, {"plan_sha256": "1" * 64},
                                    {"environment": {"synthetic": "rebuilt"}}])
def test_prepare_refuses_a_work_directory_written_under_another_fingerprint(sealed_stages, tmp_path, change):
    root, expected = sealed_stages
    work = tmp_path / "work"
    prepare(root, expected, work, synthetic_fingerprint(expected))
    before = tree(work)
    with pytest.raises(StaleWork, match="not empty"):
        prepare(root, expected, work, synthetic_fingerprint(expected, **change))
    assert tree(work) == before


def test_prepare_refuses_a_non_empty_work_directory_without_a_fingerprint(sealed_stages, tmp_path):
    root, expected = sealed_stages
    work = tmp_path / "work"
    (work / "fits" / "h1__S1__normal15").mkdir(parents=True)
    (work / "fits" / "h1__S1__normal15" / "final.json").write_text("{}")
    with pytest.raises(StaleWork, match="carries fingerprint None"):
        prepare(root, expected, work, synthetic_fingerprint(expected))
    assert tree(work) == {"fits/h1__S1__normal15/final.json": b"{}"}


def test_prepare_under_the_same_fingerprint_keeps_the_plan_and_designs_it_wrote(sealed_stages, tmp_path):
    root, expected = sealed_stages
    work, fp = tmp_path / "work", synthetic_fingerprint(expected)
    first = prepare(root, expected, work, fp)
    before = tree(work)
    assert FINGERPRINT_NAME in before and json.loads(before[FINGERPRINT_NAME])["fingerprint"] == fp.digest == first.fingerprint
    assert all(sha256_file(work / "designs" / f"{d}.npz") == v["sha256"] for d, v in first.designs.items())
    assert prepare(root, expected, work, fp).model_dump_json() == first.model_dump_json() and tree(work) == before


def test_fits_and_frequentist_results_of_another_run_are_refused_in_a_new_work_directory(sealed_stages, tmp_path):
    # the audit's case: work computed under sealed set X survives on the volume, then a run under
    # another fingerprint prepares the same design ids
    root, expected = sealed_stages
    old_fp, new_fp = synthetic_fingerprint(expected), synthetic_fingerprint(expected, run_token="corrected-run-02")
    old, new = tmp_path / "old", tmp_path / "new"
    old_plan = prepare(root, expected, old, old_fp)
    spec = next(s for s in old_plan.fits if s.fit_id == "h1__S1__normal15")
    run_frequentist(old, old_fp, ["h1__S1"], reps=5, n_boot=5)
    fit_one(old, old_plan, old_fp, spec, stub_fit_factory())
    new_plan = prepare(root, expected, new, new_fp)
    assert set(new_plan.designs) == set(old_plan.designs)
    shutil.copytree(old / "fits", new / "fits")
    shutil.copytree(old / "freq", new / "freq")
    calls = []
    with pytest.raises(StaleWork, match="final.json|full.json"):
        fit_one(new, new_plan, new_fp, spec, stub_fit_factory(calls=calls))
    assert calls == []
    with pytest.raises(StaleWork, match="freq/h1__S1.json"):
        run_frequentist(new, new_fp, ["h1__S1"], reps=5, n_boot=5)
    for phase in (lambda: run_frequentist(old, new_fp, ["h1__S1"]), lambda: load_plan(old, new_fp),
                  lambda: _assemble(root, expected, old, tmp_path / "out", new_fp)):
        with pytest.raises(StaleWork, match="carries fingerprint"):
            phase()
    assert not (tmp_path / "out").exists()


def test_assemble_refuses_a_stale_fit_attempt_frequentist_result_or_design(sealed_stages, tmp_path):
    root, expected = sealed_stages
    work, out = tmp_path / "work", tmp_path / "out"
    fp = _compute(root, expected, work, focal_mean=1.2, reps=5, n_boot=5)
    plan = load_plan(work, fp)
    final_path = work / "fits" / "h4__S1__normal15" / "final.json"
    freq_path = work / "freq" / "h4__S1.json"
    design_file = work / "designs" / "h4__S1.npz"
    good = {p: p.read_bytes() for p in (final_path, freq_path, design_file)}

    def tampered(path, mutate, match):
        rec = json.loads(good[path])
        mutate(rec)
        path.write_text(json.dumps(rec))
        with pytest.raises(StaleWork, match=match):
            _assemble(root, expected, work, out, fp)
        path.write_bytes(good[path])

    tampered(final_path, lambda r: r.update(fingerprint="0" * 64), "final.json carries")
    tampered(final_path, lambda r: r.update(design_sha256="0" * 64), "final.json carries")
    tampered(final_path, lambda r: r["attempts"][0].update(fingerprint="0" * 64), "attempt full carries")
    tampered(final_path, lambda r: r["attempts"][0].pop("design_sha256"), "attempt full carries")
    tampered(freq_path, lambda r: r.update(fingerprint="0" * 64), "h4__S1.json carries")
    tampered(freq_path, lambda r: r.update(design_sha256="0" * 64), "h4__S1.json carries")
    design, _ = load_planned_design(work, plan, "h4__S1")
    design.y[0] = 1.0 - design.y[0]
    save_design(design, design_file)
    with pytest.raises(StaleWork, match="h4__S1.npz has sha256"):
        _assemble(root, expected, work, out, fp)
    design_file.write_bytes(good[design_file])
    assert not out.exists()
    res = _assemble(root, expected, work, out, fp)
    assert res["fits"]["h4__S1__normal15"]["design_sha256"] == plan.designs["h4__S1"]["sha256"] == sha256_file(design_file)
    assert res["frequentist"]["h4__S1"]["design_sha256"] == plan.designs["h4__S1"]["sha256"]
    assert {f["fingerprint"] for f in res["fits"].values()} | {f["fingerprint"] for f in res["frequentist"].values()} == {fp.digest}


# ---- script digest ------------------------------------------------------------------------------------

def test_the_script_digest_is_over_relative_paths_of_the_runner_the_modules_the_guard_and_extra_code(tmp_path):
    package = RUNNER.parent / "stage_d"
    files = script_files(package, RUNNER)
    guard = Path(code_sha256.__code__.co_filename).resolve().parent      # where the shared guard was loaded from
    assert {k: v.resolve() for k, v in files.items()} == {
        "run_stage_d.py": RUNNER, **{f"stage_d/{f}": package / f for f in SCRIPT_FILES},
        "run_guard/v8_manifest.py": guard / "v8_manifest.py", "run_guard/v8_run_guard.py": guard / "v8_run_guard.py"}
    assert script_sha256(package, RUNNER) == code_sha256(files)
    other = tmp_path / "power_v9_fast.py"
    other.write_text("X = 1\n")
    with_extra = script_sha256(package, RUNNER, extra={"power_v9/power_v9_fast.py": other})
    assert with_extra == code_sha256({**files, "power_v9/power_v9_fast.py": other}) != script_sha256(package, RUNNER)
    other.write_text("X = 2\n")
    assert script_sha256(package, RUNNER, extra={"power_v9/power_v9_fast.py": other}) != with_extra
    assert {f.name for f in package.glob("*.py")} - set(SCRIPT_FILES) == {"__init__.py", "synthetic.py"}


# ---- a resumed run gives the result tables of a fresh run ---------------------------------------------

class KilledAfter:
    """An on_commit hook that raises on its `n`-th call: the container dies right after a checkpoint."""

    def __init__(self, n: int):
        self.n, self.calls = n, 0

    def __call__(self):
        self.calls += 1
        if self.calls == self.n:
            raise Interrupted(f"killed after commit {self.n}")


# A design file is a zip archive (np.savez_compressed) whose member headers carry the time of writing, so
# the sha256 of a design written by another `prepare` differs although its arrays are the same. Within a
# run the plan pins that sha256; across two runs it is left out of the comparison, with the clock values.
VOLATILE = ("seconds", "generated_utc", "created_utc", "design_sha256")


def _stable(obj):
    if isinstance(obj, dict):
        return {k: _stable(v) for k, v in obj.items() if k not in VOLATILE}
    if isinstance(obj, list):
        return [_stable(v) for v in obj]
    return obj


def _files(d: Path, pattern: str) -> dict[str, bytes]:
    return {str(p.relative_to(d)): p.read_bytes() for p in sorted(d.rglob(pattern)) if p.is_file()}


def _records(d: Path) -> dict[str, dict]:
    return {name: _stable(json.loads(data)) for name, data in _files(d, "*.json").items()}


def _arrays(d: Path) -> dict[str, dict[str, bytes]]:
    out = {}
    for p in sorted(d.glob("*.npz")):
        with np.load(p, allow_pickle=False) as z:
            out[p.name] = {k: f"{z[k].dtype}:{z[k].shape}:".encode() + z[k].tobytes() for k in sorted(z.files)}
    return out


def _lines(work: Path, kind: str) -> dict[str, int]:
    return {p.name: len(p.read_text().splitlines()) for p in (work / "freq_checkpoints" / kind).glob("*.jsonl")}


def test_a_run_killed_and_resumed_gives_the_result_tables_of_a_fresh_run_under_the_same_fingerprint(
        sealed_stages, tmp_path, monkeypatch):
    root, expected = sealed_stages
    fp = synthetic_fingerprint(expected)
    fresh, resumed = tmp_path / "fresh", tmp_path / "resumed"
    _compute(root, expected, fresh, focal_mean=1.2, fail_full=True)
    res_fresh = _assemble(root, expected, fresh, tmp_path / "out_fresh", fp)

    # The resumed run checkpoints every 7 null replicates and every 11 bootstrap resamples (20 and 30 in all),
    # so a null cell takes 3 commits and a bootstrap 3, and a kill can land inside either.
    monkeypatch.setattr(frequentist, "NULL_CHUNK", 7)
    monkeypatch.setattr(frequentist, "BOOTSTRAP_CHUNK", 11)
    plan = prepare(root, expected, resumed, fp)
    first = next(iter(plan.designs))
    cells = len(frequentist.null_cells(plan.designs[first]["kind"]))
    assert cells >= 6
    with pytest.raises(Interrupted):       # cell 0 complete (3 commits), cell 1 killed after 2 of its 3 commits
        run_frequentist(resumed, fp, [first], reps=20, n_boot=30, on_commit=KilledAfter(5))
    assert sorted(_lines(resumed, "null").values()) == [14, 20] and not (resumed / "freq" / f"{first}.json").exists()
    with pytest.raises(Interrupted):       # the rest of cell 1 (1), the other cells (3 each), 2 of 3 bootstrap commits
        run_frequentist(resumed, fp, [first], reps=20, n_boot=30, on_commit=KilledAfter(1 + 3 * (cells - 2) + 2))
    assert set(_lines(resumed, "null").values()) == {20} and len(_lines(resumed, "null")) == cells
    assert _lines(resumed, "bootstrap") == {f"{first}.jsonl": 22} and not (resumed / "freq" / f"{first}.json").exists()
    for kill_at in (4, 37, 90):            # then three more kills wherever they land in the remaining designs
        with pytest.raises(Interrupted):
            run_frequentist(resumed, fp, reps=20, n_boot=30, on_commit=KilledAfter(kill_at))
    assert prepare(root, expected, resumed, fp).model_dump_json() == plan.model_dump_json()
    run_frequentist(resumed, fp, reps=20, n_boot=30)

    # the reduced-model sequence of h1__S1 is killed between attempts, then every fit is run to the end
    spec = next(s for s in plan.fits if s.fit_id == "h1__S1__normal15")
    inner = _stub(spec, 1.2, True)

    def dies_on_rerun(design, prior, sampler, keep_slope, seed):
        if sampler.target_accept == 0.99:
            raise Interrupted
        return inner(design, prior, sampler, keep_slope, seed)

    with pytest.raises(Interrupted):
        fit_one(resumed, plan, fp, spec, dies_on_rerun)
    assert sorted(p.name for p in (resumed / "fits" / spec.fit_id).iterdir()) == ["full.json"]
    _fit_all(resumed, plan, fp, 1.2, fail_full=True)
    res_resumed = _assemble(root, expected, resumed, tmp_path / "out_resumed", fp)

    tables = _files(tmp_path / "out_fresh" / "tables", "*.csv")
    assert len(tables) >= 15 and _files(tmp_path / "out_resumed" / "tables", "*.csv") == tables      # byte for byte
    assert _files(resumed / "freq_checkpoints", "*.jsonl") == _files(fresh / "freq_checkpoints", "*.jsonl")
    assert _arrays(resumed / "designs") == _arrays(fresh / "designs") and len(_arrays(fresh / "designs")) == len(plan.designs)
    assert _records(resumed / "freq") == _records(fresh / "freq") and len(_records(fresh / "freq")) == len(plan.designs)
    assert _records(resumed / "fits") == _records(fresh / "fits")
    assert res_fresh["fits"]["h1__S1__normal15"]["outcome"] == res_resumed["fits"]["h1__S1__normal15"]["outcome"] == "reduced"
    for res in (res_fresh, res_resumed):
        for design in res["designs"].values():
            design.pop("sha256")
    assert _stable(res_resumed) == _stable(res_fresh)
