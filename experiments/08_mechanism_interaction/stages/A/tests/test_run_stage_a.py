import hashlib
import io
import json
import shutil
from datetime import date
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from prereg.log import append
from v8_manifest import INPUTS_COLUMNS, MANIFEST_COLUMNS, code_sha256, read_inputs, verify_output_dir
from v8_run_guard import COMMIT_ENV, RunNotAuthorized

import run_stage_a
from stage_a.build import run_stage_a as run_pure
from stage_a.flags import PublicationDatesMissing
from stage_a.loaders import UnknownClinicalStage, load_clinical_indication
from stage_a.schemas import CROSSTAB_COLUMNS, FUNNEL_COLUMNS, HYPOTHESIS_COLUMNS, OUTCOME_SELECTION_COLUMNS, \
    TRAIT_CODING_COLUMNS
from stage_a.registered import REGISTERED_COUNTS, RegisteredCountMismatch, obtained_counts
from tests.synth import DATES, WORLD_COUNTS, World, world_inputs, write_world

PREREG = Path(__file__).resolve().parents[3] / "PREREG.md"
TOKEN = "stagea-token-0001"
COMMIT = "0123456789abcdef0123456789abcdef01234567"
STAGE_A = Path(__file__).resolve().parents[1]
APPROVED = {"data": {"id": "9tzfk", "attributes": {"pending_registration_approval": False, "public": True}}}
PENDING = {"data": {"id": "9tzfk", "attributes": {"pending_registration_approval": True, "public": False}}}


def prereg_with(tmp_path, *events: str) -> Path:
    """A copy of the frozen PREREG.md with `events` appended by `prereg log`."""
    p = tmp_path / "PREREG.md"
    if not p.exists():
        shutil.copyfile(PREREG, p)
    for e in events:
        append(p, "2026-10-02", e, "nothing run")
    return p


SEALED_EXTRAS = {"INPUTS.tsv", "run_info.json", "registered_count_check.json"}
OUTPUTS = {"hypotheses.csv": HYPOTHESIS_COLUMNS, "funnel.csv": FUNNEL_COLUMNS,
           "outcome_gwas_selection.csv": OUTCOME_SELECTION_COLUMNS, "mechanism_crosstab.csv": CROSSTAB_COLUMNS,
           "outcome_trait_coding.tsv": TRAIT_CODING_COLUMNS}


def read(path):
    return pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",")


def _run(tmp_path, world: World, name: str, monkeypatch) -> tuple:
    monkeypatch.setattr(run_stage_a, "PUBLICATION_DATES", dict(DATES))
    monkeypatch.setattr(run_stage_a, "REGISTERED_COUNTS", dict(WORLD_COUNTS))
    paths = write_world(world, tmp_path / f"in_{name}")
    out = tmp_path / f"out_{name}"
    prereg = prereg_with(tmp_path) if (tmp_path / "PREREG.md").exists() else \
        prereg_with(tmp_path, f"RUN_START stage=A token={TOKEN}")
    manifest = run_stage_a.main(["--paths-json", str(paths), "--out-dir", str(out), "--run-token", TOKEN,
                                 "--prereg", str(prereg), "--repo-root", str(tmp_path)], osf_fetch=lambda: APPROVED,
                                repo_commit=COMMIT)
    return out, manifest


def test_cli_writes_the_interface_files_manifest_and_inputs(tmp_path, monkeypatch):
    out, manifest = _run(tmp_path, World(), "a", monkeypatch)
    for name, cols in OUTPUTS.items():
        assert list(read(out / name).columns) == cols
    m = pd.read_csv(manifest, sep="\t", dtype=str, keep_default_na=False)
    assert tuple(m.columns) == MANIFEST_COLUMNS
    rows = m.set_index("path")
    assert set(rows.index) == set(OUTPUTS) | SEALED_EXTRAS
    for name in OUTPUTS:
        assert rows.loc[name, "sha256"] == hashlib.sha256((out / name).read_bytes()).hexdigest()
        assert int(rows.loc[name, "rows"]) == len(read(out / name))
    assert rows["script_sha256"].nunique() == 1
    assert set(verify_output_dir(out, list(OUTPUTS))) == set(OUTPUTS) | SEALED_EXTRAS
    inputs = read_inputs(out / "INPUTS.tsv")
    assert pd.read_csv(out / "INPUTS.tsv", sep="\t", dtype=str).columns.tolist() == list(INPUTS_COLUMNS)
    assert [r.name for r in inputs] == [*run_stage_a.StageAPaths.model_fields, "stage_code"]
    for r in inputs[:-1]:
        assert not Path(r.path).is_absolute() and r.path.startswith("in_a/")
        assert r.sha256 == hashlib.sha256((tmp_path / r.path).read_bytes()).hexdigest()
    code = inputs[-1]
    assert (code.path, code.source, code.release) == ("experiments/08_mechanism_interaction/stages/A", "git commit", COMMIT)
    assert code.sha256 == rows["script_sha256"].iloc[0] == run_stage_a.stage_code_sha256()
    assert json.loads((out / "run_info.json").read_text()) == {
        "stage": "A", "run_token": TOKEN, "repo_commit": COMMIT, "code_sha256": code.sha256}


def test_file_path_matches_in_memory_pipeline(tmp_path, monkeypatch):
    out, _ = _run(tmp_path, World(), "a", monkeypatch)
    from_files = pd.read_csv(out / "hypotheses.csv", dtype=str, keep_default_na=False)
    in_memory = run_pure(world_inputs(World())).hypotheses
    expected = pd.read_csv(io.StringIO(in_memory.to_csv(index=False)), dtype=str, keep_default_na=False)
    pd.testing.assert_frame_equal(from_files, expected)


def test_output_files_identical_when_only_phase_differs_among_phase2_or_later(tmp_path, monkeypatch):
    base = World()
    swapped = World(clinical=[(d, i, "APPROVAL" if s in ("PHASE_2", "PHASE_2_3") else s) for d, i, s in base.clinical])
    a, _ = _run(tmp_path, base, "a", monkeypatch)
    b, _ = _run(tmp_path, swapped, "b", monkeypatch)
    for name in OUTPUTS:
        assert (a / name).read_bytes() == (b / name).read_bytes(), name


def test_no_phase_string_in_any_output_file(tmp_path, monkeypatch):
    out, manifest = _run(tmp_path, World(), "a", monkeypatch)
    for f in [*(out / n for n in OUTPUTS), manifest]:
        text = f.read_text()
        for w in ("PHASE_2", "PHASE_3", "PHASE_1", "APPROVAL", "maxClinicalStage", "clinicalReport"):
            assert w not in text, (f.name, w)


@pytest.mark.parametrize("events,osf", [
    ((), APPROVED),                                                  # no log entry at all
    (("Pre-stage note (A, publication dates) stage A",), APPROVED),  # a mention, not a RUN_START
    ((f"RUN_START stage=B token={TOKEN}",), APPROVED),               # the token belongs to stage B
    ((f"RUN_START stage=A token={TOKEN}",), PENDING),                # logged, registration not approved
])
def test_cli_refuses_without_a_logged_run_and_an_approved_registration(tmp_path, monkeypatch, events, osf):
    monkeypatch.setattr(run_stage_a, "PUBLICATION_DATES", dict(DATES))
    paths = write_world(World(), tmp_path / "in")
    prereg = prereg_with(tmp_path, *events)
    with pytest.raises(RunNotAuthorized):
        run_stage_a.main(["--paths-json", str(paths), "--out-dir", str(tmp_path / "out"), "--run-token", TOKEN,
                          "--prereg", str(prereg)], osf_fetch=lambda: osf)
    assert not (tmp_path / "out").exists()


def test_cli_default_prereg_is_the_frozen_plan_which_logs_no_stage_a_run(tmp_path, monkeypatch):
    monkeypatch.setattr(run_stage_a, "PUBLICATION_DATES", dict(DATES))
    paths = write_world(World(), tmp_path / "in")
    assert run_stage_a.PREREG == PREREG
    with pytest.raises(RunNotAuthorized, match="RUN_START stage=A"):
        run_stage_a.main(["--paths-json", str(paths), "--out-dir", str(tmp_path / "out"), "--run-token", TOKEN],
                         osf_fetch=lambda: APPROVED)


def test_cli_refuses_without_publication_dates(tmp_path, monkeypatch):
    monkeypatch.setattr(run_stage_a, "PUBLICATION_DATES", {"ukbppp": date(2030, 1, 1), "decode": None, "interval": None})
    paths = write_world(World(), tmp_path / "in")
    prereg = prereg_with(tmp_path, f"RUN_START stage=A token={TOKEN}")
    with pytest.raises(PublicationDatesMissing):
        run_stage_a.main(["--paths-json", str(paths), "--out-dir", str(tmp_path / "out"), "--run-token", TOKEN,
                          "--prereg", str(prereg)], osf_fetch=lambda: APPROVED)
    assert not (tmp_path / "out").exists()


# ---- the code digest, the commit and the output directory of a run ---------------------------------

def test_the_script_digest_covers_the_runner_the_package_and_the_shared_guard_by_relative_path():
    files = {str(p.relative_to(STAGE_A)): p for p in [STAGE_A / "run_stage_a.py", *(STAGE_A / "stage_a").glob("*.py")]}
    guard = STAGE_A.parent / "run_guard"
    files.update({f"run_guard/{n}": guard / n for n in ("v8_manifest.py", "v8_run_guard.py")})
    assert {"run_stage_a.py", "stage_a/instruments.py", "stage_a/build.py", "run_guard/v8_manifest.py"} <= set(files)
    assert run_stage_a.stage_code_sha256() == code_sha256(files)


def _main(tmp_path, out, monkeypatch, token=TOKEN, repo_commit=COMMIT, repo_root=None, registered=WORLD_COUNTS):
    monkeypatch.setattr(run_stage_a, "PUBLICATION_DATES", dict(DATES))
    monkeypatch.setattr(run_stage_a, "REGISTERED_COUNTS", dict(registered))
    paths = tmp_path / "in" / "paths.json"
    if not paths.exists():
        write_world(World(), tmp_path / "in")
    prereg = tmp_path / "PREREG.md"
    if not prereg.exists():
        prereg_with(tmp_path, f"RUN_START stage=A token={TOKEN}", "RUN_START stage=A token=stagea-token-0002")
    return run_stage_a.main(["--paths-json", str(paths), "--out-dir", str(out), "--run-token", token, "--prereg",
                             str(prereg), "--repo-root", str(repo_root or tmp_path)], osf_fetch=lambda: APPROVED,
                            repo_commit=repo_commit)


def snapshot(d: Path) -> dict[str, bytes]:
    return {str(p.relative_to(d)): p.read_bytes() for p in sorted(d.rglob("*")) if p.is_file()}


def test_a_run_refuses_an_output_directory_holding_files_of_no_identifiable_run(tmp_path, monkeypatch):
    out = tmp_path / "out"
    out.mkdir()
    (out / "hypotheses.csv").write_text("left over\n")
    with pytest.raises(run_stage_a.OutputDirNotThisRun, match="belongs to run token None"):
        _main(tmp_path, out, monkeypatch)
    assert snapshot(out) == {"hypotheses.csv": b"left over\n"}


def test_a_run_refuses_the_output_directory_of_another_run_token_and_leaves_it_as_it_was(tmp_path, monkeypatch):
    out = tmp_path / "out"
    _main(tmp_path, out, monkeypatch)
    before = snapshot(out)
    with pytest.raises(run_stage_a.OutputDirNotThisRun, match=f"belongs to run token '{TOKEN}', not 'stagea-token-0002'"):
        _main(tmp_path, out, monkeypatch, token="stagea-token-0002")
    assert snapshot(out) == before
    other = tmp_path / "out_2"
    _main(tmp_path, other, monkeypatch, token="stagea-token-0002")       # the other run in its own empty directory
    assert json.loads((other / "run_info.json").read_text())["run_token"] == "stagea-token-0002"
    assert (other / "hypotheses.csv").read_bytes() == before["hypotheses.csv"]


def test_a_finished_run_is_not_overwritten_by_the_same_token(tmp_path, monkeypatch):
    out = tmp_path / "out"
    _main(tmp_path, out, monkeypatch)
    before = snapshot(out)
    with pytest.raises(run_stage_a.OutputDirNotThisRun, match="already holds the finished output"):
        _main(tmp_path, out, monkeypatch)
    assert snapshot(out) == before


def test_an_unfinished_directory_of_the_same_token_is_resumed(tmp_path, monkeypatch):
    out = tmp_path / "out"
    first = snapshot(_main(tmp_path, out, monkeypatch).parent)
    (out / "MANIFEST.tsv").unlink()                      # as left by a run that died before sealing
    (out / "funnel.csv").write_text("half written")
    _main(tmp_path, out, monkeypatch)
    again = snapshot(out)
    assert set(again) == set(first)
    assert {n: again[n] for n in OUTPUTS} == {n: first[n] for n in OUTPUTS}
    assert set(verify_output_dir(out, list(OUTPUTS))) == set(OUTPUTS) | SEALED_EXTRAS


def test_an_empty_output_directory_is_accepted(tmp_path, monkeypatch):
    out = tmp_path / "out"
    out.mkdir()
    assert _main(tmp_path, out, monkeypatch).name == "MANIFEST.tsv"


@pytest.mark.parametrize("baked", ["", "not-a-commit"])
def test_a_run_in_an_image_built_from_an_unclean_tree_is_refused_before_writing(tmp_path, monkeypatch, baked):
    monkeypatch.setenv(COMMIT_ENV, baked)
    with pytest.raises(RunNotAuthorized, match="not built from a clean commit"):
        _main(tmp_path, tmp_path / "out", monkeypatch, repo_commit=None)
    assert not (tmp_path / "out").exists()


def test_a_run_records_the_commit_baked_into_its_image(tmp_path, monkeypatch):
    monkeypatch.setenv(COMMIT_ENV, "f" * 40)
    out = tmp_path / "out"
    _main(tmp_path, out, monkeypatch, repo_commit=None)
    assert json.loads((out / "run_info.json").read_text())["repo_commit"] == "f" * 40
    assert read_inputs(out / "INPUTS.tsv")[-1].release == "f" * 40


def test_a_run_outside_an_image_asks_git_and_refuses_a_directory_that_is_no_clean_repository(tmp_path, monkeypatch):
    monkeypatch.delenv(COMMIT_ENV, raising=False)
    with pytest.raises(RunNotAuthorized):                # tmp_path is not a git repository
        _main(tmp_path, tmp_path / "out", monkeypatch, repo_commit=None)
    assert not (tmp_path / "out").exists()


def test_modal_launcher_checks_the_tree_first_and_bakes_the_commit_into_the_image():
    src = (STAGE_A / "modal_stage_a.py").read_text()
    main = src.split("def main(run_token: str):")[1]
    assert main.index("clean_commit(HERE.parents[3])") < main.index("require_run(") < main.index("run_remote.remote(")
    assert "image = image.env({COMMIT_ENV: baked_commit(REPO_LOCAL)})" in src


# ---- registered outcome-blind counts ---------------------------------------------------------------

def test_the_registered_counts_are_the_numbers_of_the_plans_sample_size_section():
    text = PREREG.read_text(encoding="utf-8").split("## Sample size", 1)[1].split("**Planning power", 1)[0]
    text = " ".join(text.split())
    for sentence in ("S1 has 5,064 held-out hypotheses (4,778 gene–indication pairs, 384 genes, 493 indications)",
                     "The H1 model has 3,253: 582 abundance-aligned in 69 genes and 2,671 function-blocking in 170 "
                     "genes, with 11 genes in both classes; 1,811 are *other*",
                     "698 hypotheses in 168 genes on neurologic or psychiatric indications and 4,366 on other indications",
                     "S4 has 5,355",
                     "H1 has 582 aligned hypotheses in 69 genes and 408 blocking in 30 genes"):
        assert sentence in text, sentence
    assert REGISTERED_COUNTS == {
        "s1_hypotheses": 5064, "s1_gene_indication_pairs": 4778, "s1_genes": 384, "s1_indications": 493,
        "h1_hypotheses": 3253, "h1_aligned_hypotheses": 582, "h1_aligned_genes": 69, "h1_blocking_hypotheses": 2671,
        "h1_blocking_genes": 170, "h1_genes_in_both_classes": 11, "s1_other_hypotheses": 1811,
        "h4_neuro_psych_hypotheses": 698, "h4_neuro_psych_genes": 168, "h4_other_indication_hypotheses": 4366,
        "s4_hypotheses": 5355, "h1_blood_secreted_aligned_hypotheses": 582, "h1_blood_secreted_aligned_genes": 69,
        "h1_blood_secreted_blocking_hypotheses": 408, "h1_blood_secreted_blocking_genes": 30}
    r = REGISTERED_COUNTS      # the plan's own arithmetic
    assert r["h1_aligned_hypotheses"] + r["h1_blocking_hypotheses"] == r["h1_hypotheses"]
    assert r["h1_hypotheses"] + r["s1_other_hypotheses"] == r["s1_hypotheses"]
    assert r["h4_neuro_psych_hypotheses"] + r["h4_other_indication_hypotheses"] == r["s1_hypotheses"]
    assert run_stage_a.REGISTERED_COUNTS is REGISTERED_COUNTS            # what a real run compares against


def test_obtained_counts_on_a_hand_worked_table():
    rows = [  # gene, indication, class, in_s1, heldout, in_s4, neuro_psych, blood
        ("G1", "D1", "aligned", True, True, True, False, True),
        ("G1", "D1", "blocking", True, True, True, False, True),      # same pair, second class: G1 is in both
        ("G1", "D2", "blocking", True, True, True, True, True),
        ("G2", "D1", "blocking", True, True, True, False, False),
        ("G3", "D3", "other", True, True, True, True, False),
        ("G4", "D1", "aligned", True, False, True, False, True),      # pilot key: in the pooled set only
        ("G5", "D4", "blocking", False, True, True, False, True),     # overlap unknown: S4 only
        ("G6", "D5", "aligned", False, True, False, True, True),      # in no counted set
    ]
    hyp = pd.DataFrame(rows, columns=["gene_symbol", "indication_id", "mechanism_class", "in_s1", "heldout", "in_s4",
                                      "neuro_psych", "blood_secreted_hpa"])
    assert obtained_counts(hyp) == {
        "s1_hypotheses": 5, "s1_gene_indication_pairs": 4, "s1_genes": 3, "s1_indications": 3,
        "h1_hypotheses": 4, "h1_aligned_hypotheses": 1, "h1_aligned_genes": 1, "h1_blocking_hypotheses": 3,
        "h1_blocking_genes": 2, "h1_genes_in_both_classes": 1, "s1_other_hypotheses": 1,
        "h4_neuro_psych_hypotheses": 2, "h4_neuro_psych_genes": 2, "h4_other_indication_hypotheses": 3,
        "s4_hypotheses": 6, "h1_blood_secreted_aligned_hypotheses": 1, "h1_blood_secreted_aligned_genes": 1,
        "h1_blood_secreted_blocking_hypotheses": 2, "h1_blood_secreted_blocking_genes": 1}


def test_the_synthetic_world_gives_the_counts_written_out_by_hand():
    assert obtained_counts(run_pure(world_inputs(World())).hypotheses) == WORLD_COUNTS


def test_matching_counts_are_written_and_sealed(tmp_path, monkeypatch):
    out = tmp_path / "out"
    manifest = _main(tmp_path, out, monkeypatch)
    check = json.loads((out / "registered_count_check.json").read_text())
    assert check["all_match"] is True and check["source"].startswith("PREREG.md §Sample size")
    assert check["quantities"] == {k: {"registered": v, "obtained": v, "match": True} for k, v in WORLD_COUNTS.items()}
    assert "registered_count_check.json" in pd.read_csv(manifest, sep="\t", dtype=str)["path"].tolist()


def test_a_run_whose_counts_differ_from_the_registered_ones_is_not_sealed(tmp_path, monkeypatch):
    out = tmp_path / "out"
    with pytest.raises(RegisteredCountMismatch, match="registered count") as err:
        _main(tmp_path, out, monkeypatch, registered=REGISTERED_COUNTS)
    assert "'s1_hypotheses': (5064, 7)" in str(err.value) and "new run token" in str(err.value)
    assert not (out / "MANIFEST.tsv").exists() and not (out / "INPUTS.tsv").exists()
    check = json.loads((out / "registered_count_check.json").read_text())
    assert check["all_match"] is False
    assert check["quantities"]["s1_hypotheses"] == {"registered": 5064, "obtained": 7, "match": False}
    assert check["quantities"]["s4_hypotheses"] == {"registered": 5355, "obtained": 8, "match": False}
    assert set(check["quantities"]) == set(REGISTERED_COUNTS)
    assert (out / "hypotheses.csv").is_file()                 # the tables stay, for the investigation
    # the directory is this token's and unfinished; any other token is refused, so the rerun starts empty
    with pytest.raises(run_stage_a.OutputDirNotThisRun):
        _main(tmp_path, out, monkeypatch, token="stagea-token-0002")


@pytest.mark.parametrize("quantity", ["s1_hypotheses", "h1_genes_in_both_classes", "s4_hypotheses",
                                      "h1_blood_secreted_blocking_genes", "h4_neuro_psych_genes"])
def test_one_differing_count_is_enough_and_is_named(tmp_path, monkeypatch, quantity):
    out = tmp_path / "out"
    with pytest.raises(RegisteredCountMismatch, match=f"1 registered count.*'{quantity}'"):
        _main(tmp_path, out, monkeypatch, registered={**WORLD_COUNTS, quantity: WORLD_COUNTS[quantity] + 1})
    check = json.loads((out / "registered_count_check.json").read_text())
    assert [k for k, q in check["quantities"].items() if not q["match"]] == [quantity]
    assert not (out / "MANIFEST.tsv").exists()


def test_a_registered_quantity_without_a_counterpart_is_refused(tmp_path, monkeypatch):
    with pytest.raises(RegisteredCountMismatch, match="no counterpart"):
        _main(tmp_path, tmp_path / "out", monkeypatch, registered={**WORLD_COUNTS, "s99_hypotheses": 1})


# ---- the deCODE diagnostic in funnel.csv ------------------------------------------------------------

def test_funnel_reports_the_decode_diagnostic_and_the_lists_do_not_depend_on_it(tmp_path, monkeypatch):
    out, _ = _run(tmp_path, World(), "a", monkeypatch)
    f = read(out / "funnel.csv").set_index("step")
    diag = f[f.index.str.startswith("diagnostic:")]
    assert list(diag.index) == ["diagnostic:decode_strict_seqids_st02_symbols_differ_from_st01_gene",
                                "diagnostic:decode_strict_symbols_only_through_st02_disagreement"]
    assert list(diag["remaining"]) == [0, 0] and list(diag["excluded"]) == [0, 0]      # the default world agrees
    assert diag["reason"].str.startswith("diagnostic, changes no list").all()
    # one strict SeqId whose ST02 row names another gene: both rows count it, and GZ joins the strict list
    w = World()
    w.decode_st02 = [*w.decode_st02, {"gene": "GZ", "SeqId": "1000-1", "cis_trans": "cis"}]
    out_b, _ = _run(tmp_path, w, "b", monkeypatch)
    fb = read(out_b / "funnel.csv").set_index("step")
    assert list(fb.loc[diag.index, "remaining"]) == [1, 1]
    assert "of 2 strict deCODE SeqIds" in fb.loc[diag.index[0], "reason"]
    assert "of 3 gene symbols" in fb.loc[diag.index[1], "reason"]
    same = [s for s in f.index if not s.startswith("diagnostic:")]
    assert f.loc[same].equals(fb.loc[same])
    assert (out / "hypotheses.csv").read_bytes() == (out_b / "hypotheses.csv").read_bytes()


def test_loader_reduces_stage_to_boolean(tmp_path):
    p = tmp_path / "ci.parquet"
    stages = ["PHASE_1", "PHASE_2", "PHASE_2_3", "PHASE_3", "PREAPPROVAL", "APPROVAL", "EARLY_PHASE_1", None]
    pq.write_table(pa.Table.from_pydict({"drugId": [f"D{i}" for i in range(len(stages))], "diseaseId": ["M"] * len(stages),
                                         "maxClinicalStage": stages}), p)
    ci = load_clinical_indication(p)
    assert list(ci.columns) == ["drugId", "diseaseId", "phase2plus"]
    assert list(ci["phase2plus"]) == [False, True, True, True, True, True, False, False]


def test_loader_rejects_unregistered_stage(tmp_path):
    p = tmp_path / "ci.parquet"
    pq.write_table(pa.Table.from_pydict({"drugId": ["D"], "diseaseId": ["M"], "maxClinicalStage": ["PHASE_4"]}), p)
    with pytest.raises(UnknownClinicalStage):
        load_clinical_indication(p)


def test_paths_json_rejects_unknown_key(tmp_path):
    p = tmp_path / "paths.json"
    p.write_text(json.dumps({"not_an_input": "x"}))
    with pytest.raises(KeyError):
        run_stage_a.resolve_paths(p, tmp_path, tmp_path, tmp_path)


def test_trait_coding_file_is_what_stage_b_reads(tmp_path, monkeypatch):
    out, _ = _run(tmp_path, World(), "a", monkeypatch)
    coding = pd.read_csv(out / "outcome_trait_coding.tsv", sep="\t", dtype=str, keep_default_na=False)
    hyp = pd.read_csv(out / "hypotheses.csv", dtype=str, keep_default_na=False)
    # stage_b.units.load_trait_coding: first two columns, a boolean it lower-cases, one row per accession.
    assert list(coding.columns[:2]) == ["outcome_accession", "risk_coded"]
    assert set(coding["risk_coded"].str.lower()) == {"true"}
    assert coding["outcome_accession"].is_unique
    assert set(coding["outcome_accession"]) == set(hyp["outcome_accession"])
    by_acc = coding.set_index("outcome_accession")
    for r in hyp.itertuples():
        assert by_acc.loc[r.outcome_accession, "outcome_source"] == r.outcome_source
        assert (by_acc.loc[r.outcome_accession, "n_case"], by_acc.loc[r.outcome_accession, "n_control"]) == (
            r.outcome_n_case, r.outcome_n_control)


def test_inputs_root_resolves_the_volume_layout(tmp_path):
    p = run_stage_a.resolve_paths(None, tmp_path / "vol" / "08_mechanism_interaction", tmp_path / "exp", tmp_path / "repo")
    vol = tmp_path / "vol" / "08_mechanism_interaction"
    assert p.ot_clinical_indication.path == vol / "feasibility/v2_all_indications/inputs/ot_26.09/clinical_indication.parquet"
    assert p.ukbppp_cis_list.path == vol / "inputs/ukbppp/ukbppp_cis_pqtl_genes.json"
    assert p.karim2026_workbook.path == vol / "inputs/karim2026/v2_supp_tables.xlsx"
    assert p.decode_supplement.path == tmp_path / "exp" / "feasibility/inputs/ferkingstad2021_MOESM4_ESM.xlsx"
    assert p.frozen_candidates_v34.path == tmp_path / "repo" / "data/frozen_candidates_v34.csv"
