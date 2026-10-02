import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest
from prereg.log import append
from v8_manifest import ManifestError, code_sha256, read_inputs, verify_output_dir, write_manifest
from v8_run_guard import RunNotAuthorized

import run_stage_c
from chembl import ChemblNotVerified
from ctgov import AactNotVerified
from run_stage_c import run
from tests.conftest import HYPOTHESES
from tests.test_chembl import make_db

PREREG = Path(__file__).resolve().parents[3] / "PREREG.md"
TOKEN = "stagec-token-0001"
COMMIT = "0123456789abcdef0123456789abcdef01234567"
B_SEAL = "b" * 64   # stage C needs stage B's seal in the log and opens none of its files


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def chain(a_manifest: Path) -> tuple[str, ...]:
    """The log entries a stage C run needs, in order."""
    return (f"SEAL stage=A manifest_sha256={sha(a_manifest)}", f"SEAL stage=B manifest_sha256={B_SEAL}",
            f"RUN_START stage=C token={TOKEN}")


def prereg_with(path: Path, *events: str) -> Path:
    shutil.copyfile(PREREG, path)
    for e in events:
        append(path, "2026-10-02", e, "nothing run")
    return path


@pytest.fixture
def layout(tmp_path, e2e):
    """An inputs root laid out like the experiment directory, built from the e2e fixture."""
    root = tmp_path / "inputs_root"
    ot = root / run_stage_c.OT_DIR_REL
    ot.mkdir(parents=True)
    src = tmp_path / "ot"
    for name in ("clinical_indication", "drug_molecule"):
        (ot / f"{name}.parquet").write_bytes((src / f"{name}.parquet").read_bytes())
    cr = root / run_stage_c.CLINICAL_REPORT_REL
    cr.parent.mkdir(parents=True)
    cr.write_bytes((src / "clinical_report.parquet").read_bytes())
    pins = tmp_path / "pins.json"
    pins.write_text(json.dumps({"inputs": {f"ot_{n}": {"sha256": hashlib.sha256((ot / f"{n}.parquet").read_bytes()).hexdigest()}
                                           for n in ("clinical_indication", "drug_molecule")}}))
    hypotheses = tmp_path / "stages" / "A" / "output" / "hypotheses.csv"
    a_manifest = write_manifest(hypotheses.parent, [hypotheses], "c" * 64, [])
    prereg = prereg_with(tmp_path / "PREREG.md", *chain(a_manifest))
    chembl = tmp_path / "chembl_37"
    chembl.mkdir()
    sha = make_db(chembl / "chembl_37.db")
    return {"tmp": tmp_path, "root": root, "pins": pins, "prereg": prereg, "chembl": chembl, "chembl_sha": sha, "aact": e2e["aact"],
            "hypotheses": hypotheses, "a_manifest": a_manifest}


def _run(layout, out, chembl=True, commit=lambda: None, prereg=None):
    return run(TOKEN, layout["hypotheses"], layout["root"], layout["aact"], layout["chembl"] if chembl else None, out,
               prereg=prereg or layout["prereg"], pins=layout["pins"], commit=commit, roots=[("", layout["tmp"])],
               repo_commit=COMMIT)


def test_run_reads_the_inputs_root_a_verified_aact_archive_and_a_verified_chembl(layout, tmp_path):
    (layout["chembl"] / "VERIFIED.json").write_text(json.dumps({"chembl_37.db": layout["chembl_sha"]}))
    out = tmp_path / "out"
    manifest = pd.read_csv(_run(layout, out), sep="\t", dtype=str, keep_default_na=False)
    rows = manifest.set_index("path")
    assert rows.loc["outcomes.csv", "sha256"] == hashlib.sha256((out / "outcomes.csv").read_bytes()).hexdigest()
    assert set(verify_output_dir(out, ["outcomes.csv"])) == set(rows.index)
    inputs = {r.name: r for r in read_inputs(out / "INPUTS.tsv")}
    assert inputs["chembl_37"].sha256 == layout["chembl_sha"] and inputs["chembl_37"].path == "chembl_37/chembl_37.db"
    studies_sha = hashlib.sha256((layout["aact"] / "studies.txt").read_bytes()).hexdigest()
    assert (inputs["aact_studies"].path, inputs["aact_studies"].sha256) == ("aact_20260930/studies.txt", studies_sha)
    assert inputs["hypotheses"].path == "stages/A/output/hypotheses.csv"
    assert all(not Path(r.path).is_absolute() for r in inputs.values())
    info = json.loads((out / "run_info.json").read_text())
    assert info["clinicaltrials_gov_source"].startswith("AACT flat files 2026-09-30")
    assert info["chembl_cross_check"] == "run"
    script = rows["script_sha256"].iloc[0]
    assert set(rows["script_sha256"]) == {script} and script == run_stage_c.stage_code_sha256()
    assert (info["stage"], info["run_token"], info["repo_commit"], info["code_sha256"]) == ("C", TOKEN, COMMIT, script)
    code = inputs["stage_code"]
    assert (code.path, code.sha256, code.source, code.release) == (
        "experiments/08_mechanism_interaction/stages/C", script, "git commit", COMMIT)
    assert {"linkage_audit.csv", "partial_date_boundary.csv", "phase23_diagnostic.csv"} <= set(rows.index)
    assert (rows.loc["linkage_audit.csv", "rows"], rows.loc["partial_date_boundary.csv", "rows"],
            rows.loc["phase23_diagnostic.csv", "rows"]) == ("11", "0", "0")
    audit = pd.read_csv(out / "linkage_audit.csv", dtype=str, keep_default_na=False)
    assert list(audit.columns) == ["step", "count", "description"]
    for name in ("partial_date_boundary.csv", "phase23_diagnostic.csv"):
        assert (out / name).read_text() == "hypothesis_id,quantity,registered,alternative\n"
    assert {"post_freeze_updates.csv", "trials_snapshot/studies_linked.tsv", "INPUTS.tsv"} <= set(rows.index)
    assert rows.loc["trials_snapshot/studies_linked.tsv", "rows"] == "7"
    df = pd.read_csv(out / "outcomes.csv", dtype=str, keep_default_na=False).set_index("hypothesis_id")
    assert dict(df["status_24"]) == {h[0]: h[3] for h in HYPOTHESES}
    diag = pd.read_csv(out / "post_freeze_updates.csv", dtype=str, keep_default_na=False).set_index("hypothesis_id")
    assert list(diag.columns) == ["n_linked_trials", "n_not_in_archive", "n_first_submitted_after_freeze",
                                  "n_last_update_after_freeze"]
    assert diag.loc["h_adv", "n_last_update_after_freeze"] == "1"


def test_run_refuses_chembl_without_verified_json_before_writing(layout, tmp_path):
    with pytest.raises(ChemblNotVerified):
        _run(layout, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_run_refuses_an_unverified_aact_archive_before_writing(layout, tmp_path):
    (layout["aact"] / "studies.txt").write_text((layout["aact"] / "studies.txt").read_text() + "\n")
    with pytest.raises(AactNotVerified):
        _run(layout, tmp_path / "out", chembl=False)
    assert not (tmp_path / "out").exists()


def test_run_without_chembl_records_the_cross_check_not_run(layout, tmp_path):
    out = tmp_path / "out"
    _run(layout, out, chembl=False)
    assert json.loads((out / "run_info.json").read_text())["chembl_cross_check"].startswith("not run")
    assert "chembl_37" not in {r.name for r in read_inputs(out / "INPUTS.tsv")}
    assert set(pd.read_csv(out / "outcomes.csv", dtype=str, keep_default_na=False)["chembl_max_phase_for_ind"]) == {""}


def test_run_commits_once_the_trial_rows_are_written_and_at_the_end(layout, tmp_path):
    out = tmp_path / "out"
    seen = []
    _run(layout, out, chembl=False, commit=lambda: seen.append(((out / "trials_snapshot" / "studies_linked.tsv").exists(),
                                                                (out / "MANIFEST.tsv").exists())))
    assert seen == [(True, False), (True, True)]


@pytest.mark.parametrize("events", [
    (),
    ("Pre-stage note (C, AACT route): stage C reads the 2026-09-30 snapshot",),
    ("<SEAL_A>", "<SEAL_B>", f"stage C run started token={TOKEN}"),
    ("<SEAL_A>", "<SEAL_B>", f"RUN_START stage=B token={TOKEN}"),
    (f"RUN_START stage=C token={TOKEN}",),                                   # no predecessor seal
    ("<SEAL_B>", f"RUN_START stage=C token={TOKEN}"),                        # stage A not sealed
    ("<SEAL_A>", f"RUN_START stage=C token={TOKEN}"),                        # stage B not sealed
    ("<SEAL_A>", f"RUN_START stage=C token={TOKEN}", "<SEAL_B>"),            # stage B sealed after the start
    ("<SEAL_B>", f"RUN_START stage=C token={TOKEN}", "<SEAL_A>"),            # stage A sealed after the start
    (f"SEAL stage=A manifest_sha256={'0' * 64}", "<SEAL_B>", f"RUN_START stage=C token={TOKEN}"),   # wrong hash
    ("<SEAL_A>", f"SEAL stage=A manifest_sha256={'0' * 64}", "<SEAL_B>", f"RUN_START stage=C token={TOKEN}"),  # last seal counts
])
def test_run_refuses_without_the_sealed_chain_and_an_exact_stage_c_run_start_before_writing(layout, tmp_path, events):
    seal_a, seal_b, _ = chain(layout["a_manifest"])
    prereg = prereg_with(tmp_path / "other_PREREG.md",
                         *(e.replace("<SEAL_A>", seal_a).replace("<SEAL_B>", seal_b) for e in events))
    with pytest.raises(RunNotAuthorized):
        _run(layout, tmp_path / "out", chembl=False, prereg=prereg)
    assert not (tmp_path / "out").exists()


def test_run_passes_when_a_wrong_seal_of_a_is_followed_by_the_right_one(layout, tmp_path):
    seal_a, seal_b, start = chain(layout["a_manifest"])
    prereg = prereg_with(tmp_path / "other_PREREG.md", f"SEAL stage=A manifest_sha256={'0' * 64}", seal_a, seal_b, start)
    _run(layout, tmp_path / "out", chembl=False, prereg=prereg)
    assert (tmp_path / "out" / "MANIFEST.tsv").is_file()


def test_run_refuses_a_hypotheses_file_that_is_not_the_one_a_sealed_before_writing(layout, tmp_path):
    h = layout["hypotheses"]
    h.write_text(h.read_text().replace("P1", "P9", 1))
    with pytest.raises(ManifestError, match="hypotheses.csv sha256"):
        _run(layout, tmp_path / "out", chembl=False)
    assert not (tmp_path / "out").exists()


def test_run_refuses_a_resealed_a_manifest_the_log_does_not_name(layout, tmp_path):
    h = layout["hypotheses"]
    h.write_text(h.read_text().replace("P1", "P9", 1))
    write_manifest(h.parent, [h], "c" * 64, [])
    with pytest.raises(RunNotAuthorized, match="stage A"):
        _run(layout, tmp_path / "out", chembl=False)
    assert not (tmp_path / "out").exists()


def test_run_defaults_to_the_frozen_plan_which_logs_no_stage_c_run():
    assert run_stage_c.PREREG == PREREG
    with pytest.raises(RunNotAuthorized, match="RUN_START stage=C"):
        run_stage_c.require_run(PREREG, "C", TOKEN)


def test_the_script_digest_is_over_this_directory_and_the_shared_guard_by_relative_path():
    here = Path(run_stage_c.__file__).resolve().parent
    files = {p.name: p for p in here.glob("*.py") if not p.name.startswith("v8_")}
    assert {"run_stage_c.py", "rules.py", "ctgov.py", "ot_phase.py", "diagnostics.py", "modal_stage_c.py"} <= set(files)
    guard = Path(code_sha256.__code__.co_filename).resolve().parent
    files.update({f"run_guard/{n}": guard / n for n in ("v8_manifest.py", "v8_run_guard.py")})
    assert run_stage_c.stage_code_sha256() == code_sha256(files)


def test_modal_launcher_checks_the_tree_first_and_bakes_the_commit_into_the_image():
    src = (Path(run_stage_c.__file__).resolve().parent / "modal_stage_c.py").read_text()
    main = src.split("def main(run_token: str):")[1]
    assert main.index("clean_commit(HERE.parents[3])") < main.index("require_run(") < main.index("run_remote.remote(")
    assert "image = image.env({COMMIT_ENV: baked_commit(HERE.parents[3])})" in src
    assert "repo_commit=run_commit(None))" in src.split("def run_remote(run_token: str)")[1].split("def tests()")[0]
