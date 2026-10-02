import csv
import hashlib

import pytest

from stage_d.guard import StageInputError, sha256_file, verify_inputs
from stage_d.pipeline import prepare
from stage_d.synthetic import synthetic_fingerprint


def _rewrite_manifest(path, mutate):
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
        header = list(rows[0])
    mutate(rows)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return sha256_file(path)


def test_sealed_stages_pass(sealed_stages):
    root, expected = sealed_stages
    seals = verify_inputs(root, expected)
    assert set(seals) == {"A", "B", "C"}
    assert seals["A"].files["hypotheses.csv"] == hashlib.sha256((root / "A/output/hypotheses.csv").read_bytes()).hexdigest()


@pytest.mark.parametrize("stage", ["A", "B", "C"])
def test_missing_manifest_refuses(sealed_stages, stage):
    root, expected = sealed_stages
    (root / stage / "output" / "MANIFEST.tsv").unlink()
    with pytest.raises(StageInputError, match="not sealed"):
        verify_inputs(root, expected)


@pytest.mark.parametrize("stage,fname", [("A", "hypotheses.csv"), ("B", "evidence.csv"), ("C", "outcomes.csv")])
def test_edited_output_refuses(sealed_stages, stage, fname):
    root, expected = sealed_stages
    path = root / stage / "output" / fname
    path.write_text(path.read_text().replace("True", "False", 1))
    with pytest.raises(StageInputError, match="sha256"):
        verify_inputs(root, expected)


def test_manifest_not_matching_logged_hash_refuses(sealed_stages):
    root, expected = sealed_stages
    with pytest.raises(StageInputError, match="logged"):
        verify_inputs(root, {**expected, "B": "0" * 64})


def test_resealed_manifest_still_refused_without_new_logged_hash(sealed_stages):
    root, expected = sealed_stages
    out = root / "C" / "output"
    (out / "outcomes.csv").write_text((out / "outcomes.csv").read_text().replace("advanced", "active", 1))
    new = sha256_file(out / "outcomes.csv")
    _rewrite_manifest(out / "MANIFEST.tsv", lambda rows: [r.update(sha256=new) for r in rows if r["path"] == "outcomes.csv"])
    with pytest.raises(StageInputError, match="logged"):
        verify_inputs(root, expected)


def test_required_file_absent_from_manifest_refuses(sealed_stages):
    root, expected = sealed_stages
    new = _rewrite_manifest(root / "A" / "output" / "MANIFEST.tsv",
                            lambda rows: rows.__setitem__(slice(None), [r for r in rows if r["path"] != "funnel.csv"]))
    with pytest.raises(StageInputError, match="does not list"):
        verify_inputs(root, {**expected, "A": new})


def test_row_count_mismatch_refuses(sealed_stages):
    root, expected = sealed_stages
    new = _rewrite_manifest(root / "B" / "output" / "MANIFEST.tsv",
                            lambda rows: [r.update(rows=str(int(r["rows"]) + 1)) for r in rows])
    with pytest.raises(StageInputError, match="rows"):
        verify_inputs(root, {**expected, "B": new})


def test_broken_chain_refuses(sealed_stages):
    root, expected = sealed_stages
    out = root / "C" / "output"
    a_hyp = sha256_file(root / "A" / "output" / "hypotheses.csv")
    (out / "INPUTS.tsv").write_text((out / "INPUTS.tsv").read_text().replace(a_hyp, "f" * 64))
    new_inputs = sha256_file(out / "INPUTS.tsv")
    new = _rewrite_manifest(out / "MANIFEST.tsv",
                            lambda rows: [r.update(sha256=new_inputs) for r in rows if r["path"] == "INPUTS.tsv"])
    with pytest.raises(StageInputError, match="does not name A"):
        verify_inputs(root, {**expected, "C": new})


def test_edited_inputs_tsv_refuses_although_d_never_reads_the_inputs_it_names(sealed_stages):
    root, expected = sealed_stages
    out = root / "B" / "output"
    (out / "INPUTS.tsv").write_text((out / "INPUTS.tsv").read_text() + "x\tother.csv\t" + "e" * 64 + "\t\t\t\t\n")
    with pytest.raises(StageInputError, match="INPUTS.tsv sha256"):
        verify_inputs(root, expected)


def test_prepare_refuses_before_writing_anything(sealed_stages, tmp_path):
    root, expected = sealed_stages
    work = tmp_path / "work"
    with pytest.raises(StageInputError):
        prepare(root, {**expected, "A": "0" * 64}, work, synthetic_fingerprint({**expected, "A": "0" * 64}))
    assert not work.exists()


# ---- nothing unlisted in a sealed stage directory ---------------------------------------------------

@pytest.mark.parametrize("stage,name", [("A", "notes.txt"), ("A", "unit_plan.json"), ("A", ".gitignore"),
                                        ("B", "evidence_old.csv"), ("B", ".gitignore"), ("B", "units/u1/result.json"),
                                        ("C", "unit_plan.json"), ("C", "trials_snapshot/extra.json"), ("C", ".DS_Store")])
def test_an_unlisted_file_in_a_sealed_stage_directory_refuses(sealed_stages, stage, name):
    root, expected = sealed_stages
    path = root / stage / "output" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("not in the manifest")
    with pytest.raises(StageInputError, match=f"stage {stage}: .*does not list: .*{name.split('/')[-1]}"):
        verify_inputs(root, expected)


@pytest.mark.parametrize("stage,name", [("B", "unit_plan.json"), ("C", ".gitignore")])
def test_the_named_side_files_of_a_stage_are_accepted(sealed_stages, stage, name):
    root, expected = sealed_stages
    (root / stage / "output" / name).write_text("{}")
    assert set(verify_inputs(root, expected)) == {"A", "B", "C"}
    assert name not in verify_inputs(root, expected)[stage].files


def test_a_manifest_listing_a_path_twice_or_itself_refuses(sealed_stages):
    root, expected = sealed_stages
    manifest = root / "C" / "output" / "MANIFEST.tsv"
    twice = _rewrite_manifest(manifest, lambda rows: rows.append(dict(rows[0])))
    with pytest.raises(StageInputError, match="more than once"):
        verify_inputs(root, {**expected, "C": twice})
    itself = _rewrite_manifest(manifest, lambda rows: rows.__setitem__(-1, {**rows[0], "path": "MANIFEST.tsv"}))
    with pytest.raises(StageInputError, match="lists itself"):
        verify_inputs(root, {**expected, "C": itself})
