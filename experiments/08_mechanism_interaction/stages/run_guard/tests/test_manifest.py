import csv
import hashlib
import json
from pathlib import Path

import pytest

from v8_manifest import (INPUTS_COLUMNS, MANIFEST_COLUMNS, InputRecord, ManifestError, code_record, code_sha256,
                         count_rows, files_under, guard_files, listed_file, portable_path, read_inputs, read_manifest,
                         relative_files, verify_listed, verify_output_dir, write_manifest)
from v8_test_report import run_pytest, write_report


@pytest.fixture
def out(tmp_path) -> Path:
    d = tmp_path / "output"
    (d / "snap").mkdir(parents=True)
    (d / "table.csv").write_text('id,text\n1,"two\nlines"\n2,x\n')
    (d / "coding.tsv").write_text("a\tb\n1\t2\n")
    (d / "snap" / "info.json").write_text("{}")
    return d


INPUTS = [InputRecord("hypotheses", "experiments/08_mechanism_interaction/stages/A/output/hypotheses.csv", "a" * 64),
          InputRecord("studies", "pqtl-v8-inputs:aact_20260930/studies.txt", "b" * 64, source="AACT")]


def test_roundtrip_lists_every_output_and_inputs_tsv(out):
    m = write_manifest(out, [out / "table.csv", out / "coding.tsv", out / "snap" / "info.json"], "c" * 64, INPUTS)
    entries = {e.path: e for e in read_manifest(m)}
    assert set(entries) == {"table.csv", "coding.tsv", "snap/info.json", "INPUTS.tsv"}
    assert (entries["table.csv"].rows, entries["coding.tsv"].rows, entries["snap/info.json"].rows) == (2, 1, None)
    assert entries["INPUTS.tsv"].sha256 == hashlib.sha256((out / "INPUTS.tsv").read_bytes()).hexdigest()
    assert read_inputs(out / "INPUTS.tsv") == INPUTS
    with m.open() as f:
        assert tuple(next(csv.reader(f, delimiter="\t"))) == MANIFEST_COLUMNS
    with (out / "INPUTS.tsv").open() as f:
        assert tuple(next(csv.reader(f, delimiter="\t"))) == INPUTS_COLUMNS
    assert set(verify_output_dir(out, ["table.csv"])) == set(entries)


@pytest.mark.parametrize("bad", ["/Users/x/data.csv", "../data.csv", "pqtl-v8-inputs:/abs/studies.txt",
                                 "pqtl-v8-inputs:../x", ""])
def test_absolute_or_escaping_input_paths_are_refused(out, bad):
    with pytest.raises(ManifestError):
        write_manifest(out, [out / "table.csv"], "c" * 64, [InputRecord("x", bad, "a" * 64)])


def test_portable_path_uses_the_first_root_that_contains_it(tmp_path):
    repo, vol = tmp_path / "repo", tmp_path / "inputs"
    (repo / "data").mkdir(parents=True)
    (vol / "aact").mkdir(parents=True)
    roots = [("", repo), ("pqtl-v8-inputs:", vol), ("experiments/08_mechanism_interaction", tmp_path / "exp")]
    assert portable_path(repo / "data" / "x.csv", roots) == "data/x.csv"
    assert portable_path(vol / "aact" / "studies.txt", roots) == "pqtl-v8-inputs:aact/studies.txt"
    assert portable_path(tmp_path / "exp" / "PREREG.md", roots) == "experiments/08_mechanism_interaction/PREREG.md"
    with pytest.raises(ManifestError):
        portable_path(tmp_path / "elsewhere.csv", roots)


def test_one_changed_byte_in_an_output_or_in_inputs_tsv_is_refused(out):
    write_manifest(out, [out / "table.csv"], "c" * 64, INPUTS)
    (out / "INPUTS.tsv").write_text((out / "INPUTS.tsv").read_text().replace("aaaa", "aaab", 1))
    with pytest.raises(ManifestError, match="INPUTS.tsv sha256"):
        verify_output_dir(out, ["table.csv"])
    write_manifest(out, [out / "table.csv"], "c" * 64, INPUTS)
    (out / "table.csv").write_text((out / "table.csv").read_text().replace("x", "y"))
    with pytest.raises(ManifestError, match="table.csv sha256"):
        verify_output_dir(out, ["table.csv"])


def test_a_required_file_or_inputs_tsv_missing_from_the_manifest_is_refused(out):
    write_manifest(out, [out / "table.csv"], "c" * 64, INPUTS)
    with pytest.raises(ManifestError, match="coding.tsv"):
        verify_output_dir(out, ["table.csv", "coding.tsv"])
    lines = (out / "MANIFEST.tsv").read_text().splitlines()
    (out / "MANIFEST.tsv").write_text("\n".join(x for x in lines if not x.startswith("INPUTS.tsv")) + "\n")
    with pytest.raises(ManifestError, match="INPUTS.tsv"):
        verify_output_dir(out, ["table.csv"])


def test_an_old_style_header_is_refused(out):
    (out / "MANIFEST.tsv").write_text("path\trows\tsha256\tscript_sha256\tinputs_sha256\tutc\nx.csv\t1\ta\tb\tc\td\n")
    with pytest.raises(ManifestError, match="header"):
        read_manifest(out / "MANIFEST.tsv")


def test_count_rows_respects_quoted_newlines(out):
    assert count_rows(out / "table.csv") == 2
    assert count_rows(out / "snap" / "info.json") is None


# ---- hardening: symlinks, escaping paths, malformed hashes -----------------------------------------

def test_a_listed_file_that_is_a_symlink_is_refused(out, tmp_path):
    outside = tmp_path / "outside.csv"
    outside.write_text((out / "table.csv").read_text())
    (out / "table.csv").unlink()
    (out / "table.csv").symlink_to(outside)
    write_manifest(out, [out / "coding.tsv"], "c" * 64, INPUTS)
    lines = (out / "MANIFEST.tsv").read_text().splitlines()
    digest = hashlib.sha256(outside.read_bytes()).hexdigest()
    lines.append(f"table.csv\t2\t{digest}\t{'c' * 64}\t2026-10-01T00:00:00+00:00")
    (out / "MANIFEST.tsv").write_text("\n".join(lines) + "\n")
    with pytest.raises(ManifestError, match="symlink"):
        verify_output_dir(out, ["table.csv"])
    with pytest.raises(ManifestError, match="symlink"):
        verify_listed(out, ["table.csv"])


def test_a_listed_file_under_a_symlinked_directory_is_refused(out, tmp_path):
    write_manifest(out, [out / "snap" / "info.json"], "c" * 64, INPUTS)
    real = tmp_path / "elsewhere"
    (out / "snap").rename(real)
    (out / "snap").symlink_to(real, target_is_directory=True)
    with pytest.raises(ManifestError, match="symlink"):
        verify_output_dir(out, [])


def test_listed_file_requires_the_resolved_path_inside_the_directory(out, tmp_path):
    assert listed_file(out, "snap/info.json") == out / "snap" / "info.json"
    with pytest.raises(ManifestError):
        listed_file(out, "../outside.csv")


@pytest.mark.parametrize("column,value", [("sha256", "A" * 64), ("sha256", "a" * 63), ("sha256", ""),
                                          ("script_sha256", "c" * 65), ("script_sha256", "g" * 64)])
def test_malformed_hashes_in_a_manifest_are_refused_on_reading(out, column, value):
    write_manifest(out, [out / "table.csv"], "c" * 64, INPUTS)
    with (out / "MANIFEST.tsv").open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    rows[0][column] = value
    with (out / "MANIFEST.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(MANIFEST_COLUMNS), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with pytest.raises(ManifestError, match=column):
        read_manifest(out / "MANIFEST.tsv")


def test_a_malformed_input_hash_is_refused_on_writing_and_on_reading(out):
    with pytest.raises(ManifestError, match="64 lowercase hex"):
        write_manifest(out, [out / "table.csv"], "c" * 64, [InputRecord("x", "data/x.csv", "XYZ")])
    with pytest.raises(ManifestError, match="script_sha256"):
        write_manifest(out, [out / "table.csv"], "not-a-hash", INPUTS)
    write_manifest(out, [out / "table.csv"], "c" * 64, INPUTS)
    (out / "INPUTS.tsv").write_text((out / "INPUTS.tsv").read_text().replace("a" * 64, "A" * 64))
    with pytest.raises(ManifestError, match="64 lowercase hex"):
        read_inputs(out / "INPUTS.tsv")


def test_verify_listed_checks_only_the_named_files(out):
    write_manifest(out, [out / "table.csv", out / "coding.tsv"], "c" * 64, INPUTS)
    (out / "coding.tsv").unlink()
    assert verify_listed(out, ["table.csv"]) == {"table.csv": hashlib.sha256((out / "table.csv").read_bytes()).hexdigest()}
    with pytest.raises(ManifestError, match="missing"):
        verify_listed(out, ["coding.tsv"])
    with pytest.raises(ManifestError, match="does not list"):
        verify_listed(out, ["funnel.csv"])
    (out / "table.csv").write_text((out / "table.csv").read_text().replace("x", "y"))
    with pytest.raises(ManifestError, match="table.csv sha256"):
        verify_listed(out, ["table.csv"])


# ---- reader: one entry per path, one script digest, UTC times, well-formed rows ----------------------

ALL = ("table.csv", "coding.tsv", "snap/info.json")


def sealed(out: Path) -> Path:
    return write_manifest(out, [out / n for n in ALL], "c" * 64, INPUTS)


def rewrite(path: Path, columns, mutate) -> None:
    """Rewrite a MANIFEST.tsv / INPUTS.tsv after `mutate(rows)` changed its list of row dicts in place."""
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    mutate(rows)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def test_a_path_listed_twice_is_refused_on_reading(out):
    m = sealed(out)
    rewrite(m, MANIFEST_COLUMNS, lambda rows: rows.append(dict(rows[0])))
    with pytest.raises(ManifestError, match=r"lists \['table.csv'\] more than once"):
        read_manifest(m)
    with pytest.raises(ManifestError, match="more than once"):
        verify_output_dir(out, ["table.csv"])
    with pytest.raises(ManifestError, match="more than once"):
        verify_listed(out, ["coding.tsv"])


def test_a_manifest_with_two_script_digests_is_refused(out):
    m = sealed(out)
    assert len({e.script_sha256 for e in read_manifest(m)}) == 1
    rewrite(m, MANIFEST_COLUMNS, lambda rows: rows[1].update(script_sha256="d" * 64))
    with pytest.raises(ManifestError, match="2 different script_sha256"):
        read_manifest(m)


def test_a_manifest_that_lists_itself_is_refused(out):
    m = sealed(out)
    digest = hashlib.sha256(m.read_bytes()).hexdigest()
    rewrite(m, MANIFEST_COLUMNS, lambda rows: rows.append({**rows[0], "path": "MANIFEST.tsv", "rows": "4", "sha256": digest}))
    with pytest.raises(ManifestError, match="lists itself"):
        read_manifest(m)
    with pytest.raises(ManifestError, match="cannot list itself"):
        write_manifest(out, [out / "table.csv", out / "MANIFEST.tsv"], "c" * 64, INPUTS)


@pytest.mark.parametrize("value", ["", "yesterday", "2026-10-01", "2026-10-01T00:00:00", "2026-10-01T02:00:00+02:00",
                                   "2026-13-01T00:00:00+00:00", "1790000000"])
def test_a_time_that_is_not_utc_is_refused(out, value):
    m = sealed(out)
    rewrite(m, MANIFEST_COLUMNS, lambda rows: rows[2].update(utc_time=value))
    with pytest.raises(ManifestError, match="utc_time"):
        read_manifest(m)


@pytest.mark.parametrize("value", ["2026-10-01T00:00:00+00:00", "2026-10-01T00:00:00Z", "2026-10-01T00:00:00.250000+00:00"])
def test_a_utc_time_is_accepted(out, value):
    m = sealed(out)
    rewrite(m, MANIFEST_COLUMNS, lambda rows: [r.update(utc_time=value) for r in rows])
    assert {e.utc_time for e in read_manifest(m)} == {value}
    assert set(verify_output_dir(out, list(ALL))) == {*ALL, "INPUTS.tsv"}


@pytest.mark.parametrize("value", ["-1", "two", "1.0", " 2"])
def test_a_row_count_that_is_not_a_non_negative_integer_is_refused(out, value):
    m = sealed(out)
    rewrite(m, MANIFEST_COLUMNS, lambda rows: rows[0].update(rows=value))
    with pytest.raises(ManifestError, match="rows"):
        read_manifest(m)


@pytest.mark.parametrize("name,columns", [("MANIFEST.tsv", MANIFEST_COLUMNS), ("INPUTS.tsv", INPUTS_COLUMNS)])
@pytest.mark.parametrize("damage", [lambda line: line + "\textra", lambda line: line.rsplit("\t", 1)[0]])
def test_a_row_with_the_wrong_number_of_fields_is_refused(out, name, columns, damage):
    sealed(out)
    lines = (out / name).read_text().splitlines()
    (out / name).write_text("\n".join([lines[0], damage(lines[1]), *lines[2:]]) + "\n")
    read = read_manifest if name == "MANIFEST.tsv" else read_inputs
    with pytest.raises(ManifestError, match=f"line 2 does not have {len(columns)} tab-separated fields"):
        read(out / name)


# ---- reader: INPUTS.tsv -----------------------------------------------------------------------------

def test_input_names_are_unique_and_non_empty_on_reading(out):
    sealed(out)
    path = out / "INPUTS.tsv"
    good = path.read_bytes()
    assert [r.name for r in read_inputs(path)] == ["hypotheses", "studies"]
    rewrite(path, INPUTS_COLUMNS, lambda rows: rows[1].update(name="hypotheses"))
    with pytest.raises(ManifestError, match=r"input names \['hypotheses'\] appear more than once"):
        read_inputs(path)
    path.write_bytes(good)
    rewrite(path, INPUTS_COLUMNS, lambda rows: rows[0].update(name=""))
    with pytest.raises(ManifestError, match="has no name"):
        read_inputs(path)
    with pytest.raises(ManifestError, match="duplicate input names"):
        write_manifest(out, [out / "table.csv"], "c" * 64, [INPUTS[0], INPUTS[0]])


def test_an_inputs_tsv_that_does_not_parse_is_refused_although_its_hash_matches_the_manifest(out):
    sealed(out)
    (out / "INPUTS.tsv").write_text("name\tpath\nx\ty\n")
    digest = hashlib.sha256((out / "INPUTS.tsv").read_bytes()).hexdigest()
    rewrite(out / "MANIFEST.tsv", MANIFEST_COLUMNS,
            lambda rows: [r.update(sha256=digest, rows="1") for r in rows if r["path"] == "INPUTS.tsv"])
    with pytest.raises(ManifestError, match="INPUTS.tsv: header"):
        verify_output_dir(out, list(ALL))
    (out / "INPUTS.tsv").unlink()
    with pytest.raises(ManifestError, match="does not exist"):
        read_inputs(out / "INPUTS.tsv")


# ---- reader: MANIFEST.tsv and INPUTS.tsv are themselves plain files inside the directory -----------

@pytest.mark.parametrize("name", ["MANIFEST.tsv", "INPUTS.tsv"])
def test_a_manifest_or_inputs_file_that_is_a_symlink_is_refused(out, tmp_path, name):
    sealed(out)
    outside = tmp_path / f"outside_{name}"
    (out / name).rename(outside)
    (out / name).symlink_to(outside)
    read = read_manifest if name == "MANIFEST.tsv" else read_inputs
    with pytest.raises(ManifestError, match="symlink"):
        read(out / name)
    with pytest.raises(ManifestError, match="symlink"):
        verify_output_dir(out, list(ALL))
    with pytest.raises(ManifestError, match="symlink"):
        verify_listed(out, ["INPUTS.tsv"])


def test_an_output_directory_that_is_itself_a_symlink_is_refused_on_verifying_and_on_writing(out, tmp_path):
    sealed(out)
    alias = tmp_path / "alias"
    alias.symlink_to(out, target_is_directory=True)
    assert all((alias / n).read_bytes() == (out / n).read_bytes() for n in ("MANIFEST.tsv", *ALL))
    with pytest.raises(ManifestError, match="alias is a symlink"):
        verify_output_dir(alias, list(ALL))
    assert set(verify_output_dir(out, list(ALL))) == {*ALL, "INPUTS.tsv"}
    before = {p.name: p.read_bytes() for p in out.iterdir() if p.is_file()}
    with pytest.raises(ManifestError, match="alias is a symlink"):
        write_manifest(alias, [alias / n for n in ALL], "d" * 64, [])
    assert {p.name: p.read_bytes() for p in out.iterdir() if p.is_file()} == before


def test_a_symlinked_output_directory_is_refused_before_a_first_manifest_is_written(out, tmp_path):
    alias = tmp_path / "alias"
    alias.symlink_to(out, target_is_directory=True)
    with pytest.raises(ManifestError, match="alias is a symlink"):
        write_manifest(alias, [alias / n for n in ALL], "c" * 64, INPUTS)
    assert not (out / "INPUTS.tsv").exists() and not (out / "MANIFEST.tsv").exists()


# ---- reader: nothing unlisted in a sealed output directory ------------------------------------------

def test_files_under_lists_files_and_symlinks_without_following_them(out, tmp_path):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (elsewhere / "big.bin").write_text("x")
    (out / "linkdir").symlink_to(elsewhere, target_is_directory=True)
    (out / "snap" / "link.json").symlink_to(elsewhere / "big.bin")
    assert files_under(out) == ["coding.tsv", "linkdir", "snap/info.json", "snap/link.json", "table.csv"]


@pytest.mark.parametrize("extra", ["notes.txt", "snap/extra.json", "deep/er/file.csv", ".hidden"])
def test_an_unlisted_file_in_a_sealed_directory_is_refused(out, extra):
    sealed(out)
    assert set(verify_output_dir(out, list(ALL))) == {*ALL, "INPUTS.tsv"}
    (out / extra).parent.mkdir(parents=True, exist_ok=True)
    (out / extra).write_text("side file")
    with pytest.raises(ManifestError, match="does not list") as err:
        verify_output_dir(out, list(ALL))
    assert repr([extra]) in str(err.value)
    assert set(verify_output_dir(out, list(ALL), side_files=[extra])) == {*ALL, "INPUTS.tsv"}
    with pytest.raises(ManifestError, match="does not list"):
        verify_output_dir(out, list(ALL), side_files=["another.txt"])


def test_an_output_left_out_of_the_manifest_is_refused(out):
    write_manifest(out, [out / "table.csv"], "c" * 64, INPUTS)
    with pytest.raises(ManifestError, match=r"does not list: \['coding.tsv', 'snap/info.json'\]"):
        verify_output_dir(out, ["table.csv"])
    assert set(verify_listed(out, ["table.csv"])) == {"table.csv"}      # a partial copy is checked file by file


def test_an_unlisted_symlink_in_a_sealed_directory_is_refused(out, tmp_path):
    sealed(out)
    (tmp_path / "elsewhere").mkdir()
    (out / "more").symlink_to(tmp_path / "elsewhere", target_is_directory=True)
    with pytest.raises(ManifestError, match=r"does not list: \['more'\]"):
        verify_output_dir(out, list(ALL))


@pytest.mark.parametrize("bad", ["/abs/side.txt", "../side.txt"])
def test_a_side_file_name_must_be_relative(out, bad):
    sealed(out)
    with pytest.raises(ManifestError, match="must be relative"):
        verify_output_dir(out, list(ALL), side_files=[bad])


# ---- code digest ------------------------------------------------------------------------------------

def framed(*items: tuple[str, bytes]) -> str:
    """The digest written out by hand: per file, 8-byte path length, path, 8-byte content length, content."""
    h = hashlib.sha256()
    for rel, data in items:
        h.update(len(rel.encode()).to_bytes(8, "big") + rel.encode() + len(data).to_bytes(8, "big") + data)
    return h.hexdigest()


@pytest.fixture
def code(tmp_path) -> Path:
    root = tmp_path / "stage"
    (root / "pkg").mkdir(parents=True)
    (root / "run.py").write_bytes(b"import pkg\n")
    (root / "pkg" / "rules.py").write_bytes(b"X = 1\n")
    (root / "pkg" / "__init__.py").write_bytes(b"")
    return root


def test_code_digest_frames_each_relative_path_and_its_bytes_in_path_order(code):
    files = relative_files(code, [code / "run.py", code / "pkg" / "rules.py", code / "pkg" / "__init__.py"])
    assert list(files) == ["run.py", "pkg/rules.py", "pkg/__init__.py"]
    expected = framed(("pkg/__init__.py", b""), ("pkg/rules.py", b"X = 1\n"), ("run.py", b"import pkg\n"))
    assert code_sha256(files) == expected
    assert code_sha256(dict(reversed(list(files.items())))) == expected          # the order given does not matter


def test_code_digest_depends_on_the_directory_a_file_sits_in_not_only_its_name(code, tmp_path):
    flat = tmp_path / "flat"
    flat.mkdir()
    for name, src in (("run.py", code / "run.py"), ("rules.py", code / "pkg" / "rules.py")):
        (flat / name).write_bytes(src.read_bytes())
    nested = relative_files(code, [code / "run.py", code / "pkg" / "rules.py"])
    moved = relative_files(flat, [flat / "run.py", flat / "rules.py"])
    assert [p.name for p in nested.values()] == [p.name for p in moved.values()]   # same basenames, same bytes
    assert [p.read_bytes() for p in nested.values()] == [p.read_bytes() for p in moved.values()]
    assert code_sha256(nested) != code_sha256(moved)


def test_code_digest_keeps_two_files_with_one_basename_apart(code):
    (code / "other").mkdir()
    (code / "other" / "rules.py").write_bytes(b"X = 2\n")
    both = relative_files(code, [code / "pkg" / "rules.py", code / "other" / "rules.py"])
    before = code_sha256(both)
    swapped = {"pkg/rules.py": both["other/rules.py"], "other/rules.py": both["pkg/rules.py"]}
    assert code_sha256(swapped) != before
    (code / "other" / "rules.py").write_bytes(b"X = 3\n")
    assert code_sha256(both) != before


def test_code_digest_does_not_confuse_a_path_with_the_bytes_that_follow_it(tmp_path):
    for name, data in (("ab", b"c"), ("a", b"bc"), ("abc", b""), ("a", b""), ("b", b"c")):
        d = tmp_path / f"{name}_{len(data)}"
        d.mkdir()
        (d / name).write_bytes(data)
    one, two, three = ({n: tmp_path / f"{n}_{k}" / n} for n, k in (("ab", 1), ("a", 2), ("abc", 0)))
    split = {"a": tmp_path / "a_0" / "a", "b": tmp_path / "b_1" / "b"}
    assert len({code_sha256(x) for x in (one, two, three, split)}) == 4


def test_code_digest_is_the_same_wherever_the_tree_is_checked_out(code, tmp_path):
    names = ["run.py", "pkg/rules.py", "pkg/__init__.py"]
    copy = tmp_path / "another" / "place"
    for n in names:
        (copy / n).parent.mkdir(parents=True, exist_ok=True)
        (copy / n).write_bytes((code / n).read_bytes())
    assert code_sha256(relative_files(code, [code / n for n in names])) == \
        code_sha256(relative_files(copy, [copy / n for n in names]))


def test_relative_files_refuses_a_file_outside_the_root_or_listed_twice(code, tmp_path):
    (tmp_path / "outside.py").write_text("x")
    with pytest.raises(ManifestError, match="is not under"):
        relative_files(code, [code / "run.py", tmp_path / "outside.py"])
    with pytest.raises(ManifestError, match="listed twice"):
        relative_files(code, [code / "run.py", code / "pkg" / ".." / "run.py"])


@pytest.mark.parametrize("key", ["/abs/run.py", "../run.py", "", "./run.py", "pkg//rules.py", "pkg/"])
def test_code_digest_refuses_a_key_that_is_not_a_canonical_relative_path(code, key):
    with pytest.raises(ManifestError):
        code_sha256({key: code / "run.py"})


def test_guard_files_are_the_shared_modules_under_their_repository_names():
    files = guard_files()
    assert set(files) == {"run_guard/v8_manifest.py", "run_guard/v8_run_guard.py"}
    assert all(p.is_file() and p.name == key.split("/")[1] for key, p in files.items())
    assert "def code_sha256" in files["run_guard/v8_manifest.py"].read_text()


def test_code_record_carries_the_digest_and_the_commit(out):
    rec = code_record("C", "d" * 64, "0123456789abcdef0123456789abcdef01234567")
    assert (rec.name, rec.path, rec.sha256, rec.source, rec.release) == (
        "stage_code", "experiments/08_mechanism_interaction/stages/C", "d" * 64, "git commit",
        "0123456789abcdef0123456789abcdef01234567")
    write_manifest(out, [out / "table.csv", out / "coding.tsv", out / "snap" / "info.json"], "d" * 64, [*INPUTS, rec])
    assert read_inputs(out / "INPUTS.tsv")[-1] == rec
    for bad in ("", "0123456", "g" * 40, "A" * 40, "0" * 64):
        with pytest.raises(ManifestError, match="commit"):
            code_record("C", "d" * 64, bad)
    with pytest.raises(ManifestError, match="64 lowercase hex"):
        code_record("C", "short", "0" * 40)


# ---- test reports ------------------------------------------------------------------------------------

def test_write_report_archives_the_previous_report_instead_of_overwriting_it(tmp_path):
    path = tmp_path / "tests" / "modal_tests.json"
    write_report(path, {"summary": "1 passed"})
    first = path.read_bytes()
    write_report(path, {"summary": "2 passed"})
    (old,) = (tmp_path / "tests" / "superseded").iterdir()
    assert old.read_bytes() == first and old.name.startswith("modal_tests_") and old.suffix == ".json"
    assert json.loads(path.read_text()) == {"summary": "2 passed"}


def test_run_pytest_reports_the_counts_of_the_summary_line(tmp_path):
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "test_x.py").write_text(
        "import pytest\n\ndef test_a():\n    assert True\n\ndef test_b():\n    assert False\n\n"
        "@pytest.mark.skip\ndef test_c():\n    pass\n")
    report = run_pytest(tmp_path, [tmp_path])
    assert report["counts"] == {"failed": 1, "passed": 1, "skipped": 1} and report["returncode"] == 1
