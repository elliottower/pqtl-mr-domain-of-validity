"""The pre-analysis validation of GWAS Catalog files read without a harmonised copy (stage_b/validate.py),
on collected files written to a temporary volume root."""
import gzip
import hashlib
import json
import math
import os
import shutil

import pytest
from test_outcome_files import (MAP_HEADER, SSF_BETA_HEADER, SSF_OR_CI_HEADER, ci, exact_or_ci_row, laskar_row,
                                ssf_row, tsv)

from stage_b.author_formats import AUTHOR_FORMATS, LASKAR_HEADER, header_sha256
from stage_b.checkpoint import collect_entry
from stage_b.collect import read_record, record_path
from stage_b.fetch import VolumeFetcher
from stage_b.outcome_files import REJECT_REASONS
from stage_b.schemas import CollectError, CollectRecord, CollectTask, InputContractError, OutcomeSpec, SourceUnreadable
from stage_b.status import volume_state
from stage_b.validate import (REPORT_NAME, RULES, VALIDATION_DETAIL, VALIDATION_DIR, apply_validation, check_classification,
                              classification, failed_only_for_standard_error, gwas_catalog_records, in_scope,
                              require_validation, restorable, restore_validated, result_path, result_sha256, validate_file,
                              validation_report)


def put(root, key: str, data: bytes, name: str | None = None, layout: str = "gwas_ssf", rsid_rule: str = "column",
        build: str = "GRCh38", source: str = "gwas_catalog") -> CollectRecord:
    """A collected file and its record on the volume root."""
    name = name or f"{key}.tsv"
    raw = root / "raw" / source / name
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(data)
    record = CollectRecord(status="collected", source=source, key=key, name=name, path=raw.relative_to(root).as_posix(),
                           bytes=len(data), sha256=hashlib.sha256(data).hexdigest(), layout=layout, build=build,
                           rsid_rule=rsid_rule)
    record_path(root, source, key).parent.mkdir(parents=True, exist_ok=True)
    record_path(root, source, key).write_text(record.model_dump_json())
    return record


def mixed_rows() -> list[list[str]]:
    good = [ssf_row(pos=str(1000 + i), rsid=f"rs{i}") for i in range(10)]
    good.append(ssf_row(pos="2000", rsid="rs3"))                     # rs3 on two rows
    good.append(ssf_row(pos="2001", rsid="#NA"))                     # read, without an rsID
    bad = [ssf_row(se="0"), ssf_row(ea="N"), ssf_row(ea="A", oa="A"), ssf_row(p="1.5"), ssf_row()[:-1],
           ssf_row(effect="x"), ssf_row(chrom="")]
    return good + bad


def test_a_file_gives_aggregate_counts_by_reason_and_passes(tmp_path):
    record = put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, mixed_rows()))
    got = validate_file(tmp_path, record)
    assert (got["rows"], got["rows_kept"], got["row_width_errors"], got["parse_errors"]) == (19, 12, 1, 1)
    assert got["rejected"] == {**dict.fromkeys(REJECT_REASONS, 0), "se_nonpositive_or_nonfinite": 1, "invalid_allele": 1,
                               "equal_alleles": 1, "p_outside_0_1": 1, "row_width": 1, "unparseable_number": 1,
                               "no_position": 1}
    assert (got["rows_with_rsid"], got["rsid_recovery_rate"], got["duplicate_rsids"], got["rows_with_duplicate_rsid"]) == (
        11, pytest.approx(11 / 12), 1, 2)
    assert (got["header"], got["header_agrees"], got["passed"], got["reasons"], got["integrity"], got["ci_p_check"]) == (
        SSF_BETA_HEADER, True, True, [], None, None)
    assert (got["uncertainty_mode"], got["se_census"]) == ("native_se", {"rows": 19, "missing": 0, "present": 15,
                                                                         "not_counted": 4})
    assert json.loads(result_path(tmp_path, "GCST1").read_text()) == got
    bound = apply_validation(tmp_path, record, got)
    assert bound == record.model_copy(update={"uncertainty_mode": "native_se", "header_sha256": header_sha256(SSF_BETA_HEADER),
                                              "validation_sha256": result_sha256(got)}) == read_record(tmp_path, "gwas_catalog", "GCST1")
    assert apply_validation(tmp_path, bound, got) == bound                     # already bound: nothing moved again
    assert [CollectRecord.model_validate_json(p.read_text()) for p in
            (tmp_path / "superseded" / "collect" / "gwas_catalog").glob("GCST1__*.json")] == [record]
    text = json.dumps(got)
    assert "rs3" not in text and "0.01" not in text                             # no per-variant value
    assert not (tmp_path / "validation" / "work").exists() or not any((tmp_path / "validation" / "work").iterdir())


@pytest.mark.parametrize("good,bad,passed", [(2, 2, True), (1, 2, False), (0, 0, False)])
def test_a_file_passes_only_when_at_least_half_of_its_rows_are_read(tmp_path, good, bad, passed):
    rows = [ssf_row(pos=str(1000 + i)) for i in range(good)] + [ssf_row(ea="additional copy") for _ in range(bad)]
    got = validate_file(tmp_path, put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, rows)))
    assert got["passed"] is passed
    if good == 0:                                     # no row to take a census of: no mode, before any row is read
        assert got["reasons"] == ["uncertainty mode: no data row with a position and a valid pair of alleles to choose "
                                  "the uncertainty mode from"] and "rows" not in got
    elif not passed:
        assert got["reasons"] == [f"{good} of {good + bad} data rows read, fewer than 50%"]


def test_a_failing_file_is_made_unreadable_its_collected_record_moved_and_the_analysis_reads_it_as_unavailable(tmp_path):
    rows = [ssf_row(pos=str(1000 + i), ea="additional copy") for i in range(9)] + [ssf_row()]
    record = put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, rows))
    got = validate_file(tmp_path, record)
    new = apply_validation(tmp_path, record, got)
    assert new.status == "unreadable" and new.detail == VALIDATION_DETAIL + "1 of 10 data rows read, fewer than 50%"
    assert (new.name, new.path, new.bytes, new.sha256, new.layout) == (record.name, record.path, record.bytes, record.sha256,
                                                                         "gwas_ssf")
    assert read_record(tmp_path, "gwas_catalog", "GCST1") == new
    moved = list((tmp_path / "superseded" / "collect" / "gwas_catalog").glob("GCST1__*.json"))
    assert [CollectRecord.model_validate_json(p.read_text()) for p in moved] == [record]
    entry = collect_entry(CollectTask(source="gwas_catalog", key="GCST1"), new)
    assert (entry["status"], entry["sha256"], entry["reason"]) == ("unreadable", record.sha256, new.detail)
    fetcher = VolumeFetcher(tmp_path, "", tmp_path / "a", tmp_path / "e")
    with pytest.raises(SourceUnreadable, match="fewer than 50%"):
        fetcher.outcome_region(OutcomeSpec(accession="GCST1", source="gwas_catalog", n_case=1, n_control=1, risk_coded=True),
                               "1", 1000, 100)
    assert in_scope(new) and validate_file(tmp_path, new) == got                 # the result is kept with the record
    assert apply_validation(tmp_path, new, got) == new


@pytest.mark.parametrize("header,rule,reason", [
    (SSF_BETA_HEADER, "ukbppp_map", "the header gives the rsID rule column, the record ukbppp_map"),
    ([c for c in SSF_BETA_HEADER if c != "standard_error"], "column", "the header has no standard_error column"),
])
def test_a_header_that_disagrees_with_the_column_map_of_its_record_fails(tmp_path, header, rule, reason):
    rows = [ssf_row()[:len(header)] for _ in range(3)]
    got = validate_file(tmp_path, put(tmp_path, "GCST1", tsv(header, rows), rsid_rule=rule))
    assert (got["header_agrees"], got["passed"]) == (False, False) and reason in got["reasons"][0]
    assert "rows" not in got


def interrupted_then_resumed(tmp_path, data: bytes, header_rows: int = 4) -> list[tuple[dict, dict]]:
    """For every commit k of an uninterrupted validation: the call killed at commit k, then resumed;
    (resumed result, uninterrupted result), without their times. Where the kill falls in the read
    pass, a row written after the checkpoint is left in the shard, as a kill would leave it."""
    commits = []
    whole = validate_file(tmp_path / "whole", put(tmp_path / "whole", "GCST1", data), lambda: commits.append(1),
                          checkpoint_rows=header_rows)
    out = []
    for k in range(1, len(commits)):
        root, count = tmp_path / f"k{k}", []

        def dies_at_commit_k(k=k, count=count):
            count.append(1)
            if len(count) == k:
                raise KeyboardInterrupt

        record = put(root, "GCST1", data)
        with pytest.raises(KeyboardInterrupt):
            validate_file(root, record, dies_at_commit_k, checkpoint_rows=header_rows)
        shard = list((root / "validation" / "work").glob("*/variants.tsv"))
        if shard and not result_path(root, "GCST1").exists():
            with shard[0].open("ab") as fh:
                fh.write(b"rs999\t1\t5\tA\tG\n")
        resumed = validate_file(root, record, checkpoint_rows=header_rows)
        out.append(tuple({key: v for key, v in r.items() if key != "utc"} for r in (resumed, whole)))
    return out


def test_an_interrupted_validation_resumes_in_either_pass_and_gives_the_uninterrupted_result(tmp_path):
    pairs = interrupted_then_resumed(tmp_path, tsv(SSF_BETA_HEADER, mixed_rows() * 3))
    assert len(pairs) >= 12                                   # census and read pass each checkpoint every 4 of 57 rows
    for resumed, whole in pairs:
        assert resumed == whole


def test_the_ukbppp_map_rule_reports_rsid_recovery_through_the_map_of_each_chromosome(tmp_path):
    header = [c for c in SSF_BETA_HEADER if c != "rsid"]
    rows = [ssf_row(pos=str(1000 + k))[:-1] for k in range(4)] + [ssf_row(chrom="2", pos="7")[:-1]]
    record = put(tmp_path, "GCST1", tsv(header, rows), rsid_rule="ukbppp_map")
    map_rows = [["x", "A", "G", "rs100", 9, 1000], ["x", "G", "A", "rs101", 9, 1001],
                ["x", "A", "G", "rs102", 9, 1002], ["x", "A", "G", "rs902", 9, 1002]]   # two rsIDs at 1002: none taken
    put(tmp_path, "1", gzip.compress(tsv(MAP_HEADER, map_rows), mtime=0), name="map_chr1.tsv.gz", layout="",
        rsid_rule="", source="ukbppp_rsid_map")
    got = validate_file(tmp_path, record)
    assert (got["rows_kept"], got["rows_with_rsid"], got["rsid_recovery_rate"]) == (5, 2, pytest.approx(0.4))
    assert got["rows_on_a_chromosome_without_a_collected_rsid_map"] == 1 and got["passed"] is True


def laskar_file(rows) -> bytes:
    return gzip.compress(tsv(list(LASKAR_HEADER), rows, eol="\r\n"), mtime=0)


@pytest.mark.parametrize("se_for_ci,passed", [(None, True), (0.09, False)])
def test_the_laskar_file_passes_only_under_the_or_ci_se_integrity_rule(tmp_path, se_for_ci, passed):
    rows = [ci(1.0 + 0.01 * i, 0.05, se_for_ci=se_for_ci) for i in range(200)] + [laskar_row(odds="NA")]
    record = put(tmp_path, "GCST008226", laskar_file(rows), name=AUTHOR_FORMATS["GCST008226"].data_file, layout="author",
                 build="GRCh37")
    got = validate_file(tmp_path, record)
    integrity = got["integrity"]
    assert (integrity["rows_checked"], integrity["rows_not_checkable"], integrity["passed"]) == (200, 1, passed)
    assert got["passed"] is passed and got["rows_kept"] == 200
    new = apply_validation(tmp_path, record, got)
    assert new.status == ("collected" if passed else "unreadable")
    if not passed:
        assert new.detail == VALIDATION_DETAIL + "OR/CI/SE integrity: 0 of 200 checked rows agree, fewer than 99%"


def test_the_report_holds_aggregate_diagnostics_and_lists_what_was_not_validated(tmp_path):
    ok = put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, mixed_rows()))
    bad = put(tmp_path, "GCST2", tsv(SSF_BETA_HEADER, [ssf_row(ea="additional copy")] * 3))
    results = {"GCST1": validate_file(tmp_path, ok), "GCST2": validate_file(tmp_path, bad), "GCST3": "OSError: gone"}
    indexed = CollectRecord(status="remote_indexed", source="gwas_catalog", key="GCST4", name="GCST4.h.tsv.gz")
    at_collect = CollectRecord(status="unreadable", source="gwas_catalog", key="GCST5", detail="no other_allele column")
    absent = CollectRecord(status="absent", source="gwas_catalog", key="GCST6")
    report = validation_report([ok, bad, indexed, at_collect, absent], results)
    assert (report["files_validated"], report["passed"], report["failed"], report["errors"]) == (2, 1, ["GCST2"], ["GCST3"])
    assert (report["unreadable_at_collect"], report["absent"], report["not_validated_remote_indexed"]) == (
        ["GCST5"], ["GCST6"], ["GCST4"])
    assert "tabix" in report["not_validated_remote_indexed_reason"]
    assert report["rules"] == RULES and set(report["files"]) == {"GCST1", "GCST2", "GCST3"}
    assert not any(in_scope(r) for r in (indexed, at_collect, absent))


def test_a_result_is_reused_while_the_record_names_the_same_bytes_and_redone_for_other_bytes(tmp_path):
    record = put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, mixed_rows()))
    first = validate_file(tmp_path, record)
    shutil.copy(tmp_path / record.path, tmp_path / "kept.tsv")
    (tmp_path / record.path).unlink()
    assert validate_file(tmp_path, record) == first                              # not read again
    other = put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, mixed_rows()[:12]))
    again = validate_file(tmp_path, other)
    assert (again["rows"], again["sha256"]) == (12, other.sha256) and math.isclose(again["rsid_recovery_rate"], 11 / 12)
    assert len(list((tmp_path / "validation" / "superseded").glob("GCST1__*.json"))) == 1


def made_unreadable(root, key: str = "GCST1", se: str = "0", ea: str = "A") -> tuple[CollectRecord, CollectRecord, dict]:
    """A file whose validation failed for its standard error alone (a zero SE on every row), or, with
    `ea`, for another reason."""
    rows = [ssf_row(pos=str(1000 + i), se=se, ea=ea) for i in range(3)]
    record = put(root, key, tsv(SSF_BETA_HEADER, rows))
    got = validate_file(root, record)
    assert got["passed"] is False
    return record, apply_validation(root, record, got), got


def all_files(root) -> dict:
    return {p.relative_to(root).as_posix(): p.read_bytes() for p in root.rglob("*") if p.is_file()}


def test_restore_moves_the_collected_record_back_archives_the_rest_and_the_next_validation_reads_the_file_again(tmp_path):
    record, unreadable, failed = made_unreadable(tmp_path)
    before = all_files(tmp_path)
    out = restore_validated(tmp_path, restorable(tmp_path, "GCST1"))
    assert read_record(tmp_path, "gwas_catalog", "GCST1") == record
    after = all_files(tmp_path)
    assert sorted(before.values()) == sorted(after.values())                 # every byte kept: moved, none deleted
    assert after[out["moved"]["unreadable_record"]] == before[record_path(tmp_path, "gwas_catalog", "GCST1")
                                                             .relative_to(tmp_path).as_posix()]
    assert json.loads(after[out["moved"]["validation_result"]]) == failed and not result_path(tmp_path, "GCST1").exists()
    (tmp_path / record.path).unlink()                                         # a reused result would not need the file
    with pytest.raises(CollectError, match="is missing"):
        validate_file(tmp_path, record)


@pytest.mark.parametrize("damage", ["bytes", "collected", "not_validation", "two_superseded", "not_the_standard_error"])
def test_restore_is_refused_and_moves_nothing_unless_the_bytes_match_and_one_validation_unreadable_record_is_undone(
        tmp_path, damage):
    record, unreadable, _ = made_unreadable(tmp_path, ea="additional copy" if damage == "not_the_standard_error" else "A")
    if damage == "bytes":
        (tmp_path / record.path).write_bytes(b"chromosome\n")
    elif damage == "collected":
        restore_validated(tmp_path, restorable(tmp_path, "GCST1"))
    elif damage == "not_validation":
        record_path(tmp_path, "gwas_catalog", "GCST1").write_text(
            unreadable.model_copy(update={"detail": "no other_allele column"}).model_dump_json())
    elif damage == "two_superseded":
        extra = tmp_path / "superseded" / "collect" / "gwas_catalog" / "GCST1__20990101T000000Z.json"
        extra.write_text(record.model_dump_json())
    before = all_files(tmp_path)
    with pytest.raises(CollectError):
        restorable(tmp_path, "GCST1")
    assert all_files(tmp_path) == before


def test_an_interrupted_restore_is_finished_by_the_next_one(tmp_path):
    record, _, _ = made_unreadable(tmp_path)
    plan = restorable(tmp_path, "GCST1")
    os.replace(record_path(tmp_path, "gwas_catalog", "GCST1"),              # stopped after archiving the unreadable record
               tmp_path / "superseded" / "collect" / "gwas_catalog" / "GCST1__20990101T000000Z__unreadable.json")
    restore_validated(tmp_path, restorable(tmp_path, "GCST1"))
    assert read_record(tmp_path, "gwas_catalog", "GCST1") == record and plan["current"].status == "unreadable"


HARMONISED_HEADER = ["hm_variant_id", "hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta",
                     "hm_odds_ratio", "hm_effect_allele_frequency", "p_value", "standard_error"]


def harmonised_row(i: int, beta: str = "0.1", rsid: str | None = None, chrom: str = "1") -> list[str]:
    return ["x", rsid if rsid is not None else f"rs{i}", chrom, str(1000 + i), "A", "G", beta, "NA", "0.2", "1e-4", "0.1"]


def test_a_harmonised_file_is_validated_by_its_own_reader_which_applies_no_allele_or_se_rule(tmp_path):
    rows = [harmonised_row(i) for i in range(6)] + [harmonised_row(6, rsid="rs0"), harmonised_row(7, rsid=".")]
    rows += [harmonised_row(8)[:-1], harmonised_row(9, beta="x"), harmonised_row(10, chrom="NA")]
    rows[1][4] = "additional copy"                          # an allele the GWAS-SSF readers reject; this reader keeps it
    rows[2][10] = "0"                                       # a zero standard error, kept likewise
    record = put(tmp_path, "GCST1", gzip.compress(tsv(HARMONISED_HEADER, rows), mtime=0), name="GCST1.h.tsv.gz", layout="",
                 rsid_rule="")
    assert in_scope(record)
    got = validate_file(tmp_path, record)
    assert (got["rows"], got["rows_kept"], got["row_width_errors"], got["parse_errors"], got["rejected"]["no_position"]) == (
        11, 8, 1, 1, 1)
    assert (got["rows_with_rsid"], got["duplicate_rsids"], got["passed"]) == (7, 1, True)
    assert sum(got["rejected"].values()) == 3


def test_a_harmonised_file_that_fails_is_made_unreadable_through_the_same_superseded_record(tmp_path):
    header = [c for c in HARMONISED_HEADER if c != "standard_error"]
    record = put(tmp_path, "GCST1", gzip.compress(tsv(header, [harmonised_row(0)[:-1]]), mtime=0), name="GCST1.h.tsv.gz",
                 layout="", rsid_rule="")
    got = validate_file(tmp_path, record)
    assert (got["header_agrees"], got["passed"]) == (False, False) and "standard_error" in got["reasons"][0]
    new = apply_validation(tmp_path, record, got)
    assert new.status == "unreadable" and new.detail.startswith(VALIDATION_DETAIL) and new.sha256 == record.sha256
    assert len(list((tmp_path / "superseded" / "collect" / "gwas_catalog").glob("GCST1__*.json"))) == 1
    fetcher = VolumeFetcher(tmp_path, "", tmp_path / "a", tmp_path / "e")
    with pytest.raises(SourceUnreadable, match="pre-analysis validation"):
        fetcher.outcome_region(OutcomeSpec(accession="GCST1", source="gwas_catalog", n_case=1, n_control=1, risk_coded=True),
                               "1", 1000, 100)


def test_the_spawn_guard_requires_a_report_covering_every_current_gwas_catalog_whole_file(tmp_path):
    ok = put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, mixed_rows()))
    bad = put(tmp_path, "GCST2", tsv(SSF_BETA_HEADER, [ssf_row(ea="additional copy")] * 3))
    with pytest.raises(InputContractError, match="no pre-analysis validation report"):
        require_validation(None, [ok, bad])
    results = {"GCST1": validate_file(tmp_path, ok), "GCST2": validate_file(tmp_path, bad)}
    report = validation_report([ok, bad], results)
    with pytest.raises(InputContractError, match=r"GCST1: the record's uncertainty mode \(none\) is not bound"):
        require_validation(report, [ok, bad])                                   # validated, but not yet applied
    ok = apply_validation(tmp_path, ok, results["GCST1"])
    with pytest.raises(InputContractError, match="GCST2: the record is collected and the report says the file failed"):
        require_validation(report, [ok, bad])
    bad = apply_validation(tmp_path, bad, results["GCST2"])
    require_validation(report, [ok, bad])
    for field, value in (("uncertainty_mode", "pvalue_coloc"), ("header_sha256", "0" * 64), ("validation_sha256", "0" * 64)):
        with pytest.raises(InputContractError, match="GCST1: the record's uncertainty mode"):
            require_validation(report, [ok.model_copy(update={field: value}), bad])
    with pytest.raises(InputContractError, match="GCST1: its result was made under other rules"):
        require_validation({**report, "files": {**report["files"], "GCST1": {**results["GCST1"], "rules": {}}}}, [ok, bad])
    with pytest.raises(InputContractError, match="GCST1: its validation ended in an error"):
        require_validation({**report, "files": {**report["files"], "GCST1": "OSError: gone"}}, [ok, bad])
    with pytest.raises(InputContractError, match="made under other rules"):
        require_validation({**report, "rules": {**RULES, "min_kept_fraction_percent": 40}}, [ok, bad])
    recollected = put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, mixed_rows()[:12]))      # stale: other bytes since the report
    with pytest.raises(InputContractError, match="GCST1: the report validated other bytes"):
        require_validation(report, [recollected, bad])
    (tmp_path / VALIDATION_DIR / REPORT_NAME).write_text(json.dumps(report))
    state = volume_state(tmp_path)
    assert state["outcome_validation"] == report
    assert {r["key"]: r["sha256"] for r in state["gwas_catalog_records"]} == {"GCST1": recollected.sha256, "GCST2": bad.sha256}


# ---- uncertainty modes: chosen at validation, bound to the record, restored files read again --------------

def or_ci_file(n: int = 60, p_scale: float = 1.0) -> bytes:
    """An odds-ratio GWAS-SSF file with `#NA` standard error on every row and 95% limits; `p_scale`
    multiplies every published p (1: each p is the Wald p of its interval)."""
    rows = []
    for i in range(n):
        row = exact_or_ci_row(0.02 * (i - n / 2), 0.03)
        row[1], row[7] = str(1000 + i), repr(min(float(row[7]) * p_scale, 0.999))
        row[8] = f"rs{i}"
        rows.append(row)
    return tsv(SSF_OR_CI_HEADER, rows)


def pvalue_file(n: int = 60, se: str = "NA") -> bytes:
    """A beta GWAS-SSF file (the beta a z-score) with p-values, frequencies and `NA` standard error."""
    return tsv(SSF_BETA_HEADER, [ssf_row(pos=str(1000 + i), effect=f"{(i - n / 2) / 10:.1f}", se=se, p=f"{0.5 / (1 + i):.4g}",
                                         rsid=f"rs{i}", eaf=f"{0.05 + 0.01 * i:.2f}") for i in range(n)])


OLD_RULES = {k: RULES[k] for k in ("min_kept_fraction_percent", "header", "integrity")}   # the rules before the modes


def failed_under_the_old_rules(root, key: str, data: bytes) -> tuple[CollectRecord, CollectRecord]:
    """A file as the volume holds it now: collected, its result under the earlier rules (every row
    rejected for its standard error), and its record made unreadable by that result."""
    record = put(root, key, data)
    rows = data.decode().count("\n") - 1
    old = {"key": key, "sha256": record.sha256, "header_agrees": True, "rows": rows, "rows_kept": 0, "integrity": None,
           "rejected": {"se_nonpositive_or_nonfinite": rows}, "passed": False, "rules": OLD_RULES,
           "reasons": [f"0 of {rows} data rows read, fewer than 50%"]}
    result_path(root, key).parent.mkdir(parents=True, exist_ok=True)
    result_path(root, key).write_text(json.dumps(old))
    return record, apply_validation(root, record, old)


def test_files_failing_only_for_a_missing_se_are_restored_validated_in_their_mode_and_bound(tmp_path):
    collected = {k: failed_under_the_old_rules(tmp_path, k, d)[0] for k, d in (("GCST1", or_ci_file()), ("GCST2", pvalue_file()))}
    cnv = put(tmp_path, "GCST3", tsv(SSF_BETA_HEADER, [ssf_row(pos=str(1000 + i), ea="additional copy") for i in range(4)]))
    apply_validation(tmp_path, cnv, validate_file(tmp_path, cnv))
    with pytest.raises(CollectError, match="did not fail for the standard error alone"):
        restorable(tmp_path, "GCST3")
    for key in collected:
        restore_validated(tmp_path, restorable(tmp_path, key))
        assert read_record(tmp_path, "gwas_catalog", key) == collected[key]
    records = [r for r in gwas_catalog_records(tmp_path) if in_scope(r)]
    results = {r.key: validate_file(tmp_path, r) for r in records}
    assert {k: (r["passed"], r["uncertainty_mode"]) for k, r in results.items()} == {
        "GCST1": (True, "or_ci_derived_se"), "GCST2": (True, "pvalue_coloc"), "GCST3": (False, "")}
    assert results["GCST1"]["se_census"] == {"rows": 60, "missing": 60, "present": 0, "not_counted": 0}
    check = results["GCST1"]["ci_p_check"]
    assert (check["rows_checked"], check["rows_agree"], check["passed"]) == (60, 60, True) and results["GCST2"]["ci_p_check"] is None
    for r in records:
        apply_validation(tmp_path, r, results[r.key])
    report = validation_report(gwas_catalog_records(tmp_path), results)
    after = gwas_catalog_records(tmp_path)
    require_validation(report, after)
    bound = {r.key: r for r in after}
    for key, mode in (("GCST1", "or_ci_derived_se"), ("GCST2", "pvalue_coloc")):
        assert (bound[key].status, bound[key].uncertainty_mode, bound[key].sha256) == ("collected", mode, collected[key].sha256)
        assert bound[key].validation_sha256 == result_sha256(report["files"][key])
        entry = collect_entry(CollectTask(source="gwas_catalog", key=key), bound[key])
        assert entry["uncertainty_mode"] == mode and entry != collect_entry(CollectTask(source="gwas_catalog", key=key),
                                                                            collected[key])
    assert bound["GCST3"].status == "unreadable"
    fetcher = VolumeFetcher(tmp_path, "", tmp_path / "a", tmp_path / "e")
    spec = OutcomeSpec(accession="GCST1", source="gwas_catalog", n_case=1, n_control=1, risk_coded=True)
    df = fetcher.outcome_region(spec, "1", 1030, 100)
    assert fetcher.outcome_mode(spec) == "or_ci_derived_se" and len(df) == 60 and df["se"].to_numpy() == pytest.approx(0.03)
    spec2 = spec.model_copy(update={"accession": "GCST2"})
    df2 = fetcher.outcome_region(spec2, "1", 1030, 100)
    assert fetcher.outcome_mode(spec2) == "pvalue_coloc" and len(df2) == 60 and df2["se"].isna().all()
    assert fetcher.outcome_mode(spec.model_copy(update={"source": "finngen"})) == "native_se"


def test_a_file_whose_se_is_missing_on_some_rows_and_present_on_others_fails_closed(tmp_path):
    rows = [ssf_row(pos=str(1000 + i), se="NA" if i % 3 else "0.01", rsid=f"rs{i}") for i in range(30)]
    record = put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, rows))
    got = validate_file(tmp_path, record)
    assert (got["passed"], got["uncertainty_mode"], got["se_census"]["missing"], got["se_census"]["present"]) == (False, "", 20, 10)
    assert got["reasons"] == ["uncertainty mode: standard_error is missing in 20 and present in 10 rows; the uncertainty mode "
                              "is chosen for the whole file, never row by row"] and "rows_kept" not in got
    new = apply_validation(tmp_path, record, got)
    assert new.status == "unreadable" and new.uncertainty_mode == "" and "never row by row" in new.detail
    assert not failed_only_for_standard_error(got)


def test_an_or_ci_file_whose_p_values_disagree_with_its_intervals_fails_the_ci_versus_p_check(tmp_path):
    record = put(tmp_path, "GCST1", or_ci_file(p_scale=3.0))           # every p 0.48 log10 above its interval's
    got = validate_file(tmp_path, record)
    assert (got["uncertainty_mode"], got["rows_kept"], got["passed"]) == ("or_ci_derived_se", 60, False)
    assert got["ci_p_check"]["rows_agree"] < 0.95 * got["ci_p_check"]["rows_checked"]
    assert got["reasons"] == [f"CI-versus-p check: {got['ci_p_check']['rows_agree']} of {got['ci_p_check']['rows_checked']} "
                              "checkable rows agree, fewer than 95%"]
    assert apply_validation(tmp_path, record, got).status == "unreadable"


def test_a_pvalue_coloc_file_counts_its_mode_specific_rejections(tmp_path):
    rows = [ssf_row(pos=str(1000 + i), se="NA", rsid=f"rs{i}") for i in range(20)]
    rows += [ssf_row(pos="2000", se="NA", p="0"), ssf_row(pos="2001", se="NA", p="#NA"), ssf_row(pos="2002", se="NA", eaf="1"),
             ssf_row(pos="2003", se="NA", eaf="#NA"), ssf_row(pos="2004", se="NA", eaf="1.5"), ssf_row(pos="2005", se="NA", effect="NA")]
    got = validate_file(tmp_path, put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, rows)))
    assert (got["uncertainty_mode"], got["rows"], got["rows_kept"], got["passed"]) == ("pvalue_coloc", 26, 20, True)
    assert {k: v for k, v in got["rejected"].items() if v} == {
        "p_zero_pvalue_coloc": 1, "p_missing_pvalue_coloc": 1, "eaf_not_inside_0_1_pvalue_coloc": 2, "eaf_outside_0_1": 1,
        "beta_nonfinite": 1}


def test_a_file_made_unreadable_under_other_rules_is_read_again_and_its_record_left_as_it_is(tmp_path):
    record, unreadable, failed = made_unreadable(tmp_path, ea="additional copy")
    stale = {**failed, "rules": OLD_RULES}
    result_path(tmp_path, "GCST1").write_text(json.dumps(stale))
    again = validate_file(tmp_path, unreadable)                         # the bytes are checked against the record and read
    assert again["rules"] == RULES and again["passed"] is False
    assert apply_validation(tmp_path, unreadable, again) == unreadable == read_record(tmp_path, "gwas_catalog", "GCST1")
    require_validation(validation_report([unreadable], {"GCST1": again}), [unreadable])
    archived = [json.loads(p.read_text()) for p in (tmp_path / "validation" / "superseded").glob("GCST1__*.json")]
    assert archived == [stale]                                          # the earlier result is archived, not deleted


def test_the_classification_check_compares_the_files_by_mode_failed_and_error(tmp_path):
    report = {"files": {"GCST1": {"passed": True, "uncertainty_mode": "native_se"},
                        "GCST2": {"passed": True, "uncertainty_mode": "or_ci_derived_se"},
                        "GCST3": {"passed": True, "uncertainty_mode": "pvalue_coloc"},
                        "GCST4": {"passed": False, "uncertainty_mode": ""}, "GCST5": "OSError: gone"}}
    assert classification(report) == {"native_se": ["GCST1"], "or_ci_derived_se": ["GCST2"], "pvalue_coloc": ["GCST3"],
                                      "failed": ["GCST4"], "error": ["GCST5"]}
    expected = {"native_se": ["GCST1"], "or_ci_derived_se": ["GCST2"], "pvalue_coloc": ["GCST3"], "failed": ["GCST4"]}
    got = check_classification(report, expected)
    assert (got["agrees"], got["mismatches"]) == (False, {"error": {"expected_not_found": [], "found_not_expected": ["GCST5"]}})
    report["files"]["GCST5"] = {"passed": True, "uncertainty_mode": "native_se"}
    assert check_classification(report, {**expected, "native_se": ["GCST1", "GCST5"]})["agrees"] is True
    assert check_classification(report, {**expected, "native_se": ["GCST1", "GCST5", "GCST9"]})["mismatches"] == {
        "native_se": {"expected_not_found": ["GCST9"], "found_not_expected": []}}
    with pytest.raises(InputContractError, match="names classes"):
        check_classification(report, {"unreadable": ["GCST4"]})
