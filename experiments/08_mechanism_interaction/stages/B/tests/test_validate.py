"""The pre-analysis validation of GWAS Catalog files read without a harmonised copy (stage_b/validate.py),
on collected files written to a temporary volume root."""
import gzip
import hashlib
import json
import math
import shutil

import pytest
from test_outcome_files import MAP_HEADER, SSF_BETA_HEADER, ci, laskar_row, ssf_row, tsv

from stage_b.author_formats import AUTHOR_FORMATS, LASKAR_HEADER
from stage_b.checkpoint import collect_entry
from stage_b.collect import read_record, record_path
from stage_b.fetch import VolumeFetcher
from stage_b.outcome_files import REJECT_REASONS
from stage_b.schemas import CollectRecord, CollectTask, InputContractError, OutcomeSpec, SourceUnreadable
from stage_b.status import volume_state
from stage_b.validate import (REPORT_NAME, RULES, VALIDATION_DETAIL, VALIDATION_DIR, apply_validation, in_scope,
                              require_validation, result_path, validate_file, validation_report)


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
    assert (got["header"], got["header_agrees"], got["passed"], got["reasons"], got["integrity"]) == (
        SSF_BETA_HEADER, True, True, [], None)
    assert json.loads(result_path(tmp_path, "GCST1").read_text()) == got
    assert apply_validation(tmp_path, record, got) == record == read_record(tmp_path, "gwas_catalog", "GCST1")
    text = json.dumps(got)
    assert "rs3" not in text and "0.01" not in text                             # no per-variant value
    assert not (tmp_path / "validation" / "work").exists() or not any((tmp_path / "validation" / "work").iterdir())


@pytest.mark.parametrize("good,bad,passed", [(2, 2, True), (1, 2, False), (0, 0, False)])
def test_a_file_passes_only_when_at_least_half_of_its_rows_are_read(tmp_path, good, bad, passed):
    rows = [ssf_row(pos=str(1000 + i)) for i in range(good)] + [ssf_row(ea="additional copy") for _ in range(bad)]
    got = validate_file(tmp_path, put(tmp_path, "GCST1", tsv(SSF_BETA_HEADER, rows)))
    assert got["passed"] is passed
    if not passed:
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


def test_an_interrupted_validation_resumes_from_its_checkpoint_and_gives_the_uninterrupted_result(tmp_path):
    data = tsv(SSF_BETA_HEADER, mixed_rows() * 3)
    whole = validate_file(tmp_path / "a", put(tmp_path / "a", "GCST1", data), checkpoint_rows=4)
    record = put(tmp_path / "b", "GCST1", data)
    commits = []

    def dies_at_the_third_commit():
        commits.append(1)
        if len(commits) == 3:
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        validate_file(tmp_path / "b", record, dies_at_the_third_commit, checkpoint_rows=4)
    state = json.loads(next((tmp_path / "b" / "validation" / "work").glob("*/state.json")).read_text())
    assert state["rows"] == 12 and not result_path(tmp_path / "b", "GCST1").exists()
    with (next((tmp_path / "b" / "validation" / "work").glob("*/variants.tsv"))).open("ab") as fh:
        fh.write(b"rs999\t1\t5\tA\tG\n")                       # a row written after the checkpoint, before the kill
    resumed = validate_file(tmp_path / "b", record, checkpoint_rows=4)
    assert {k: v for k, v in resumed.items() if k != "utc"} == {k: v for k, v in whole.items() if k != "utc"}


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
    with pytest.raises(InputContractError, match="GCST2: the record is collected and the report says the file failed"):
        require_validation(report, [ok, bad])                                   # validated, but not yet applied
    bad = apply_validation(tmp_path, bad, results["GCST2"])
    require_validation(report, [ok, bad])
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
