"""GWAS Catalog outcome files without a harmonised copy (stage_b/outcome_files.py, author_formats.py),
end to end against the local HTTP server: the study directory is resolved, the file collected with
its MD5 checked, and read back as the canonical outcome table on the build its record names."""
import gzip
import hashlib
import math

import numpy as np
import pandas as pd
import pytest
from fake_remote import FakeRemote
from test_pipeline import COLLECT, COMMIT, H4, TOOLS, FakeFetcher, StubBackend, fp, unit

from stage_b import remote as remote_module
from stage_b.assemble import build_evidence, write_outputs
from stage_b.author_formats import AUTHOR_FORMATS, LASKAR_HEADER, SPARK_HEADER, header_sha256
from stage_b.checkpoint import collect_entry
from stage_b.collect import collect_one, collect_tasks, read_record, record_path, supersede_absent
from stage_b.fetch import Endpoints, RemoteSources, VolumeFetcher, gwas_catalog_study_dir
from stage_b.outcome_files import (MODE_REJECT_REASONS, REJECT_REASONS, AuthorReader, CiPCheck, OrCiCheck, SsfReader,
                                   choose_uncertainty_mode, ci_z_for_level, half_unit, log10_two_sided_p, p_value,
                                   parse_md5sum, parse_meta_yaml, ssf_columns)
from stage_b.pipeline import DirStore, process_unit
from stage_b.schemas import (CollectError, CollectRecord, CollectTask, HypothesisInput, OutcomeSpec, SourceAbsent,
                             SourceUnreadable)
from v8_manifest import sha256_file

SSF38 = "GCST90000101"
SSF37 = "GCST90000102"
CENTER, HALF = 1_005_000, 3_000


@pytest.fixture(autouse=True)
def no_waiting(monkeypatch):
    monkeypatch.setattr(remote_module, "WAIT_S", 0.0)


@pytest.fixture
def remote():
    with FakeRemote() as r:
        yield r


def md5(b: bytes) -> str:
    return hashlib.md5(b, usedforsecurity=False).hexdigest()


def tsv(header: list[str], rows: list[list], eol: str = "\n") -> bytes:
    return "".join("\t".join(str(x) for x in line) + eol for line in [header, *rows]).encode()


def meta(name: str, data: bytes, build: str = "GRCh38", coordinate: str | None = "1-based", declared: str | None = None) -> bytes:
    lines = ["# Study meta-data", "gwas_id: X", "author_notes: 'genome_assembly: NCBI36 is a note, not a key'",
             f"genome_assembly: {build}"]
    lines += [f"coordinate_system: {coordinate}"] if coordinate is not None else []
    lines += ["samples:", "  - sample_size: 1000", f"data_file_name: {name}", "file_type: GWAS-SSF v1.0",
              f"data_file_md5sum: {declared if declared is not None else md5(data)}", "is_harmonised: false"]
    return ("\n".join(lines) + "\n").encode()


def study(remote: FakeRemote, accession: str, files: dict[str, bytes]) -> None:
    """A GWAS Catalog study directory: a listing of `files` and the files themselves."""
    base = gwas_catalog_study_dir("/gwascat", accession)
    remote.files[base] = "".join(f'<a href="{n}">{n}</a>\n' for n in ["harmonised/", *files]).encode()
    for name, data in files.items():
        remote.files[base + name] = data


def ssf_study(remote: FakeRemote, accession: str, data: bytes, gz: bool = True, md5sum: str | None = None, **meta_kw) -> str:
    name = f"{accession}.tsv" + (".gz" if gz else "")
    body = gzip.compress(data, mtime=0) if gz else data
    meta_bytes = meta(name, body, **meta_kw)
    listed = md5sum if md5sum is not None else f"{md5(body)}  {name}\n{md5(meta_bytes)}  {name}-meta.yaml\n"
    study(remote, accession, {name: body, f"{name}-meta.yaml": meta_bytes, "md5sum.txt": listed.encode()})
    return name


def sources(remote: FakeRemote, tmp_path) -> RemoteSources:
    return RemoteSources(tmp_path / "cache", "", lambda: None, Endpoints(gwascat_ftp=f"{remote.base}/gwascat"))


def fetcher(remote: FakeRemote, root) -> VolumeFetcher:
    return VolumeFetcher(root, "", root / "unused_annotation", root / "unused_excluded",
                         Endpoints(gwascat_ftp=f"{remote.base}/gwascat"))


def collect(remote, tmp_path, accession: str) -> CollectRecord:
    return collect_one(CollectTask(source="gwas_catalog", key=accession), sources(remote, tmp_path), tmp_path / "vol")


def spec(accession: str) -> OutcomeSpec:
    return OutcomeSpec(accession=accession, source="gwas_catalog", n_case=1000, n_control=9000, risk_coded=True)


SSF_OR_HEADER = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "odds_ratio", "standard_error",
                 "effect_allele_frequency", "p_value", "rsid", "ci_upper", "ci_lower", "n"]


def or_rows() -> list[list]:
    rows = [[1, 1_000_000 + 1000 * i, "G", "A", round(math.exp(0.01 * i), 6), 0.02, 0.3, 0.5, f"rs{i}", 1, 1, 5000]
            for i in range(10)]
    rows[4][8] = "#NA"                                  # no rsID for this variant
    rows[5][2] = "additional copy"                      # not an allele: the row is dropped
    rows[6][6] = "#NA"                                  # no frequency
    rows[7][0] = 2                                      # another chromosome
    return rows


# ---- GWAS-SSF: collected whole, MD5 checked, read on the build its metadata names ---------------------

def test_a_grch38_gwas_ssf_file_with_an_odds_ratio_only_is_collected_and_read_as_log_or(remote, tmp_path):
    name = ssf_study(remote, SSF38, tsv(SSF_OR_HEADER, or_rows()))
    record = collect(remote, tmp_path, SSF38)
    data = remote.files[gwas_catalog_study_dir("/gwascat", SSF38) + name]
    assert (record.status, record.name, record.md5, record.layout, record.build, record.position_offset, record.rsid_rule) == (
        "collected", name, md5(data), "gwas_ssf", "GRCh38", 0, "column")
    assert record.sha256 == hashlib.sha256(data).hexdigest() and (tmp_path / "vol" / record.path).read_bytes() == data
    f = fetcher(remote, tmp_path / "vol")
    assert f.outcome_build(spec(SSF38)) == "GRCh38"
    mark = len(remote.log)
    df = f.outcome_region(spec(SSF38), "1", CENTER, HALF)
    assert remote.log[mark:] == []                                           # read from the volume
    assert list(df["pos"]) == [1_002_000, 1_003_000, 1_004_000, 1_006_000, 1_008_000]
    assert list(df["rsid"]) == ["rs2", "rs3", "", "rs6", "rs8"]
    assert df["beta"].to_numpy() == pytest.approx([0.02, 0.03, 0.04, 0.06, 0.08], abs=1e-6)
    assert list(df["se"]) == [0.02] * 5 and list(df["n"]) == [5000.0] * 5 and list(df["p"]) == [0.5] * 5
    assert df.loc[df["rsid"] == "rs6", "eaf"].isna().all() and (df.loc[df["rsid"] != "rs6", "eaf"] == 0.3).all()
    assert set(df["ea"]) == {"G"} and set(df["oa"]) == {"A"}


def test_a_grch37_zero_based_file_is_read_on_grch37_with_one_added_to_each_position(remote, tmp_path):
    header = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error",
              "effect_allele_frequency", "p_value", "rs_id"]
    edge = CENTER + HALF                                                     # 0-based: 1-based CENTER + HALF + 1
    rows = [[1, edge - 1, "C", "T", 0.1, 0.01, 0.2, 0.01, "rs1"], [1, edge, "C", "T", 0.2, 0.01, 0.2, 0.01, "rs2"],
            [23, CENTER, "C", "T", 0.3, 0.01, 0.2, 0.01, "rs3"]]
    ssf_study(remote, SSF37, tsv(header, rows), gz=False, build="GRCh37", coordinate="0-based")
    record = collect(remote, tmp_path, SSF37)
    assert (record.layout, record.build, record.position_offset, record.rsid_rule) == ("gwas_ssf", "GRCh37", 1, "column")
    f = fetcher(remote, tmp_path / "vol")
    assert f.outcome_build(spec(SSF37)) == "GRCh37"
    df = f.outcome_region(spec(SSF37), "1", CENTER, HALF)
    assert list(df["rsid"]) == ["rs1"] and list(df["pos"]) == [edge] and list(df["beta"]) == [0.1]
    assert list(f.outcome_region(spec(SSF37), "X", CENTER + 1, 0)["rsid"]) == ["rs3"]     # chromosome 23 is X


def test_the_unit_takes_the_outcome_window_on_the_build_the_outcome_file_is_on(tmp_path):
    class Recording(FakeFetcher):
        def __init__(self):
            super().__init__()
            self.centers: dict[str, int] = {}

        def positions(self, sentinel):
            return {"GRCh37": sentinel.pos - 12_345, "GRCh38": sentinel.pos}

        def outcome_build(self, spec):
            return "GRCh37" if spec.accession == "G37" else "GRCh38"

        def outcome_region(self, spec, chrom, center, half_width):
            self.centers[spec.accession] = center
            return super().outcome_region(spec, chrom, center, half_width)

    u = unit(("G37", "G38")).model_copy(update={"outcomes": tuple(spec(a) for a in ("G37", "G38"))})
    rec = Recording()
    res = process_unit(u, rec, StubBackend(H4), DirStore(tmp_path / "u"), fp(u), TOOLS, COLLECT)
    pos = u.sentinel.pos
    assert rec.centers == {"G37": pos - 12_345, "G38": pos}
    assert res["outcomes"]["G37"]["coloc_run"] and res["positions"]["GRCh37"] == pos - 12_345


MAP_HEADER = ["ID", "REF", "ALT", "rsid", "POS19", "POS38"]


@pytest.mark.parametrize("build,pos_column", [("GRCh38", 5), ("GRCh37", 4)])
def test_a_file_without_an_rsid_column_takes_rsids_from_the_ukbppp_map_on_its_own_build(remote, tmp_path, build, pos_column):
    header = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error",
              "effect_allele_frequency", "p_value"]
    rows = [[1, 1_004_000 + k, "G", "A", 0.1, 0.01, 0.3, 0.5] for k in range(4)]
    ssf_study(remote, SSF38, tsv(header, rows), build=build)
    record = collect(remote, tmp_path, SSF38)
    assert (record.rsid_rule, record.build) == ("ukbppp_map", build)
    other = 9_000_000                                                       # a position on the other build
    def row(k, ref, alt, rs):
        p = 1_004_000 + k
        both = (p, other) if pos_column == 4 else (other, p)
        return [f"1:{both[0]}:{ref}:{alt}:imp:v1", ref, alt, rs, *both]
    map_rows = [row(0, "A", "G", "rs100"),                                   # the alleles as a pair, either order
                row(1, "C", "T", "rs101"),                                   # other alleles: no match
                row(2, "G", "A", "rs102"), row(2, "G", "A", "rs902"),        # two rsIDs: none is taken
                row(3, "G", "A", "rs103")]
    root = tmp_path / "vol"
    raw = root / "raw" / "ukbppp_rsid_map" / "olink_rsid_map_chr1.tsv.gz"
    raw.parent.mkdir(parents=True)
    raw.write_bytes(gzip.compress(tsv(MAP_HEADER, map_rows), mtime=0))
    record_path(root, "ukbppp_rsid_map", "1").parent.mkdir(parents=True)
    record_path(root, "ukbppp_rsid_map", "1").write_text(CollectRecord(
        status="collected", source="ukbppp_rsid_map", key="1", name=raw.name, path=raw.relative_to(root).as_posix(),
        bytes=raw.stat().st_size, sha256=sha256_file(raw)).model_dump_json())
    df = fetcher(remote, root).outcome_region(spec(SSF38), "1", CENTER, HALF)
    assert list(df["rsid"]) == ["rs100", "", "", "rs103"]


def test_a_downloaded_file_whose_md5_is_not_the_md5sum_txt_one_is_refused_and_not_recorded(remote, tmp_path):
    data = gzip.compress(tsv(SSF_OR_HEADER, or_rows()), mtime=0)
    wrong = "0" * 32
    ssf_study(remote, SSF38, tsv(SSF_OR_HEADER, or_rows()), declared=wrong,
              md5sum=f"{wrong}  {SSF38}.tsv.gz\n")
    assert md5(data) != wrong
    with pytest.raises(CollectError, match="MD5 other than the one the source declares"):
        collect(remote, tmp_path, SSF38)
    assert not record_path(tmp_path / "vol", "gwas_catalog", SSF38).exists()


def test_a_meta_yaml_rewritten_after_md5sum_txt_is_read_and_the_difference_is_noted_on_the_record(remote, tmp_path):
    name = ssf_study(remote, SSF38, tsv(SSF_OR_HEADER, or_rows()))
    base = gwas_catalog_study_dir("/gwascat", SSF38)
    meta_bytes = remote.files[base + name + "-meta.yaml"]
    remote.files[base + "md5sum.txt"] = remote.files[base + "md5sum.txt"].replace(md5(meta_bytes).encode(), b"f" * 32)
    record = collect(remote, tmp_path, SSF38)
    assert record.status == "collected" and record.meta_sha256 == hashlib.sha256(meta_bytes).hexdigest()
    assert record.detail == f"md5sum.txt lists another MD5 for {name}-meta.yaml than that of the file served"
    clean = ssf_study(remote, SSF37, tsv(SSF_OR_HEADER, or_rows()))
    assert clean and collect(remote, tmp_path, SSF37).detail == ""


UNREADABLE = {
    "build not GRCh37/38": ({"build": "NCBI36"}, SSF_OR_HEADER, "genome_assembly 'NCBI36'"),
    "another coordinate system": ({"coordinate": "2-based"}, SSF_OR_HEADER, "coordinate_system '2-based'"),
    "no other_allele": ({}, [c for c in SSF_OR_HEADER if c != "other_allele"], "no other_allele column"),
    "no effect": ({}, ["chromosome", "base_pair_location", "effect_allele", "other_allele", "hazard_ratio", "standard_error",
                       "effect_allele_frequency", "p_value"], "neither a beta nor an odds_ratio"),
    "both effects, column 4 neither": ({}, ["chromosome", "base_pair_location", "effect_allele", "other_allele", "z",
                                            "standard_error", "beta", "odds_ratio"], "column 4 is neither"),
    "md5s disagree": ({"declared": "e" * 32}, SSF_OR_HEADER, "md5sum.txt disagree"),
}


@pytest.mark.parametrize("case", list(UNREADABLE))
def test_a_file_no_fixed_rule_reads_is_recorded_unreadable_with_its_reason_and_reads_as_unavailable(remote, tmp_path, case):
    meta_kw, header, reason = UNREADABLE[case]
    rows = [[1, 1_004_000] + ["A"] * (len(header) - 2)]
    ssf_study(remote, SSF38, tsv(header, rows), **meta_kw)
    record = collect(remote, tmp_path, SSF38)
    assert record.status == "unreadable" and reason in record.detail and record.path == "" and record.sha256 == ""
    assert not (tmp_path / "vol" / "raw").exists()                          # nothing downloaded
    with pytest.raises(SourceUnreadable, match=reason.replace("(", r"\(").replace(")", r"\)")) as err:
        fetcher(remote, tmp_path / "vol").outcome_region(spec(SSF38), "1", CENTER, HALF)
    assert isinstance(err.value, SourceAbsent)                               # the pipeline's unavailable consequence


def test_unreadable_and_unharmonised_records_enter_the_collect_digest_with_how_they_are_read():
    task = CollectTask(source="gwas_catalog", key="GCST1")
    ssf = CollectRecord(status="collected", source="gwas_catalog", key="GCST1", name="GCST1.tsv.gz", path="raw/gwas_catalog/x",
                        bytes=5, sha256="a" * 64, layout="gwas_ssf", build="GRCh37", position_offset=1, rsid_rule="ukbppp_map",
                        meta_sha256="b" * 64)
    assert collect_entry(task, ssf) == {"source": "gwas_catalog", "key": "GCST1", "status": "collected", "name": "GCST1.tsv.gz",
                                        "bytes": 5, "sha256": "a" * 64, "layout": "gwas_ssf", "build": "GRCh37",
                                        "position_offset": 1, "rsid_rule": "ukbppp_map", "meta_sha256": "b" * 64}
    assert collect_entry(task, ssf.model_copy(update={"build": "GRCh38"})) != collect_entry(task, ssf)
    harmonised = ssf.model_copy(update={"layout": "", "build": "", "position_offset": 0, "rsid_rule": "", "meta_sha256": ""})
    assert set(collect_entry(task, harmonised)) == {"source", "key", "status", "name", "bytes", "sha256"}
    bad = CollectRecord(status="unreadable", source="gwas_catalog", key="GCST1", detail="GWAS Catalog GCST1: no beta")
    assert collect_entry(task, bad) == {"source": "gwas_catalog", "key": "GCST1", "status": "unreadable", "name": "",
                                        "bytes": None, "sha256": "", "reason": "GWAS Catalog GCST1: no beta"}


def test_a_unit_with_a_gwas_catalog_outcome_reads_the_rsid_map_of_its_chromosome_once():
    u = unit(("F_ok",)).model_copy(update={"outcomes": (spec("GCST1"),)})
    assert [(t.source, t.key) for t in collect_tasks([u])] == [("decode", "1_1"), ("decode_smp", "1_1"), ("gwas_catalog", "GCST1"),
                                                               ("ukbppp_rsid_map", "1")]
    assert ("ukbppp_rsid_map", "1") in [(t.source, t.key) for t in collect_tasks([unit(("F_ok",))])]   # read for LD by every unit


# ---- the re-collection of accessions recorded absent under the harmonised-only rule -----------------

def test_only_an_absent_gwas_catalog_record_is_superseded_and_it_is_moved_not_deleted(remote, tmp_path):
    root = tmp_path / "vol"
    old = CollectRecord(status="absent", source="gwas_catalog", key=SSF38, detail=f"GWAS Catalog {SSF38}: 0 harmonised files")
    record_path(root, "gwas_catalog", SSF38).parent.mkdir(parents=True)
    record_path(root, "gwas_catalog", SSF38).write_text(old.model_dump_json())
    ssf_study(remote, SSF38, tsv(SSF_OR_HEADER, or_rows()))
    task = CollectTask(source="gwas_catalog", key=SSF38)
    assert collect_one(task, sources(remote, tmp_path), root) == old         # a record is never redone in place
    moved = supersede_absent(root, task)
    assert moved is not None and moved.parent == root / "superseded" / "collect" / "gwas_catalog"
    assert CollectRecord.model_validate_json(moved.read_text()) == old
    new = collect_one(task, sources(remote, tmp_path), root)
    assert new.status == "collected" and read_record(root, "gwas_catalog", SSF38) == new
    assert supersede_absent(root, task) is None and read_record(root, "gwas_catalog", SSF38) == new   # a collected record stays
    decode = CollectTask(source="decode", key="1_1")
    record_path(root, "decode", "1_1").parent.mkdir(parents=True)
    record_path(root, "decode", "1_1").write_text(CollectRecord(status="absent", source="decode", key="1_1").model_dump_json())
    assert supersede_absent(root, decode) is None and record_path(root, "decode", "1_1").is_file()


# ---- reviewed author formats ------------------------------------------------------------------------

def laskar_rows() -> list[list]:
    return [["rs1", 1, 1_004_000, "A", "G", "G", round(math.exp(0.2), 6), 1, 1, 0.05, 0.001, 900, 100, 0.3, 0.35],
            ["rs2", 1, 1_004_100, "A", "G", "A", round(math.exp(-0.1), 6), 1, 1, 0.04, 0.2, 900, 100, 0.3, 0.35],
            ["rs3", 1, 1_004_200, "A", "G", "T", 1.1, 1, 1, 0.04, 0.2, 900, 100, 0.3, 0.35],      # effect allele is neither
            ["rs4", 1, 1_004_300, "C", "C", "C", 1.1, 1, 1, 0.04, 0.2, 900, 100, 0.3, 0.35],      # no other allele
            ["rs5", 1, 1_104_300, "A", "G", "G", 1.1, 1, 1, 0.04, 0.2, 900, 100, 0.3, 0.35]]      # outside the window


def test_the_laskar_author_format_is_read_by_its_reviewed_map(remote, tmp_path):
    fmt = AUTHOR_FORMATS["GCST008226"]
    data = gzip.compress(tsv(list(LASKAR_HEADER), laskar_rows(), eol="\r\n"), mtime=0)     # the deposit's CRLF line ends
    study(remote, "GCST008226", {fmt.data_file: data, "read-me_Laskar_31231134.txt": b"## readme\n"})
    record = collect(remote, tmp_path, "GCST008226")
    assert (record.status, record.name, record.layout, record.build, record.rsid_rule, record.md5) == (
        "collected", fmt.data_file, "author", "GRCh37", "column", "")
    f = fetcher(remote, tmp_path / "vol")
    assert f.outcome_build(spec("GCST008226")) == "GRCh37"
    df = f.outcome_region(spec("GCST008226"), "1", CENTER, HALF)
    assert df[["rsid", "ea", "oa"]].values.tolist() == [["rs1", "G", "A"], ["rs2", "A", "G"]]
    assert df["beta"].to_numpy() == pytest.approx([0.2, -0.1], abs=1e-6)
    assert list(df["se"]) == [0.05, 0.04] and list(df["p"]) == [0.001, 0.2] and list(df["n"]) == [1000.0, 1000.0]
    assert df["eaf"].isna().all()                                            # cases and controls only: no overall frequency


def test_an_author_file_whose_header_is_not_the_reviewed_one_is_unreadable(remote, tmp_path):
    fmt = AUTHOR_FORMATS["GCST008226"]
    header = ["Allele_2" if c == "Allele_1" else "Allele_1" if c == "Allele_2" else c for c in LASKAR_HEADER]
    study(remote, "GCST008226", {fmt.data_file: gzip.compress(tsv(header, laskar_rows()), mtime=0)})
    record = collect(remote, tmp_path, "GCST008226")
    assert record.status == "unreadable" and "not the one its reviewed author-format map was written for" in record.detail


def test_the_spark_asd_deposit_is_unreadable_for_what_its_documents_leave_unstated(remote, tmp_path):
    fmt = AUTHOR_FORMATS["GCST010514"]
    study(remote, "GCST010514", {fmt.data_file: tsv(list(SPARK_HEADER), [[1, 1_004_000, "rs1", "a", "g", 0.1, 0.01, 0.5,
                                                                          "++", 1000]]), "md5sum.txt": b""})
    record = collect(remote, tmp_path, "GCST010514")
    assert record.status == "unreadable" and "the effect allele" in record.detail and "genome build" in record.detail


def test_each_author_format_names_the_header_it_was_reviewed_against():
    for fmt in AUTHOR_FORMATS.values():
        assert header_sha256(fmt.header) == fmt.header_sha256
        assert all(len(sha) == 64 for sha in fmt.documents.values())
        if fmt.readable:
            mapped = [fmt.rsid, fmt.chrom, fmt.pos, fmt.effect_allele, *fmt.allele_pair, fmt.odds_ratio or fmt.beta, fmt.se,
                      *fmt.n_sum, *([fmt.p] if fmt.p else [])]
            assert set(mapped) <= set(fmt.header) and fmt.build in ("GRCh37", "GRCh38")
        else:
            assert fmt.reason


# ---- the column rules on the header lines of the 66 deposits with a GWAS-SSF-named file ----------------

REAL_HEADERS = {   # header line (tab-separated) -> (effect column, rsID column) or the start of the unreadable reason
    "chromosome base_pair_location reference_allele alternative_allele odds_ratio OR_se p_value n_samples effects":
        "the header has no effect_allele, other_allele, standard_error column",
    "chromosome variant_id base_pair_location effect_allele other_allele N effect_allele_frequency T SE_T P_noSPA beta "
    "standard_error p_value CONVERGE": ("beta", None),
    "chromosome variant_id base_pair_location effect_allele odds_ratio standard_error ci_lower ci_upper p_value":
        "the header has no other_allele column",
    "variant_id p_value chromosome base_pair_location effect_allele other_allele beta standard_error I2": ("beta", None),
    "chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rs_id":
        ("beta", "rs_id"),
    "chromosome base_pair_location effect_allele other_allele odds_ratio standard_error effect_allele_frequency p_value rs_id "
    "ci_lower ci_upper n": ("odds_ratio", "rs_id"),
    "chromosome base_pair_location effect_allele other_allele odds_ratio standard_error effect_allele_frequency p_value rsid "
    "ci_upper ci_lower alt n case_af num_cases control_af num_controls r2 q_pval i2 direction": ("odds_ratio", "rsid"),
    "Gene Z_score pval": "the header has no chromosome, base_pair_location, effect_allele, other_allele, standard_error column",
    "chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value":
        ("beta", None),
    "chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rsid logP "
    "n_case N_total isq_het p_het key rsid2": ("beta", "rsid"),
    "chromosome base_pair_location effect_allele other_allele odds_ratio standard_error effect_allele_frequency p_value "
    "ci_upper ci_lower rsid n beta se_beta": ("odds_ratio", "rsid"),
}


@pytest.mark.parametrize("line", list(REAL_HEADERS))
def test_the_column_rules_on_each_header_layout_the_deposits_use(line):
    expected = REAL_HEADERS[line]
    if isinstance(expected, str):
        with pytest.raises(SourceUnreadable, match=expected):
            ssf_columns(line.split(), "GWAS Catalog X")
    else:
        cols = ssf_columns(line.split(), "GWAS Catalog X")
        assert (cols.effect, cols.rsid) == expected


def test_md5sum_parsing_keeps_names_with_one_md5():
    assert parse_md5sum("A" * 32 + "  README copy.md\n" + "b" * 32 + " *x.tsv.gz\r\n" + "c" * 32 + "  d\n" + "e" * 32 + "  d\n") == {
        "README copy.md": "a" * 32, "x.tsv.gz": "b" * 32}


# ---- -meta.yaml: PyYAML's SafeLoader, duplicate keys refused, the keys read checked ---------------------

META_TEXT = ("# Study meta-data\n"
             "genome_assembly: GRCh37   # a trailing comment\n"
             "author_notes: 'genome_assembly: NCBI36 is a note, not a key'\n"
             "data_file_name: \"x.tsv: gz\"\n"
             "coordinate_system: 1-based\n"
             "samples:\n  - genome_assembly: NCBI36\n    sample_size: 3\n")


@pytest.mark.parametrize("encoding", ["lf", "crlf", "bom"])
def test_the_meta_yaml_reads_top_level_keys_through_comments_quoted_colons_crlf_and_a_bom(encoding):
    text = META_TEXT.replace("\n", "\r\n") if encoding == "crlf" else META_TEXT
    raw = (b"\xef\xbb\xbf" if encoding == "bom" else b"") + text.encode()
    meta = parse_meta_yaml(raw, "GWAS Catalog X")
    assert (meta.genome_assembly, meta.data_file_name, meta.coordinate_system, meta.data_file_md5sum) == (
        "GRCh37", "x.tsv: gz", "1-based", None)


def test_a_plain_scalar_is_read_as_the_string_it_is_written_as_so_an_all_digit_md5_keeps_its_digits():
    raw = (b"genome_assembly: GRCh38\ndata_file_name: 2024\ncoordinate_system:\n"
           b"data_file_md5sum: 01234567890123456789012345678901\n")
    meta = parse_meta_yaml(raw, "GWAS Catalog X")
    assert (meta.data_file_md5sum, meta.data_file_name, meta.coordinate_system) == (
        "01234567890123456789012345678901", "2024", None)


@pytest.mark.parametrize("raw,reason", [
    (b"genome_assembly: GRCh37\ngenome_assembly: GRCh38\ndata_file_name: x\n", "'genome_assembly' is given twice"),
    (b"genome_assembly: GRCh37\ndata_file_name: x\nsamples:\n  n: 1\n  n: 2\n", "'n' is given twice"),
    ("genome_assembly: GRCh37\ndata_file_name: café\n".encode("latin-1"), "not UTF-8"),
    ("genome_assembly: GRCh37\ndata_file_name: x\n".encode("utf-16"), "not UTF-8"),
    (b"genome_assembly: [GRCh37\ndata_file_name: x\n", "does not parse as YAML"),
    (b"- genome_assembly: GRCh37\n", "not a mapping"),
    (b"data_file_name: x\n", "no string value for genome_assembly"),
    (b"genome_assembly: GRCh37\ndata_file_name: x\ndata_file_md5sum: [a, b]\n", "no string value for data_file_md5sum"),
    (b"genome_assembly:\n  build: GRCh37\ndata_file_name: x\n", "no string value for genome_assembly"),
])
def test_a_meta_yaml_that_is_not_utf8_does_not_parse_repeats_a_key_or_lacks_a_string_is_refused(raw, reason):
    with pytest.raises(SourceUnreadable, match=reason):
        parse_meta_yaml(raw, "GWAS Catalog X")


def test_a_study_whose_meta_yaml_repeats_a_key_is_recorded_unreadable_before_any_download(remote, tmp_path):
    name = ssf_study(remote, SSF38, tsv(SSF_OR_HEADER, or_rows()))
    base = gwas_catalog_study_dir("/gwascat", SSF38)
    remote.files[base + name + "-meta.yaml"] += b"genome_assembly: GRCh37\n"
    record = collect(remote, tmp_path, SSF38)
    assert record.status == "unreadable" and "'genome_assembly' is given twice" in record.detail
    assert not (tmp_path / "vol" / "raw").exists()


# ---- rows: alleles, whitespace, p-values and the reasons a row is not read, in both readers --------------

SSF_BETA_HEADER = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error",
                   "effect_allele_frequency", "p_value", "rsid"]
SSF_OR_ONLY_HEADER = [c if c != "beta" else "odds_ratio" for c in SSF_BETA_HEADER]


def ssf_row(ea="A", oa="G", effect="0.1", se="0.01", p="0.5", rsid="rs1", chrom="1", pos="1000", eaf="0.3") -> list[str]:
    return [chrom, pos, ea, oa, effect, se, eaf, p, rsid]


def laskar_row(ea="G", a1="A", a2="G", odds="1.2", se="0.05", p="0.01", rsid="rs1", lo="1.1", hi="1.3") -> list[str]:
    return [rsid, "1", "1000", a1, a2, ea, odds, lo, hi, se, p, "900", "100", "0.3", "0.35"]


def ssf(header=SSF_BETA_HEADER) -> SsfReader:
    return SsfReader(header, 0, "GWAS Catalog X")


def laskar() -> AuthorReader:
    return AuthorReader(AUTHOR_FORMATS["GCST008226"], list(LASKAR_HEADER))


ALLELES = [      # (effect allele, other allele) -> the alleles kept, or the reason the row is not read
    (("A", "G"), ("A", "G")),
    (("N", "G"), "invalid_allele"),
    (("A", "N"), "invalid_allele"),
    (("N", "D"), "invalid_allele"),
    (("", "G"), "invalid_allele"),
    (("a", "g"), ("A", "G")),
    ((" a ", "g "), ("A", "G")),                  # surrounding whitespace removed before the allele rule
    (("A", "A"), "equal_alleles"),
    (("ac", "AC "), "equal_alleles"),
]


@pytest.mark.parametrize("alleles,expected", ALLELES)
def test_the_gwas_ssf_reader_keeps_a_row_only_with_two_valid_different_alleles(alleles, expected):
    got = ssf().row(ssf_row(*alleles))
    assert got == expected if isinstance(expected, str) else (got["ea"], got["oa"]) == expected


@pytest.mark.parametrize("alleles,expected", ALLELES)
def test_the_author_reader_keeps_a_row_only_with_two_valid_different_alleles(alleles, expected):
    effect, other = alleles
    got = laskar().row(laskar_row(ea=effect, a1=other, a2=effect))      # the pair holds both; the other is not the effect
    assert got == expected if isinstance(expected, str) else (got["ea"], got["oa"]) == expected


def test_the_author_reader_drops_a_row_whose_effect_allele_is_neither_allele_of_the_pair():
    assert laskar().row(laskar_row(ea="T", a1="A", a2="G")) == "effect_allele_not_in_pair"
    assert laskar().row(laskar_row(ea="g", a1=" a", a2="G "))["oa"] == "A"


@pytest.mark.parametrize("text,expected", [("rs12", "rs12"), (" rs12 ", "rs12"), ("rs12x", ""), (" 1:1000:A:G", ""),
                                           ("#NA", "")])
def test_an_rsid_is_validated_after_its_surrounding_whitespace_is_removed_in_both_readers(text, expected):
    assert ssf().row(ssf_row(rsid=text))["rsid"] == expected
    assert laskar().row(laskar_row(rsid=text))["rsid"] == expected


NAN = float("nan")
P_VALUES = [     # (p column, text) -> p, or the reason the row is not read
    ("p_value", "0.5", 0.5), ("p_value", "0", 0.0), ("p_value", "1", 1.0), ("p_value", "3e-300", 3e-300),
    ("p_value", "#NA", NAN), ("p_value", "", NAN), ("p_value", "1.5", "p_outside_0_1"), ("p_value", "-0.1", "p_outside_0_1"),
    ("p_value", "inf", "p_outside_0_1"),
    ("neg_log_10_p_value", "2", 0.01), ("neg_log_10_p_value", "0", 1.0), ("neg_log_10_p_value", "NA", NAN),
    ("neg_log_10_p_value", "400", 0.0), ("neg_log_10_p_value", "1e6", 0.0), ("neg_log_10_p_value", "inf", 0.0),
    ("neg_log_10_p_value", "-1", "neg_log10_p_negative"), ("neg_log_10_p_value", "-inf", "neg_log10_p_negative"),
]


@pytest.mark.parametrize("column,text,expected", P_VALUES)
def test_an_ordinary_and_a_neg_log10_p_value_are_each_read_on_their_own_branch(column, text, expected):
    header = [column if c == "p_value" else c for c in SSF_BETA_HEADER]
    got = ssf(header).row(ssf_row(p=text))
    if isinstance(expected, str):
        assert got == expected
    elif math.isnan(expected):
        assert math.isnan(got["p"])
    else:
        assert got["p"] == pytest.approx(expected, rel=1e-12, abs=0.0)
    assert p_value(0.5, False) == 0.5 and p_value(0.5, True) == pytest.approx(10 ** -0.5)


REJECTED = [     # (reader, header, row) -> the reason
    ("ssf", SSF_BETA_HEADER, ssf_row(se="0"), "se_nonpositive_or_nonfinite"),
    ("ssf", SSF_BETA_HEADER, ssf_row(se="-0.01"), "se_nonpositive_or_nonfinite"),
    ("ssf", SSF_BETA_HEADER, ssf_row(se="inf"), "se_nonpositive_or_nonfinite"),
    ("ssf", SSF_BETA_HEADER, ssf_row(se="#NA"), "se_nonpositive_or_nonfinite"),
    ("ssf", SSF_OR_ONLY_HEADER, ssf_row(effect="0"), "or_nonpositive_or_nonfinite"),
    ("ssf", SSF_OR_ONLY_HEADER, ssf_row(effect="-1.2"), "or_nonpositive_or_nonfinite"),
    ("ssf", SSF_OR_ONLY_HEADER, ssf_row(effect="inf"), "or_nonpositive_or_nonfinite"),
    ("ssf", SSF_OR_ONLY_HEADER, ssf_row(effect="NA"), "or_nonpositive_or_nonfinite"),
    ("ssf", SSF_BETA_HEADER, ssf_row(effect="NA"), "beta_nonfinite"),
    ("ssf", SSF_BETA_HEADER, ssf_row(effect="inf"), "beta_nonfinite"),
    ("ssf", SSF_BETA_HEADER, ssf_row(effect="-inf", se="0"), "beta_nonfinite"),
    ("ssf", SSF_BETA_HEADER, ssf_row(eaf="1.2"), "eaf_outside_0_1"),
    ("ssf", SSF_BETA_HEADER, ssf_row(eaf="-0.01"), "eaf_outside_0_1"),
    ("ssf", SSF_BETA_HEADER, ssf_row(eaf="inf"), "eaf_outside_0_1"),
    ("ssf", SSF_BETA_HEADER, ssf_row(effect="0.1x"), "unparseable_number"),
    ("ssf", SSF_BETA_HEADER, ssf_row(pos="12a"), "unparseable_number"),
    ("ssf", SSF_BETA_HEADER, ssf_row(chrom="#NA"), "no_position"),
    ("ssf", SSF_BETA_HEADER, ssf_row()[:-1], "row_width"),
    ("author", LASKAR_HEADER, laskar_row(se="0"), "se_nonpositive_or_nonfinite"),
    ("author", LASKAR_HEADER, laskar_row(se="NA"), "se_nonpositive_or_nonfinite"),
    ("author", LASKAR_HEADER, laskar_row(odds="0"), "or_nonpositive_or_nonfinite"),
    ("author", LASKAR_HEADER, laskar_row(odds="-inf"), "or_nonpositive_or_nonfinite"),
    ("author", LASKAR_HEADER, laskar_row(p="1.01"), "p_outside_0_1"),
    ("author", LASKAR_HEADER, laskar_row(ea="N", a1="N", a2="N"), "invalid_allele"),
    ("author", LASKAR_HEADER, laskar_row(ea="A", a1="A", a2="A"), "equal_alleles"),
    ("author", LASKAR_HEADER, laskar_row() + ["extra"], "row_width"),
]


@pytest.mark.parametrize("kind,header,row,reason", REJECTED)
def test_a_row_is_not_read_for_the_first_reason_that_applies(kind, header, row, reason):
    reader = ssf(header) if kind == "ssf" else laskar()
    assert reader.row(row) == reason


def test_a_row_with_a_finite_positive_odds_ratio_reads_as_its_log_and_a_beta_row_keeps_its_sign():
    assert ssf(SSF_OR_ONLY_HEADER).row(ssf_row(effect="2.0"))["beta"] == pytest.approx(math.log(2.0))
    assert ssf().row(ssf_row(effect="-0.3"))["beta"] == -0.3
    assert laskar().row(laskar_row(odds="0.5"))["beta"] == pytest.approx(math.log(0.5))


def test_the_window_reader_skips_rejected_rows_and_rows_outside_the_window_and_keeps_the_rest():
    rows = [ssf_row(pos="1000"), ssf_row(pos="1001", se="0"), ssf_row(pos="1002")[:-1], ssf_row(pos="9000"),
            ssf_row(pos="1003", ea="a", rsid=" rs9 ")]
    df = ssf().window(rows, "1", 1000, 10)
    assert list(df["pos"]) == [1000, 1003] and list(df["rsid"]) == ["rs1", "rs9"] and list(df["ea"]) == ["A", "A"]


# ---- the Laskar OR / CI / SE integrity rule ------------------------------------------------------------------

def ci(odds: float, se: float, digits: int = 6, se_for_ci: float | None = None) -> list[str]:
    """A Laskar row whose limits are exp(ln OR -/+ 1.96 se_for_ci), printed with `digits` decimals."""
    s = se if se_for_ci is None else se_for_ci
    lo, hi = math.exp(math.log(odds) - 1.96 * s), math.exp(math.log(odds) + 1.96 * s)
    return laskar_row(odds=f"{odds:.{digits}f}", se=f"{se:.{digits}f}", lo=f"{lo:.{digits}f}", hi=f"{hi:.{digits}f}")


def run_check(rows) -> dict:
    check = OrCiCheck(AUTHOR_FORMATS["GCST008226"], list(LASKAR_HEADER))
    for r in rows:
        check.add(r)
    return check.result()


def consistent(k: int) -> list[list[str]]:
    return [ci(0.5 + 0.01 * i, 0.02 + 0.0001 * i) for i in range(k)]


def test_limits_equal_to_exp_log_or_plus_minus_1_96_se_pass():
    got = run_check(consistent(1000))
    assert (got["rows_checked"], got["rows_agree"], got["rows_disagree"], got["passed"]) == (1000, 1000, 0, True)


def test_a_standard_error_on_another_scale_than_log_or_fails():
    rows = [ci(1.5 + 0.01 * i, 0.03 * (1.5 + 0.01 * i), se_for_ci=0.03) for i in range(500)]   # SE of OR, not of ln(OR)
    got = run_check(rows)
    assert (got["rows_checked"], got["rows_agree"], got["passed"]) == (500, 0, False)


@pytest.mark.parametrize("bad,passed", [(10, True), (11, False)])
def test_the_file_passes_only_when_at_least_99_percent_of_checked_rows_agree(bad, passed):
    rows = consistent(1000 - bad) + [ci(1.2, 0.05, se_for_ci=0.08) for _ in range(bad)]
    got = run_check(rows)
    assert (got["rows_checked"], got["rows_agree"], got["passed"]) == (1000, 1000 - bad, passed)


def test_the_tolerance_is_the_printed_rounding_where_that_is_looser_than_one_percent():
    # OR 0.01, SE 0.1: the limits are 0.00822 and 0.01216, 18% and 22% from what two decimals print
    assert run_check([laskar_row(odds="0.01", se="0.1", lo="0.01", hi="0.01")])["passed"] is True
    assert run_check([laskar_row(odds="0.0100", se="0.1000", lo="0.0082", hi="0.0122")])["passed"] is True
    assert run_check([laskar_row(odds="0.0100", se="0.1000", lo="0.0100", hi="0.0100")])["passed"] is False
    assert (half_unit("1.23"), half_unit("1.2e-05"), half_unit("3")) == (pytest.approx(0.005), pytest.approx(5e-7), 0.5)


def test_rows_without_finite_positive_or_ci_and_se_are_not_checked_and_no_checked_row_fails():
    rows = [laskar_row(odds="NA"), laskar_row(odds="0"), laskar_row(se="0"), laskar_row(lo="#NA"), laskar_row()[:-1]]
    got = run_check(rows)
    assert (got["rows_checked"], got["rows_not_checkable"], got["passed"]) == (0, 5, False)


def test_a_header_without_the_confidence_limit_columns_makes_the_study_unreadable():
    fmt = AUTHOR_FORMATS["GCST008226"]
    with pytest.raises(SourceUnreadable, match="no confidence-limit column ci_lower, ci_upper"):
        OrCiCheck(fmt, [c for c in LASKAR_HEADER if c not in ("ci_lower", "ci_upper")])
    with pytest.raises(SourceUnreadable, match=r"no confidence-limit column \(none mapped\)"):
        OrCiCheck(fmt.model_copy(update={"ci_lower": ""}), list(LASKAR_HEADER))


# ---- an unreadable outcome: the record stays `unreadable`, the pipeline applies the unavailable consequence ----

def test_an_unreadable_outcome_keeps_its_unreadable_record_while_its_hypothesis_takes_the_unavailable_consequence(
        remote, tmp_path):
    header = [c for c in SSF_OR_HEADER if c != "other_allele"]
    ssf_study(remote, SSF38, tsv(header, [[1, 1_004_000] + ["A"] * (len(header) - 2)]))
    root = tmp_path / "vol"
    record = collect(remote, tmp_path, SSF38)
    assert record.status == "unreadable" and "no other_allele column" in record.detail
    volume = fetcher(remote, root)

    class Mixed(FakeFetcher):
        def outcome_region(self, spec, chrom, center, half_width):
            if spec.source == "gwas_catalog":
                return volume.outcome_region(spec, chrom, center, half_width)
            return super().outcome_region(spec, chrom, center, half_width)

        def outcome_build(self, spec):
            return volume.outcome_build(spec) if spec.source == "gwas_catalog" else super().outcome_build(spec)

    u = unit(("F_ok",))
    u = u.model_copy(update={"outcomes": (*u.outcomes, spec(SSF38))})
    store = DirStore(tmp_path / "u")
    res = process_unit(u, Mixed(), StubBackend(H4), store, fp(u), TOOLS, COLLECT)
    out = res["outcomes"][SSF38]
    assert (out["coloc_run"], out["not_run_reason"]) == (False, "outcome_file_unavailable")
    assert "no other_allele column" in out["detail"]
    assert store.json(f"outcome__{SSF38}.meta.json")["status"] == "unavailable"
    assert read_record(root, "gwas_catalog", SSF38).status == "unreadable"              # the provenance record
    assert collect_entry(CollectTask(source="gwas_catalog", key=SSF38), record)["status"] == "unreadable"
    h = HypothesisInput(hypothesis_id="h1", gene_symbol="G", gene_ensembl="ENSG1", direction="decrease",
                        instrument_source="decode", instrument_assay_id="1_1", platform="SomaScan", outcome_accession=SSF38,
                        outcome_source="gwas_catalog", outcome_n_case=1000, outcome_n_control=9000)
    row = build_evidence([h], {"h1": u.unit_key}, {u.unit_key: ("decode", "1_1", "ENSG1")}, {u.unit_key: res}, {})[0]
    assert (row.coloc_run, row.not_run_reason, row.evidence_state) == (False, "outcome_file_unavailable", "inconclusive")
    src = tmp_path / "in.csv"
    src.write_text("x")
    write_outputs(tmp_path / "out", [row], [], [src], {"hypotheses": src}, [("", tmp_path)], script_root=tmp_path,
                  run_token="stageb-token-0001", repo_commit=COMMIT, tools=TOOLS, collected=[record])
    published = pd.read_csv(tmp_path / "out" / "collected_files.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert published[["key", "status"]].values.tolist() == [[SSF38, "unreadable"]]
    assert "no other_allele column" in published.loc[0, "detail"]


# ---- uncertainty modes: native_se, or_ci_derived_se, pvalue_coloc --------------------------------------------

SSF_OR_CI_HEADER = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "odds_ratio", "standard_error",
                    "effect_allele_frequency", "p_value", "rsid", "ci_upper", "ci_lower"]


def or_ci_row(odds="1.2", lower="1.1", upper="1.3", p="0.01", eaf="0.3", se="#NA", pos="1000") -> list[str]:
    return ["1", pos, "A", "G", odds, se, eaf, p, "rs1", upper, lower]


def or_ci(header=SSF_OR_CI_HEADER) -> SsfReader:
    return SsfReader(header, 0, "GWAS Catalog X", "or_ci_derived_se")


def pvalue(header=SSF_BETA_HEADER) -> SsfReader:
    return SsfReader(header, 0, "GWAS Catalog X", "pvalue_coloc")


def exact_or_ci_row(beta: float, se: float, p: float | None = None, digits: int | None = None) -> list[str]:
    """A row whose limits are exp(beta -/+ 1.96 se), at full precision or printed with `digits` significant digits."""
    fmt = (lambda x: repr(x)) if digits is None else (lambda x: f"{x:.{digits}g}")
    p = math.erfc(abs(beta / se) / math.sqrt(2)) if p is None else p
    return or_ci_row(odds=fmt(math.exp(beta)), lower=fmt(math.exp(beta - 1.96 * se)), upper=fmt(math.exp(beta + 1.96 * se)),
                     p=repr(p))


def test_each_mode_has_its_own_reasons_and_every_reason_is_listed_once():
    assert len(set(REJECT_REASONS)) == len(REJECT_REASONS)
    own = [r for reasons in MODE_REJECT_REASONS.values() for r in reasons]
    assert len(own) == len(set(own)) and set(own) <= set(REJECT_REASONS)
    assert ci_z_for_level(0.95) == 1.96 and ci_z_for_level(0.90) == pytest.approx(1.6448536269514722, rel=1e-12)
    with pytest.raises(ValueError):
        ci_z_for_level(95)


def test_the_ci_route_gives_the_se_of_ln_or_from_the_width_of_the_log_interval_whatever_its_centre():
    got = or_ci().row(exact_or_ci_row(math.log(2.0), 0.1))
    assert got["beta"] == pytest.approx(math.log(2.0), rel=1e-15) and got["se"] == pytest.approx(0.1, rel=1e-12)
    # symmetric on the natural scale, not on the log scale: SE from ln U - ln L only, beta from OR only
    got = or_ci().row(or_ci_row(odds="2.0", lower="1.8", upper="2.2"))
    assert got["se"] == pytest.approx((math.log(2.2) - math.log(1.8)) / (2 * 1.96), rel=1e-15)
    assert got["beta"] == pytest.approx(math.log(2.0), rel=1e-15)
    assert or_ci().row(or_ci_row(se="0.05"))["se"] == pytest.approx((math.log(1.3) - math.log(1.1)) / 3.92, rel=1e-15)


def test_ci_derived_se_recovers_the_se_within_the_rounding_of_six_printed_digits_over_many_rows():
    rng = np.random.default_rng()
    beta, se = rng.normal(0, 0.3, 20_000), rng.uniform(0.005, 0.5, 20_000)
    derived = np.array([or_ci().row(exact_or_ci_row(b, s, digits=6))["se"] for b, s in zip(beta, se)])
    # each limit carries a relative rounding of at most 5e-6, so ln U - ln L is off by at most 1e-5
    assert np.max(np.abs(derived - se)) <= 1e-5 / 3.92 * 1.0001


@pytest.mark.parametrize("odds,lower,upper,expected", [
    ("1.295", "1.10", "1.29", "read"),                  # 0.005 above U: within half a unit of OR (0.0005) and of U (0.005)
    ("1.30", "1.10", "1.28", "or_outside_ci"),          # 0.02 above U: beyond 0.005 + 0.005
    ("0.95", "0.955", "1.10", "read"),                  # 0.005 below L: within 0.005 + 0.0005
    ("0.95", "0.97", "1.10", "or_outside_ci"),
    ("1.2000", "1.2001", "1.3", "read"),                # 1e-4 below L: within 5e-5 + 5e-5
    ("1.2000", "1.2002", "1.3", "or_outside_ci"),       # 2e-4 below L: beyond 1e-4
])
def test_an_odds_ratio_outside_its_interval_is_read_only_within_the_rounding_of_the_printed_decimals(odds, lower, upper, expected):
    got = or_ci().row(or_ci_row(odds=odds, lower=lower, upper=upper))
    assert (got if isinstance(got, str) else "read") == expected


@pytest.mark.parametrize("lower,upper,reason", [
    ("#NA", "1.3", "ci_nonpositive_or_nonfinite"), ("1.1", "NA", "ci_nonpositive_or_nonfinite"),
    ("0", "1.3", "ci_nonpositive_or_nonfinite"), ("-0.1", "1.3", "ci_nonpositive_or_nonfinite"),
    ("1.1", "inf", "ci_nonpositive_or_nonfinite"), ("1.1", "-inf", "ci_nonpositive_or_nonfinite"),
    ("1.3", "1.1", "ci_not_increasing"), ("1.2", "1.2", "ci_not_increasing"), ("1.1x", "1.3", "unparseable_number"),
])
def test_malformed_reversed_and_nonpositive_limits_are_rejected_each_with_its_reason(lower, upper, reason):
    assert or_ci().row(or_ci_row(odds="1.2", lower=lower, upper=upper)) == reason
    if reason != "unparseable_number":                                    # a field that is not a number comes first
        assert or_ci().row(or_ci_row(odds="0", lower=lower, upper=upper)) == "or_nonpositive_or_nonfinite"


NEG_LOG10 = [c if c != "p_value" else "neg_log_10_p_value" for c in SSF_BETA_HEADER]


def test_p_zero_is_read_in_native_se_and_in_or_ci_derived_se_and_rejected_in_pvalue_coloc():
    assert ssf().row(ssf_row(p="0"))["p"] == 0.0
    assert or_ci().row(or_ci_row(p="0"))["p"] == 0.0
    assert pvalue().row(ssf_row(p="0", se="NA")) == "p_zero_pvalue_coloc"
    assert pvalue(NEG_LOG10).row(ssf_row(p="400", se="NA")) == "p_zero_pvalue_coloc"      # 10^-400 is 0 as a double
    assert pvalue().row(ssf_row(p="#NA", se="NA")) == "p_missing_pvalue_coloc"
    assert math.isnan(ssf().row(ssf_row(p="#NA"))["p"]) and math.isnan(or_ci().row(or_ci_row(p="#NA"))["p"])
    assert pvalue().row(ssf_row(p="1", se="NA"))["p"] == 1.0 and pvalue().row(ssf_row(p="1.5", se="NA")) == "p_outside_0_1"


@pytest.mark.parametrize("eaf,native,pvalue_mode", [
    ("0.3", 0.3, 0.3), ("0", 0.0, "eaf_not_inside_0_1_pvalue_coloc"), ("1", 1.0, "eaf_not_inside_0_1_pvalue_coloc"),
    ("#NA", "nan", "eaf_not_inside_0_1_pvalue_coloc"), ("1.0001", "eaf_outside_0_1", "eaf_outside_0_1"),
    ("-0.0001", "eaf_outside_0_1", "eaf_outside_0_1"),
])
def test_effect_allele_frequency_bounds_in_each_mode(eaf, native, pvalue_mode):
    for got, want in ((ssf().row(ssf_row(eaf=eaf)), native), (or_ci().row(or_ci_row(eaf=eaf)), native),
                      (pvalue().row(ssf_row(eaf=eaf, se="NA")), pvalue_mode)):
        if want == "nan":
            assert math.isnan(got["eaf"])
        elif isinstance(want, str):
            assert got == want
        else:
            assert got["eaf"] == want


def test_a_pvalue_coloc_row_keeps_beta_for_the_direction_reads_no_se_and_ignores_the_se_field():
    for token in ("NA", "#NA", "0.01", "x"):
        got = pvalue().row(ssf_row(effect="-3.2", se=token, p="0.0014"))
        assert (got["beta"], got["p"], got["eaf"]) == (-3.2, 0.0014, 0.3) and math.isnan(got["se"])
    assert pvalue().row(ssf_row(effect="NA", se="NA")) == "beta_nonfinite"


def test_a_reader_refuses_a_mode_its_header_or_format_cannot_give():
    with pytest.raises(SourceUnreadable, match="or_ci_derived_se needs an odds_ratio effect with ci_lower and ci_upper"):
        or_ci(SSF_BETA_HEADER)
    with pytest.raises(SourceUnreadable, match="or_ci_derived_se needs"):
        or_ci([c for c in SSF_OR_CI_HEADER if c != "ci_lower"])
    with pytest.raises(SourceUnreadable, match="pvalue_coloc needs a beta effect and a p-value column"):
        pvalue(SSF_OR_CI_HEADER)
    with pytest.raises(SourceUnreadable, match="pvalue_coloc needs"):
        pvalue([c for c in SSF_BETA_HEADER if c != "p_value"])
    for mode in ("or_ci_derived_se", "pvalue_coloc"):
        with pytest.raises(SourceUnreadable, match="read in native_se only"):
            AuthorReader(AUTHOR_FORMATS["GCST008226"], list(LASKAR_HEADER), mode)


def test_the_se_census_counts_only_rows_with_a_position_and_a_valid_pair_of_alleles():
    states = [ssf().se_state(r) for r in (ssf_row(se="#NA"), ssf_row(se="NA"), ssf_row(se=""), ssf_row(se=" #NA "),
                                          ssf_row(se="0.01"), ssf_row(se="0"), ssf_row(se="x"), ssf_row(ea="N", se="#NA"),
                                          ssf_row(chrom="NA", se="#NA"), ssf_row()[:-1])]
    assert states == ["missing", "missing", "missing", "missing", "present", "present", "present", "", "", ""]
    assert [laskar().se_state(r) for r in (laskar_row(se="NA"), laskar_row(se="0.05"), laskar_row(ea="T", a1="A", a2="G"))] == [
        "missing", "present", ""]


OR_CI_COLS = ssf_columns(SSF_OR_CI_HEADER, "X")
BETA_COLS = ssf_columns(SSF_BETA_HEADER, "X")


@pytest.mark.parametrize("cols,missing,present,mode,reason", [
    (BETA_COLS, 0, 10, "native_se", ""), (OR_CI_COLS, 0, 10, "native_se", ""), (None, 0, 10, "native_se", ""),
    (OR_CI_COLS, 10, 0, "or_ci_derived_se", ""), (BETA_COLS, 10, 0, "pvalue_coloc", ""),
    (OR_CI_COLS, 9, 1, None, "missing in 9 and present in 1 rows"), (BETA_COLS, 1, 9, None, "missing in 1 and present in 9 rows"),
    (None, 10, 0, None, "missing in every row and the header gives no other route"),
    (ssf_columns(SSF_OR_ONLY_HEADER, "X"), 10, 0, None, "no other route"),                      # an odds ratio without limits
    (ssf_columns(SSF_BETA_HEADER + ["ci_lower", "ci_upper"], "X"), 10, 0, None, "no other route"),   # a beta with limits
    (ssf_columns([c for c in SSF_BETA_HEADER if c != "effect_allele_frequency"], "X"), 10, 0, None, "no other route"),
    (BETA_COLS, 0, 0, None, "no data row with a position and a valid pair of alleles"),
])
def test_the_mode_is_chosen_once_per_file_and_a_file_mixing_missing_and_present_se_fails_closed(cols, missing, present, mode,
                                                                                                    reason):
    got, why = choose_uncertainty_mode(cols, missing, present)
    assert got == mode and reason in why and (why == "") == (mode is not None)


def ci_p(rows) -> dict:
    reader = or_ci()
    check = CiPCheck(reader)
    for r in rows:
        got = reader.row(r)
        assert not isinstance(got, str), got
        check.add(got, r)
    return check.result()


def test_the_ci_versus_p_check_agrees_where_the_published_p_is_the_wald_p_of_the_interval():
    rng = np.random.default_rng()
    # |z| <= 0.2 * 5 / 0.05 = 20 almost surely: every published p is a positive double below 1
    rows = [exact_or_ci_row(b, s) for b, s in zip(rng.normal(0, 0.2, 2000), rng.uniform(0.05, 0.2, 2000))]
    got = ci_p(rows)
    assert (got["rows_checked"], got["rows_agree"], got["passed"]) == (2000, 2000, True)
    assert set(got) == {"rule", "ci_level", "z", "log10_tolerance", "min_agreement_percent", "rows_checked", "rows_agree",
                        "rows_disagree", "rows_not_checkable", "passed"}                    # counts and pass/fail only


@pytest.mark.parametrize("bad,passed", [(5, True), (6, False)])
def test_the_ci_versus_p_check_passes_only_when_at_least_95_percent_of_checkable_rows_agree(bad, passed):
    good = [exact_or_ci_row(0.1 + 0.001 * i, 0.05) for i in range(100 - bad)]
    off = [exact_or_ci_row(0.1, 0.05, p=10 * math.erfc(2 / math.sqrt(2)) / 1.5) for _ in range(bad)]   # 0.82 log10 off
    got = ci_p(good + off)
    assert (got["rows_checked"], got["rows_agree"], got["rows_disagree"], got["passed"]) == (100, 100 - bad, bad, passed)


def test_the_ci_versus_p_check_allows_the_printed_rounding_of_p_where_it_is_looser_than_0_1():
    beta = 1.8119 * 0.1                                       # z = 1.8119: p = 0.0700, 0.155 below 0.1 on the log10 scale
    assert ci_p([or_ci_row(odds=repr(math.exp(beta)), lower=repr(math.exp(beta - 0.196)), upper=repr(math.exp(beta + 0.196)),
                           p="0.1")])["rows_agree"] == 1      # "0.1" may stand for anything in [0.05, 0.15]
    assert ci_p([or_ci_row(odds=repr(math.exp(beta)), lower=repr(math.exp(beta - 0.196)), upper=repr(math.exp(beta + 0.196)),
                           p="0.10")])["rows_agree"] == 0     # "0.10": [0.095, 0.105]


def test_rows_with_p_missing_zero_or_one_are_not_checkable_and_no_checkable_row_fails_the_file():
    got = ci_p([or_ci_row(p="#NA"), or_ci_row(p="0"), or_ci_row(p="1")])
    assert (got["rows_checked"], got["rows_not_checkable"], got["passed"]) == (0, 3, False)
    with pytest.raises(ValueError, match="applies to or_ci_derived_se"):
        CiPCheck(ssf())


def test_log10_of_the_two_sided_p_is_exact_where_erfc_is_a_double_and_its_asymptotic_series_takes_over(monkeypatch):
    zs = (1.0, 5.0, 20.0, 28.0, 35.0)
    exact = [math.log10(math.erfc(z / math.sqrt(2))) for z in zs]
    assert [log10_two_sided_p(z) for z in zs] == pytest.approx(exact, rel=1e-12)
    assert log10_two_sided_p(-5.0) == log10_two_sided_p(5.0)
    monkeypatch.setattr(math, "erfc", lambda x: 0.0)           # force the series, and compare it where erfc is exact
    assert [log10_two_sided_p(z) for z in zs[2:]] == pytest.approx(exact[2:], abs=1e-6)
    beyond = [log10_two_sided_p(z) for z in (40.0, 60.0, 200.0)]
    assert all(math.isfinite(v) for v in beyond) and beyond == sorted(beyond, reverse=True)
