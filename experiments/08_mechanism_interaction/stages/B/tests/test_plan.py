"""`plan` on synthetic copies of the tables it reads: the three supplementary workbooks in their
own sheet layouts, the OpenGWAS listing, the deCODE folder listing, and a stage A output."""
import hashlib
import json
import re
from pathlib import Path

import openpyxl
import pandas as pd
import pytest
from v8_manifest import sha256_file

from stage_b.collect import collect_tasks
from stage_b.plan import OPTIONAL_TABLES, PLAN_TABLES, decode_listing_files, make_plan, pinned_tables
from stage_b.schemas import HYPOTHESIS_INPUT_COLUMNS, InputContractError, InstrumentUnit, SourceFile
from stage_b.sentinels import DECODE_ST02, INTERVAL_ST4, UKBPPP_ST9, sheet_header

ETAG = "0123456789abcdef0123456789abcdef"
LISTING = {"dlTokenValidDays": "30", "directoryName": "", "files": [
    {"Key": "1_1_G1_Protein_one.txt.gz", "Size": 950_000_001, "ETag": f'"{ETAG}"'},
    {"Key": "2_2_G2_a.txt.gz", "Size": 5, "ETag": '"aa"'}, {"Key": "2_2_G2_b.txt.gz", "Size": 6, "ETag": '"bb"'},
    {"Key": "READMEproteomics.txt", "Size": 1966, "ETag": '"cc"'}]}


def workbook(path: Path, sheet: str, header_row: int, rows: list[list]) -> Path:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sheet
    for _ in range(header_row - 1):
        ws.append(["supplementary table title"])
    for row in rows:
        ws.append(row)
    wb.save(path)
    return path


@pytest.fixture
def tables(tmp_path) -> dict[str, Path]:
    root = tmp_path / "exp"
    paths = {name: root / rel for name, rel in PLAN_TABLES.items()}
    for p in paths.values():
        p.parent.mkdir(parents=True, exist_ok=True)
    workbook(paths["ukbppp_st9"], "ST9", 5, [
        ["UKBPPP ProteinID", "rsID", "CHROM", "GENPOS (hg38)", "log10(p) (discovery)", "cis/trans"],
        ["G1:P1:OID1:v1", "rs11", 1, 1000, 80.0, "cis"], ["G1:P1:OID1:v1", "rs12", 1, 2000, 300.0, "trans"],
        ["G4:P4:OID4:v1", "rs41", 4, 4000, 30.0, "cis"]])
    workbook(paths["decode_st02"], "ST02", 3, [
        ["SeqId", "variant", "chr\n(var.)", "pos\n(var.)", "cis/\ntrans", "Rank\n(cond.\nsign.)", "-Log10(P)\n(adj.)"],
        ["1_1", "rs21", "chr1", 1100, "cis", 1, 90.0], ["1_1", "rs22", "chr1", 1200, "cis", 2, 200.0],
        ["2_2", "rs23", "chr2", 2100, "cis", 1, 50.0], ["3_3", "rs24", "chr3", 3100, "cis", 1, 40.0]])
    workbook(paths["interval_st4"], "ST4 - pQTL summary", 5, [
        ["SOMAmer ID", "Target fullname", "Sentinel variant*", "Chr", "Pos", "cis/ trans", "Meta-analysis", None, None],
        [None, None, None, None, None, None, "beta", "SE", "p"],
        ["G5.5.5.3", "Prot five", "rs51", 5, 5000, "cis", 0.5, 0.01, 1e-60]])
    paths["opengwas_gwasinfo"].write_text(json.dumps([{"id": "prot-a-5", "trait": "Prot five"}, {"id": "ieu-a-1", "trait": "x"}]))
    paths["decode_listing"].write_text(json.dumps(LISTING))
    return paths


def hypothesis(hid: str, source: str, assay: str, acc: str, outcome_source: str, **extra) -> dict:
    row = {c: "" for c in HYPOTHESIS_INPUT_COLUMNS}
    row.update({"hypothesis_id": hid, "gene_symbol": f"G{hid}", "gene_ensembl": f"ENSG{assay}", "direction": "decrease",
                "instrument_source": source, "instrument_assay_id": assay, "platform": "Olink" if source == "ukbppp" else "SomaScan",
                "outcome_accession": acc, "outcome_source": outcome_source, "outcome_n_case": 100, "outcome_n_control": 900, **extra})
    return row


@pytest.fixture
def stage_a(tmp_path) -> tuple[Path, Path]:
    out = tmp_path / "A"
    out.mkdir()
    rows = [hypothesis("h1", "decode", "1_1", "GCST1", "gwas_catalog"), hypothesis("h2", "decode", "1_1", "FINNGEN_R12_X", "finngen"),
            hypothesis("h3", "decode", "2_2", "GCST1", "gwas_catalog"), hypothesis("h4", "decode", "3_3", "GCST1", "gwas_catalog"),
            hypothesis("h5", "ukbppp", "OID1", "GCST2", "gwas_catalog"), hypothesis("h6", "interval", "G5.5.5.3", "ieu-b-1", "opengwas")]
    pd.DataFrame(rows).assign(indication_id="EFO_1").to_csv(out / "hypotheses.csv", index=False)
    pd.DataFrame({"outcome_accession": ["GCST1", "GCST2", "FINNGEN_R12_X", "ieu-b-1"], "risk_coded": [True] * 4}).to_csv(
        out / "outcome_trait_coding.tsv", sep="\t", index=False)
    return out / "hypotheses.csv", out / "outcome_trait_coding.tsv"


def test_plan_locates_each_instrument_by_identity_and_records_every_file_it_read(tables, stage_a):
    units_text, plan = make_plan(*stage_a, tables)
    units = {u.unit_key: u for u in map(InstrumentUnit.model_validate_json, units_text.splitlines())}
    decode_key, ukb_key, interval_key = "decode__1_1__ENSG1_1", "ukbppp__OID1__ENSGOID1", "interval__G5.5.5.3__ENSGG5.5.5.3"
    assert set(units) == {decode_key, ukb_key, interval_key}
    decode = units[decode_key]
    assert decode.pqtl_locator == "1_1_G1_Protein_one.txt.gz" and decode.smp_listing is None
    assert decode.pqtl_listing == SourceFile(name="1_1_G1_Protein_one.txt.gz", size=950_000_001, etag=ETAG)
    assert (decode.sentinel.rsid, [o.accession for o in decode.outcomes]) == ("rs21", ["FINNGEN_R12_X", "GCST1"])
    assert (units[ukb_key].pqtl_locator, units[ukb_key].sentinel.rsid) == ("OID1", "rs11")
    assert (units[interval_key].pqtl_locator, units[interval_key].sentinel.build) == ("prot-a-5", "GRCh37")
    assert "http" not in units_text and "token" not in units_text.lower()
    assert plan["hypothesis_unit"] == {"h1": decode_key, "h2": decode_key, "h5": ukb_key, "h6": interval_key}
    assert plan["unit_ids"] == {decode_key: ["decode", "1_1", "ENSG1_1"], ukb_key: ["ukbppp", "OID1", "ENSGOID1"],
                                interval_key: ["interval", "G5.5.5.3", "ENSGG5.5.5.3"]}
    assert set(plan["unresolved"]) == {"h3", "h4"}                    # two files for the SeqId; no file for the SeqId
    assert all(r.startswith("regional_file_unavailable: the deCODE folder listing names no single file") for r in plan["unresolved"].values())
    assert plan["tables"] == {name: sha256_file(path) for name, path in tables.items()} and set(plan["tables"]) == set(PLAN_TABLES)
    assert plan["a_outputs"] == {"hypotheses.csv": sha256_file(stage_a[0]), "outcome_trait_coding.tsv": sha256_file(stage_a[1])}
    assert plan["units_sha256"] == hashlib.sha256(units_text.encode()).hexdigest()
    assert (plan["units"], plan["units_by_source"]) == (3, {"decode": 1, "interval": 1, "ukbppp": 1})
    assert plan["collect_tasks_by_source"] == {"decode": 1, "gwas_catalog": 2, "ukbppp": 1, "ukbppp_rsid_map": 1}
    assert [(t.source, t.key, t.name) for t in collect_tasks(units.values())][0] == ("decode", "1_1", "1_1_G1_Protein_one.txt.gz")
    assert make_plan(*stage_a, tables)[0] == units_text               # the same inputs give the same units


def test_plan_gives_an_assay_that_serves_two_genes_one_unit_per_gene_and_an_assay_placeholder_none(tables, tmp_path):
    out = tmp_path / "A2"
    out.mkdir()
    complex_a = {"gene_symbol": "LA", "gene_ensembl": "ENSG_A"}
    complex_c = {"gene_symbol": "LC", "gene_ensembl": "ENSG_C"}
    rows = [hypothesis("a1", "interval", "G5.5.5.3", "FINNGEN_R12_X", "finngen", **complex_a),
            hypothesis("c1", "interval", "G5.5.5.3", "FINNGEN_R12_X", "finngen", **complex_c),
            hypothesis("c2", "interval", "G5.5.5.3", "GCST1", "gwas_catalog", **complex_c),
            hypothesis("d1", "decode", "1_1", "GCST1", "gwas_catalog", gene_ensembl="ENSG_D"),
            hypothesis("e1", "decode", "1_1", "GCST2", "gwas_catalog", gene_ensembl="ENSG_E"),
            hypothesis("n1", "interval", "epigraphdb:no_assay", "GCST1", "gwas_catalog", gene_ensembl="ENSG_N1"),
            hypothesis("n2", "interval", "epigraphdb:no_assay", "GCST1", "gwas_catalog", gene_ensembl="ENSG_N2")]
    pd.DataFrame(rows).to_csv(out / "hypotheses.csv", index=False)
    pd.DataFrame({"outcome_accession": ["GCST1", "GCST2", "FINNGEN_R12_X"], "risk_coded": [True] * 3}).to_csv(
        out / "outcome_trait_coding.tsv", sep="\t", index=False)
    units_text, plan = make_plan(out / "hypotheses.csv", out / "outcome_trait_coding.tsv", tables)
    units = {u.unit_key: u for u in map(InstrumentUnit.model_validate_json, units_text.splitlines())}
    assert set(units) == {"interval__G5.5.5.3__ENSG_A", "interval__G5.5.5.3__ENSG_C", "decode__1_1__ENSG_D", "decode__1_1__ENSG_E"}
    a, c = units["interval__G5.5.5.3__ENSG_A"], units["interval__G5.5.5.3__ENSG_C"]
    assert (a.gene_ensembl, a.gene_symbol, [o.accession for o in a.outcomes]) == ("ENSG_A", "LA", ["FINNGEN_R12_X"])
    assert (c.gene_ensembl, c.gene_symbol, [o.accession for o in c.outcomes]) == ("ENSG_C", "LC", ["FINNGEN_R12_X", "GCST1"])
    assert (a.sentinel, a.pqtl_locator) == (c.sentinel, c.pqtl_locator) == (a.sentinel, "prot-a-5") and a.sentinel.rsid == "rs51"
    d, e = units["decode__1_1__ENSG_D"], units["decode__1_1__ENSG_E"]
    assert (d.sentinel, d.pqtl_locator, d.pqtl_listing) == (e.sentinel, e.pqtl_locator, e.pqtl_listing)
    assert plan["hypothesis_unit"] == {"a1": "interval__G5.5.5.3__ENSG_A", "c1": "interval__G5.5.5.3__ENSG_C",
                                       "c2": "interval__G5.5.5.3__ENSG_C", "d1": "decode__1_1__ENSG_D", "e1": "decode__1_1__ENSG_E"}
    assert plan["hypothesis_source_units"] == {"d1": {"decode": "decode__1_1__ENSG_D"}, "e1": {"decode": "decode__1_1__ENSG_E"}}
    assert plan["unit_ids"]["interval__G5.5.5.3__ENSG_C"] == ["interval", "G5.5.5.3", "ENSG_C"]
    assert (plan["units"], plan["units_by_source"]) == (4, {"decode": 2, "interval": 2})
    # the two deCODE units read one file: one task, and one task per GWAS Catalog accession; a unit with a GWAS
    # Catalog outcome also reads the UKB-PPP rsID map of its sentinel's chromosome
    assert plan["collect_tasks_by_source"] == {"decode": 1, "gwas_catalog": 2, "ukbppp_rsid_map": 2}
    maps = sorted({("ukbppp_rsid_map", u.sentinel.chrom) for u in (c, d, e)})
    assert len(maps) == 2 and [(t.source, t.key) for t in collect_tasks(units.values())] == [
        ("decode", "1_1"), ("gwas_catalog", "GCST1"), ("gwas_catalog", "GCST2"), *maps]
    # the EpiGraphDB placeholder names no assay: no sentinel, no unit, whatever the number of genes carrying it
    assert set(plan["unresolved"]) == {"n1", "n2"}
    assert all(r == "regional_file_unavailable: no interval sentinel for assays ['epigraphdb:no_assay']"
               for r in plan["unresolved"].values())
    assert all(re.fullmatch(r"[A-Za-z0-9_.-]+", k) for k in units)


def test_the_smp_listing_is_used_only_when_given(tables, stage_a, tmp_path):
    smp = tmp_path / "smp.json"
    smp.write_text(json.dumps({"files": [{"Key": "1_1_G1_Protein_one.txt.gz", "Size": 7, "ETag": '"dd"'}]}))
    units_text, plan = make_plan(*stage_a, {**tables, "decode_smp_listing": smp})
    decode = next(u for u in map(InstrumentUnit.model_validate_json, units_text.splitlines()) if u.source == "decode")
    assert decode.smp_listing == SourceFile(name="1_1_G1_Protein_one.txt.gz", size=7, etag="dd")
    assert plan["collect_tasks_by_source"]["decode_smp"] == 1 and plan["tables"]["decode_smp_listing"] == sha256_file(smp)
    with pytest.raises(InputContractError, match="plan needs the tables"):
        make_plan(*stage_a, {k: v for k, v in tables.items() if k != "decode_listing"})


def test_a_listing_names_one_file_per_seqid_or_none():
    files = decode_listing_files(LISTING)
    assert set(files) == {"1_1"} and files["1_1"].etag == ETAG


def test_tables_are_read_only_when_they_are_the_pinned_files(tables, tmp_path):
    root = tmp_path / "exp"
    manifest = tmp_path / "modal_inputs_manifest.json"

    def pin(extra: dict[str, str] | None = None) -> None:
        files = [{"path": rel, "sha256": sha256_file(root / rel)} for rel in PLAN_TABLES.values()]
        manifest.write_text(json.dumps({"files": files + [{"path": k, "sha256": v} for k, v in (extra or {}).items()]}))

    pin()
    assert pinned_tables(manifest, root) == tables
    smp = root / OPTIONAL_TABLES["decode_smp_listing"]
    smp.write_text(json.dumps({"files": []}))
    pin({OPTIONAL_TABLES["decode_smp_listing"]: sha256_file(smp)})
    assert pinned_tables(manifest, root) == {**tables, "decode_smp_listing": smp}
    tables["decode_listing"].write_text(json.dumps({**LISTING, "files": LISTING["files"][:1]}))
    with pytest.raises(InputContractError, match="differs from the pin"):
        pinned_tables(manifest, root)
    manifest.write_text(json.dumps({"files": [{"path": PLAN_TABLES["ukbppp_st9"], "sha256": "0" * 64}]}))
    with pytest.raises(InputContractError, match="pins no sha256"):
        pinned_tables(manifest, root)


def test_sheet_headers_are_read_without_the_rows(tables):
    for name, (sheet, row, columns) in (("ukbppp_st9", UKBPPP_ST9), ("decode_st02", DECODE_ST02), ("interval_st4", INTERVAL_ST4)):
        assert set(columns) <= set(sheet_header(tables[name], sheet, row))
    assert sheet_header(tables["decode_st02"], "ST02", 3)[2] == "chr (var.)"
