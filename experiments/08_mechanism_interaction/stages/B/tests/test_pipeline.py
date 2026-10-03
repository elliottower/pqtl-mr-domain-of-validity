"""Unit pipeline and assembly on a synthetic region, with a fake fetcher and a stub coloc backend.

The stub returns fixed posteriors, so these tests check wiring: which variants reach coloc,
how results become evidence rows, the checkpoint/resume behavior, and the output contract."""
import json
import platform
import re
import shutil
from importlib.metadata import PackageNotFoundError
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from helpers import table
from v8_manifest import code_sha256, read_inputs, verify_output_dir
from v8_run_guard import PLAN_SHA256

from stage_b.assemble import build_evidence, collect_unit_dir, regional_rows, sha256_file, write_outputs
from stage_b.checkpoint import (FINGERPRINT_NAME, PINNED_SOURCES, PY_PACKAGES, R_VERSIONS, REMOTE_RELEASE, collect_digest,
                                collect_entry, common_tools, package_files, package_sha256, source_pins, tool_versions,
                                unit_fingerprint, unit_tools, verified_source_pins)
from stage_b.coloc_backend import AbfResult, ColocTask, SusieResult
from stage_b.pipeline import DirStore, process_unit, reset_unavailable
from stage_b.schemas import (EVIDENCE_COLUMNS, OUTCOME_BUILD, PRIMARY_P1, PRIMARY_P2, PRIMARY_P12, WINDOW_PRIMARY, CollectError,
                             CollectRecord, CollectTask, HypothesisInput, InputContractError, InstrumentUnit, OutcomeSpec,
                             RetryableSourceError, Sentinel, SourceAbsent, SourceFile, StaleCheckpointError)

INTERFACES = Path(__file__).resolve().parents[2] / "INTERFACES.md"
PACKAGE = Path(__file__).resolve().parents[1] / "stage_b"
PINS = dict(zip(PINNED_SOURCES, ("a" * 64, "b" * 64)))
CODE = package_sha256(PACKAGE)
TOOLS = {"python": "3.12.0", "python:numpy": "2.1.3", "R": "4.4.3", "R:coloc": "5.2.3", "R:susieR": "0.12.35",
         "R:jsonlite": "2.0.0", "bcftools": "1.21", "bcftools:htslib": "1.21", "tabix:htslib": "1.21"}
COMMIT = "0123456789abcdef0123456789abcdef01234567"
COLLECT = "f" * 64          # the collect digest of a unit, where the test collects no file


def fp(u: InstrumentUnit, pins=None, code: str = CODE, plan: str = PLAN_SHA256, tools=None, collect: str = COLLECT) -> str:
    return unit_fingerprint(u, PINS if pins is None else pins, code, plan, TOOLS if tools is None else tools, collect)


def run(u: InstrumentUnit, fetcher, backend, store: DirStore) -> dict:
    """process_unit under the unit's own fingerprint, as modal_stage_b.run_unit calls it."""
    return process_unit(u, fetcher, backend, store, fp(u), TOOLS, COLLECT)
N = 120
RSIDS = [f"rs{i}" for i in range(N)]
POS = [1_000_000 + 1000 * i for i in range(N)]
SENTINEL = Sentinel(source="decode", assay_id="1_1", rsid="rs60", chrom="1", pos=POS[60], build="GRCh38", neg_log10_p=50)


def pqtl_table(sign: float = 1.0):
    return table([{"rsid": r, "pos": p, "ea": "G", "oa": "A", "eaf": 0.3, "beta": sign * 0.2, "n": 35000.0}
                  for r, p in zip(RSIDS, POS)])


def outcome_table(k: int):
    return table([{"rsid": r, "pos": p, "ea": "G", "oa": "A", "eaf": 0.3, "beta": 0.05} for r, p in zip(RSIDS[:k], POS[:k])])


class StubBackend:
    def __init__(self, h4: dict[str, float], lead: str = "rs60", susie_cs: int = 2):
        self.h4, self.lead, self.susie_cs = h4, lead, susie_cs
        self.tasks: list[ColocTask] = []

    def run(self, tasks):
        self.tasks += tasks
        out = {}
        for t in tasks:
            if t.method == "susie":
                out[t.id] = SusieResult(n_cs1=self.susie_cs, n_cs2=1, n_pairs=2, max_pp_h4=0.97)
                continue
            h4 = next(v for k, v in self.h4.items() if t.id.startswith(k))
            lead = self.lead if self.lead in t.d1.snp else t.d1.snp[0]
            snp = {s: (1.0 if s == lead else 0.0) for s in t.d1.snp}
            r = (1 - h4) / 4
            out[t.id] = AbfResult(pp=(r, r, r, r, h4), nsnps=len(t.d1.snp), snp_pp_h4=snp)
        return out


class FakeFetcher:
    def __init__(self, fail_outcomes=(), qtl_error: Exception | None = None):
        self.fail_outcomes = set(fail_outcomes)
        self.qtl_error = qtl_error
        self.calls: dict[str, int] = {}

    def _count(self, k):
        self.calls[k] = self.calls.get(k, 0) + 1

    def positions(self, sentinel):
        return {"GRCh37": None, "GRCh38": sentinel.pos}

    def pqtl_region(self, unit, chrom, center, half_width, smp=False):
        self._count("pqtl_smp" if smp else "pqtl")
        return pqtl_table(-1.0 if smp else 1.0)

    def outcome_region(self, spec, chrom, center, half_width):
        self._count(spec.accession)
        if spec.accession in self.fail_outcomes:
            raise SourceAbsent(f"{spec.accession} gone")
        return outcome_table({"F_ok": N, "F_small": 30}.get(spec.accession, N))

    def outcome_build(self, spec):
        return OUTCOME_BUILD[spec.source]

    def ld_panel(self, chrom, center, half_width):
        self._count("ld")
        rng = np.random.default_rng()
        d = rng.integers(0, 3, size=(N, 500)).astype(float)
        d[61] = d[60]
        meta = pd.DataFrame({"rsid": RSIDS, "ref": "A", "alt": "G", "pos": POS})
        return meta, d

    def vep(self, rsids, build):
        self._count("vep")
        return [{"id": "rs60", "transcript_consequences": [{"gene_id": "ENSG1", "consequence_terms": ["missense_variant"]}]}]

    def qtl_regions(self, gene, chrom, center, half_width):
        self._count("qtl")
        if self.qtl_error is not None:
            raise self.qtl_error
        def q(tid):
            t = outcome_table(N)
            t["molecular_trait_id"], t["gene_id"] = tid, "ENSG1"
            return t
        return {t: {"ge": q("ENSG1"), "leafcutter": q(f"intron_{t}")} for t in ("liver", "blood", "thyroid")}


def unit(outcomes=("F_ok", "F_small", "F_missing")) -> InstrumentUnit:
    return InstrumentUnit(unit_key="decode__1_1__ENSG1", source="decode", assay_id="1_1", gene_symbol="G", gene_ensembl="ENSG1",
                          platform="SomaScan", sentinel=SENTINEL, pqtl_locator="1_1_G_G.txt.gz",
                          pqtl_listing=SourceFile(name="1_1_G_G.txt.gz", size=950_000_000, etag="a" * 32),
                          smp_listing=SourceFile(name="1_1_G_G.txt.gz", size=940_000_000, etag="b" * 32),
                          outcomes=tuple(OutcomeSpec(accession=a, source="finngen", n_case=1000, n_control=9000,
                                                     risk_coded=True) for a in outcomes))


def hyp(hid: str, direction: str, acc: str, assay: str = "1_1") -> HypothesisInput:
    return HypothesisInput(hypothesis_id=hid, gene_symbol="G", gene_ensembl="ENSG1", direction=direction,
                           instrument_source="decode", instrument_assay_id=assay, platform="SomaScan",
                           outcome_accession=acc, outcome_source="finngen", outcome_n_case=1000, outcome_n_control=9000)


H4 = {"primary": 0.9, "s15a": 0.8, "s15b": 0.95, "s15c": 0.85, "s16": 0.9, "sqtl::": 0.9, "eqtl::": 0.1}


def test_unit_to_evidence_rows(tmp_path):
    backend = StubBackend(H4)
    res = run(unit(), FakeFetcher(fail_outcomes={"F_missing"}), backend, DirStore(tmp_path / "u"))
    hyps = [hyp("h1", "decrease", "F_ok"), hyp("h2", "increase", "F_ok"), hyp("h3", "ambiguous", "F_ok"),
            hyp("h4", "decrease", "F_small"), hyp("h5", "decrease", "F_missing"), hyp("h6", "decrease", "F_ok", "9_9")]
    st29 = {("soma", "1_1"): {"N platforms tested": 2, "cis pQTL on both and high correlation (> 0.5)": "Y",
                              "PAV olink": "N", "PAV soma": "N"}}
    rows = {r.hypothesis_id: r for r in build_evidence(
        hyps, {h: "decode__1_1__ENSG1" for h in ("h1", "h2", "h3", "h4", "h5")}, {"decode__1_1__ENSG1": ("decode", "1_1", "ENSG1")},
        {"decode__1_1__ENSG1": res}, st29)}
    h1 = rows["h1"]
    assert (h1.coloc_run, h1.evidence_state, h1.S, h1.E, h1.lead_variant, h1.genetic_direction) == \
        (True, "supportive", 1, pytest.approx(0.9), "rs60", 1)
    assert (h1.s15a_pp_h4, h1.s15b_pp_h4, h1.s15c_pp_h4, h1.s15d_pp_h4, h1.s15g_pp_h4) == (0.8, 0.95, 0.85, 0.9, 0.97)
    assert (h1.n_shared, h1.frac_pqtl_retained, h1.low_coverage, h1.s15f_low_coverage_excluded) == (N, 1.0, False, False)
    assert (h1.protein_altering, h1.platform_concordant, h1.splicing_candidate) == (True, "concordant", True)
    assert h1.s16_evidence_state == "contradictory"          # SMP release flips the pQTL sign in this fixture
    assert (rows["h2"].evidence_state, rows["h2"].E) == ("contradictory", pytest.approx(-0.9))
    assert (rows["h3"].evidence_state, rows["h3"].E) == ("inconclusive", 0.0)
    h4 = rows["h4"]
    assert (h4.coloc_run, h4.not_run_reason, h4.n_shared, h4.low_coverage, h4.pp_h4) == \
        (False, "fewer_than_50_shared", 30, True, "")
    assert h4.frac_pqtl_retained == pytest.approx(30 / N)
    assert (rows["h5"].not_run_reason, rows["h5"].evidence_state) == ("outcome_file_unavailable", "inconclusive")
    assert (rows["h6"].not_run_reason, rows["h6"].platform_concordant) == ("regional_file_unavailable", "untested")
    primary = next(t for t in backend.tasks if t.id == "primary")
    assert len(primary.d1.snp) == N and primary.p12 == 5e-6 and primary.d2.s == pytest.approx(0.1)
    assert {t.p12 for t in backend.tasks if t.id in ("s15a", "s15b")} == {1e-6, 1e-5}
    # the registered settings: p1 = p2 = 1e-4, p12 = 5e-6, +/-500 kb, coloc.abf (coloc.susie for S15g only)
    assert (PRIMARY_P1, PRIMARY_P2, PRIMARY_P12, WINDOW_PRIMARY) == (1e-4, 1e-4, 5e-6, 500_000)
    assert {(t.p1, t.p2) for t in backend.tasks} == {(1e-4, 1e-4)}
    assert (primary.method, {t.method for t in backend.tasks if t.id == "s15g"}) == ("abf", {"susie"})
    assert {t.method for t in backend.tasks if t.id != "s15g"} == {"abf"}


def test_primary_window_is_500kb_and_s15c_is_1mb(tmp_path):
    wide = [f"rs{i}" for i in range(200)]
    far = table([{"rsid": r, "pos": SENTINEL.pos + (i - 100) * 9000, "ea": "G", "oa": "A"} for i, r in enumerate(wide)])

    class WideFetcher(FakeFetcher):
        def pqtl_region(self, unit, chrom, center, half_width, smp=False):
            return far

        def outcome_region(self, spec, chrom, center, half_width):
            return far

    backend = StubBackend(H4)
    run(unit(("F_ok",)), WideFetcher(), backend, DirStore(tmp_path / "u"))
    by_id = {t.id: t for t in backend.tasks}
    inside = [r for i, r in enumerate(wide) if abs((i - 100) * 9000) <= 500_000]
    assert sorted(by_id["primary"].d1.snp) == sorted(inside)
    assert len(inside) == 111
    assert len(by_id["s15c"].d1.snp) == 200


def test_s15g_falls_back_to_abf_without_two_credible_sets(tmp_path):
    res = run(unit(("F_ok",)), FakeFetcher(), StubBackend(H4, susie_cs=1), DirStore(tmp_path / "u"))
    rec = res["outcomes"]["F_ok"]
    assert rec["s15g_pp_h4"] == 0.9 and rec["s15g_method"].startswith("coloc.abf")


def test_checkpoints_resume_without_refetching(tmp_path):
    store = DirStore(tmp_path / "u")
    crashing = FakeFetcher(qtl_error=RuntimeError("container killed"))
    with pytest.raises(RuntimeError):
        run(unit(("F_ok",)), crashing, StubBackend(H4), store)
    assert store.has("pqtl.meta.json") and store.has("coloc__F_ok.json") and not store.has("result.json")
    again = FakeFetcher()
    res = run(unit(("F_ok",)), again, StubBackend(H4), store)
    assert again.calls == {"qtl": 1, "pqtl_smp": 1}
    assert res["splicing"]["splicing_candidate"] is True
    third = FakeFetcher()
    assert run(unit(("F_ok",)), third, StubBackend(H4), store) == res and third.calls == {}


def test_failed_eqtl_catalogue_leaves_splicing_missing(tmp_path):
    res = run(unit(("F_ok",)), FakeFetcher(qtl_error=SourceAbsent("gone")), StubBackend(H4), DirStore(tmp_path / "u"))
    assert res["splicing"] == {"query_failed": True, "detail": "gone", "splicing_candidate": ""}


def test_reset_unavailable_retries_only_the_failed_region(tmp_path):
    store = DirStore(tmp_path / "u")
    run(unit(("F_ok", "F_missing")), FakeFetcher(fail_outcomes={"F_missing"}), StubBackend(H4), store)
    removed = reset_unavailable(store)
    assert set(removed) == {"outcome__F_missing.meta.json", "result.json"}
    f = FakeFetcher()
    res = run(unit(("F_ok", "F_missing")), f, StubBackend(H4), store)
    assert f.calls.get("F_missing") == 1 and "F_ok" not in f.calls and "pqtl" not in f.calls
    assert res["outcomes"]["F_missing"]["coloc_run"] is True


def test_regional_unavailable_unit(tmp_path):
    class NoPqtl(FakeFetcher):
        def pqtl_region(self, unit, chrom, center, half_width, smp=False):
            raise SourceAbsent("deCODE 1_1: HTTP 404")

    res = run(unit(("F_ok",)), NoPqtl(), StubBackend(H4), DirStore(tmp_path / "u"))
    row = build_evidence([hyp("h1", "decrease", "F_ok")], {"h1": "decode__1_1__ENSG1"}, {"decode__1_1__ENSG1": ("decode", "1_1", "ENSG1")},
                         {"decode__1_1__ENSG1": res}, {})[0]
    assert (row.coloc_run, row.not_run_reason, row.S, row.E, row.s16_evidence_state) == \
        (False, "regional_file_unavailable", 0, 0.0, "")


def test_evidence_columns_match_interfaces_md():
    text = INTERFACES.read_text().split("## Stage B")[1].split("## Stage C")[0]
    cols = []
    for line in text.splitlines():
        m = re.match(r"^\| ([A-Za-z0-9_, …]+) \|", line)
        if not m or m.group(1).strip() == "column":
            continue
        for part in [p.strip() for p in m.group(1).split(",")]:
            if part == "pp_h0 … pp_h4":
                cols += [f"pp_h{i}" for i in range(5)]
            elif part == "s15a_pp_h4 … s15g_pp_h4":
                cols += [f"s15{c}_pp_h4" for c in "abcdefg"]
            else:
                cols.append(part)
    assert cols == EVIDENCE_COLUMNS


def test_outputs_and_manifest(tmp_path):
    res = run(unit(("F_ok",)), FakeFetcher(), StubBackend(H4), DirStore(tmp_path / "u"))
    rows = build_evidence([hyp("h1", "decrease", "F_ok")], {"h1": "decode__1_1__ENSG1"}, {"decode__1_1__ENSG1": ("decode", "1_1", "ENSG1")},
                          {"decode__1_1__ENSG1": res}, {})
    src = tmp_path / "in.csv"
    src.write_text("x")
    manifest = write_outputs(tmp_path / "out", rows, [], [src], {"hypotheses": src}, [("", tmp_path)],
                             script_root=tmp_path, run_token="stageb-token-0001", repo_commit=COMMIT, tools=TOOLS)
    ev = pd.read_csv(tmp_path / "out" / "evidence.csv", dtype=str, keep_default_na=False)
    assert list(ev.columns) == EVIDENCE_COLUMNS and ev.loc[0, "evidence_state"] == "supportive"
    m = pd.read_csv(manifest, sep="\t", dtype=str)
    assert m.set_index("path").loc["evidence.csv", "sha256"] == sha256_file(tmp_path / "out" / "evidence.csv")
    assert m.set_index("path").loc["evidence.csv", "rows"] == "1"
    assert set(verify_output_dir(tmp_path / "out", ["evidence.csv"])) == {
        "evidence.csv", "regional_manifest.tsv", "collected_files.tsv", "run_info.json", "INPUTS.tsv"}
    # the output directory holds those tables and MANIFEST.tsv and nothing else: no regional extract, no whole file
    assert sorted(p.name for p in (tmp_path / "out").rglob("*")) == [
        "INPUTS.tsv", "MANIFEST.tsv", "collected_files.tsv", "evidence.csv", "regional_manifest.tsv", "run_info.json"]
    assert sorted(p.name for p in (tmp_path / "u").iterdir() if p.name.endswith((".tsv.gz", ".npz"))) == [
        "ld.npz", "outcome__F_ok.tsv.gz", "pqtl.tsv.gz", "pqtl_smp.tsv.gz"]          # the extracts stay in the unit directory
    inp, code = read_inputs(tmp_path / "out" / "INPUTS.tsv")
    assert (inp.name, inp.path, inp.sha256) == ("hypotheses", "in.csv", sha256_file(src))
    script = m["script_sha256"].iloc[0]
    assert (code.name, code.path, code.sha256, code.release) == (
        "stage_code", "experiments/08_mechanism_interaction/stages/B", script, COMMIT)
    guard = Path(__file__).resolve().parents[2] / "run_guard"
    assert script == code_sha256({"in.csv": src, "run_guard/v8_manifest.py": guard / "v8_manifest.py",
                                  "run_guard/v8_run_guard.py": guard / "v8_run_guard.py"})
    assert json.loads((tmp_path / "out" / "run_info.json").read_text()) == {
        "stage": "B", "run_token": "stageb-token-0001", "repo_commit": COMMIT, "code_sha256": script, "tools": TOOLS}


def test_s17_sentinel_p_is_the_outcome_p_at_the_sentinel(tmp_path):
    res = run(unit(), FakeFetcher(fail_outcomes={"F_missing"}), StubBackend(H4), DirStore(tmp_path / "u"))
    assert res["outcomes"]["F_ok"]["s17_sentinel_p"] == pytest.approx(1e-3)
    assert res["outcomes"]["F_small"]["s17_sentinel_p"] is None      # rs60 is not among the 30 outcome variants
    hyps = [hyp("h1", "decrease", "F_ok"), hyp("h4", "decrease", "F_small"), hyp("h5", "decrease", "F_missing")]
    rows = {r.hypothesis_id: r for r in build_evidence(hyps, {h.hypothesis_id: "decode__1_1__ENSG1" for h in hyps},
                                                       {"decode__1_1__ENSG1": ("decode", "1_1", "ENSG1")}, {"decode__1_1__ENSG1": res}, {})}
    assert (rows["h1"].s17_sentinel_p, rows["h4"].s17_sentinel_p, rows["h5"].s17_sentinel_p) == (pytest.approx(1e-3), "", "")


def ukb_unit() -> InstrumentUnit:
    s = SENTINEL.model_copy(update={"source": "ukbppp", "assay_id": "OID1"})
    return unit(("F_ok",)).model_copy(update={"unit_key": "ukbppp__OID1__ENSG1", "source": "ukbppp", "assay_id": "OID1",
                                              "platform": "Olink", "sentinel": s, "pqtl_locator": "OID1",
                                              "pqtl_listing": None, "smp_listing": None})


def test_per_source_states_apply_the_primary_rule_with_each_instrument(tmp_path):
    dec = run(unit(("F_ok",)), FakeFetcher(), StubBackend(H4), DirStore(tmp_path / "d"))
    ukb = run(ukb_unit(), FakeFetcher(), StubBackend({**H4, "primary": 0.5}), DirStore(tmp_path / "k"))
    results = {"decode__1_1__ENSG1": dec, "ukbppp__OID1__ENSG1": ukb}
    ids = {"decode__1_1__ENSG1": ("decode", "1_1", "ENSG1"), "ukbppp__OID1__ENSG1": ("ukbppp", "OID1", "ENSG1")}
    hyps = [hyp("both", "decrease", "F_ok"), hyp("opp", "increase", "F_ok"), hyp("nofile", "decrease", "F_ok"),
            hyp("decode_only", "decrease", "F_ok")]
    source_units = {"both": {"decode": "decode__1_1__ENSG1", "ukbppp": "ukbppp__OID1__ENSG1"},
                    "opp": {"decode": "decode__1_1__ENSG1", "ukbppp": "ukbppp__OID1__ENSG1"},
                    "nofile": {"decode": "decode__1_1__ENSG1", "ukbppp": ""}, "decode_only": {"decode": "decode__1_1__ENSG1"}}
    rows = {r.hypothesis_id: r for r in build_evidence(hyps, {h.hypothesis_id: "decode__1_1__ENSG1" for h in hyps}, ids,
                                                       results, {}, source_units)}
    for r in rows.values():
        assert r.evidence_state_decode == r.evidence_state          # the selected source reproduces the primary state
    assert (rows["both"].evidence_state, rows["both"].evidence_state_ukbppp) == ("supportive", "inconclusive")
    assert (rows["opp"].evidence_state, rows["opp"].evidence_state_ukbppp) == ("contradictory", "inconclusive")
    assert rows["nofile"].evidence_state_ukbppp == "inconclusive"
    assert rows["decode_only"].evidence_state_ukbppp == ""
    ukb_supportive = run(ukb_unit(), FakeFetcher(), StubBackend(H4), DirStore(tmp_path / "k2"))
    row = build_evidence([hyps[0]], {"both": "decode__1_1__ENSG1"}, ids, {**results, "ukbppp__OID1__ENSG1": ukb_supportive}, {},
                         {"both": source_units["both"]})[0]
    assert row.evidence_state_ukbppp == "supportive"


def test_two_units_of_one_assay_run_the_gene_specific_steps_each_for_its_own_gene(tmp_path):
    class Recording(FakeFetcher):                         # VEP and the eQTL Catalogue hold the gene ENSG1 only
        def qtl_regions(self, gene, chrom, center, half_width):
            self.genes = [*getattr(self, "genes", []), gene]
            return super().qtl_regions(gene, chrom, center, half_width)

    one = unit(("F_ok",))
    two = one.model_copy(update={"unit_key": "decode__1_1__ENSG2", "gene_symbol": "G2", "gene_ensembl": "ENSG2"})
    assert fp(one) != fp(two)
    fetchers = {one.unit_key: Recording(), two.unit_key: Recording()}
    results = {u.unit_key: run(u, fetchers[u.unit_key], StubBackend(H4), DirStore(tmp_path / u.unit_key)) for u in (one, two)}
    r1, r2 = results[one.unit_key], results[two.unit_key]
    assert (fetchers[one.unit_key].genes, fetchers[two.unit_key].genes) == (["ENSG1"], ["ENSG2"])
    assert r1["outcomes"]["F_ok"]["pp"] == r2["outcomes"]["F_ok"]["pp"] and (r1["unit_key"], r2["unit_key"]) == (one.unit_key, two.unit_key)
    assert (r1["vep"]["rs60"]["hit"], r2["vep"]["rs60"]["hit"]) == (True, False)
    assert (r1["splicing"]["splicing_candidate"], r2["splicing"]["splicing_candidate"]) == (True, "")
    assert {t["lead_event"] for t in r2["splicing"]["tissues"].values()} == {None}
    h1, h2 = hyp("h1", "decrease", "F_ok"), hyp("h2", "decrease", "F_ok").model_copy(update={"gene_ensembl": "ENSG2"})
    ids = {one.unit_key: ("decode", "1_1", "ENSG1"), two.unit_key: ("decode", "1_1", "ENSG2")}
    rows = build_evidence([h1, h2], {"h1": one.unit_key, "h2": two.unit_key}, ids, results, {},
                          {"h1": {"decode": one.unit_key}, "h2": {"decode": two.unit_key}})
    assert [(r.evidence_state, r.protein_altering, r.splicing_candidate, r.evidence_state_decode) for r in rows] == [
        ("supportive", True, True, "supportive"), ("supportive", False, "", "supportive")]
    with pytest.raises(InputContractError, match="is an instrument for ENSG1"):       # the other gene's unit is refused
        build_evidence([h2], {"h2": one.unit_key}, ids, results, {})
    with pytest.raises(InputContractError, match="is an instrument for ENSG1"):
        build_evidence([h2], {"h2": two.unit_key}, ids, results, {}, {"h2": {"decode": one.unit_key}})


# ---- checkpoints bound to the run fingerprint ----------------------------------------------------

def files(root: Path) -> dict[str, bytes]:
    return {p.name: p.read_bytes() for p in sorted(root.iterdir())}


def test_the_fingerprint_changes_with_the_unit_the_pins_the_code_and_the_plan():
    u = unit(("F_ok",))
    same = InstrumentUnit.model_validate_json(u.model_dump_json())
    assert fp(u) == fp(same) and len(fp(u)) == 64
    variants = {
        "another outcome list": fp(unit(("F_ok", "F_small"))),
        "another outcome record": fp(u.model_copy(update={"outcomes": (u.outcomes[0].model_copy(update={"n_case": 1001}),)})),
        "another file name": fp(u.model_copy(update={"pqtl_locator": "1_1_G_H.txt.gz"})),
        "another listed size": fp(u.model_copy(update={"pqtl_listing": u.pqtl_listing.model_copy(update={"size": 1})})),
        "another listed ETag": fp(u.model_copy(update={"pqtl_listing": u.pqtl_listing.model_copy(update={"etag": "c" * 32})})),
        "no SMP file": fp(u.model_copy(update={"smp_listing": None})),
        "another sentinel": fp(u.model_copy(update={"sentinel": SENTINEL.model_copy(update={"pos": SENTINEL.pos + 1})})),
        "another pinned file": fp(u, pins={**PINS, PINNED_SOURCES[0]: "c" * 64}),
        "other code": fp(u, code="d" * 64),
        "another plan": fp(u, plan="e" * 64),
        "another coloc": fp(u, tools={**TOOLS, "R:coloc": "5.2.4"}),
        "another R": fp(u, tools={**TOOLS, "R": "4.5.0"}),
        "another htslib": fp(u, tools={**TOOLS, "bcftools:htslib": "1.22"}),
        "another numpy": fp(u, tools={**TOOLS, "python:numpy": "2.1.4"}),
        "a tool more": fp(u, tools={**TOOLS, "python:scipy": "1.0"}),
        "other collected bytes": fp(u, collect="e" * 64),
    }
    assert fp(u, tools=dict(reversed(list(TOOLS.items())))) == fp(u)          # the order of the record does not matter
    with pytest.raises(InputContractError, match="needs the tool versions"):
        fp(u, tools={})
    for digest in ("", "f" * 63, "F" * 64):
        with pytest.raises(InputContractError, match="needs the collect digest"):
            fp(u, collect=digest)
    assert len({fp(u), *variants.values()}) == len(variants) + 1, variants


def test_a_finished_unit_is_not_returned_to_a_run_with_another_unit_record(tmp_path):
    store = DirStore(tmp_path / "u")
    run(unit(("F_ok",)), FakeFetcher(), StubBackend(H4), store)
    before = files(store.root)
    grown, fetcher = unit(("F_ok", "F_small")), FakeFetcher()
    assert grown.unit_key == unit(("F_ok",)).unit_key          # same directory name, another unit record
    with pytest.raises(StaleCheckpointError, match="not resumed"):
        run(grown, fetcher, StubBackend(H4), store)
    assert fetcher.calls == {} and files(store.root) == before


@pytest.mark.parametrize("change", [{"pins": {**PINS, PINNED_SOURCES[1]: "c" * 64}}, {"code": "d" * 64}, {"plan": "e" * 64},
                                    {"tools": {**TOOLS, "R:susieR": "0.12.36"}}, {"tools": {**TOOLS, "bcftools": "1.22"}},
                                    {"collect": "e" * 64}])
def test_a_half_finished_unit_is_not_resumed_under_other_pins_code_plan_tool_versions_or_collected_bytes(tmp_path, change):
    store, u = DirStore(tmp_path / "u"), unit(("F_ok",))
    with pytest.raises(RuntimeError, match="container killed"):
        run(u, FakeFetcher(qtl_error=RuntimeError("container killed")), StubBackend(H4), store)
    assert store.has("coloc__F_ok.json") and not store.has("result.json")
    before, fetcher = files(store.root), FakeFetcher()
    with pytest.raises(StaleCheckpointError):
        process_unit(u, fetcher, StubBackend(H4), store, fp(u, **change), change.get("tools", TOOLS),
                     change.get("collect", COLLECT))
    assert fetcher.calls == {} and files(store.root) == before
    assert run(u, FakeFetcher(), StubBackend(H4), store)["fingerprint"] == fp(u)   # the run that wrote it still resumes


def test_a_checkpoint_directory_written_without_a_fingerprint_is_refused(tmp_path):
    store, u = DirStore(tmp_path / "u"), unit(("F_ok",))
    run(u, FakeFetcher(), StubBackend(H4), store)
    (store.root / FINGERPRINT_NAME).unlink()
    fetcher = FakeFetcher()
    with pytest.raises(StaleCheckpointError, match=f"no {FINGERPRINT_NAME}"):
        run(u, fetcher, StubBackend(H4), store)
    assert fetcher.calls == {}


def test_a_result_carrying_another_fingerprint_is_refused_by_the_unit_and_by_assembly(tmp_path):
    store, u = DirStore(tmp_path / "u"), unit(("F_ok",))
    res = run(u, FakeFetcher(), StubBackend(H4), store)
    assert res["fingerprint"] == fp(u)
    assert json.loads((store.root / FINGERPRINT_NAME).read_text()) == {"fingerprint": fp(u), "tools": TOOLS,
                                                                       "collect_sha256": COLLECT}
    assert res["collect_sha256"] == COLLECT and collect_unit_dir(u, store.root, fp(u), COLLECT)[0] == res
    with pytest.raises(StaleCheckpointError):
        collect_unit_dir(u, store.root, fp(u, code="d" * 64), COLLECT)
    with pytest.raises(StaleCheckpointError):          # the collect records at assembly are not the ones the unit ran on
        collect_unit_dir(u, store.root, fp(u, collect="e" * 64), "e" * 64)
    with pytest.raises(StaleCheckpointError):          # the fingerprint alone is not enough: the digest must be the unit's too
        collect_unit_dir(u, store.root, fp(u), "e" * 64)
    store.put_json("result.json", {**res, "collect_sha256": "e" * 64})
    with pytest.raises(StaleCheckpointError):
        run(u, FakeFetcher(), StubBackend(H4), store)
    with pytest.raises(StaleCheckpointError):
        collect_unit_dir(u, store.root, fp(u), COLLECT)
    store.put_json("result.json", {**res, "fingerprint": "0" * 64})
    with pytest.raises(StaleCheckpointError):
        run(u, FakeFetcher(), StubBackend(H4), store)
    with pytest.raises(StaleCheckpointError):
        collect_unit_dir(u, store.root, fp(u), COLLECT)


def test_source_pins_are_read_from_the_inputs_manifest_and_checked_against_the_files(tmp_path):
    root = tmp_path / "exp"
    (root / "inputs" / "decode").mkdir(parents=True)
    for rel, text in zip(PINNED_SOURCES, ("annotated", "excluded")):
        (root / rel).write_text(text)
    manifest = tmp_path / "modal_inputs_manifest.json"
    manifest.write_text(json.dumps({"files": [{"path": rel, "sha256": sha256_file(root / rel)} for rel in PINNED_SOURCES]
                                    + [{"path": "inputs/ukbppp/other.xlsx", "sha256": "9" * 64}]}))
    pins = {rel: sha256_file(root / rel) for rel in PINNED_SOURCES}
    assert source_pins(manifest) == pins and verified_source_pins(manifest, root) == pins
    (root / PINNED_SOURCES[1]).write_text("excluded, edited")
    with pytest.raises(InputContractError, match="differs from the pin"):
        verified_source_pins(manifest, root)
    manifest.write_text(json.dumps({"files": [{"path": PINNED_SOURCES[0], "sha256": pins[PINNED_SOURCES[0]]}]}))
    with pytest.raises(InputContractError, match="pins no sha256"):
        source_pins(manifest)


# ---- tool versions in the fingerprint ------------------------------------------------------------

VERSION_OUTPUT = {
    "Rscript": "4.4.3\n5.2.3\n0.12.35\n2.0.0\n",
    "bcftools": "bcftools 1.21\nUsing htslib 1.21\nCopyright (C) 2024 Genome Research Ltd.\n",
    "tabix": "tabix (htslib) 1.21\nCopyright (C) 2024 Genome Research Ltd.\n",
}


def fake_tools(outputs=None, missing=()):
    outputs = {**VERSION_OUTPUT, **(outputs or {})}
    calls = []

    def run_cmd(cmd):
        calls.append(cmd)
        if cmd[0] in missing:
            raise InputContractError(f"{cmd[0]} could not be run for its version")
        return outputs[cmd[0]]

    return run_cmd, calls


def test_tool_versions_are_read_from_the_tools_themselves():
    run_cmd, calls = fake_tools()
    got = tool_versions(run_cmd, lambda p: f"{len(p)}.0")
    assert got == {"python": got["python"], **{f"python:{p}": f"{len(p)}.0" for p in PY_PACKAGES}, "R": "4.4.3",
                   "R:coloc": "5.2.3", "R:susieR": "0.12.35", "R:jsonlite": "2.0.0", "bcftools": "1.21",
                   "bcftools:htslib": "1.21", "tabix:htslib": "1.21"}
    assert re.fullmatch(r"3\.\d+\.\d+\S*", got["python"])
    assert calls == [["Rscript", "-e", R_VERSIONS], ["bcftools", "--version"], ["tabix", "--version"]]
    for name in ("getRversion()", 'packageVersion("coloc")', 'packageVersion("susieR")', 'packageVersion("jsonlite")'):
        assert name in R_VERSIONS
    assert R_VERSIONS.index("coloc") < R_VERSIONS.index("susieR") < R_VERSIONS.index("jsonlite")


def test_tool_versions_follow_what_the_tools_report():
    run_cmd, _ = fake_tools({"Rscript": "4.5.0\n6.0.0\n0.14.2\n2.0.0\n",
                             "bcftools": "bcftools 1.22\nUsing htslib 1.22.1\n", "tabix": "tabix (htslib) 1.20\n"})
    got = tool_versions(run_cmd, lambda p: "9")
    assert (got["R"], got["R:coloc"], got["R:susieR"], got["bcftools"], got["bcftools:htslib"], got["tabix:htslib"]) == (
        "4.5.0", "6.0.0", "0.14.2", "1.22", "1.22.1", "1.20")


@pytest.mark.parametrize("outputs,match", [
    ({"Rscript": "4.4.3\n5.2.3\n"}, "Rscript reported"),
    ({"bcftools": "bcftools 1.21\n"}, "bcftools --version"),
    ({"bcftools": "samtools 1.21\nUsing htslib 1.21\n"}, "bcftools --version"),
    ({"tabix": ""}, "tabix --version"),
    ({"tabix": "Version: 1.21\n"}, "tabix --version"),
])
def test_tool_versions_refuse_output_they_cannot_read(outputs, match):
    with pytest.raises(InputContractError, match=match):
        tool_versions(fake_tools(outputs)[0], lambda p: "1")


@pytest.mark.parametrize("tool", ["Rscript", "bcftools", "tabix"])
def test_a_missing_tool_is_an_error_not_an_omitted_version(tool):
    with pytest.raises(InputContractError, match=tool):
        tool_versions(fake_tools(missing={tool})[0], lambda p: "1")


def test_a_python_package_that_is_not_installed_is_an_error():
    def absent(p):
        raise PackageNotFoundError(p)

    with pytest.raises(InputContractError, match="numpy is not installed"):
        tool_versions(fake_tools()[0], absent)


IN_IMAGE = all(shutil.which(t) for t in ("Rscript", "bcftools", "tabix"))


@pytest.mark.skipif(not IN_IMAGE, reason="needs R, bcftools and tabix: runs in the stage B image")
def test_in_the_image_the_tools_report_the_versions_the_image_pins():
    src = (PACKAGE.parent / "modal_stage_b.py").read_text()
    pinned = dict(re.findall(r'"([A-Za-z0-9_-]+)==([^"]+)"', src.split("PY_PACKAGES = [")[1].split("]")[0]))
    r_version = re.search(r'R_VERSION = "([^"]+)"', src).group(1)
    htslib = re.search(r'HTSLIB_VERSION = "([^"]+)"', src).group(1)
    r_packages = dict(re.findall(r'"(\w+)": "([^"]+)"', src.split("R_PACKAGES = {")[1].split("}")[0]))
    got = tool_versions()
    assert got == {"python": platform.python_version(), **{f"python:{p}": pinned[p] for p in PY_PACKAGES}, "R": r_version,
                   "R:coloc": r_packages["coloc"], "R:susieR": r_packages["susieR"], "R:jsonlite": r_packages["jsonlite"],
                   "bcftools": htslib, "bcftools:htslib": htslib, "tabix:htslib": htslib}
    assert len(fp(unit(("F_ok",)), tools=got)) == 64 and fp(unit(("F_ok",)), tools=got) != fp(unit(("F_ok",)))


def test_the_recorded_python_packages_are_the_ones_the_image_pins():
    src = (PACKAGE.parent / "modal_stage_b.py").read_text()
    pins = src.split("PY_PACKAGES = [")[1].split("]")[0]
    assert sorted(re.findall(r'"([A-Za-z0-9_-]+)==', pins)) == sorted(PY_PACKAGES)
    assert "tools = tool_versions()" in src
    for tool in ("coloc", "susieR", "jsonlite"):
        assert f'"{tool}"' in src


def test_assembly_recomputes_the_fingerprint_from_the_tool_versions_the_unit_recorded(tmp_path):
    store, u = DirStore(tmp_path / "u"), unit(("F_ok",))
    res = run(u, FakeFetcher(), StubBackend(H4), store)
    tools = unit_tools(store.root)
    assert tools == TOOLS
    assert collect_unit_dir(u, store.root, unit_fingerprint(u, PINS, CODE, PLAN_SHA256, tools, COLLECT), COLLECT)[0] == res
    # a tool version edited in the unit directory no longer matches the fingerprint the unit carries
    record = json.loads((store.root / FINGERPRINT_NAME).read_text())
    record["tools"]["R:coloc"] = "5.2.4"
    (store.root / FINGERPRINT_NAME).write_text(json.dumps(record))
    with pytest.raises(StaleCheckpointError):
        collect_unit_dir(u, store.root, unit_fingerprint(u, PINS, CODE, PLAN_SHA256, unit_tools(store.root), COLLECT), COLLECT)
    (store.root / FINGERPRINT_NAME).write_text(json.dumps({"fingerprint": fp(u)}))
    with pytest.raises(StaleCheckpointError, match="records no tool versions"):
        unit_tools(store.root)
    with pytest.raises(StaleCheckpointError, match="records no tool versions"):
        unit_tools(tmp_path / "absent")


def test_assembly_refuses_units_run_under_different_tool_versions():
    assert common_tools({"u1": TOOLS, "u2": dict(reversed(list(TOOLS.items())))}) == TOOLS
    with pytest.raises(StaleCheckpointError, match="2 different sets of tool versions"):
        common_tools({"u1": TOOLS, "u2": {**TOOLS, "R:coloc": "5.2.4"}, "u3": TOOLS})
    with pytest.raises(StaleCheckpointError):
        common_tools({})


def test_the_code_digest_covers_every_module_the_r_script_and_the_shared_guard_by_relative_path():
    files = package_files(PACKAGE)
    assert {"stage_b/pipeline.py", "stage_b/checkpoint.py", "stage_b/coloc_run.R", "run_guard/v8_manifest.py",
            "run_guard/v8_run_guard.py"} <= set(files)
    assert all(k.startswith(("stage_b/", "run_guard/")) for k in files)
    assert {k for k in files if k.startswith("stage_b/") and k.endswith(".py")} == {
        f"stage_b/{p.name}" for p in PACKAGE.glob("*.py")}
    assert CODE == code_sha256(files)


# ---- a resumed unit gives the files a fresh run gives ---------------------------------------------

class SteadyFetcher(FakeFetcher):
    """FakeFetcher with a fixed LD panel, so two runs of a unit see the same reference data."""

    def ld_panel(self, chrom, center, half_width):
        self._count("ld")
        d = ((np.arange(N)[:, None] * 7 + np.arange(500)[None, :] * (np.arange(N)[:, None] % 5 + 1)) % 3).astype(float)
        d[61] = d[60]
        return pd.DataFrame({"rsid": RSIDS, "ref": "A", "alt": "G", "pos": POS}), d


class Killed(RuntimeError):
    pass


class DiesAt(SteadyFetcher):
    """A fetcher whose container is killed the first time step `at` is reached."""

    def __init__(self, at: str):
        super().__init__(fail_outcomes={"F_missing"})
        self.at = at

    def _count(self, k):
        if k == self.at:
            raise Killed(k)
        super()._count(k)


TIMED = re.compile(r"\.meta\.json$")


def unit_files(root: Path) -> dict[str, bytes]:
    """Every file of a unit directory, byte for byte, except the two places that record when a step
    ran: the `retrieved_utc` of a regional meta file is dropped, and ld.npz (a zip archive, whose
    member headers carry the time of writing) is compared by the arrays it holds."""
    out = {}
    for p in sorted(root.iterdir()):
        if TIMED.search(p.name):
            meta = json.loads(p.read_text())
            assert set(meta) == {"status", "variants", "sha256", "window", "detail", "retrieved_utc"}
            meta.pop("retrieved_utc")
            out[p.name] = json.dumps(meta, sort_keys=True).encode()
        elif p.suffix == ".npz":
            with np.load(p, allow_pickle=False) as z:
                out[p.name] = b"".join(f"{k}:{z[k].dtype}:{z[k].shape}:".encode() + z[k].tobytes() for k in sorted(z.files))
        else:
            out[p.name] = p.read_bytes()
    return out


def tables(tmp: Path, u: InstrumentUnit, unit_dir: Path) -> dict[str, bytes]:
    """evidence.csv and regional_manifest.tsv as assembly writes them from one unit directory."""
    res, metas = collect_unit_dir(u, unit_dir, fp(u), COLLECT)
    hyps = [hyp("h1", "decrease", "F_ok"), hyp("h2", "increase", "F_ok"), hyp("h4", "decrease", "F_small"),
            hyp("h5", "decrease", "F_missing")]
    rows = build_evidence(hyps, {h.hypothesis_id: u.unit_key for h in hyps}, {u.unit_key: (u.source, u.assay_id, u.gene_ensembl)},
                          {u.unit_key: res}, {}, {"h1": {"decode": u.unit_key}})
    for m in metas:
        m["retrieved_utc"] = "fixed"
    src = tmp / "in.csv"
    src.write_text("x")
    write_outputs(tmp / "out", rows, regional_rows({u.unit_key: metas}), [src], {"hypotheses": src}, [("", tmp)],
                  script_root=tmp, run_token="stageb-token-0001", repo_commit=COMMIT, tools=TOOLS)
    return {n: (tmp / "out" / n).read_bytes() for n in ("evidence.csv", "regional_manifest.tsv")}


@pytest.mark.parametrize("at", ["pqtl", "ld", "F_ok", "F_small", "vep", "qtl", "pqtl_smp"])
def test_a_unit_killed_at_any_step_and_resumed_gives_the_files_and_tables_of_a_fresh_run(tmp_path, at):
    u = unit()
    fresh = DirStore(tmp_path / "fresh" / "u")
    done = run(u, SteadyFetcher(fail_outcomes={"F_missing"}), StubBackend(H4), fresh)
    resumed = DirStore(tmp_path / "resumed" / "u")
    with pytest.raises(Killed, match=at):
        run(u, DiesAt(at), StubBackend(H4), resumed)
    assert not resumed.has("result.json")
    assert run(u, SteadyFetcher(fail_outcomes={"F_missing"}), StubBackend(H4), resumed) == done
    assert unit_files(resumed.root) == unit_files(fresh.root)
    assert (resumed.root / "result.json").read_bytes() == (fresh.root / "result.json").read_bytes()
    assert tables(tmp_path / "resumed", u, resumed.root) == tables(tmp_path / "fresh", u, fresh.root)


def test_a_unit_killed_twice_and_resumed_twice_gives_the_tables_of_a_fresh_run(tmp_path):
    u = unit()
    fresh = DirStore(tmp_path / "fresh" / "u")
    run(u, SteadyFetcher(fail_outcomes={"F_missing"}), StubBackend(H4), fresh)
    resumed = DirStore(tmp_path / "resumed" / "u")
    for at in ("F_small", "pqtl_smp"):
        with pytest.raises(Killed, match=at):
            run(u, DiesAt(at), StubBackend(H4), resumed)
    run(u, SteadyFetcher(fail_outcomes={"F_missing"}), StubBackend(H4), resumed)
    assert unit_files(resumed.root) == unit_files(fresh.root)
    assert tables(tmp_path / "resumed", u, resumed.root) == tables(tmp_path / "fresh", u, fresh.root)


# ---- only a definitive absence is recorded; every other fault leaves the unit unfinished ----------

class Faulty(SteadyFetcher):
    """A fetcher that raises `error` the first time step `at` is reached."""

    def __init__(self, at: str, error: Exception):
        super().__init__(fail_outcomes={"F_missing"})
        self.at, self.error = at, error

    def _count(self, k):
        if k == self.at:
            raise self.error
        super()._count(k)


FAULTS = [RetryableSourceError("auth", "deCODE 1_1: HTTP 403"), RetryableSourceError("rate_limit", "OpenGWAS x: HTTP 429"),
          RetryableSourceError("server", "x: HTTP 503"), RetryableSourceError("timeout", "x: ReadTimeout"),
          RetryableSourceError("connection", "x: ChunkedEncodingError"), RetryableSourceError("corrupt", "x: EOFError"),
          CollectError("no collect record for decode 1_1")]


@pytest.mark.parametrize("error", FAULTS, ids=lambda e: getattr(e, "kind", type(e).__name__))
@pytest.mark.parametrize("at", ["pqtl", "ld", "F_ok", "vep", "qtl", "pqtl_smp"])
def test_a_fault_that_is_not_an_absence_records_nothing_and_the_resumed_unit_equals_a_fresh_one(tmp_path, at, error):
    u = unit()
    fresh = DirStore(tmp_path / "fresh" / "u")
    done = run(u, SteadyFetcher(fail_outcomes={"F_missing"}), StubBackend(H4), fresh)
    store = DirStore(tmp_path / "faulty" / "u")
    with pytest.raises(type(error)) as raised:
        run(u, Faulty(at, error), StubBackend(H4), store)
    assert raised.value is error and not store.has("result.json")
    step = {"pqtl": "pqtl.meta.json", "ld": "ld.npz", "F_ok": "outcome__F_ok.meta.json", "vep": "vep__rs60.json",
            "qtl": "splicing.json", "pqtl_smp": "pqtl_smp.meta.json"}[at]
    assert not store.has(step)                                   # the step that failed left no file at all
    unavailable = sorted(p.name for p in store.root.glob("*.meta.json") if json.loads(p.read_text())["status"] != "ok")
    assert unavailable in ([], ["outcome__F_missing.meta.json"])  # only the outcome the source does not hold
    assert run(u, SteadyFetcher(fail_outcomes={"F_missing"}), StubBackend(H4), store) == done
    assert unit_files(store.root) == unit_files(fresh.root)


def test_a_unit_record_cannot_hold_a_url_and_the_fingerprint_reads_the_record_only():
    u = unit(("F_ok",))
    for update in ({"pqtl_locator": "https://download.example/folder/token/1_1_G_G.txt.gz"},
                   {"pqtl_locator": "1_1_G_G.txt.gz?token=abc"},
                   {"pqtl_listing": {"name": "folder/token/1_1_G_G.txt.gz"}}, {"smp_listing": {"name": "https://x/y.txt.gz"}}):
        with pytest.raises(ValueError, match="not a URL"):
            InstrumentUnit(**{**u.model_dump(), **update})
    record = json.dumps(u.model_dump(mode="json"))
    assert "http" not in record and "token" not in record
    assert set(InstrumentUnit.model_fields) == {"unit_key", "source", "assay_id", "gene_symbol", "gene_ensembl", "platform",
                                                "sentinel", "pqtl_locator", "pqtl_listing", "smp_listing", "outcomes"}


def test_collected_files_are_listed_in_the_outputs(tmp_path):
    res = run(unit(("F_ok",)), FakeFetcher(), StubBackend(H4), DirStore(tmp_path / "u"))
    rows = build_evidence([hyp("h1", "decrease", "F_ok")], {"h1": "decode__1_1__ENSG1"}, {"decode__1_1__ENSG1": ("decode", "1_1", "ENSG1")},
                          {"decode__1_1__ENSG1": res}, {})
    src = tmp_path / "in.csv"
    src.write_text("x")
    records = [CollectRecord(status="collected", source="decode", key="1_1", name="1_1_G_G.txt.gz", path="raw/decode/1_1_G_G.txt.gz",
                             bytes=5, sha256="d" * 64, etag="a" * 32,
                             source_url="https://download.decode.is/s3/download?token=<token>&file=1_1_G_G.txt.gz",
                             utc="2026-10-02T00:00:00+00:00"),
               CollectRecord(status="absent", source="gwas_catalog", key="GCST1", detail="GWAS Catalog GCST1 listing: HTTP 404"),
               CollectRecord(status="remote_indexed", source="gwas_catalog", key="GCST2", name="2-GCST2-EFO_2.h.tsv.gz",
                             source_url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST000001-GCST001000/"
                                        "GCST2/harmonised/2-GCST2-EFO_2.h.tsv.gz"),
               CollectRecord(status="collected", source="ukbppp", key="OID1", name="G_P_OID1_v1.tar", path="raw/ukbppp/G_P_OID1_v1.tar",
                             bytes=7, sha256="e" * 64, md5="b" * 32, source_url="synapse:syn900001"),
               CollectRecord(status="collected", source="decode_smp", key="1_1", name="1_1_G_G.txt.gz",
                             path="raw/decode_smp/1_1_G_G.txt.gz", bytes=6, sha256="f" * 64,
                             source_url="http://127.0.0.1:8765/decode_smp/s3/download?token=<token>&file=1_1_G_G.txt.gz")]
    write_outputs(tmp_path / "out", rows, [], [src], {"hypotheses": src}, [("", tmp_path)], script_root=tmp_path,
                  run_token="stageb-token-0001", repo_commit=COMMIT, tools=TOOLS, collected=records)
    text = (tmp_path / "out" / "collected_files.tsv").read_text()
    table = pd.read_csv(tmp_path / "out" / "collected_files.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert table[["source", "key", "status", "bytes", "sha256", "source_url"]].values.tolist() == [
        ["decode", "1_1", "collected", "5", "d" * 64, "download.decode.is"],
        ["decode_smp", "1_1", "collected", "6", "f" * 64, "127.0.0.1"],
        ["gwas_catalog", "GCST1", "absent", "", "", ""],
        ["gwas_catalog", "GCST2", "remote_indexed", "", "", "ftp.ebi.ac.uk"],
        ["ukbppp", "OID1", "collected", "7", "e" * 64, "synapse"]]
    # the table publishes the host only: no scheme, port, path, query, token placeholder or Synapse entity
    assert not any(part in text for part in ("://", "/s3/", "token", "?", "8765", "harmonised", "syn900001"))
    assert records[0].source_url.startswith("https://download.decode.is/s3/download?")   # the record keeps its own
    assert "collected_files.tsv" in verify_output_dir(tmp_path / "out", ["evidence.csv"])


def test_a_variant_ensembl_does_not_know_leaves_protein_altering_missing_unless_st29_flags_it(tmp_path):
    class NoVep(FakeFetcher):
        def vep(self, rsids, build):
            raise SourceAbsent("VEP GRCh38: HTTP 400, the id is not known to Ensembl")

    res = run(unit(("F_ok",)), NoVep(), StubBackend(H4), DirStore(tmp_path / "u"))
    assert res["vep"]["rs60"]["hit"] is None and res["outcomes"]["F_ok"]["coloc_run"] is True
    args = ([hyp("h1", "decrease", "F_ok")], {"h1": "decode__1_1__ENSG1"}, {"decode__1_1__ENSG1": ("decode", "1_1", "ENSG1")}, {"decode__1_1__ENSG1": res})
    assert build_evidence(*args, {})[0].protein_altering == ""
    flagged = {("soma", "1_1"): {"N platforms tested": 1, "cis pQTL on both and high correlation (> 0.5)": "N",
                                 "PAV olink": "N", "PAV soma": "Y"}}
    assert build_evidence(*args, flagged)[0].protein_altering is True


# ---- the collect digest: the fingerprint covers the bytes the collect phase recorded -----------------

def collected(source: str, key: str, name: str, sha: str, **extra) -> CollectRecord:
    return CollectRecord(status="collected", source=source, key=key, name=name, path=f"raw/{source}/{name}", bytes=950_000_000,
                         sha256=sha, **extra)


def catalog_unit() -> InstrumentUnit:
    outcomes = tuple(OutcomeSpec(accession=a, source="gwas_catalog", n_case=1000, n_control=9000, risk_coded=True)
                     for a in ("GCST1", "GCST2", "GCST3"))
    return unit().model_copy(update={"outcomes": outcomes})


def unit_records(u: InstrumentUnit, decode_sha: str = "1" * 64) -> dict[tuple[str, str], CollectRecord]:
    """One record per whole file of `catalog_unit` (its GWAS Catalog outcomes add the rsID map of its
    chromosome): three collected, two absent, one queried by region."""
    return {("decode", "1_1"): collected("decode", "1_1", u.pqtl_listing.name, decode_sha, etag=u.pqtl_listing.etag),
            ("decode_smp", "1_1"): CollectRecord(status="absent", source="decode_smp", key="1_1", detail="decode_smp 1_1: HTTP 404"),
            ("gwas_catalog", "GCST1"): collected("gwas_catalog", "GCST1", "1-GCST1-EFO_1.h.tsv.gz", "2" * 64),
            ("gwas_catalog", "GCST2"): CollectRecord(status="remote_indexed", source="gwas_catalog", key="GCST2",
                                                     name="2-GCST2-EFO_2.h.tsv.gz"),
            ("gwas_catalog", "GCST3"): CollectRecord(status="absent", source="gwas_catalog", key="GCST3",
                                                     detail="GWAS Catalog GCST3 listing: HTTP 404"),
            ("ukbppp_rsid_map", "1"): collected("ukbppp_rsid_map", "1", "olink_rsid_map_chr1.tsv.gz", "5" * 64)}


def digest(u: InstrumentUnit, records: dict[tuple[str, str], CollectRecord]) -> str:
    return collect_digest(u, lambda source, key: records[(source, key)])


def test_the_collect_digest_covers_identity_status_size_and_sha256_and_no_address():
    u = catalog_unit()
    records = unit_records(u)
    entries = {(s, k): collect_entry(CollectTask(source=s, key=k, name=u.pqtl_listing.name if s.startswith("decode") else ""), r)
               for (s, k), r in records.items()}
    assert entries[("decode", "1_1")] == {"source": "decode", "key": "1_1", "status": "collected", "name": "1_1_G_G.txt.gz",
                                          "bytes": 950_000_000, "sha256": "1" * 64}
    assert entries[("decode_smp", "1_1")] == {"source": "decode_smp", "key": "1_1", "status": "absent",
                                              "name": "1_1_G_G.txt.gz", "bytes": None, "sha256": ""}
    assert entries[("gwas_catalog", "GCST2")] == {"source": "gwas_catalog", "key": "GCST2", "status": "remote_indexed",
                                                  "name": "2-GCST2-EFO_2.h.tsv.gz", "bytes": None, "sha256": "",
                                                  "release": REMOTE_RELEASE["gwas_catalog"]}
    assert entries[("gwas_catalog", "GCST3")] == {"source": "gwas_catalog", "key": "GCST3", "status": "absent", "name": "",
                                                  "bytes": None, "sha256": ""}
    assert REMOTE_RELEASE == {"gwas_catalog": "r2026-09-13"}
    base = digest(u, records)
    assert len(base) == 64 and base == digest(u, dict(reversed(list(records.items()))))
    # what is not identity, status, size or sha256 leaves the digest as it was: an address, a time, a detail
    decode = records[("decode", "1_1")]
    moved = decode.model_copy(update={"source_url": "https://download.example/s3/download?token=<token>&file=x",
                                      "utc": "2026-10-09T00:00:00+00:00", "last_modified": "another day", "detail": "note"})
    assert digest(u, {**records, ("decode", "1_1"): moved}) == base
    changes = {
        "other bytes, same name size and ETag": decode.model_copy(update={"sha256": "3" * 64}),
        "another size": decode.model_copy(update={"bytes": 950_000_001}),
        "absent where it was collected": CollectRecord(status="absent", source="decode", key="1_1"),
    }
    digests = {name: digest(u, {**records, ("decode", "1_1"): r}) for name, r in changes.items()}
    digests["a file now collected that was absent"] = digest(u, {**records, ("gwas_catalog", "GCST3"): collected(
        "gwas_catalog", "GCST3", "3-GCST3-EFO_3.h.tsv.gz", "4" * 64)})
    digests["an indexed file under another name"] = digest(u, {**records, ("gwas_catalog", "GCST2"): records[
        ("gwas_catalog", "GCST2")].model_copy(update={"name": "2-GCST2-EFO_9.h.tsv.gz"})})
    assert len({base, *digests.values()}) == len(digests) + 1, digests
    for bad, match in ((decode.model_copy(update={"sha256": ""}), "lacks its name, size or sha256"),
                       (decode.model_copy(update={"bytes": None}), "lacks its name, size or sha256"),
                       (decode.model_copy(update={"name": "1_1_G_H.txt.gz"}), "the pinned listing"),
                       (records[("gwas_catalog", "GCST1")], "was given for decode 1_1")):
        with pytest.raises(CollectError, match=match):
            digest(u, {**records, ("decode", "1_1"): bad})
    with pytest.raises(KeyError):                                  # a file without a record gives no digest
        digest(u, {k: r for k, r in records.items() if k != ("gwas_catalog", "GCST1")})
    assert "http" not in json.dumps(list(entries.values())) and "token" not in json.dumps(list(entries.values()))


def test_a_checkpoint_is_refused_when_the_file_was_collected_again_with_other_bytes_under_the_same_name_size_and_etag(tmp_path):
    u = catalog_unit()
    first = unit_records(u)
    again = {**first, ("decode", "1_1"): first[("decode", "1_1")].model_copy(update={"sha256": "9" * 64})}
    same = {f: getattr(first[("decode", "1_1")], f) for f in ("name", "bytes", "etag")}
    assert same == {f: getattr(again[("decode", "1_1")], f) for f in same} and same["etag"] == u.pqtl_listing.etag
    old, new = digest(u, first), digest(u, again)
    assert old != new and fp(u, collect=old) != fp(u, collect=new)
    store = DirStore(tmp_path / "u")
    done = process_unit(u, FakeFetcher(fail_outcomes={"GCST3"}), StubBackend(H4), store, fp(u, collect=old), TOOLS, old)
    assert done["collect_sha256"] == old and json.loads((store.root / FINGERPRINT_NAME).read_text())["collect_sha256"] == old
    before, fetcher = files(store.root), FakeFetcher()
    with pytest.raises(StaleCheckpointError, match="not resumed"):     # the unit record, and so name, size and ETag, is unchanged
        process_unit(u, fetcher, StubBackend(H4), store, fp(u, collect=new), TOOLS, new)
    assert fetcher.calls == {} and files(store.root) == before
    with pytest.raises(StaleCheckpointError):                          # and assembly, which recomputes the digest, refuses it
        collect_unit_dir(u, store.root, fp(u, collect=new), new)
    assert collect_unit_dir(u, store.root, fp(u, collect=old), old)[0] == done
    assert process_unit(u, FakeFetcher(), StubBackend(H4), store, fp(u, collect=old), TOOLS, old) == done


# ---- S16 without a pinned SMP-normalized listing: nothing fetched, the column empty ---------------------

def test_without_an_smp_listing_s16_fetches_nothing_and_its_column_is_empty_for_every_row(tmp_path):
    decode = unit(("F_ok", "F_small", "F_missing")).model_copy(update={"smp_listing": None})
    fetcher, backend = FakeFetcher(fail_outcomes={"F_missing"}), StubBackend(H4)
    res = run(decode, fetcher, backend, DirStore(tmp_path / "d"))
    assert "pqtl_smp" not in fetcher.calls and res["s16"] == {} and not any(t.id == "s16" for t in backend.tasks)
    assert not list((tmp_path / "d").glob("pqtl_smp*")) and not list((tmp_path / "d").glob("s16__*"))
    ukb = run(ukb_unit(), FakeFetcher(), StubBackend(H4), DirStore(tmp_path / "k"))
    hyps = [hyp("h1", "decrease", "F_ok"), hyp("h4", "decrease", "F_small"), hyp("h5", "decrease", "F_missing"),
            hyp("h6", "decrease", "F_ok", "9_9"),
            hyp("k1", "decrease", "F_ok", "OID1").model_copy(update={"instrument_source": "ukbppp", "platform": "Olink"})]
    rows = build_evidence(hyps, {"h1": "decode__1_1__ENSG1", "h4": "decode__1_1__ENSG1", "h5": "decode__1_1__ENSG1", "k1": "ukbppp__OID1__ENSG1"},
                          {"decode__1_1__ENSG1": ("decode", "1_1", "ENSG1"), "ukbppp__OID1__ENSG1": ("ukbppp", "OID1", "ENSG1")},
                          {"decode__1_1__ENSG1": res, "ukbppp__OID1__ENSG1": ukb}, {})
    assert [r.s16_evidence_state for r in rows] == [""] * 5
    assert [r.evidence_state for r in rows] == ["supportive", "inconclusive", "inconclusive", "inconclusive", "supportive"]


def test_s16_state_is_set_only_for_a_decode_instrument_with_an_smp_colocalization(tmp_path):
    dec = run(unit(("F_ok",)), FakeFetcher(), StubBackend(H4), DirStore(tmp_path / "d"))
    ukb = run(ukb_unit(), FakeFetcher(), StubBackend(H4), DirStore(tmp_path / "k"))
    hyps = [hyp("d1", "decrease", "F_ok"),
            hyp("k1", "decrease", "F_ok", "OID1").model_copy(update={"instrument_source": "ukbppp", "platform": "Olink"})]
    rows = build_evidence(hyps, {"d1": "decode__1_1__ENSG1", "k1": "ukbppp__OID1__ENSG1"},
                          {"decode__1_1__ENSG1": ("decode", "1_1", "ENSG1"), "ukbppp__OID1__ENSG1": ("ukbppp", "OID1", "ENSG1")},
                          {"decode__1_1__ENSG1": dec, "ukbppp__OID1__ENSG1": ukb}, {})
    assert [(r.evidence_state, r.s16_evidence_state) for r in rows] == [("supportive", "contradictory"), ("supportive", "")]
