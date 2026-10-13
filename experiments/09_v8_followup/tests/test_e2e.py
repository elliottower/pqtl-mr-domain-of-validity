"""End to end on V8's synthetic stage tables (stage_d.synthetic), in the analysis image: prepare,
prepare_q4, Q1 point and bootstrap, one Bayesian fit at smoke sampler settings, and the fields
assemble reads. No study data."""
import json

import numpy as np
import pandas as pd
import pytest
from openpyxl import Workbook

from followup.common import write_json
from followup.q1 import MODELS, run_bootstrap, summarize_bootstrap
from followup.run import FITS, load_frame, prepare, prepare_q4, q1_point
from stage_d.bayes import fit
from stage_d.constants import SamplerSettings
from stage_d.designs import load_design
from stage_d.evidence import evidence_score, evidence_state
from stage_d.fitting import run_fit
from stage_d.synthetic import make_tables, write_stage_outputs

SHARED = ["coloc_run", "not_run_reason", "pp_h0", "pp_h1", "pp_h2", "pp_h3", "pp_h4", "n_shared", "lead_variant",
          "genetic_direction"]


def realistic_tables():
    h, ev, out = make_tables(seed=11, n_genes=70, hyps_per_gene=8)
    g = h["gene_ensembl"].str.extract(r"(\d+)$")[0].astype(int)
    h["gene_ensembl"] = [f"ENSG{x:011d}" for x in g]
    h["instrument_assay_id"] = np.where(h["instrument_source"] == "ukbppp", [f"OID{x:05d}" for x in g],
                                        np.where(h["instrument_source"] == "decode", [f"{x}_1" for x in g], "IV"))
    # one colocalization per unique genetic unit, as in stage B
    j = h[["hypothesis_id", "instrument_source", "instrument_assay_id", "gene_ensembl", "outcome_accession", "direction"]].merge(
        ev, on="hypothesis_id")
    unit = j.groupby(["instrument_source", "instrument_assay_id", "gene_ensembl", "outcome_accession"], sort=False).ngroup()
    lead = j.groupby(unit).head(1).set_index(unit[j.groupby(unit).head(1).index])
    for c in SHARED:
        j[c] = lead[c].reindex(unit).to_numpy()
    j["evidence_state"] = evidence_state(j["pp_h4"], j["genetic_direction"], j["direction"])
    j["S"] = (j["evidence_state"] == "supportive").astype(int)
    j["E"] = evidence_score(j["pp_h4"], j["genetic_direction"], j["direction"])
    for src in ("ukbppp", "decode"):
        own = j["instrument_source"] == src
        j[f"evidence_state_{src}"] = np.where(own, j["evidence_state"], j[f"evidence_state_{src}"])
    return h, j[ev.columns], out


def spreadsheets(tmp, genes):
    sun, ferk = tmp / "sun.xlsx", tmp / "ferk.xlsx"
    wb = Workbook()
    ws = wb.active
    ws.title = "ST9"
    for _ in range(4):
        ws.append([])
    ws.append(["Variant ID", "UKBPPP ProteinID", "log10(p) (discovery)", "cis/trans"])
    for x in genes:
        ws.append(["v", f"G{x}:P1:OID{x:05d}:v1", 20.0 + x, "cis"])
        ws.append(["v", f"G{x}:P1:OID{x:05d}:v1", 999.0, "trans"])
    wb.save(sun)
    wb = Workbook()
    ws = wb.active
    ws.title = "ST08"
    ws.append(["SUPPLEMENTARY TABLE 8"])
    ws.append([])
    ws.append(["pQTL_ID", "SeqId", "-Log10(P)\n(unadj.)"])
    for x in genes:
        ws.append(["p", f"{x}_1", 15.0 + x])
    wb.save(ferk)
    return sun, ferk


@pytest.fixture(scope="module")
def run_root(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("e2e")
    h, ev, out = realistic_tables()
    seals = write_stage_outputs(tmp / "stages", h, ev, out)
    genes = sorted(set(h["gene_ensembl"].str[4:].astype(int)))
    sun, ferk = spreadsheets(tmp, genes)
    root = tmp / "vol"
    prepare(tmp / "stages", seals, sun, ferk, root)
    rng = np.random.default_rng(2)
    gs = sorted(h["gene_ensembl"].unique())
    inds = sorted(h["indication_id"].unique())
    evid = pd.DataFrame({"studyLocusId": [f"L{i}" for i in range(300)], "targetId": rng.choice(gs, 300),
                         "diseaseId": rng.choice(inds + ["CHILD"], 300), "resourceScore": rng.uniform(0.05, 1, 300),
                         "publicationDate": rng.choice(["2008-01-01", "2016-05-05", None], 300)})
    evid["score"] = evid["resourceScore"]
    disease = pd.DataFrame({"id": inds + ["CHILD"],
                            "descendants": [np.array(["CHILD"]) if i == inds[0] else np.array([], dtype=object)
                                            for i in inds] + [None]})
    prepare_q4(evid, disease, root)
    return root


def test_prepare_writes_every_result_file_with_the_declaration(run_root):
    res = run_root / "results"
    for name in ("q1_descriptive.json", "q3_balanced_accuracy.json", "q3_population.json", "q5_counts.json",
                 "q7_event_counts.json", "q4_tables.json"):
        assert "post hoc and exploratory" in json.loads((res / name).read_text())["declaration"]
    assert (res / "q2_pp_h4_ge_0.5.csv").exists() and (res / "q1_class_share_by_neff_quintile.png").exists()
    df = load_frame(run_root)
    assert df["assay_cis_signal_strength"].notna().sum() > 0
    assert set(df["coloc_class"]) <= {"not_run", "shared", "distinct", "outcome_only", "no_outcome_signal",
                                      "unresolved_h0h1", "unresolved_h2", "unresolved_h3", "unresolved_h4"}


def test_every_fit_has_a_design_on_disk(run_root):
    for spec in FITS:
        assert (run_root / "work" / "designs" / f"{spec.design_id}.npz").exists()


def test_q1_models_fit_bootstrap_and_summarize(run_root, tmp_path):
    for model in MODELS:
        pop, spec, point = q1_point(run_root, model)
        if spec is None:
            continue
        path = tmp_path / f"{model}.jsonl"
        run_bootstrap(pop, spec, model, path, reps=20)
        out = summarize_bootstrap(point, path, reps=20)
        write_json(tmp_path / f"{model}.json", out)
        assert out["replicates"] == 20


def test_one_fit_runs_and_carries_what_assemble_reads(run_root):
    smoke = SamplerSettings(chains=2, draws=150, tune=150, target_accept=0.9)
    spec = next(s for s in FITS if s.fit_id == "q4_h1_G").model_copy(update={"primary": smoke, "rerun": smoke})
    path = run_root / "work" / "designs" / f"{spec.design_id}.npz"
    final = run_fit(spec, load_design(path), run_root / "work" / "fits", fit, fingerprint="e2e", design_sha256="x")
    att = final["attempts"][-1]
    for key in ("focal", "focal_posterior_sd", "focal_prior_sd", "prior_dominated", "coefficients",
                "stratum_odds_ratios", "marginal"):
        assert key in att["summary"]
    assert "G" in att["summary"]["coefficients"] and "GxA" in att["summary"]["coefficients"]
