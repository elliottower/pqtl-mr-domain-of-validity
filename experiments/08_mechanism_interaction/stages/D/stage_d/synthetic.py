"""SYNTHETIC stage A/B/C outputs for tests and the Modal smoke run. Every value is simulated;
no study data is read. Gene ids are `SYNTH_G###`, indications `SYNTH_I###`, programs `SYNTH_P###`.
"""
import hashlib
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

from stage_d.constants import PLAN_SHA256
from stage_d.evidence import evidence_score, evidence_state
from stage_d.fingerprint import RunFingerprint
from v8_manifest import InputRecord, write_manifest

SYNTHETIC_SCRIPT_SHA256 = hashlib.sha256(b"synthetic").hexdigest()
SYNTHETIC_TOKEN = "synthetic-token-0001"


def synthetic_fingerprint(seals: dict[str, str], **changes) -> RunFingerprint:
    """SYNTHETIC: a run fingerprint over the given seals; `changes` replace single components."""
    return RunFingerprint(**{"plan_sha256": PLAN_SHA256, "seals": seals, "script_sha256": SYNTHETIC_SCRIPT_SHA256,
                             "pins": {"synthetic": "test"}, "run_token": SYNTHETIC_TOKEN,
                             "environment": {"synthetic": "test"}, **changes})


def _hid(*parts) -> str:
    return hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


def stage_a_s12(frame: pd.DataFrame) -> list[bool]:
    """Stage A's S12 for the synthetic held-out S1 rows, written out here rather than taken from
    stage_d.sets.s12_select, so the join's agreement check compares two implementations: per gene,
    visit rows by (specificity_rank, indication_id, lowest program id); keep a row when none of its
    programs was kept before in that gene."""
    kept = [False] * len(frame)
    rows = sorted(range(len(frame)), key=lambda i: (frame["gene_ensembl"].iloc[i], int(frame["specificity_rank"].iloc[i]),
                                                     frame["indication_id"].iloc[i], min(frame["programs"].iloc[i])))
    used: dict[str, set[str]] = {}
    for i in rows:
        seen = used.setdefault(frame["gene_ensembl"].iloc[i], set())
        if seen.isdisjoint(frame["programs"].iloc[i]):
            kept[i] = True
            seen.update(frame["programs"].iloc[i])
    return kept


def make_tables(seed: int = 1, n_genes: int = 90, n_indications: int = 40, hyps_per_gene: int = 8,
                p_coloc: float = 0.8, or_blocking: float = 1.3, or_ratio: float = 3.0,
                base_advance: float = 0.3) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    gene_class = rng.choice(["aligned", "blocking", "other"], size=n_genes, p=[0.3, 0.5, 0.2])
    crossover = set(rng.choice(n_genes, size=6, replace=False).tolist())
    hyp, seen = [], set()
    for g in range(n_genes):
        pool = [f"SYNTH_P{g:03d}_{j}" for j in range(3)] + [f"SYNTH_PSHARED_{rng.integers(0, 10)}"]
        for _ in range(hyps_per_gene):
            ind = int(rng.integers(0, n_indications))
            cls = str(gene_class[g])
            if g in crossover and rng.random() < 0.5:
                cls = "aligned" if cls != "aligned" else "blocking"
            direction = "decrease" if cls != "other" else str(rng.choice(["decrease", "increase", "ambiguous"]))
            key = (f"SYNTH_ENSG{g:05d}", f"SYNTH_I{ind:03d}", direction, cls)
            if key in seen:
                continue
            seen.add(key)
            progs = sorted(set(rng.choice(pool, size=int(rng.integers(1, 3)), replace=False).tolist()))
            hyp.append({"hypothesis_id": _hid(*key), "gene_symbol": f"SYNTH_G{g:03d}", "gene_ensembl": key[0],
                        "indication_id": key[1], "indication_name": f"synthetic indication {ind}",
                        "direction": direction, "mechanism_class": cls, "class_reason": f"synthetic rule {cls}",
                        "drug_program_ids": ";".join(progs), "g": g, "ind": ind})
    h = pd.DataFrame(hyp)
    n = len(h)
    tier = rng.choice([1, 2, 3], size=n, p=[0.85, 0.1, 0.05])
    overlap = rng.choice(["no", "unknown", "yes"], size=n, p=[0.8, 0.15, 0.05])
    heldout = rng.random(n) < 0.9
    strict = rng.random(n) < 0.97
    restored = (rng.random(n) < 0.02) & ~strict
    in_s1 = (tier == 1) & (overlap == "no") & strict & ~restored
    ind_neuro = rng.random(n_indications) < 0.2
    ind_psych_only = ind_neuro & (rng.random(n_indications) < 0.3)
    h = h.assign(
        eligible=True, heldout=heldout, specificity_rank=h["ind"] % 7,
        instrument_source=rng.choice(["ukbppp", "decode", "interval"], size=n, p=[0.6, 0.3, 0.1]),
        instrument_assay_id=[f"SYNTH_A{g}" for g in h["g"]],
        platform=None, outcome_accession=[f"SYNTH_GCST{i:04d}" for i in h["ind"]], outcome_source="synthetic",
        outcome_tier=tier, overlap=overlap, outcome_n_case=rng.integers(500, 20000, n),
        outcome_n_control=rng.integers(5000, 400000, n),
        blood_secreted_hpa=(h["mechanism_class"] == "aligned") | (rng.random(n) < 0.3),
        secreted_uniprot=rng.random(n) < 0.4,
        neuro_psych=ind_neuro[h["ind"]], neuro_only=ind_neuro[h["ind"]] & ~ind_psych_only[h["ind"]],
        psych_only=ind_psych_only[h["ind"]], oncology=(h["ind"] % 5 == 0),
        pilot_indication=h["ind"] < 3, karim_launched=rng.random(n) < 0.02,
        single_protein_row=rng.random(n) < 0.9, conflicted_row_restored=restored, strict_instrument=strict,
        in_s1=in_s1, in_s4=in_s1 | ((tier == 1) & (overlap == "unknown") & strict & ~restored),
        in_s10=in_s1 | ((tier > 1) & (overlap == "no") & strict & ~restored), in_s12=False,
        subtype_restricted=tier == 2, phenotype_broader=tier == 3, sample_overlap=overlap == "yes",
        overlap_unknown=overlap == "unknown",
        pre_pqtl_publication_date=np.where(h["g"] % 2 == 0, "2021-06-01", "2023-10-01"),
    )
    h["platform"] = np.where(h["instrument_source"] == "ukbppp", "Olink", "SomaScan")
    h["outcome_neff"] = 4 / (1 / h["outcome_n_case"] + 1 / h["outcome_n_control"])
    h.loc[h.index[:3], "outcome_neff"] = np.nan
    h["programs"] = h["drug_program_ids"].str.split(";")
    base = (h["in_s1"] & h["heldout"]).to_numpy()
    h.loc[base, "in_s12"] = stage_a_s12(h.loc[base])
    h["variant"], h["in_s8"], h["in_s9"], h["in_s21"] = "primary", h["in_s1"], h["in_s1"], h["in_s1"]
    h["conflicted_row_restored"] = False   # restored rows exist only as variant s21 rows (below)

    coloc = rng.random(n) < p_coloc
    gd = rng.choice([-1, 1], size=n, p=[0.3, 0.7])
    pp = np.where(coloc, rng.beta(0.6, 0.6, n), np.nan)
    state = evidence_state(pp, gd, h["direction"])
    e = evidence_score(pp, gd, h["direction"])
    jitter = lambda: np.clip(pp + rng.normal(0, 0.05, n), 0, 1)  # noqa: E731
    ev = pd.DataFrame({
        "hypothesis_id": h["hypothesis_id"], "coloc_run": coloc,
        "not_run_reason": np.where(coloc, "", "fewer_than_50_shared"),
        "pp_h0": np.where(coloc, 0.0, np.nan), "pp_h1": np.where(coloc, 0.0, np.nan),
        "pp_h2": np.where(coloc, 0.0, np.nan), "pp_h3": np.where(coloc, 1 - pp, np.nan), "pp_h4": pp,
        "n_shared": np.where(coloc, rng.integers(50, 3000, n), -1), "frac_pqtl_retained": rng.uniform(0.3, 1, n),
        "frac_outcome_retained": rng.uniform(0.3, 1, n), "sentinel_or_proxy_retained": rng.random(n) < 0.9,
        "low_coverage": rng.random(n) < 0.1, "lead_variant": [f"rs{i}" for i in range(n)],
        "genetic_direction": gd, "evidence_state": state, "S": (state == "supportive").astype(int), "E": e,
        "protein_altering": rng.random(n) < 0.1,
        "platform_concordant": rng.choice(["concordant", "discordant", "untested"], size=n),
        "splicing_candidate": np.where(rng.random(n) < 0.1, None, rng.random(n) < 0.05),
        **{f"s15{x}_pp_h4": jitter() for x in "abcdefg"},
        "s15f_low_coverage_excluded": False,
        "s16_evidence_state": np.where(rng.random(n) < 0.1, state, None),
        "s17_sentinel_p": rng.uniform(0, 1, n),
    })
    # Table 12: the selected UKB-PPP or deCODE source repeats the primary state; the other source has
    # an instrument for about half the hypotheses, agreeing with the primary state 70% of the time.
    src = h["instrument_source"].to_numpy()
    other = np.where(rng.random(n) < 0.7, state, rng.choice(["supportive", "contradictory", "inconclusive"], size=n))
    has_other = rng.random(n) < 0.5
    for s in ("ukbppp", "decode"):
        ev[f"evidence_state_{s}"] = np.where(src == s, state, np.where((src != "interval") & has_other, other, None))
    ev["n_shared"] = ev["n_shared"].astype(object).where(coloc, None)

    S = ev["S"].to_numpy()
    aligned = (h["mechanism_class"] == "aligned").to_numpy()
    lo = (np.log(base_advance / (1 - base_advance)) + rng.normal(0, 0.5, n_genes)[h["g"]]
          + S * (np.log(or_blocking) + aligned * np.log(or_ratio)))
    adv = rng.random(n) < 1 / (1 + np.exp(-lo))
    status = np.where(adv, "advanced", rng.choice(["no_observed_advancement", "active", "business_only", "undated"],
                                                  size=n, p=[0.8, 0.1, 0.05, 0.05]))
    # no retrievable phase (PREREG §Missing data): half of the undated rows, chosen by gene parity
    status = np.where((status == "undated") & (h["g"].to_numpy() % 2 == 0), "no_phase", status)

    def adv_col(st):
        return np.where(st == "advanced", "1", np.where(st == "no_observed_advancement", "0", ""))

    def mature(months):
        st = status.copy()
        flip = (st == "active") & (rng.random(n) < months / 60)
        st[flip] = "no_observed_advancement"
        return st

    p2 = [date(2012, 1, 1) + timedelta(days=int(d)) for d in rng.integers(0, 4000, n)]
    t3 = rng.exponential(900, n)
    dated = rng.random(n) < 0.85
    stop = np.where(status == "no_observed_advancement",
                    rng.choice(["efficacy", "safety", "unclassified", "completed_no_successor"], size=n), "")
    stop = np.where(status == "business_only", "business", stop)
    out = pd.DataFrame({
        "hypothesis_id": h["hypothesis_id"], "status_24": status, "advanced_24": adv_col(status),
        **{f"status_{m}": mature(m) for m in (12, 36, 48)},
        "stop_code": stop, "why_stopped_texts": np.where(stop == "", "", "synthetic reason"),
        "efficacy_coded": stop == "efficacy",
        "phase2_start_date": [d.isoformat() if k else "" for d, k in zip(p2, dated)],
        "phase3_start_date": [(d + timedelta(days=int(t))).isoformat() if a and k else ""
                              for d, t, a, k in zip(p2, t3, adv, dated)],
        "time_to_phase3_days": np.where(dated, np.where(adv, t3, rng.uniform(100, 3000, n)).round(), np.nan),
        "phase3_event": np.where(dated, adv.astype(int), -1),
        "phase3_success": np.where(adv, (rng.random(n) < 0.5).astype(int).astype(str), ""),
        "approved": adv & (rng.random(n) < 0.3),
        "chembl_max_phase_for_ind": np.where(adv, 3.0, 2.0),
        "earliest_phase2_start": [d.isoformat() if k else "" for d, k in zip(p2, dated)],
        "last_phase2_end_date": [(d + timedelta(days=700)).isoformat() for d in p2],
    })
    for m in (12, 36, 48):
        out[f"advanced_{m}"] = adv_col(out[f"status_{m}"].to_numpy())
    out["phase3_event"] = out["phase3_event"].astype(object).where(dated, "")

    # S19 arm, drawn after every other primary-row draw: most aligned rows on blood-secreted targets are neutralizing biologics (the rest
    # oligonucleotides or degraders, in no arm); most blocking rows on blood-secreted targets are
    # small-molecule blockers.
    arm_draw = rng.random(n)
    h["s19_arm"] = np.where(h["blood_secreted_hpa"] & (h["mechanism_class"] == "aligned") & (arm_draw < 0.8),
                            "neutralizing_biologic",
                            np.where(h["blood_secreted_hpa"] & (h["mechanism_class"] == "blocking") & (arm_draw < 0.9),
                                     "small_molecule_blocker", ""))
    # S8 / S9 re-formed: primary rows are shared except where the set's table changes the class; a
    # changed hypothesis is written again as a variant-only row with the new class. Its evidence and
    # outcome are its source row's (neither depends on the class). Appended last, so the random
    # draws of every primary row are those of the generator without variants.
    moved = {"s8": (h["in_s1"] & (h["mechanism_class"] == "other") & (h["direction"] == "decrease")
                    & (h["g"] % 3 == 0), "blocking"),
             "s9": (h["in_s1"] & (h["mechanism_class"] == "aligned") & ~h["secreted_uniprot"], "other"),
             "s21": (h["in_s1"] & (h["mechanism_class"] == "other") & (h["direction"] == "ambiguous"), "blocking")}
    hs, evs, outs = [h], [ev], [out]
    for v, (mask, cls) in moved.items():
        h.loc[mask, f"in_{v}"] = False
        x = h.loc[mask].copy()
        x["mechanism_class"], x["variant"], x["s19_arm"] = cls, v, ""
        if v == "s21":
            x["direction"], x["conflicted_row_restored"] = "decrease", True
        x["hypothesis_id"] = [_hid(e, i, d, cls, v) for e, i, d in zip(x["gene_ensembl"], x["indication_id"], x["direction"])]
        x[["in_s1", "in_s4", "in_s10", "in_s12", "in_s8", "in_s9", "in_s21"]] = False
        x[f"in_{v}"] = True
        hs.append(x)
        e = ev.loc[mask.to_numpy()].assign(hypothesis_id=x["hypothesis_id"].to_numpy())
        # the state depends on the intervention direction, which S21 can change
        st = evidence_state(e["pp_h4"], e["genetic_direction"], x["direction"])
        e = e.assign(evidence_state=st, S=(st == "supportive").astype(int),
                     E=evidence_score(e["pp_h4"], e["genetic_direction"], x["direction"]))
        for src in ("ukbppp", "decode"):
            own = (x["instrument_source"] == src).to_numpy()
            e[f"evidence_state_{src}"] = np.where(own, st, e[f"evidence_state_{src}"].to_numpy())
        evs.append(e)
        outs.append(out.loc[mask.to_numpy()].assign(hypothesis_id=x["hypothesis_id"].to_numpy()))
    h, ev, out = (pd.concat(f, ignore_index=True) for f in (hs, evs, outs))
    return h.drop(columns=["g", "ind", "programs"]), ev, out


def write_stage_outputs(root: Path, hyp: pd.DataFrame, ev: pd.DataFrame, out: pd.DataFrame) -> dict[str, str]:
    """Write synthetic A/B/C outputs with sealed MANIFEST.tsv files; return the manifest sha256s."""
    a = root / "A" / "output"
    a.mkdir(parents=True, exist_ok=True)
    hyp.to_csv(a / "hypotheses.csv", index=False)
    pd.DataFrame({"step": ["universe", "eligible"], "remaining": [len(hyp) * 2, len(hyp)],
                  "excluded": [0, len(hyp)], "reason": ["", "synthetic"]}).to_csv(a / "funnel.csv", index=False)
    pd.DataFrame({"indication_id": ["SYNTH_I000"], "accession": ["SYNTH_GCST0000"], "tier": [1],
                  "selected": [True]}).to_csv(a / "outcome_gwas_selection.csv", index=False)
    pd.DataFrame({"molecule_type": ["Antibody"], "action_type": ["INHIBITOR"], "localization": ["blood"],
                  "class": ["aligned"], "n": [1]}).to_csv(a / "mechanism_crosstab.csv", index=False)
    write_manifest(a, [a / n for n in ("hypotheses.csv", "funnel.csv", "outcome_gwas_selection.csv",
                                         "mechanism_crosstab.csv")], SYNTHETIC_SCRIPT_SHA256, [])
    from_a = [InputRecord(name="hypotheses", path="experiments/08_mechanism_interaction/stages/A/output/hypotheses.csv",
                          sha256=hashlib.sha256((a / "hypotheses.csv").read_bytes()).hexdigest())]
    for stage, name, frame in (("B", "evidence.csv", ev), ("C", "outcomes.csv", out)):
        d = root / stage / "output"
        d.mkdir(parents=True, exist_ok=True)
        frame.to_csv(d / name, index=False)
        write_manifest(d, [d / name], SYNTHETIC_SCRIPT_SHA256, from_a)
    return {s: hashlib.sha256((root / s / "output" / "MANIFEST.tsv").read_bytes()).hexdigest() for s in "ABC"}
