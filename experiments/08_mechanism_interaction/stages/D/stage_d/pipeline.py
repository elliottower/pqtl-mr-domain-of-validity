"""Stage D orchestration without PyMC: prepare (guard, join, sets, gates, designs, fit list),
frequentist checks, and assemble (decisions, results.json, tables, figures, MANIFEST.tsv).
The Bayesian fits between prepare and assemble run in fitting.py on Linux.

Work files (designs, checkpoints, fit attempts) go to `work_dir`; `assemble` writes
`results.json` first, then `tables/`, `figures/`, and `MANIFEST.tsv` last.

The work directory is bound to one run fingerprint (stage_d/fingerprint.py). `prepare` claims it
and refuses a non-empty directory under another fingerprint; `run_frequentist`, the fits and
`assemble` refuse unless it carries this run's. The plan records each design file's sha256; a fit
or frequentist result is accepted only with this run's fingerprint and that sha256, and `assemble`
re-hashes every design file and checks every fit attempt, final fit and frequentist result.
"""
import json
import math
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import BaseModel
from v8_manifest import (INPUTS_NAME, MANIFEST_NAME, InputRecord, code_record, code_sha256, guard_files, relative_files,
                         sha256_file)
from v8_manifest import write_manifest as write_v8_manifest

from stage_d.constants import (FOREST_MAIN_SETS, FREEZE_COMMIT, N_BOOTSTRAP, N_NULL_REPS, OSF_REGISTRATION,
                               PLAN_SHA256, PROB_THRESHOLD, SEED, SamplerSettings)
from stage_d.cox import CoxFailure, cox_two_way
from stage_d.decisions import ModelEvidence, against_h1, decide_h1, decide_h4, h1_wording
from stage_d.descriptive import (stratum_2x2_table, table1_funnel, table3_by_class, table4_other_reasons,
                                 table5_evidence, table6_missingness, table7_outcomes, table11_per_indication,
                                 table12_platform, table13_neuro, table15_karim)
from stage_d.designs import (COVARIATES, Design, DesignError, build_design, component_labels, load_design, save_design,
                             within_gene_estimable)
from stage_d.figures import forest_plot
from stage_d.fingerprint import RunFingerprint, StaleWork, check_record, claim, require
from stage_d.fitspec import FitSpec
from stage_d.frequentist import frequentist_check
from stage_d.gates import StratumGate, h1_gate, h3_gate, h4_gate
from stage_d.guard import verify_inputs
from stage_d.join import derive, join_stages, load_stage_tables
from stage_d.posterior import PRIMARY_PRIOR, PRIORS
from stage_d.sets import SET_ORDER, form_sets

ALL_PRIORS = tuple(PRIORS)
NAMED_MODELS = (  # design id prefix, kind, y column, role
    ("h1", "h1", "y", "H1"), ("h4", "h4", "y", "H4"), ("h2", "h2", "y", "H2"), ("h3", "h3", "y", "H3"),
    ("h1_E", "h1_E", "y", "continuous_H1"), ("h2_E", "h2_E", "y", "continuous_H2"),
    ("within", "within", "y", "within_gene"),
    ("h1_p3s", "h1", "phase3_success", "phase3_success_H1"), ("h2_p3s", "h2", "phase3_success", "phase3_success_H2"),
)
SCRIPT_FILES = ("bayes.py", "constants.py", "cox.py", "decisions.py", "descriptive.py", "designs.py", "environment.py",
                "evidence.py", "figures.py", "fingerprint.py", "fitspec.py", "fitting.py", "frequentist.py", "gates.py",
                "guard.py", "join.py", "pipeline.py", "posterior.py", "schemas.py", "sets.py")


class IncompleteRun(RuntimeError):
    pass


class Plan(BaseModel):
    created_utc: str
    fingerprint: str
    seals: dict
    derived: dict
    sets: list[dict]
    gates: dict
    set_gates: dict
    designs: dict
    skipped: dict[str, str]
    fits: list[FitSpec]


def load_joined(stages_root: Path, expected: dict[str, str]):
    seals = verify_inputs(stages_root, expected)
    hyp, ev, out = load_stage_tables(stages_root)
    df, info = derive(join_stages(hyp, ev, out))
    return seals, df, info, form_sets(df)


def _y_varies(frame: pd.DataFrame, y_col: str) -> bool:
    y = frame[y_col].dropna()
    return len(y) > 1 and y.nunique() == 2


def design_path(work_dir: Path, design_id: str) -> Path:
    return work_dir / "designs" / f"{design_id}.npz"


def load_planned_design(work_dir: Path, plan: Plan, design_id: str) -> tuple[Design, str]:
    """The design and its file's sha256, after checking the file is the one the plan recorded."""
    path = design_path(work_dir, design_id)
    digest = sha256_file(path)
    if digest != plan.designs[design_id]["sha256"]:
        raise StaleWork(f"{path} has sha256 {digest}, plan.json recorded {plan.designs[design_id]['sha256']}")
    return load_design(path), digest


def prepare(stages_root: Path, expected: dict[str, str], work_dir: Path, fp: RunFingerprint,
            samplers: tuple[SamplerSettings, SamplerSettings] | None = None) -> Plan:
    """Verify the sealed inputs, claim `work_dir` for `fp`, then write the designs and plan.json.
    A work directory already holding this run's plan.json is not rewritten (the design files are
    re-hashed against it), so results computed from those designs stay valid."""
    if fp.seals != expected:
        raise ValueError(f"the fingerprint names seals {fp.seals}, the run guard returned {expected}")
    seals, df, info, sets = load_joined(stages_root, expected)
    claim(work_dir, fp)
    if (work_dir / "plan.json").is_file():
        plan = load_plan(work_dir, fp)
        for design_id in plan.designs:
            load_planned_design(work_dir, plan, design_id)
        return plan
    s1 = sets["S1"].frame
    gates = {"h1": h1_gate(s1), "h4": h4_gate(s1), "h3": h3_gate(s1)}
    within_ok, within_counts = within_gene_estimable(s1.loc[s1["y"].notna()])
    designs, skipped, fits = {}, {}, []
    sampler_kw = {"primary": samplers[0], "rerun": samplers[1]} if samplers else {}

    def add(design_id, kind, frame, y_col, set_id, role, priors):
        try:
            d = build_design(design_id, kind, frame, y_col)
        except DesignError as err:
            skipped[design_id] = str(err)
            return
        if len(np.unique(d.y)) < 2:
            skipped[design_id] = "outcome does not vary among analysed rows"
            return
        save_design(d, design_path(work_dir, design_id))
        designs[design_id] = {"sha256": sha256_file(design_path(work_dir, design_id)),
                              "kind": kind, "set_id": set_id, "role": role, "y": y_col, "n": d.n,
                              "columns": list(d.columns), "dropped_covariates": list(d.dropped_covariates),
                              "genes": d.n_gene, "indications": d.n_indication, "programs": d.W.shape[1],
                              "components": int(len(np.unique(d.component)))}
        for prior in priors:
            fits.append(FitSpec(fit_id=f"{design_id}__{prior}", design_id=design_id, kind=kind, set_id=set_id,
                                prior=prior, role=role, **sampler_kw))

    s_varies = s1["S"].nunique() == 2
    for prefix, kind, y_col, role in NAMED_MODELS:
        design_id = f"{prefix}__S1"
        base = "h4" if kind == "h4" else "h1" if kind in ("h1", "h1_E", "within") else "h2"
        if base == "h1" and not gates["h1"].structurally_estimable:
            skipped[design_id] = "H1 structurally non-estimable: " + "; ".join(gates["h1"].reasons)
        elif base == "h4" and not gates["h4"].structurally_estimable:
            skipped[design_id] = "H4 structurally non-estimable: " + "; ".join(gates["h4"].reasons)
        elif base == "h2" and kind != "h3" and not s_varies:
            skipped[design_id] = "S does not vary in S1"
        elif kind == "h3" and not gates["h3"].passes:
            skipped[design_id] = (f"H3 gate: {gates['h3'].n_supportive} supportive, "
                                  f"{gates['h3'].n_contradictory} contradictory (< 20 each)")
        elif kind == "within" and not within_ok:
            skipped[design_id] = f"within-gene S:A not estimable: {within_counts}"
        elif y_col != "y" and not _y_varies(s1, y_col):
            skipped[design_id] = f"{y_col} does not vary"
        else:
            add(design_id, kind, s1, y_col, "S1", role, ALL_PRIORS)

    set_gates = {}
    for set_id in SET_ORDER[1:]:
        aset = sets[set_id]
        if not aset.info.formed:
            continue
        for model in aset.info.models:
            gate = h1_gate(aset.frame) if model == "h1" else h4_gate(aset.frame)
            set_gates[f"{model}__{set_id}"] = gate.model_dump()
            if gate.structurally_estimable:
                add(f"{model}__{set_id}", model, aset.frame, "y", set_id, "set_row", (PRIMARY_PRIOR,))
            else:
                skipped[f"{model}__{set_id}"] = "structurally non-estimable: " + "; ".join(gate.reasons)

    plan = Plan(created_utc=datetime.now(timezone.utc).isoformat(), fingerprint=fp.digest,
                seals={k: v.model_dump() for k, v in seals.items()},
                derived=info.model_dump(), sets=[sets[s].info.model_dump() for s in SET_ORDER],
                gates={**{k: v.model_dump() for k, v in gates.items()}, "within": {"estimable": within_ok, **within_counts}},
                set_gates=set_gates, designs=designs, skipped=skipped, fits=fits)
    (work_dir / "plan.json").write_text(plan.model_dump_json(indent=2))
    return plan


def load_plan(work_dir: Path, fp: RunFingerprint) -> Plan:
    """plan.json of a work directory bound to `fp`; any other directory or plan raises StaleWork."""
    require(work_dir, fp)
    plan = Plan.model_validate_json((work_dir / "plan.json").read_text())
    if plan.fingerprint != fp.digest:
        raise StaleWork(f"{work_dir / 'plan.json'} carries fingerprint {plan.fingerprint}, this run is {fp.digest}")
    return plan


def run_frequentist(work_dir: Path, fp: RunFingerprint, design_ids=None, reps: int = N_NULL_REPS,
                    n_boot: int = N_BOOTSTRAP, on_commit=None) -> dict:
    plan = load_plan(work_dir, fp)
    out = {}
    for design_id in design_ids or plan.designs:
        d, design_sha = load_planned_design(work_dir, plan, design_id)
        path = work_dir / "freq" / f"{design_id}.json"
        if path.exists():
            out[design_id] = json.loads(path.read_text())
            check_record(out[design_id], str(path), fp.digest, design_sha)
            continue
        res = frequentist_check(d, work_dir / "freq_checkpoints", reps, n_boot, on_commit)
        res.update({"null_reps_per_cell": reps, "bootstrap_resamples": n_boot, "fingerprint": fp.digest,
                    "design_sha256": design_sha})
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=2, default=float))
        if on_commit:
            on_commit()
        out[design_id] = res
    return out


# ----------------------------------------------------------------------------- assemble


def reported_attempt(final: dict | None) -> dict | None:
    """The passing attempt, or the last attempt run when none passed."""
    if not final or not final["attempts"]:
        return None
    return next((a for a in final["attempts"] if a["passed"]), final["attempts"][-1])


def model_evidence(gate: StratumGate, final: dict | None) -> ModelEvidence:
    if final is None:
        return ModelEvidence(structurally_estimable=gate.structurally_estimable, reliable=gate.reliable,
                             fit_outcome="not_fitted")
    att = reported_attempt(final)
    s = att["summary"]
    return ModelEvidence(structurally_estimable=gate.structurally_estimable, reliable=gate.reliable,
                         fit_outcome=final["outcome"], prior_dominated=s["prior_dominated"],
                         pr_predicted=s["focal"]["pr_predicted"], pr_opposite=s["focal"]["pr_opposite"])


def cox_models(s1_pre: pd.DataFrame) -> dict:
    dated = s1_pre.loc[s1_pre["time_to_phase3_days"].notna()]
    out = {"excluded_no_dated_phase2": int(s1_pre["time_to_phase3_days"].isna().sum()), "n_dated": len(dated)}
    for name, rows, terms in (("H1", dated.loc[dated["cls"].isin(["aligned", "blocking"])], ("S", "A", "SxA")),
                              ("H2", dated, ("S",))):
        if rows.empty or rows["phase3_event"].sum() == 0:
            out[name] = {"fitted": False, "reason": "no rows or no events"}
            continue
        A = (rows["cls"] == "aligned").to_numpy(dtype=float)
        S = rows["S"].to_numpy(dtype=float)
        cols = {"S": S, "A": A, "SxA": S * A}
        X, names = [cols[t] for t in terms], list(terms)
        for col, label in COVARIATES:
            v = rows[col].to_numpy(dtype=float)
            if np.ptp(v) > 0:
                X.append(v)
                names.append(label)
        comp = component_labels(list(rows["gene_ensembl"]), list(rows["programs"]))
        _, ind = np.unique(rows["indication_id"].to_numpy(), return_inverse=True)
        try:
            coefs = cox_two_way(rows["time_to_phase3_days"].to_numpy(dtype=float),
                                rows["phase3_event"].to_numpy(dtype=float), np.column_stack(X), comp, ind, names)
            out[name] = {"fitted": True, "n": len(rows), "events": int(rows["phase3_event"].sum()), "coefficients": coefs}
        except (CoxFailure, np.linalg.LinAlgError) as err:
            out[name] = {"fitted": False, "reason": str(err)}
    return out


def _focal_row(final: dict | None) -> dict:
    att = reported_attempt(final)
    if att is None:
        return {"fit_outcome": "not_fitted"}
    s = att["summary"]
    return {"fit_outcome": final["outcome"], "attempt": att["attempt"], "prior_dominated": s["prior_dominated"],
            "focal": s["focal"], "stratum_odds_ratios": s["stratum_odds_ratios"]}


def set_table(plan: Plan, finals: dict, freq: dict) -> list[dict]:
    rows = []
    for info in plan.sets:
        for model in info["models"]:
            design_id = f"{model}__{info['set_id']}"
            fit_id = f"{design_id}__{PRIMARY_PRIOR}"
            f = freq.get(design_id, {})
            status = ("not formed: " + info["reason"] if not info["formed"] else
                      "skipped: " + plan.skipped[design_id] if design_id in plan.skipped else "fitted")
            rows.append({"set_id": info["set_id"], "model": model, "definition": info["definition"], "status": status,
                         "n_pre_outcome": info["n_pre_outcome"], "n_analysed": info["n_analysed"],
                         "gate": plan.set_gates.get(design_id) or (plan.gates.get(model) if info["set_id"] == "S1" else None),
                         **_focal_row(finals.get(fit_id)),
                         "frequentist": {k: f.get(k) for k in ("estimate", "se", "z", "p_one_sided_unadjusted",
                                                               "p_one_sided_size_corrected", "reject_size_corrected")}
                         | {"crit_z": f.get("size_correction", {}).get("crit_z"), "bootstrap": f.get("bootstrap")}})
    return rows


def table9(finals: dict, freq: dict) -> dict[str, list[dict]]:
    coefs, marginal, ppc, prior_pc, sens = [], [], [], [], []
    for fit_id, final in finals.items():
        spec = final["spec"]
        for att in final["attempts"]:
            s = att["summary"]
            key = {"fit_id": fit_id, "role": spec["role"], "set_id": spec["set_id"], "prior": spec["prior"],
                   "attempt": att["attempt"], "passed": att["passed"], **att["diagnostics"]}
            for name, v in s["coefficients"].items():
                coefs.append({**key, "term": name, **v})
            if s["marginal"]:
                for stratum, pr in s["marginal"]["probabilities"].items():
                    for ev, v in pr.items():
                        marginal.append({**key, "quantity": f"P(advance | {ev}, {stratum})", **v})
                for stratum, v in s["marginal"]["risk_difference"].items():
                    marginal.append({**key, "quantity": f"RD({stratum})", **v})
                marginal.append({**key, "quantity": s["marginal"]["definition"], **s["marginal"]["risk_difference_interaction"]})
            for g, v in s["posterior_predictive"].items():
                ppc.append({**key, "group": g, **v})
            for g, v in s["prior_predictive"].items():
                prior_pc.append({**key, "group": g, **v})
            sens.append({**key, "focal": s["focal"]["quantity"], **{k: v for k, v in s["focal"].items() if k != "quantity"},
                         "prior_dominated": s["prior_dominated"]})
    freq_rows = []
    for design_id, f in freq.items():
        if not f.get("fitted"):
            freq_rows.append({"design_id": design_id, "fitted": False, "reason": f.get("reason")})
            continue
        for term, b in f["coefficients"].items():
            freq_rows.append({"design_id": design_id, "term": term, "estimate": b, "se_two_way": f["se_all"][term]})
    return {"table09_coefficients_bayes": coefs, "table09_coefficients_frequentist": freq_rows,
            "table09_marginal": marginal, "table09_posterior_predictive": ppc,
            "table09_prior_predictive": prior_pc, "table09_prior_sensitivity": sens}


def _clean(obj):
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return float(obj) if math.isfinite(obj) else None
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, (pd.Timestamp, datetime)):
        return obj.isoformat()
    return obj


def script_files(package_dir: Path, runner: Path, extra: Mapping[str, Path] | None = None) -> dict[str, Path]:
    """The runner (`run_stage_d.py`) and the stage_d modules of SCRIPT_FILES by their path relative
    to the stage D directory, the shared guard modules (`run_guard/...`), and `extra` (relative
    path -> file, for code outside the stage directory)."""
    return {**relative_files(package_dir.parent, [runner, *(package_dir / f for f in SCRIPT_FILES)]), **guard_files(),
            **(extra or {})}


def script_sha256(package_dir: Path, runner: Path, extra: Mapping[str, Path] | None = None) -> str:
    """One hash over `script_files`, each by relative path, length and bytes (v8_manifest.code_sha256)."""
    return code_sha256(script_files(package_dir, runner, extra))


def write_manifest(out_dir: Path, script_digest: str, seals: dict, repo_commit: str) -> Path:
    """Shared format (stages/run_guard/v8_manifest.py): every file under `out_dir`, and INPUTS.tsv
    naming the sealed A, B and C MANIFEST.tsv files D read and, as `stage_code`, the script digest
    and the repository commit of the run."""
    outputs = sorted(x for x in out_dir.rglob("*") if x.is_file() and x.name not in (MANIFEST_NAME, INPUTS_NAME))
    inputs = [InputRecord(name=f"stage_{s}_manifest", path=f"experiments/08_mechanism_interaction/stages/{s}/output/"
                          f"{MANIFEST_NAME}", sha256=v["manifest_sha256"]) for s, v in sorted(seals.items())]
    inputs.append(code_record("D", script_digest, repo_commit))
    return write_v8_manifest(out_dir, outputs, script_digest, inputs)


def assemble(stages_root: Path, expected: dict[str, str], work_dir: Path, out_dir: Path, *, fp: RunFingerprint,
             runner_path: Path, package_versions: dict[str, str], repo_commit: str, n_boot: int = N_BOOTSTRAP) -> dict:
    """`repo_commit` is v8_run_guard.run_commit of the caller; it is written to results.json
    (`meta.repo_commit`) and to the `stage_code` row of INPUTS.tsv."""
    plan = load_plan(work_dir, fp)
    seals, df, info, sets = load_joined(stages_root, expected)
    if {k: v.manifest_sha256 for k, v in seals.items()} != {k: v["manifest_sha256"] for k, v in plan.seals.items()}:
        raise IncompleteRun("stage manifests changed between prepare and assemble")
    design_sha = {design_id: load_planned_design(work_dir, plan, design_id)[1] for design_id in plan.designs}
    finals = {}
    for spec in plan.fits:
        p = work_dir / "fits" / spec.fit_id / "final.json"
        if not p.exists():
            raise IncompleteRun(f"fit {spec.fit_id} has no final.json")
        finals[spec.fit_id] = json.loads(p.read_text())
        for rec, what in [(finals[spec.fit_id], str(p))] + [(a, f"{p} attempt {a['attempt']}")
                                                           for a in finals[spec.fit_id]["attempts"]]:
            check_record(rec, what, fp.digest, design_sha[spec.design_id])
    freq = {}
    for design_id in plan.designs:
        p = work_dir / "freq" / f"{design_id}.json"
        if not p.exists():
            raise IncompleteRun(f"frequentist check {design_id} missing")
        freq[design_id] = json.loads(p.read_text())
        check_record(freq[design_id], str(p), fp.digest, design_sha[design_id])

    gates = {k: StratumGate.model_validate(plan.gates[k]) for k in ("h1", "h4")}
    primary = lambda prefix: finals.get(f"{prefix}__S1__{PRIMARY_PRIOR}")  # noqa: E731
    ev_h1 = model_evidence(gates["h1"], primary("h1"))
    ev_h4 = model_evidence(gates["h4"], primary("h4"))
    d_h1 = decide_h1(ev_h1)
    d_h4 = decide_h4(d_h1, ev_h4)
    h2_att = reported_attempt(primary("h2"))
    h2_focal = h2_att["summary"]["focal"] if h2_att else None
    within_att = reported_attempt(primary("within"))
    h1_att = reported_attempt(primary("h1"))
    within = {"estimable": plan.gates["within"]["estimable"], "counts": plan.gates["within"],
              "focal": within_att["summary"]["focal"] if within_att else None}
    if within_att and h1_att:
        within["sign_opposite_to_primary"] = bool(np.sign(within_att["summary"]["focal"]["median"])
                                                  != np.sign(h1_att["summary"]["focal"]["median"]))

    s1 = sets["S1"].frame
    s1_an = s1.loc[s1["y"].notna()]
    a_out = stages_root / "A" / "output"
    pilot = sets["S2"].frame.loc[~sets["S2"].frame["heldout"]]
    tables = {
        "table01_funnel": table1_funnel(pd.read_csv(a_out / "funnel.csv"), s1),
        "table02_outcome_gwas_selection": pd.read_csv(a_out / "outcome_gwas_selection.csv", dtype=str,
                                                      keep_default_na=False).to_dict(orient="records"),
        "table03_by_class": table3_by_class(s1),
        "table05_evidence": table5_evidence(s1),
        "table06_missingness": table6_missingness(s1),
        "table07_outcomes": table7_outcomes(s1),
        "table08_stratum_2x2": stratum_2x2_table(s1_an, n_boot=n_boot),
        "table11_per_indication": table11_per_indication(s1),
        "table12_platform": table12_platform(s1),
        "table13_neuro_psych": table13_neuro(s1_an, n_boot=n_boot),
        "table14_pilot_2x2": stratum_2x2_table(pilot, n_boot=n_boot, seed_tag=14),
        "table15_karim_targets": table15_karim(df),
    }
    tables["table04_mechanism"] = table4_other_reasons(pd.read_csv(a_out / "mechanism_crosstab.csv", dtype=str,
                                                                   keep_default_na=False), s1)
    set_rows = set_table(plan, finals, freq)
    tables["table10_analysis_sets"] = set_rows
    tables.update(table9(finals, freq))

    results = {
        "meta": {"generated_utc": datetime.now(timezone.utc).isoformat(), "plan_sha256": PLAN_SHA256,
                 "freeze_commit": FREEZE_COMMIT, "registration": OSF_REGISTRATION, "seed": SEED,
                 "package_versions": package_versions, "prob_threshold": PROB_THRESHOLD,
                 "script_sha256": script_sha256(Path(__file__).parent, runner_path), "repo_commit": repo_commit,
                 "run_token": fp.run_token,
                 "run_fingerprint": fp.record()},
        "inputs": plan.seals, "join": info.model_dump(), "sets": plan.sets, "gates": plan.gates,
        "set_gates": plan.set_gates, "designs": plan.designs, "skipped": plan.skipped,
        "decisions": {"H1": d_h1.model_dump(), "H4": d_h4.model_dump(), "against_H1": against_h1(ev_h1),
                      "H1_evidence": ev_h1.model_dump(), "H4_evidence": ev_h4.model_dump(),
                      "H1_wording_given_H2": h1_wording(d_h1, h2_focal["pr_predicted"] if h2_focal else None)},
        "H2": {"focal": h2_focal, "confirmatory": False},
        "within_gene": within,
        "fits": finals, "frequentist": freq,
        "secondary": {"cox_time_to_phase3": cox_models(s1)},
        "descriptive": tables,
    }
    results = _clean(results)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "results.json").write_text(json.dumps(results, indent=2))

    tdir = out_dir / "tables"
    tdir.mkdir(exist_ok=True)
    for name, obj in results["descriptive"].items():
        parts = {name: obj} if isinstance(obj, list) else {f"{name}__{k}": v for k, v in obj.items() if isinstance(v, list)}
        for fname, recs in parts.items():
            (pd.json_normalize(recs) if recs else pd.DataFrame({"note": []})).to_csv(tdir / f"{fname}.csv", index=False)
    main = [r for r in set_rows if r["set_id"] in FOREST_MAIN_SETS]
    supp = [r for r in set_rows if r["set_id"] not in FOREST_MAIN_SETS]
    forest_plot(main, out_dir / "figures" / "forest_main", "Stratum odds ratios, main analysis sets")
    forest_plot(supp, out_dir / "figures" / "forest_supplementary", "Stratum odds ratios, remaining analysis sets")
    write_manifest(out_dir, results["meta"]["script_sha256"], plan.seals, repo_commit)
    return results
