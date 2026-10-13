"""Q1: where strict evidence is lost (PREREG §Statistical models, Q1).

Descriptive (primary): colocalization classes and B, C distributions by stratum, at the hypothesis
and unique-unit levels. Models (secondary): fractional logistic regressions of B (all completed runs)
and C (completed runs with B >= 0.01), each with and without assay_cis_signal_strength, marginally
standardized, with a 10,000-replicate gene bootstrap checkpointed to JSONL.

Implementation choices the plan leaves open, fixed here: spline knots (10th, 50th, 90th percentiles)
and the standardization percentiles are computed once on the model's analysis population and held
fixed in every replicate; factor levels are those of the analysis population, and a replicate in which
a level is absent (an all-zero column) is a failed replicate under the rank rule, not a stop.
"""
import json
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from followup.common import BOOT_FAIL_MAX, BOOT_REPS, SEED, FollowupError

FACTORS = ("mechanism_class", "instrument_source", "platform", "oncology", "neuro_psych")
SPLINES = ("log10_neff", "assay_cis_signal_strength")
STRATA = ("mechanism_class", "instrument_source", "platform", "oncology", "neuro_psych", "neff_quintile",
          "strength_quintile")
SPLINE_PCTS = (10, 25, 50, 75, 90)
MODELS = {"B_with_strength": ("B", True), "B_without_strength": ("B", False),
          "C_with_strength": ("C", True), "C_without_strength": ("C", False)}
MIN_C_RUNS = 50
CHECKPOINT_EVERY = 250


# ---- descriptive ----------------------------------------------------------------------------------

def _dist(x: pd.Series) -> dict:
    x = x.dropna()
    if x.empty:
        return {"n": 0}
    return {"n": int(len(x)), "q25": float(x.quantile(0.25)), "median": float(x.median()), "q75": float(x.quantile(0.75))}


def describe(df: pd.DataFrame) -> dict:
    """Class counts and B, C distributions per stratum value; units counted once per value they take."""
    out = {}
    for level_name in ("hypothesis", "unit"):
        tables = {}
        for s in STRATA:
            f = df if level_name == "hypothesis" else df.drop_duplicates(["unit_id", s])
            key = f[s].astype(object).where(f[s].notna(), "missing")
            rows = {}
            for value, g in f.groupby(key, sort=True):
                rows[str(value)] = {"n": int(len(g)), "classes": g["coloc_class"].value_counts().to_dict(),
                                    "B": _dist(g["B"]), "C": _dist(g["C"])}
            tables[s] = rows
        overall = df if level_name == "hypothesis" else df.drop_duplicates("unit_id")
        out[level_name] = {"overall": {"n": int(len(overall)), "classes": overall["coloc_class"].value_counts().to_dict(),
                                       "B": _dist(overall["B"]), "C": _dist(overall["C"])},
                           "by_stratum": tables}
    return out


def class_share_figure(df: pd.DataFrame, path_stem: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    key = df["neff_quintile"].astype(object).where(df["neff_quintile"].notna(), "missing").astype(str)
    tab = pd.crosstab(key, df["coloc_class"], normalize="index")
    fig, ax = plt.subplots(figsize=(7, 4))
    tab.plot(kind="bar", stacked=True, ax=ax, width=0.8)
    ax.set_xlabel("outcome-GWAS effective N quintile")
    ax.set_ylabel("share of hypotheses")
    ax.legend(frameon=False, fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(path_stem.with_suffix(f".{ext}"), dpi=200)
    plt.close(fig)


# ---- design ---------------------------------------------------------------------------------------

def rcs_term(x: np.ndarray, knots: tuple[float, float, float]) -> np.ndarray:
    """The nonlinear term of a restricted cubic spline with three knots (Harrell), scaled by (t3 - t1)^2."""
    t1, t2, t3 = knots
    p = lambda v: np.maximum(v, 0.0) ** 3  # noqa: E731
    return (p(x - t1) - p(x - t2) * (t3 - t1) / (t3 - t2) + p(x - t3) * (t2 - t1) / (t3 - t2)) / (t3 - t1) ** 2


@dataclass(frozen=True)
class Spec:
    """What is fixed on the analysis population: response, factor levels, spline knots, percentiles."""
    response: str
    factors: tuple[str, ...]
    levels: dict
    splines: tuple[str, ...]
    knots: dict
    pcts: dict
    drop_platform: bool


def make_spec(pop: pd.DataFrame, response: str, with_strength: bool) -> Spec:
    splines = SPLINES if with_strength else SPLINES[:1]
    levels = {f: sorted(pop[f].astype(str).unique()) for f in FACTORS}
    knots = {s: tuple(float(v) for v in np.percentile(pop[s].astype(float), (10, 50, 90))) for s in splines}
    pcts = {s: {p: float(np.percentile(pop[s].astype(float), p)) for p in SPLINE_PCTS} for s in splines}
    spec = Spec(response, FACTORS, levels, splines, knots, pcts, drop_platform=False)
    X, _ = design(pop, spec)
    if np.linalg.matrix_rank(X) == X.shape[1]:
        return spec
    spec = Spec(response, FACTORS, levels, splines, knots, pcts, drop_platform=True)
    X, names = design(pop, spec)
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise FollowupError(f"{response}: design rank deficient after dropping platform ({names})")
    return spec


def design(frame: pd.DataFrame, spec: Spec) -> tuple[np.ndarray, list[str]]:
    cols, names = [np.ones(len(frame))], ["intercept"]
    for f in spec.factors:
        if f == "platform" and spec.drop_platform:
            continue
        v = frame[f].astype(str).to_numpy()
        for level in spec.levels[f][1:]:
            cols.append((v == level).astype(float))
            names.append(f"{f}[{level}]")
    for s in spec.splines:
        x = frame[s].astype(float).to_numpy()
        cols += [x, rcs_term(x, spec.knots[s])]
        names += [s, f"{s}'"]
    return np.column_stack(cols), names


def fit(frame: pd.DataFrame, spec: Spec):
    """Fractional logit; None when the design is rank deficient, the fit raises or does not converge."""
    X, _ = design(frame, spec)
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            res = sm.GLM(frame[spec.response].astype(float).to_numpy(), X, family=sm.families.Binomial()).fit()
        except (np.linalg.LinAlgError, ValueError, FloatingPointError, sm.tools.sm_exceptions.PerfectSeparationError):
            return None
    if not res.converged or not np.all(np.isfinite(res.params)):
        return None
    return res


def standardize(res, frame: pd.DataFrame, spec: Spec) -> dict[str, float]:
    """Marginal standardization: every row set to each level (or percentile), predicted, averaged;
    contrasts are each level minus the first level (or the 10th percentile)."""
    out = {}
    for f in spec.factors:
        if f == "platform" and spec.drop_platform:
            continue
        means = {}
        for level in spec.levels[f]:
            X, _ = design(frame.assign(**{f: level}), spec)
            means[level] = float(res.predict(X).mean())
        for level, m in means.items():
            out[f"{f}={level}"] = m
            out[f"{f}={level}-{spec.levels[f][0]}"] = m - means[spec.levels[f][0]]
    for s in spec.splines:
        means = {}
        for p, value in spec.pcts[s].items():
            X, _ = design(frame.assign(**{s: value}), spec)
            means[p] = float(res.predict(X).mean())
        for p, m in means.items():
            out[f"{s}@p{p}"] = m
            out[f"{s}@p{p}-p{SPLINE_PCTS[0]}"] = m - means[SPLINE_PCTS[0]]
    return out


# ---- populations and bootstrap -------------------------------------------------------------------

def population(df: pd.DataFrame, model: str) -> pd.DataFrame:
    response, with_strength = MODELS[model]
    pop = df.loc[df["coloc_run"].astype(bool)]
    if response == "C":
        pop = pop.loc[pop["B"] >= 0.01]
    if with_strength:
        pop = pop.loc[pop["assay_cis_signal_strength"].notna()]
    return pop.reset_index(drop=True)


def replicate(pop: pd.DataFrame, genes: np.ndarray, rng: np.random.Generator) -> pd.DataFrame:
    """Genes drawn with replacement, as many as there are genes; every hypothesis of a drawn gene,
    once per draw."""
    drawn = genes[rng.integers(0, len(genes), len(genes))]
    by_gene = pop.groupby("gene_ensembl").indices
    idx = np.concatenate([by_gene[g] for g in drawn])
    return pop.iloc[idx].reset_index(drop=True)


def boot_rng(model: str, rep: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(list(MODELS).index(model), rep)))


def run_bootstrap(pop: pd.DataFrame, spec: Spec, model: str, path: Path, on_commit=None, reps: int = BOOT_REPS) -> None:
    """Append one JSON line per replicate; resumes from the replicates already on disk."""
    done = len(path.read_text().splitlines()) if path.exists() else 0
    genes = np.array(sorted(pop["gene_ensembl"].unique()))
    buf = []
    for rep in range(done, reps):
        sample = replicate(pop, genes, boot_rng(model, rep))
        res = fit(sample, spec)
        buf.append({"rep": rep, "ok": res is not None, "est": standardize(res, sample, spec) if res is not None else None})
        if len(buf) == CHECKPOINT_EVERY or rep == reps - 1:
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("a") as fh:
                fh.write("".join(json.dumps(r) + "\n" for r in buf))
            buf = []
            if on_commit:
                on_commit()


def summarize_bootstrap(point: dict[str, float], path: Path, reps: int = BOOT_REPS) -> dict:
    lines = [json.loads(line) for line in path.read_text().splitlines()]
    if len(lines) != reps or [r["rep"] for r in lines] != list(range(reps)):
        raise FollowupError(f"{path}: {len(lines)} replicates on disk, expected {reps} in order")
    ok = [r["est"] for r in lines if r["ok"]]
    failed = reps - len(ok)
    available = failed / reps <= BOOT_FAIL_MAX
    out = {"replicates": reps, "failed": failed, "intervals_available": available, "estimates": {}}
    for k, v in point.items():
        vals = np.array([r[k] for r in ok if k in r], dtype=float)
        out["estimates"][k] = {"point": v, **({"ci95": [float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))]}
                                              if available and len(vals) else {"ci95": None})}
    return out
