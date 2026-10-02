"""Forest plots of the posterior stratum odds ratios by analysis set (PREREG §Other planned
analysis): the main plot shows S1, S2, S3, S4, S6, S7, S8 and S13; the supplementary plot the
remaining sets."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter, NullFormatter  # noqa: E402

COLORS = {"aligned": "#1b6ca8", "blocking": "#c0392b", "neuro_psych": "#6c3483", "other_indications": "#7f8c8d"}


def forest_plot(rows: list[dict], path: Path, title: str) -> bool:
    """rows: {"set_id", "stratum_odds_ratios": {stratum: {or_median, or_q05, or_q95}} | None}.
    Sets without a fit are listed with their status in place of an interval. Returns False when
    no row has an estimate."""
    entries = []
    for r in rows:
        ors = r.get("stratum_odds_ratios")
        if not ors:
            entries.append((f"{r['set_id']} ({r.get('status', 'not fitted')})", None, None))
            continue
        for stratum, v in ors.items():
            entries.append((f"{r['set_id']} {stratum}", stratum, v))
    if not any(v for _, _, v in entries):
        return False
    fig, ax = plt.subplots(figsize=(6.5, 0.32 * len(entries) + 1.2))
    for i, (label, stratum, v) in enumerate(entries):
        y = len(entries) - 1 - i
        if v is None:
            continue
        ax.errorbar(v["or_median"], y, xerr=[[v["or_median"] - v["or_q05"]], [v["or_q95"] - v["or_median"]]],
                    fmt="o", color=COLORS.get(stratum, "black"), ms=4, capsize=2)
    ax.axvline(1.0, color="0.6", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(range(len(entries)))
    ax.set_yticklabels([e[0] for e in reversed(entries)], fontsize=8)
    ax.set_xlabel("odds ratio for supportive evidence (posterior median, 90% credible interval)")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".png"), dpi=200)
    plt.close(fig)
    return True
