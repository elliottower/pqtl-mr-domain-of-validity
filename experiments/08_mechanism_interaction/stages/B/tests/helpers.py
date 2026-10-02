import pandas as pd

from stage_b.schemas import VARIANT_COLUMNS


def table(rows: list[dict]) -> pd.DataFrame:
    """Canonical regional table from partial rows (defaults fill the rest)."""
    base = {"rsid": "", "chrom": "1", "pos": 0, "ea": "A", "oa": "G", "eaf": 0.3, "beta": 0.1, "se": 0.01,
            "p": 1e-3, "n": 10000.0}
    df = pd.DataFrame([{**base, **r} for r in rows], columns=VARIANT_COLUMNS)
    df["pos"] = df["pos"].astype("int64")
    for c in ("eaf", "beta", "se", "p", "n"):
        df[c] = df[c].astype(float)
    return df
