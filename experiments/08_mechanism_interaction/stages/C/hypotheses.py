"""The only stage A read stage C makes: three columns of `A/output/hypotheses.csv`.

PREREG §Explanation of foreknowledge: the outcome script (C) never reads stage B output.
INTERFACES.md: C reads `hypothesis_id`, `drug_program_ids` and `indication_id` of
`hypotheses.csv`. Any path under a stage B directory, and any column beyond these three, is
refused.
"""
from pathlib import Path

import pandas as pd

from models import BlindingViolation, StageCError

ALLOWED_COLUMNS = ("hypothesis_id", "indication_id", "drug_program_ids")


def assert_not_stage_b(path: Path) -> None:
    parts = path.resolve().parts
    for i in range(len(parts) - 1):
        if parts[i] == "stages" and parts[i + 1] == "B":
            raise BlindingViolation(f"stage C may not read {path}")


def read_hypotheses(path: Path) -> pd.DataFrame:
    assert_not_stage_b(path)
    if path.name != "hypotheses.csv" or path.parent.name != "output" or path.parent.parent.name != "A":
        raise BlindingViolation(f"stage C reads only A/output/hypotheses.csv, not {path}")
    df = pd.read_csv(path, usecols=list(ALLOWED_COLUMNS), dtype=str, keep_default_na=False)
    if df["hypothesis_id"].duplicated().any():
        raise StageCError("duplicate hypothesis_id in hypotheses.csv")
    if (df["drug_program_ids"] == "").any():
        raise StageCError("a hypothesis has no drug programs")
    df["programs"] = [tuple(sorted(set(s.split(";")))) for s in df["drug_program_ids"]]
    return df[["hypothesis_id", "indication_id", "programs"]]
