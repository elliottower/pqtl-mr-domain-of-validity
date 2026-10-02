import ast
from pathlib import Path

import pandas as pd
import pytest

from hypotheses import ALLOWED_COLUMNS, assert_not_stage_b, read_hypotheses
from models import BlindingViolation

STAGE_C = Path(__file__).resolve().parents[1]
# Stage C's own modules. The Modal image also bakes the shared guard modules of stages/run_guard
# (v8_*.py) into this directory; they are not stage C code and have their own suite.
MODULES = sorted(p for p in STAGE_C.glob("*.py") if not p.name.startswith("v8_"))
ALLOWED_IMPORTS = {"argparse", "calendar", "collections", "csv", "dataclasses", "datetime", "enum",
                   "hashlib", "io", "json", "pathlib", "re", "shutil", "sqlite3", "urllib", "zipfile",
                   "pandas", "pyarrow", "pydantic", "requests", "v8_manifest", "v8_run_guard",
                   "aact_download", "chembl", "ctgov", "diagnostics", "hypotheses", "models", "ot_phase", "rules",
                   "modal", "run_stage_c", "v8_test_report"}  # the thin wrapper modal_stage_c.py and its `tests` function
FORBIDDEN_SUBSTRINGS = ("B/output", "stages/B", "evidence", "coloc", "pp_h", "genetic_direction")


def docstring_nodes(tree: ast.AST) -> set[int]:
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.ClassDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant):
                out.add(id(first.value))
    return out


def test_modules_exist():
    assert {p.stem for p in MODULES} >= {"models", "rules", "ctgov", "aact_download", "hypotheses", "ot_phase",
                                         "chembl", "run_stage_c"}


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_module_imports_nothing_outside_the_allowed_set(path):
    tree = ast.parse(path.read_text())
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            names.add((node.module or "").split(".")[0])
    assert names <= ALLOWED_IMPORTS, names - ALLOWED_IMPORTS
    assert "sys" not in names


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_module_code_names_no_stage_b_path_or_evidence_field(path):
    tree = ast.parse(path.read_text())
    docs = docstring_nodes(tree)
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docs:
            low = node.value.lower()
            assert not any(s.lower() in low for s in FORBIDDEN_SUBSTRINGS), (path.name, node.value)
        if isinstance(node, ast.Name):
            assert "evidence" not in node.id.lower() and "coloc" not in node.id.lower()


@pytest.mark.parametrize("path", MODULES, ids=lambda p: p.name)
def test_imports_are_at_module_top(path):
    tree = ast.parse(path.read_text())
    top = {id(n) for n in tree.body}
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            assert id(node) in top, f"{path.name}:{node.lineno}"


def test_stage_b_paths_are_refused(tmp_path):
    b = tmp_path / "stages" / "B" / "output"
    b.mkdir(parents=True)
    with pytest.raises(BlindingViolation):
        assert_not_stage_b(b / "evidence.csv")
    with pytest.raises(BlindingViolation):
        read_hypotheses(b / "hypotheses.csv")


def test_only_a_output_hypotheses_csv_is_read(tmp_path):
    other = tmp_path / "stages" / "A" / "output" / "funnel.csv"
    other.parent.mkdir(parents=True)
    other.write_text("x\n1\n")
    with pytest.raises(BlindingViolation):
        read_hypotheses(other)


def test_read_hypotheses_keeps_three_columns_only(tmp_path):
    p = tmp_path / "stages" / "A" / "output" / "hypotheses.csv"
    p.parent.mkdir(parents=True)
    pd.DataFrame({"hypothesis_id": ["h1"], "gene_symbol": ["PCSK9"], "indication_id": ["EFO_1"],
                  "mechanism_class": ["aligned"], "drug_program_ids": ["CHEMBL2;CHEMBL1;CHEMBL2"],
                  "S": ["1"]}).to_csv(p, index=False)
    df = read_hypotheses(p)
    assert list(df.columns) == ["hypothesis_id", "indication_id", "programs"]
    assert df.at[0, "programs"] == ("CHEMBL1", "CHEMBL2")
    assert set(ALLOWED_COLUMNS) == {"hypothesis_id", "indication_id", "drug_program_ids"}


def test_run_succeeds_beside_an_unreadable_stage_b_file(tmp_path, e2e):
    b = tmp_path / "stages" / "B" / "output" / "evidence.csv"
    assert b.exists()
    with pytest.raises(OSError):
        b.read_bytes()
    assert e2e["rows"]
