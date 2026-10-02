"""The dry run's scenario (synthetic_unit.scenario) as a test: the synthetic sources served by
fake_remote.py, the real collect and analyze code, tabix, bcftools and the R backend, with faults
injected. Runs where htslib and R coloc exist (the stage B image)."""
import shutil

from fake_remote import FakeRemote
from synthetic_unit import build_world, scenario
from test_coloc import R_LIBS, needs_r
from test_pipeline import TOOLS, fp

import pytest

from stage_b.coloc_backend import RscriptColoc
from stage_b.pipeline import DirStore, process_unit

needs_htslib = pytest.mark.skipif(not all(shutil.which(t) for t in ("tabix", "bgzip", "bcftools")),
                                  reason="needs tabix, bgzip and bcftools: runs in the stage B image")


@needs_r
@needs_htslib
def test_collect_then_analyze_gives_the_known_answers_whatever_the_links_tokens_and_connections_do(tmp_path):
    def analyze(unit, fetcher, units_root):
        return process_unit(unit, fetcher, RscriptColoc(r_libs=R_LIBS), DirStore(units_root / unit.unit_key), fp(unit), TOOLS)

    with FakeRemote() as remote:
        report = scenario(build_world(remote, tmp_path / "work"), tmp_path / "vol", analyze)
    assert report["passed"], {k: v for k, v in report["checks"].items() if not v}
    assert report["events"]["expired_folder_token"]["kind"] == "auth"
    assert report["events"]["missing_collect_record"]["error_class"] == "CollectError"
