import hashlib
import json
import shutil
from pathlib import Path

import pytest
from prereg.log import append, log_problems
from prereg.plan import plan_of, sha256_of

import v8_run_guard as g
from v8_run_guard import PLAN_SHA256, RunNotAuthorized, check_osf_approved, parse_log, require_run

PREREG = Path(__file__).resolve().parents[3] / "PREREG.md"
TOKEN = "a1b2c3d4e5f60718"
APPROVED = {"data": {"id": "9tzfk", "type": "registrations",
                     "attributes": {"pending_registration_approval": False, "public": True}}}


@pytest.fixture
def plan(tmp_path) -> Path:
    p = tmp_path / "PREREG.md"
    shutil.copyfile(PREREG, p)
    return p


def log(p: Path, event: str, access: str = "nothing run") -> None:
    append(p, "2026-10-02", event, access)


def approved() -> dict:
    return APPROVED


def test_current_prereg_plan_hashes_to_the_frozen_value_with_preregs_own_function():
    text = PREREG.read_text(encoding="utf-8")
    assert sha256_of(plan_of(text)) == PLAN_SHA256
    g.check_plan(text)


def test_current_prereg_authorizes_no_stage_although_it_mentions_stage_b():
    text = PREREG.read_text(encoding="utf-8")
    assert "Pre-stage note (B" in text
    for stage in "ABCD":
        with pytest.raises(RunNotAuthorized, match="RUN_START"):
            require_run(PREREG, stage, TOKEN, osf_fetch=approved,
                        manifests={s: PREREG for s in "ABC"} if stage == "D" else None)


def test_log_lines_parse_into_date_event_and_access():
    entries = parse_log(PREREG.read_text(encoding="utf-8"))
    first, last = entries[0], entries[-1]
    assert (first.date, first.event, first.access, first.chained) == ("2026-09-29", "created", "nothing run", False)
    assert last.event.startswith("Pre-stage note (B, eQTL Catalogue access)") and last.chained
    assert last.access == "nothing run"


@pytest.mark.parametrize("event", [
    "Pre-stage note (B, eQTL route) token=a1b2c3d4e5f60718",
    "stage B started, token a1b2c3d4e5f60718",
    "RUN_START stage=B token=a1b2c3d4e5f60718 (second attempt)",
    "run_start stage=B token=a1b2c3d4e5f60718",
    "RUN_START stage=B  token=a1b2c3d4e5f60718",
    "note: RUN_START stage=B token=a1b2c3d4e5f60718",
])
@pytest.mark.parametrize("stage", ["B", "C"])
def test_a_mention_is_not_a_run_start(plan, event, stage):
    log(plan, event.replace("stage=B", f"stage={stage}").replace("stage B", f"stage {stage}"))
    assert not log_problems(plan.read_text())
    with pytest.raises(RunNotAuthorized):
        require_run(plan, stage, TOKEN)


@pytest.mark.parametrize("stage", ["B", "C"])
def test_an_exact_run_start_line_after_the_predecessor_seals_passes(plan, manifests, stage):
    for s in g.PREDECESSORS[stage]:
        log(plan, f"SEAL stage={s} manifest_sha256={sha(manifests[s])}")
    log(plan, f"RUN_START stage={stage} token={TOKEN}")
    auth = require_run(plan, stage, TOKEN, manifests={"A": manifests["A"]})
    assert (auth.stage, auth.token) == (stage, TOKEN)
    assert auth.run_start_index == len(parse_log(plan.read_text()))
    assert auth.seals == {s: sha(manifests[s]) for s in g.PREDECESSORS[stage]}


def test_the_token_must_match(plan):
    log(plan, f"RUN_START stage=B token={TOKEN}")
    with pytest.raises(RunNotAuthorized, match="no entry"):
        require_run(plan, "B", "ffffffffffffffff")


def test_a_token_logged_for_another_stage_is_refused(plan):
    log(plan, f"RUN_START stage=B token={TOKEN}")
    with pytest.raises(RunNotAuthorized, match="logged for stage B"):
        require_run(plan, "C", TOKEN)


def test_a_token_logged_twice_is_refused(plan):
    log(plan, f"RUN_START stage=B token={TOKEN}")
    log(plan, f"RUN_START stage=B token={TOKEN}")
    with pytest.raises(RunNotAuthorized, match="unique"):
        require_run(plan, "B", TOKEN)


@pytest.mark.parametrize("token", ["", "short", "has space12", "semi;colon12"])
def test_malformed_tokens_are_refused(plan, token):
    with pytest.raises(RunNotAuthorized):
        require_run(plan, "B", token)


def test_a_hand_written_unchained_line_is_refused(plan):
    text = plan.read_text()
    head, fence, tail = text.rpartition("```")
    plan.write_text(head + f"2026-10-02  RUN_START stage=B token={TOKEN}  nothing run\n" + fence + tail)
    with pytest.raises(RunNotAuthorized):
        require_run(plan, "B", TOKEN)


def test_an_edited_plan_is_refused_even_with_a_run_start(plan):
    log(plan, f"RUN_START stage=B token={TOKEN}")
    plan.write_text(plan.read_text().replace("PP.H4 ≥ 0.80", "PP.H4 ≥ 0.70", 1))
    with pytest.raises(RunNotAuthorized, match="hashes to"):
        require_run(plan, "B", TOKEN)


def test_an_edited_log_entry_is_refused(plan):
    log(plan, f"RUN_START stage=B token={TOKEN}")
    plan.write_text(plan.read_text().replace("osf draft 6abd", "osf draft 7abd", 1))
    with pytest.raises(RunNotAuthorized, match="prereg check"):
        require_run(plan, "B", TOKEN)


# ---- stage A: OSF approval ----------------------------------------------------------------------

def test_stage_a_passes_with_run_start_and_approved_registration(plan):
    log(plan, f"RUN_START stage=A token={TOKEN}")
    assert require_run(plan, "A", TOKEN, osf_fetch=approved).stage == "A"


@pytest.mark.parametrize("attrs", [
    {"pending_registration_approval": True, "public": False},
    {"pending_registration_approval": False, "public": False},
    {"pending_registration_approval": True, "public": True},
    {"public": True},
    {"pending_registration_approval": "false", "public": True},
])
def test_stage_a_refuses_an_unapproved_registration(plan, attrs):
    log(plan, f"RUN_START stage=A token={TOKEN}")
    doc = {"data": {"id": "9tzfk", "attributes": attrs}}
    with pytest.raises(RunNotAuthorized, match="OSF"):
        require_run(plan, "A", TOKEN, osf_fetch=lambda: doc)


def test_stage_a_refuses_another_registration_id():
    with pytest.raises(RunNotAuthorized, match="not 9tzfk"):
        check_osf_approved({"data": {"id": "xxxxx", "attributes": APPROVED["data"]["attributes"]}})


def test_stage_a_refuses_without_run_start_even_when_approved(plan):
    with pytest.raises(RunNotAuthorized, match="RUN_START"):
        require_run(plan, "A", TOKEN, osf_fetch=approved)


class FakeResponse:
    def __init__(self, body: bytes):
        self.body = body

    def read(self):
        return self.body

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


@pytest.mark.parametrize("pat", [None, "secret-pat-value"])
def test_default_osf_fetch_queries_the_registration_endpoint(plan, monkeypatch, pat):
    seen = []

    def fake_urlopen(req, timeout):
        seen.append(req)
        return FakeResponse(json.dumps(APPROVED).encode())

    monkeypatch.setattr(g, "urlopen", fake_urlopen)
    if pat is None:
        monkeypatch.delenv("OSF_PAT", raising=False)
    else:
        monkeypatch.setenv("OSF_PAT", pat)
    log(plan, f"RUN_START stage=A token={TOKEN}")
    require_run(plan, "A", TOKEN)
    (req,) = seen
    assert req.full_url == "https://api.osf.io/v2/registrations/9tzfk/"
    headers = {k.lower(): v for k, v in req.header_items()}
    assert ("authorization" in headers) == (pat is not None)
    assert "@" not in " ".join(headers.values())


def test_default_osf_fetch_turns_an_http_error_into_a_refusal(plan, monkeypatch):
    def fail(req, timeout):
        raise OSError("HTTP Error 401: Unauthorized")

    monkeypatch.setattr(g, "urlopen", fail)
    log(plan, f"RUN_START stage=A token={TOKEN}")
    with pytest.raises(RunNotAuthorized, match="could not be read"):
        require_run(plan, "A", TOKEN)


# ---- stage D: seals ------------------------------------------------------------------------------

@pytest.fixture
def manifests(tmp_path) -> dict[str, Path]:
    out = {}
    for s in "ABC":
        p = tmp_path / s / "output" / "MANIFEST.tsv"
        p.parent.mkdir(parents=True)
        p.write_text(f"path\trows\nfile_{s}.csv\t{ord(s)}\n")
        out[s] = p
    return out


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def seal_all(plan: Path, manifests: dict[str, Path]) -> None:
    for s in "ABC":
        log(plan, f"SEAL stage={s} manifest_sha256={sha(manifests[s])}", "results not opened")


def test_stage_d_passes_with_all_three_seals_before_its_run_start(plan, manifests):
    seal_all(plan, manifests)
    log(plan, f"RUN_START stage=D token={TOKEN}", "results not opened")
    auth = require_run(plan, "D", TOKEN, manifests=manifests)
    assert auth.seals == {s: sha(manifests[s]) for s in "ABC"}


@pytest.mark.parametrize("missing", ["A", "B", "C"])
def test_stage_d_refuses_a_missing_seal(plan, manifests, missing):
    for s in "ABC":
        if s != missing:
            log(plan, f"SEAL stage={s} manifest_sha256={sha(manifests[s])}")
    log(plan, f"stage {missing} sealed, manifest_sha256={sha(manifests[missing])}")
    log(plan, f"RUN_START stage=D token={TOKEN}")
    with pytest.raises(RunNotAuthorized, match=f"SEAL stage={missing}"):
        require_run(plan, "D", TOKEN, manifests=manifests)


def test_stage_d_refuses_a_seal_logged_after_its_run_start(plan, manifests):
    log(plan, f"SEAL stage=A manifest_sha256={sha(manifests['A'])}")
    log(plan, f"SEAL stage=B manifest_sha256={sha(manifests['B'])}")
    log(plan, f"RUN_START stage=D token={TOKEN}")
    log(plan, f"SEAL stage=C manifest_sha256={sha(manifests['C'])}")
    with pytest.raises(RunNotAuthorized, match="SEAL stage=C"):
        require_run(plan, "D", TOKEN, manifests=manifests)


def test_stage_d_refuses_a_manifest_that_is_not_the_sealed_one(plan, manifests):
    seal_all(plan, manifests)
    log(plan, f"RUN_START stage=D token={TOKEN}")
    manifests["B"].write_text(manifests["B"].read_text() + "extra.csv\t1\n")
    with pytest.raises(RunNotAuthorized, match="stage B"):
        require_run(plan, "D", TOKEN, manifests=manifests)


def test_stage_d_uses_the_last_seal_of_a_stage(plan, manifests):
    log(plan, f"SEAL stage=A manifest_sha256={'0' * 64}")
    seal_all(plan, manifests)
    log(plan, f"RUN_START stage=D token={TOKEN}")
    require_run(plan, "D", TOKEN, manifests=manifests)
    log(plan, f"SEAL stage=A manifest_sha256={'0' * 64}")
    log(plan, "RUN_START stage=D token=second-run-01")
    with pytest.raises(RunNotAuthorized, match="stage A"):
        require_run(plan, "D", "second-run-01", manifests=manifests)


# ---- stages B and C: predecessor seals (A -> B -> C -> D) ---------------------------------------

def a_only(manifests: dict[str, Path]) -> dict[str, Path]:
    return {"A": manifests["A"]}


@pytest.mark.parametrize("stage", ["B", "C"])
def test_b_and_c_refuse_without_a_seal_of_a(plan, manifests, stage):
    if stage == "C":
        log(plan, f"SEAL stage=B manifest_sha256={sha(manifests['B'])}")
    log(plan, f"stage A sealed, manifest_sha256={sha(manifests['A'])}")
    log(plan, f"RUN_START stage={stage} token={TOKEN}")
    with pytest.raises(RunNotAuthorized, match="SEAL stage=A"):
        require_run(plan, stage, TOKEN, manifests=a_only(manifests))


def test_c_refuses_without_a_seal_of_b(plan, manifests):
    log(plan, f"SEAL stage=A manifest_sha256={sha(manifests['A'])}")
    log(plan, f"RUN_START stage=C token={TOKEN}")
    with pytest.raises(RunNotAuthorized, match="SEAL stage=B"):
        require_run(plan, "C", TOKEN, manifests=a_only(manifests))


@pytest.mark.parametrize("stage,late", [("B", "A"), ("C", "A"), ("C", "B")])
def test_a_predecessor_seal_logged_after_the_run_start_is_refused(plan, manifests, stage, late):
    for s in g.PREDECESSORS[stage]:
        if s != late:
            log(plan, f"SEAL stage={s} manifest_sha256={sha(manifests[s])}")
    log(plan, f"RUN_START stage={stage} token={TOKEN}")
    log(plan, f"SEAL stage={late} manifest_sha256={sha(manifests[late])}")
    with pytest.raises(RunNotAuthorized, match=f"SEAL stage={late}"):
        require_run(plan, stage, TOKEN, manifests=a_only(manifests))


@pytest.mark.parametrize("stage", ["B", "C"])
def test_b_and_c_refuse_an_a_manifest_that_is_not_the_sealed_one(plan, manifests, stage):
    for s in g.PREDECESSORS[stage]:
        log(plan, f"SEAL stage={s} manifest_sha256={sha(manifests[s])}")
    log(plan, f"RUN_START stage={stage} token={TOKEN}")
    manifests["A"].write_text(manifests["A"].read_text() + "extra.csv\t1\n")
    with pytest.raises(RunNotAuthorized, match="stage A"):
        require_run(plan, stage, TOKEN, manifests=a_only(manifests))


@pytest.mark.parametrize("stage", ["B", "C"])
def test_the_last_seal_before_the_start_counts(plan, manifests, stage):
    log(plan, f"SEAL stage=A manifest_sha256={'0' * 64}")
    for s in g.PREDECESSORS[stage]:
        log(plan, f"SEAL stage={s} manifest_sha256={sha(manifests[s])}")
    log(plan, f"RUN_START stage={stage} token={TOKEN}")
    require_run(plan, stage, TOKEN, manifests=a_only(manifests))
    log(plan, f"SEAL stage=A manifest_sha256={'0' * 64}")
    log(plan, f"RUN_START stage={stage} token=second-run-01")
    with pytest.raises(RunNotAuthorized, match="stage A"):
        require_run(plan, stage, "second-run-01", manifests=a_only(manifests))


@pytest.mark.parametrize("stage,given", [("B", {}), ("B", {"A", "B"}), ("C", {}), ("C", {"A", "B"}),
                                         ("D", {"A", "B"}), ("A", {"A"})])
def test_a_stage_is_given_exactly_the_manifests_it_reads(plan, manifests, stage, given):
    for s in g.PREDECESSORS[stage]:
        log(plan, f"SEAL stage={s} manifest_sha256={sha(manifests[s])}")
    log(plan, f"RUN_START stage={stage} token={TOKEN}")
    with pytest.raises(ValueError, match="needs the MANIFEST.tsv"):
        require_run(plan, stage, TOKEN, manifests={s: manifests[s] for s in given}, osf_fetch=approved)


def test_the_ordered_chain_authorizes_each_stage_in_turn(plan, manifests):
    log(plan, "RUN_START stage=A token=token-a-0001")
    assert require_run(plan, "A", "token-a-0001", osf_fetch=approved).seals == {}
    log(plan, f"SEAL stage=A manifest_sha256={sha(manifests['A'])}")
    log(plan, "RUN_START stage=B token=token-b-0001")
    assert require_run(plan, "B", "token-b-0001", manifests=a_only(manifests)).seals == {"A": sha(manifests["A"])}
    log(plan, f"SEAL stage=B manifest_sha256={sha(manifests['B'])}")
    log(plan, "RUN_START stage=C token=token-c-0001")
    assert set(require_run(plan, "C", "token-c-0001", manifests=a_only(manifests)).seals) == {"A", "B"}
    log(plan, f"SEAL stage=C manifest_sha256={sha(manifests['C'])}")
    log(plan, "RUN_START stage=D token=token-d-0001")
    assert require_run(plan, "D", "token-d-0001", manifests=manifests).seals == {s: sha(manifests[s]) for s in "ABC"}
    assert not log_problems(plan.read_text())


# ---- commit identity (git is always a fake here) ------------------------------------------------

HEAD = "0123456789abcdef0123456789abcdef01234567"
STATUS = ("status", "--porcelain", "--", "experiments/08_mechanism_interaction/stages",
          "experiments/08_mechanism_interaction/PREREG.md")


class FakeGit:
    def __init__(self, status: str = "", head: str = HEAD + "\n", fail: bool = False):
        self.status, self.head, self.fail, self.calls = status, head, fail, []

    def __call__(self, repo: Path, *args: str) -> str:
        self.calls.append((repo, args))
        if self.fail:
            raise RunNotAuthorized("git is not installed")
        return self.status if args[0] == "status" else self.head


def test_clean_commit_asks_git_for_the_status_of_the_stage_code_and_the_plan_then_for_head(tmp_path):
    git = FakeGit()
    assert g.clean_commit(tmp_path, git) == HEAD
    assert git.calls == [(tmp_path, STATUS), (tmp_path, ("rev-parse", "HEAD"))]


@pytest.mark.parametrize("status", [
    " M experiments/08_mechanism_interaction/PREREG.md\n",
    "?? experiments/08_mechanism_interaction/stages/\n",
    "M  experiments/08_mechanism_interaction/stages/A/run_stage_a.py\n"
    " D experiments/08_mechanism_interaction/stages/C/rules.py\n",
    "A  experiments/08_mechanism_interaction/stages/B/output/unit_plan.json\n",
])
def test_clean_commit_refuses_a_modified_staged_or_untracked_path_without_reading_head(tmp_path, status):
    git = FakeGit(status=status)
    with pytest.raises(RunNotAuthorized, match=r"uncommitted path\(s\)") as err:
        g.clean_commit(tmp_path, git)
    assert status.splitlines()[0].strip() in str(err.value)
    assert [args[0] for _, args in git.calls] == ["status"]


@pytest.mark.parametrize("head", ["", "HEAD\n", "0123456\n", HEAD.upper() + "\n", HEAD + "0\n"])
def test_clean_commit_refuses_a_head_that_is_not_a_full_sha(tmp_path, head):
    with pytest.raises(RunNotAuthorized, match="40-character commit sha"):
        g.clean_commit(tmp_path, FakeGit(head=head))


def test_baked_commit_is_the_clean_commit_or_empty_and_never_raises(tmp_path):
    assert g.baked_commit(tmp_path, FakeGit()) == HEAD
    assert g.baked_commit(tmp_path, FakeGit(status=" M experiments/08_mechanism_interaction/PREREG.md\n")) == ""
    assert g.baked_commit(tmp_path, FakeGit(fail=True)) == ""


def test_run_commit_takes_the_baked_value_and_does_not_ask_git(tmp_path):
    git = FakeGit(status="?? experiments/08_mechanism_interaction/stages/\n")
    assert g.run_commit(tmp_path, {g.COMMIT_ENV: HEAD}, git) == HEAD
    assert g.run_commit(None, {g.COMMIT_ENV: HEAD}, git) == HEAD
    assert git.calls == []


@pytest.mark.parametrize("baked", ["", "dirty", HEAD[:12]])
def test_run_commit_refuses_an_image_built_from_an_unclean_tree_even_beside_a_clean_repository(tmp_path, baked):
    git = FakeGit()
    with pytest.raises(RunNotAuthorized, match="not built from a clean commit"):
        g.run_commit(tmp_path, {g.COMMIT_ENV: baked}, git)
    assert git.calls == []


def test_run_commit_without_a_baked_value_is_the_clean_commit_of_the_repository(tmp_path):
    assert g.run_commit(tmp_path, {}, FakeGit()) == HEAD
    with pytest.raises(RunNotAuthorized, match="uncommitted"):
        g.run_commit(tmp_path, {}, FakeGit(status=" M experiments/08_mechanism_interaction/PREREG.md\n"))
    with pytest.raises(RunNotAuthorized, match="no commit identity"):
        g.run_commit(None, {}, FakeGit())


def test_run_commit_reads_the_process_environment_by_default(tmp_path, monkeypatch):
    monkeypatch.setenv(g.COMMIT_ENV, HEAD)
    assert g.run_commit(None) == HEAD
    monkeypatch.setenv(g.COMMIT_ENV, "")
    with pytest.raises(RunNotAuthorized, match="not built from a clean commit"):
        g.run_commit(tmp_path, git=FakeGit())


class Completed:
    def __init__(self, returncode: int, stdout: str = "", stderr: str = ""):
        self.returncode, self.stdout, self.stderr = returncode, stdout, stderr


def test_git_output_runs_git_in_the_repository_and_returns_stdout(tmp_path, monkeypatch):
    seen = []

    def fake_run(cmd, **kwargs):
        seen.append(cmd)
        return Completed(0, stdout="" if "status" in cmd else HEAD + "\n")

    monkeypatch.setattr(g.subprocess, "run", fake_run)
    assert g.git_output(tmp_path, "rev-parse", "HEAD") == HEAD + "\n"
    assert seen == [["git", "-C", str(tmp_path), "rev-parse", "HEAD"]]
    assert g.clean_commit(tmp_path) == HEAD          # the default git: an empty status and this HEAD
    assert seen[1:] == [["git", "-C", str(tmp_path), *STATUS], ["git", "-C", str(tmp_path), "rev-parse", "HEAD"]]


def test_git_output_turns_a_failing_or_missing_git_into_a_refusal(tmp_path, monkeypatch):
    monkeypatch.setattr(g.subprocess, "run", lambda cmd, **kw: Completed(128, stderr="fatal: not a git repository"))
    with pytest.raises(RunNotAuthorized, match="not a git repository"):
        g.git_output(tmp_path, "status")

    def missing(cmd, **kw):
        raise FileNotFoundError("git")

    monkeypatch.setattr(g.subprocess, "run", missing)
    with pytest.raises(RunNotAuthorized, match="could not run"):
        g.clean_commit(tmp_path)
    assert g.baked_commit(tmp_path) == ""
