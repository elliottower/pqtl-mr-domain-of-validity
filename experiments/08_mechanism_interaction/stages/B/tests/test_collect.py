"""The collect phase and the error classification against a local HTTP server (fake_remote.py):
resumable downloads, one record per file, absence against every other fault, and the analyze-side
reads of collected files. The server is the synthetic one of synthetic_unit.py without its
tabix-indexed files, so nothing here needs htslib."""
import ast
import gzip
import json
import subprocess
import zlib
from pathlib import Path

import pytest
import requests
from fake_remote import FakeRemote, FakeSynapse
from synthetic_unit import (CHROM, DECODE_KEY, DISTINCT, ENSEMBL_EXPECTED, ENSEMBL_REPLIES, EXCLUDED_INDEX, GENE, GRCH37,
                            N_VARIANTS, NO_DIR, NO_FILE, OID, SECOND_GENE, SENTINEL_INDEX, SHARED, TAR_NAME, TOKENS, build_world,
                            ensembl_outcomes, rsid, same_size_other_bytes, second_gene_unit)
from test_pipeline import COLLECT, H4, SENTINEL, TOOLS, FakeFetcher, StubBackend, fp

from stage_b import collect as collect_module
from stage_b import remote as remote_module
from stage_b.assemble import collect_unit_dir
from stage_b.checkpoint import FINGERPRINT_NAME, volume_collect_digest
from stage_b.collect import PURGED_NAME, collect_one, collect_tasks, collected_file, download, purge_raw, read_record, record_path
from stage_b.fetch import (DECODE_FILE_URL, ENSEMBL_UNKNOWN_VARIATION, ENSEMBL_UNKNOWN_VEP_ID, RemoteFile, RemoteSources,
                           ensembl_unknown_id, https_path)
from stage_b.pipeline import DirStore, process_unit
from stage_b.remote import attempt, check_status, classify, http
from stage_b.schemas import (WINDOW_WIDE, CollectError, CollectTask, InputContractError, OutcomeSpec, RetryableSourceError,
                             SourceAbsent, StaleCheckpointError)
from stage_b.status import error_of, marked, status_report, volume_state

HERE = Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def no_waiting(monkeypatch):
    monkeypatch.setattr(remote_module, "WAIT_S", 0.0)
    monkeypatch.setattr(collect_module, "CHUNK_BYTES", 512)


@pytest.fixture
def world(tmp_path):
    with FakeRemote() as remote:
        yield build_world(remote, tmp_path / "work", indexed=False)


def task(world, source: str, key: str) -> CollectTask:
    return next(t for t in collect_tasks(world.units) if (t.source, t.key) == (source, key))


def collect(world, root, source: str, key: str, **kw):
    return collect_one(task(world, source, key), world.sources(root / "cache", **kw), root, checkpoint_bytes=2048)


# ---- classification ---------------------------------------------------------------------------------

@pytest.mark.parametrize("code,kind", [(401, "auth"), (403, "auth"), (429, "rate_limit"), (500, "server"), (502, "server"),
                                       (503, "server"), (400, "protocol"), (416, "protocol"), (302, "protocol")])
def test_a_status_that_is_not_404_or_410_is_retryable_with_its_class(code, kind):
    with pytest.raises(RetryableSourceError) as err:
        check_status(code, "source x")
    assert err.value.kind == kind and str(err.value) == f"[{kind}] source x: HTTP {code}"


@pytest.mark.parametrize("code", [404, 410])
def test_only_404_and_410_are_a_definitive_absence(code):
    with pytest.raises(SourceAbsent, match=f"HTTP {code}"):
        check_status(code, "source x")
    assert check_status(200, "source x") is None and check_status(206, "source x") is None


@pytest.mark.parametrize("error,expected", [
    (requests.ConnectionError("https://h/folder/SECRET/f"), "connection"), (requests.ReadTimeout("https://h/SECRET"), "timeout"),
    (requests.exceptions.ChunkedEncodingError("broken"), "connection"), (ConnectionResetError("reset"), "connection"),
    (TimeoutError("slow"), "timeout"), (subprocess.TimeoutExpired("tabix", 5), "timeout"),
    (subprocess.CalledProcessError(1, ["tabix", "https://h/SECRET"]), "connection"),
    (EOFError("Compressed file ended before the end-of-stream marker was reached"), "corrupt"),
    (zlib.error("invalid block"), "corrupt"), (gzip.BadGzipFile("not gzip"), "corrupt"),
])
def test_transport_and_stream_faults_are_retryable_and_their_text_never_reaches_the_message(error, expected):
    got = classify(error, "source x")
    assert isinstance(got, RetryableSourceError) and got.kind == expected
    assert "SECRET" not in str(got) and "https" not in str(got)


def test_attempt_retries_a_retryable_fault_and_raises_auth_and_absence_at_once():
    calls = []

    def flaky(errors):
        def call():
            calls.append(1)
            if errors:
                raise errors.pop(0)
            return "ok"
        return call

    assert attempt(flaky([requests.ConnectionError("x"), RetryableSourceError("server", "x: HTTP 503")]), "x") == "ok"
    assert len(calls) == 3
    for error, n in ((RetryableSourceError("auth", "x: HTTP 403"), 1), (SourceAbsent("x: HTTP 404"), 1)):
        calls.clear()
        with pytest.raises(type(error)):
            attempt(flaky([error] * 9), "x")
        assert len(calls) == n
    calls.clear()
    with pytest.raises(RetryableSourceError) as err:
        attempt(flaky([RetryableSourceError("rate_limit", "x: HTTP 429")] * 9), "x")
    assert err.value.kind == "rate_limit" and len(calls) == remote_module.ATTEMPTS
    with pytest.raises(KeyError):                         # a programming error is not a source fault
        attempt(flaky([KeyError("bug")]), "x")


# ---- collect: one resumable download per file --------------------------------------------------------

def test_tasks_are_one_per_file_whatever_the_number_of_units_that_read_it(world):
    tasks = [(t.source, t.key) for t in collect_tasks(world.units)]
    assert tasks == [("decode", "0_0"), ("decode_smp", "0_0"), ("gwas_catalog", DISTINCT), ("gwas_catalog", "GCST90000002"),
                     ("gwas_catalog", NO_DIR), ("gwas_catalog", NO_FILE), ("ukbppp", OID), ("ukbppp_rsid_map", CHROM)]
    readers = [u.unit_key for u in world.units if any(o.accession == DISTINCT for o in u.outcomes)]
    assert len(readers) == 2 and tasks.count(("gwas_catalog", DISTINCT)) == 1
    decode = task(world, "decode", "0_0")
    assert (decode.name, decode.size, decode.etag) == (DECODE_KEY, len(world.remote.files[f"/decode/{DECODE_KEY}"]),
                                                       world.remote.etag(f"/decode/{DECODE_KEY}"))


def test_a_file_the_units_of_two_genes_read_is_downloaded_once_and_each_unit_keeps_its_own_checkpoint(world, tmp_path):
    root, path = tmp_path / "vol", f"/decode/{DECODE_KEY}"
    first = decode_unit(world)
    second = second_gene_unit(first)
    assert (first.gene_ensembl, second.gene_ensembl) == (GENE, SECOND_GENE)
    assert (first.unit_key, second.unit_key) == (f"decode__0_0__{GENE}", f"decode__0_0__{SECOND_GENE}")
    assert (first.sentinel, first.pqtl_locator, first.pqtl_listing) == (second.sentinel, second.pqtl_locator, second.pqtl_listing)
    tasks = collect_tasks([first, second])
    assert tasks == collect_tasks([first]) and [(t.source, t.key) for t in tasks].count(("decode", "0_0")) == 1
    for t in tasks:
        collect_one(t, world.sources(root / "cache"), root, checkpoint_bytes=2048)
    assert len([p for _m, p, _r in world.remote.log if p == path]) == 1          # one download for the two units
    fetcher, mark = world.fetcher(root), len(world.remote.log)
    regions = [fetcher.pqtl_region(u, CHROM, 50_000_000, WINDOW_WIDE) for u in (first, second)]
    assert world.remote.log[mark:] == [] and regions[0].equals(regions[1]) and len(regions[0]) == N_VARIANTS - 1
    digests = {u.unit_key: volume_collect_digest(root, u) for u in (first, second)}
    assert fp(first, collect=digests[first.unit_key]) != fp(second, collect=digests[second.unit_key])
    assert fp(first) != fp(second)                                               # the unit record alone tells them apart
    stores = {u.unit_key: DirStore(root / "units" / u.unit_key) for u in (first, second)}
    stores[first.unit_key].bind(fp(first), TOOLS, COLLECT)
    with pytest.raises(StaleCheckpointError):                                    # one gene's checkpoint is not the other's
        DirStore(root / "units" / first.unit_key).bind(fp(second), TOOLS, COLLECT)
    stores[second.unit_key].bind(fp(second), TOOLS, COLLECT)
    assert sorted(p.name for p in (root / "units").iterdir()) == sorted(stores)


def test_a_collected_file_is_the_source_file_and_its_record_holds_no_token(world, tmp_path):
    root = tmp_path / "vol"
    record = collect(world, root, "decode", "0_0")
    data = world.remote.files[f"/decode/{DECODE_KEY}"]
    assert (root / record.path).read_bytes() == data and record.path == f"raw/decode/{DECODE_KEY}"
    assert (record.status, record.bytes, record.etag, record.last_modified) == (
        "collected", len(data), world.remote.etag(f"/decode/{DECODE_KEY}"), "Wed, 30 Sep 2026 00:00:00 GMT")
    assert record.source_url == f"{world.remote.base}/decode/s3/download?token=<token>&file={DECODE_KEY}"
    on_disk = record_path(root, "decode", "0_0").read_text()
    assert TOKENS["decode"] not in on_disk and "signed" not in on_disk and read_record(root, "decode", "0_0") == record
    assert [p for _m, p, _r in world.remote.log] == ["/decode/s3/download", f"/decode/{DECODE_KEY}"]   # the endpoint, then its redirect
    assert collected_file(root, record) == root / record.path
    assert not list((root / "raw" / "decode").glob("*.part*"))


def test_a_stream_broken_in_every_attempt_raises_and_the_next_call_resumes_by_range(world, tmp_path):
    root, path = tmp_path / "vol", f"/decode/{DECODE_KEY}"
    size = len(world.remote.files[path])
    world.remote.faults[path] = [("break", size // 8)] * remote_module.ATTEMPTS
    with pytest.raises(RetryableSourceError) as err:
        collect(world, root, "decode", "0_0")
    part = root / "raw" / "decode" / f"{DECODE_KEY}.part"
    held = part.stat().st_size
    assert err.value.kind == "connection" and 0 < held < size
    assert part.read_bytes() == world.remote.files[path][:held]
    assert not record_path(root, "decode", "0_0").exists()            # nothing recorded, and not "unavailable"
    mark = len(world.remote.log)
    record = collect(world, root, "decode", "0_0")
    assert [rng for _m, p, rng in world.remote.log[mark:] if p == path] == [f"bytes={held}-"]
    assert (root / record.path).read_bytes() == world.remote.files[path] and record.status == "collected"


def test_an_expired_folder_token_raises_auth_and_a_reissued_link_collects_the_same_file(world, tmp_path):
    root = tmp_path / "vol"
    first = collect(world, root, "decode", "0_0")
    world.remote.tokens["decode"] = "folder-token-two"                 # the old link no longer works
    other = tmp_path / "other"
    with pytest.raises(RetryableSourceError) as err:
        collect(world, other, "decode", "0_0")                         # still the first token
    assert err.value.kind == "auth" and not (other / "collect").exists() and not (other / "raw").exists()
    assert len([m for m, p, _ in world.remote.log if p == "/decode/s3/download"]) == 2       # no retry on a refused credential
    second = collect(world, other, "decode", "0_0", decode_token="folder-token-two")
    assert second.model_dump(exclude={"utc"}) == first.model_dump(exclude={"utc"})
    mark = len(world.remote.log)
    assert collect(world, root, "decode", "0_0", decode_token="folder-token-two") == first   # recorded once, never refetched
    assert world.remote.log[mark:] == []


@pytest.mark.parametrize("fault,kind", [(("status", 429), "rate_limit"), (("status", 500), "server"), (("status", 503), "server"),
                                        (("status", 401), "auth"), (("status", 403), "auth")])
def test_a_refused_or_failing_source_leaves_no_record(world, tmp_path, fault, kind):
    root, path = tmp_path / "vol", f"/decode/{DECODE_KEY}"
    world.remote.faults[path] = [fault] * remote_module.ATTEMPTS
    with pytest.raises(RetryableSourceError) as err:
        collect(world, root, "decode", "0_0")
    assert err.value.kind == kind and not (root / "collect").exists()
    world.remote.faults[path] = []
    assert collect(world, root, "decode", "0_0").status == "collected"


def test_absence_is_recorded_for_a_404_a_file_the_listing_lacks_and_a_study_directory_without_summary_statistics(world, tmp_path):
    root = tmp_path / "vol"
    no_dir, no_file = collect(world, root, "gwas_catalog", NO_DIR), collect(world, root, "gwas_catalog", NO_FILE)
    assert (no_dir.status, no_dir.detail) == ("absent", f"GWAS Catalog {NO_DIR} listing: HTTP 404")
    assert (no_file.status, no_file.detail) == (
        "absent", f"GWAS Catalog {NO_FILE}: no harmonised file and no summary-statistics file in the study directory")
    del world.remote.files[f"/decode/{DECODE_KEY}"]
    gone = collect(world, root, "decode", "0_0")
    assert (gone.status, gone.detail) == ("absent", "decode 0_0: HTTP 404")
    unlisted = CollectTask(source="ukbppp", key="OID99999")
    record = collect_one(unlisted, world.sources(root / "cache"), root)
    assert record.status == "absent" and "lists 0 files" in record.detail
    with pytest.raises(CollectError, match="has no file on the volume"):
        collected_file(root, gone)


def test_a_truncated_gzip_is_corrupt_and_is_neither_kept_nor_recorded(world, tmp_path):
    root, path = tmp_path / "vol", f"/decode/{DECODE_KEY}"
    whole = world.remote.files[path]
    world.remote.files[path] = whole[: len(whole) // 2]               # the source itself serves half a stream
    unpinned = task(world, "decode", "0_0").model_copy(update={"size": None, "etag": ""})
    with pytest.raises(RetryableSourceError) as err:
        collect_one(unpinned, world.sources(root / "cache"), root)
    assert err.value.kind == "corrupt" and not (root / "collect").exists()
    assert not list((root / "raw" / "decode").iterdir())


def test_a_file_that_is_not_the_listed_one_is_refused(world, tmp_path):
    root, path = tmp_path / "vol", f"/decode/{DECODE_KEY}"
    world.remote.files[path] = gzip.compress(b"Chrom\tPos\n", mtime=0)   # another file under the listed name
    with pytest.raises(CollectError, match="pinned listing"):
        collect(world, root, "decode", "0_0")
    assert not (root / "collect").exists()


def test_a_download_that_cannot_be_resumed_starts_over(world, tmp_path):
    path, part = f"/decode/{DECODE_KEY}", tmp_path / "f.part"
    data = world.remote.files[path]
    url = f"{world.remote.base}/decode/s3/download?token={TOKENS['decode']}&file={DECODE_KEY}"

    def open_at(offset: int):
        return http("GET", url, "decode 0_0", headers={"Range": f"bytes={offset}-"} if offset else {}, stream=True)

    world.remote.faults[path] = [("break", len(data) // 3)]
    with pytest.raises(requests.RequestException):
        download(open_at, part, "decode 0_0", lambda: None)
    assert 0 < part.stat().st_size < len(data)
    world.remote.files[path] = gzip.compress(b"another version of the file\n" * 400, mtime=0)   # the source file changed
    found = download(open_at, part, "decode 0_0", lambda: None)
    assert part.read_bytes() == world.remote.files[path] and found["etag"] == world.remote.etag(path)
    assert download(open_at, part, "decode 0_0", lambda: None) == found          # complete: no request is made
    assert [rng for _m, p, rng in world.remote.log if p == path][-2:] != ["", ""]


def test_a_tar_is_resumed_through_a_fresh_link_and_checked_against_the_declared_md5(world, tmp_path):
    root, path = tmp_path / "vol", f"/synapse/{TAR_NAME}"
    synapse = FakeSynapse(world.remote, world.folders)
    world.remote.faults[path] = [("break", len(world.remote.files[path]) // 2)]
    record = collect(world, root, "ukbppp", OID, synapse=synapse)
    assert (root / record.path).read_bytes() == world.remote.files[path]
    assert (record.name, record.source_url, record.md5) == (TAR_NAME, "synapse:syn900001", synapse.file("syn900001")["md5"])
    assert "sig" not in record_path(root, "ukbppp", OID).read_text()
    assert len([p for _m, p, _r in world.remote.log if p == path]) == 2 and synapse.calls >= 3   # each connection asked for a link
    refused = FakeSynapse(world.remote, world.folders, refuse=RetryableSourceError("auth", "Synapse login: refused"))
    with pytest.raises(RetryableSourceError) as err:
        collect_one(task(world, "ukbppp_rsid_map", CHROM), world.sources(tmp_path / "cache2", synapse=refused), root)
    assert err.value.kind == "auth" and not record_path(root, "ukbppp_rsid_map", CHROM).exists()


# ---- analyze: whole files come from the volume, verified, or the unit stops ----------------------------

def decode_unit(world):
    return next(u for u in world.units if u.source == "decode")


def test_a_missing_collect_record_is_an_error_and_not_an_absence(world, tmp_path):
    fetcher = world.fetcher(tmp_path / "vol")
    with pytest.raises(CollectError, match="no collect record for decode 0_0"):
        fetcher.pqtl_region(decode_unit(world), CHROM, 50_000_000, WINDOW_WIDE)
    with pytest.raises(CollectError, match="no collect record for gwas_catalog"):
        fetcher.outcome_region(OutcomeSpec(accession=DISTINCT, source="gwas_catalog", n_case=1, n_control=1, risk_coded=True),
                               CHROM, 50_000_000, WINDOW_WIDE)


def test_whole_files_are_read_from_the_volume_with_no_request_and_only_after_they_match_their_record(world, tmp_path):
    root = tmp_path / "vol"
    for source, key in (("decode", "0_0"), ("decode_smp", "0_0"), ("ukbppp", OID), ("ukbppp_rsid_map", CHROM),
                        ("gwas_catalog", DISTINCT), ("gwas_catalog", NO_DIR)):
        collect(world, root, source, key)
    world.remote.tokens.update({"decode": "folder-token-two", "decode_smp": "smp-two"})       # links re-issued after collect
    world.remote.signatures.clear()
    fetcher, mark = world.fetcher(root), len(world.remote.log)
    decode = fetcher.pqtl_region(decode_unit(world), CHROM, 50_000_000, WINDOW_WIDE)
    smp = fetcher.pqtl_region(decode_unit(world), CHROM, 50_000_000, WINDOW_WIDE, smp=True)
    ukb = fetcher.pqtl_region(next(u for u in world.units if u.source == "ukbppp"), CHROM, 50_000_000, WINDOW_WIDE)
    spec = OutcomeSpec(accession=DISTINCT, source="gwas_catalog", n_case=1, n_control=1, risk_coded=True)
    outcome = fetcher.outcome_region(spec, CHROM, 50_000_000, WINDOW_WIDE)
    assert world.remote.log[mark:] == []
    assert len(decode) == len(smp) == N_VARIANTS - 1 and rsid(EXCLUDED_INDEX) not in set(decode["rsid"])
    assert len(ukb) == len(outcome) == N_VARIANTS and list(ukb["rsid"]) == [rsid(i) for i in range(N_VARIANTS)]
    lead = decode.loc[decode["p"].idxmin()]
    assert (lead["rsid"], lead["ea"], lead["oa"]) == (rsid(SENTINEL_INDEX), "G", "A") and lead["eaf"] == pytest.approx(0.3, abs=0.01)
    assert smp["beta"].to_numpy() == pytest.approx(0.9 * decode["beta"].to_numpy())
    assert ukb["beta"].to_numpy() == pytest.approx(decode.set_index("rsid")["beta"].reindex(ukb["rsid"]).fillna(ukb.set_index("rsid")["beta"]).to_numpy())
    with pytest.raises(SourceAbsent, match="HTTP 404"):
        fetcher.outcome_region(spec.model_copy(update={"accession": NO_DIR}), CHROM, 50_000_000, WINDOW_WIDE)
    raw = root / "raw" / "decode" / DECODE_KEY
    data = bytearray(raw.read_bytes())
    data[len(data) // 2] ^= 0x01                                         # one changed bit, the same size
    raw.write_bytes(bytes(data))
    with pytest.raises(CollectError, match="sha256"):
        fetcher.pqtl_region(decode_unit(world), CHROM, 50_000_000, WINDOW_WIDE)


def test_regional_queries_raise_on_an_expired_token_and_record_absence_only_for_404(world, tmp_path):
    fetcher = world.fetcher(tmp_path / "vol")
    spec = OutcomeSpec(accession=GRCH37, source="opengwas", n_case=1, n_control=1, risk_coded=True)
    assert len(fetcher.outcome_region(spec, CHROM, 50_000_000, WINDOW_WIDE)) == N_VARIANTS
    with pytest.raises(RetryableSourceError) as err:
        world.fetcher(tmp_path / "vol", opengwas_token="opengwas-token-expired").outcome_region(spec, CHROM, 50_000_000, WINDOW_WIDE)
    assert err.value.kind == "auth" and len([p for _m, p, _r in world.remote.log if p == "/opengwas/associations"]) == 2
    world.remote.faults["/opengwas/associations"] = [("status", 429)] * remote_module.ATTEMPTS
    with pytest.raises(RetryableSourceError) as err:
        fetcher.outcome_region(spec, CHROM, 50_000_000, WINDOW_WIDE)
    assert err.value.kind == "rate_limit"
    world.remote.faults["/opengwas/associations"] = [("status", 404)]
    with pytest.raises(SourceAbsent):
        fetcher.outcome_region(spec, CHROM, 50_000_000, WINDOW_WIDE)
    with pytest.raises(SourceAbsent, match="HTTP 404"):                 # a FinnGen endpoint that does not exist
        fetcher.outcome_region(OutcomeSpec(accession="FINNGEN_R12_NO_SUCH", source="finngen", n_case=1, n_control=1,
                                           risk_coded=True), CHROM, 50_000_000, WINDOW_WIDE)
    world.remote.faults["/finngen/finngen_R12_SYN_SHARED.gz"] = [("status", 503)] * remote_module.ATTEMPTS
    with pytest.raises(RetryableSourceError) as err:
        fetcher.outcome_region(OutcomeSpec(accession=SHARED, source="finngen", n_case=1, n_control=1, risk_coded=True),
                               CHROM, 50_000_000, WINDOW_WIDE)
    assert err.value.kind == "server"


# ---- Ensembl: a 400 is an absence only as its unknown-identifier message naming the requested rsID ------

def test_ensembl_knows_the_sentinel_and_reports_an_rsid_it_does_not_hold_as_unknown(world, tmp_path):
    fetcher, sentinel = world.fetcher(tmp_path / "vol"), decode_unit(world).sentinel
    assert fetcher.positions(sentinel) == {"GRCh38": sentinel.pos, "GRCh37": sentinel.pos - 12_345}
    unknown = sentinel.model_copy(update={"rsid": "rs1"})        # the server answers 400 {"error": "rs1 not found for human"}
    assert fetcher.positions(unknown) == {"GRCh38": sentinel.pos, "GRCh37": None}
    assert fetcher.vep([sentinel.rsid], "GRCh38")[0]["id"] == sentinel.rsid


def test_ensembl_error_replies_are_an_absence_only_for_the_unknown_id_message_naming_the_requested_rsid(world, tmp_path):
    sentinel = decode_unit(world).sentinel
    got = ensembl_outcomes(world, world.fetcher(tmp_path / "vol"), sentinel)
    assert set(ENSEMBL_REPLIES) >= {"unknown_id", "backend_error", "endpoint_not_found", "malformed_json", "another_rsid"}
    assert got == {"position_lookup": ENSEMBL_EXPECTED, "vep": ENSEMBL_EXPECTED}
    assert [name for name, outcome in ENSEMBL_EXPECTED.items() if outcome == "absent"] == ["unknown_id"]
    assert set(ENSEMBL_EXPECTED.values()) == {"absent", "raised:protocol", "raised:server"}
    assert fetcher_still_answers(world, tmp_path, sentinel)          # the routes were put back


def fetcher_still_answers(world, tmp_path, sentinel) -> bool:
    fetcher = world.fetcher(tmp_path / "vol")
    return fetcher.positions(sentinel)["GRCh37"] == sentinel.pos - 12_345 and bool(fetcher.vep([sentinel.rsid], "GRCh38"))


GENUINE = {ENSEMBL_UNKNOWN_VARIATION: "rs123 not found for human", ENSEMBL_UNKNOWN_VEP_ID: "No variant found with ID 'rs123'"}


@pytest.mark.parametrize("pattern", list(GENUINE), ids=["position_lookup", "vep"])
def test_the_ensembl_matcher_needs_the_error_object_the_whole_message_and_a_requested_rsid(pattern):
    message = GENUINE[pattern]

    def unknown(body, requested=("rs123",)) -> bool:
        return ensembl_unknown_id(body if isinstance(body, bytes) else json.dumps(body).encode(), pattern, set(requested))

    assert unknown({"error": message})
    assert unknown({"error": message}, requested=("rs5", "rs123"))                  # one of the ids of a batch
    assert not unknown({"error": message}, requested=("rs1234",))                   # another rsID, of which rs123 is a prefix
    assert not unknown({"error": message}, requested=())
    for body in ({"error": "Unknown backend error while querying the variation database"}, {"error": "endpoint not found"},
                 {"error": "page not found. Please check your uri and refer to our documentation https://rest.ensembl.org/"},
                 json.dumps({"error": message}).encode()[:-1], b"", b"<html>rs123 not found for human</html>",
                 {"error": message.replace("rs123", "rs124")}, {"error": message + " "}, {"error": " " + message},
                 {"error": message + "; Unknown backend error"}, {"error": message.lower() if message[0] == "N" else message.upper()},
                 {"error": message, "status": 500}, {"message": message}, {"error": [message]}, {"error": None}, [message],
                 message, {"error": "ID 'rs123' not found"}, {"error": "rs123 unknown"}, {"error": "unknown"},
                 {"error": "rs123 not found"}, {"error": "No variant found"}, {"error": "No mappings found for variant 'rs123'"},
                 {"error": "POST message too large. You have submitted 201 elements but a limit of 200 is in place. "
                           "Request smaller regions or lists of IDs"},
                 {"error": ' Cannot find "ids" key in your POST. Please check the format of your message against the documentation'}):
        assert not unknown(body), body
    other = next(m for p, m in GENUINE.items() if p is not pattern)
    assert not unknown({"error": other})                                            # the other endpoint's message
    assert not unknown({"error": "rs123 not found for homo_sapiens"})               # the species is the URL's: `human`


def test_an_ensembl_status_other_than_the_unknown_id_400_raises_with_its_class(world, tmp_path):
    fetcher, remote = world.fetcher(tmp_path / "vol"), world.remote
    sentinel = decode_unit(world).sentinel
    path = f"/ensembl/GRCh37/variation/human/{sentinel.rsid}"
    for status, kind in ((404, "protocol"), (410, "protocol"), (400, "protocol"), (429, "rate_limit"), (500, "server"), (403, "auth")):
        remote.faults[path] = [("status", status)] * remote_module.ATTEMPTS
        with pytest.raises(RetryableSourceError) as err:
            fetcher.positions(sentinel)
        assert err.value.kind == kind and str(err.value) == f"[{kind}] Ensembl GRCh37 {sentinel.rsid}: HTTP {status}"
        remote.faults[path] = []
    vep = ("POST", "/ensembl/GRCh38/vep/human/id")
    remote.json_routes[vep] = lambda q, b, h: (503, {"error": f"No variant found with ID '{sentinel.rsid}'"})   # the words alone are not a 400
    with pytest.raises(RetryableSourceError) as err:
        fetcher.vep([sentinel.rsid], "GRCh38")
    assert err.value.kind == "server"
    remote.json_routes[vep] = lambda q, b, h: (400, {"error": "x" * 5000})
    with pytest.raises(RetryableSourceError) as err:
        fetcher.vep([sentinel.rsid], "GRCh38")
    assert err.value.kind == "protocol"


def test_vep_batches_ensembl_does_not_know_add_nothing_and_only_all_of_them_unknown_is_an_absence(world, tmp_path):
    fetcher, remote = world.fetcher(tmp_path / "vol"), world.remote
    known = [rsid(i) for i in range(N_VARIANTS)]
    asked: list[list[str]] = []

    def vep(q, b, h):
        """As Ensembl: records for the ids it knows; 400 naming the first id only when it knows none."""
        ids = json.loads(b)["ids"]
        asked.append(ids)
        found = [{"id": i, "transcript_consequences": []} for i in ids if i in known]
        return (200, found) if found else (400, {"error": f"No variant found with ID '{ids[0]}'"})

    remote.json_routes[("POST", "/ensembl/GRCh38/vep/human/id")] = vep
    strangers = [f"rs{i}" for i in range(1, 251)]                    # none is in the synthetic panel
    got = fetcher.vep(known[:200] + strangers[:200] + known[200:] + ["not_an_rsid"], "GRCh38")
    assert [len(batch) for batch in asked] == [200, 200, N_VARIANTS - 200]      # the second batch is answered 400
    assert [r["id"] for r in got] == known                                       # and the batches around it are kept
    assert [r["id"] for r in fetcher.vep(known[:3] + strangers[:3], "GRCh38")] == known[:3]
    with pytest.raises(SourceAbsent, match="VEP GRCh38: HTTP 400, the id is not known to Ensembl"):
        fetcher.vep(strangers, "GRCh38")                                         # two batches, both unknown
    assert fetcher.vep(["not_an_rsid"], "GRCh38") == []
    remote.json_routes[("POST", "/ensembl/GRCh38/vep/human/id")] = lambda q, b, h: (400, {"error": "No variant found with ID 'rs999'"})
    with pytest.raises(RetryableSourceError) as err:                              # the reply names an id that was not asked
        fetcher.vep(strangers[:5], "GRCh38")
    assert err.value.kind == "protocol"


# ---- a source that keeps failing is never turned into "unavailable" ----------------------------------------

CALLS_FOR_FIFTY = -(-50 // remote_module.ATTEMPTS)          # calls whose attempts add up to at least 50 responses


def test_fifty_consecutive_server_errors_on_a_file_leave_no_record_no_file_and_nothing_absent(world, tmp_path):
    root, path = tmp_path / "vol", f"/decode/{DECODE_KEY}"
    world.remote.faults[path] = [("status", 503)] * (CALLS_FOR_FIFTY * remote_module.ATTEMPTS)
    for _ in range(CALLS_FOR_FIFTY):                             # as Modal's retries and later launches: the same call again
        with pytest.raises(RetryableSourceError) as err:
            marked(root, "collect", "decode__0_0", "decode/0_0", lambda: collect(world, root, "decode", "0_0"), lambda: None)
        assert err.value.kind == "server"
        assert not (root / "collect").exists() and not (root / "raw").exists()
    assert world.remote.faults[path] == []
    assert len([p for _m, p, _r in world.remote.log if p == path]) == CALLS_FOR_FIFTY * remote_module.ATTEMPTS >= 50
    state = volume_state(root)
    assert state["records"] == [] and state["errors"]["collect"]["decode/0_0"]["kind"] == "server"
    failed = {"call_id": "fc-1", "state": "failed", **error_of(err.value)}
    report = status_report(collect_tasks(world.units), world.units, state, {"collect": {"decode/0_0": failed}}, {"deployed": True})
    assert report["files"]["decode/0_0"]["state"] == "failed" and report["files_by_source"]["decode"] == {"failed": 1}
    assert collect(world, root, "decode", "0_0").status == "collected"   # the source held the file all along


def test_fifty_consecutive_server_errors_in_a_unit_leave_no_step_no_result_and_assembly_refuses(world, tmp_path):
    root = tmp_path / "vol"
    unit = decode_unit(world)
    path = f"/ensembl/GRCh37/variation/human/{unit.sentinel.rsid}"
    world.remote.faults[path] = [("status", 502)] * (CALLS_FOR_FIFTY * remote_module.ATTEMPTS)
    store, fetcher = DirStore(root / "units" / unit.unit_key), world.fetcher(root)
    for _ in range(CALLS_FOR_FIFTY):
        with pytest.raises(RetryableSourceError) as err:
            marked(root, "units", unit.unit_key, unit.unit_key,
                   lambda: process_unit(unit, fetcher, StubBackend(H4), store, fp(unit), TOOLS, COLLECT), lambda: None)
        assert err.value.kind == "server"
        assert sorted(p.name for p in store.root.iterdir()) == [FINGERPRINT_NAME]    # no step, no meta file, no result
    assert len([p for _m, p, _r in world.remote.log if p == path]) == CALLS_FOR_FIFTY * remote_module.ATTEMPTS >= 50
    state = volume_state(root)
    assert state["units_done"] == [] and state["errors"]["units"][unit.unit_key]["kind"] == "server"
    report = status_report([], world.units, state, {}, {"deployed": True})
    assert report["units"][unit.unit_key]["state"] == "not_started" and "done" not in report["units_by_source"]["decode"]
    with pytest.raises(InputContractError, match="has no result.json; stage B is not complete"):
        collect_unit_dir(unit, store.root, fp(unit), COLLECT)


def test_only_a_definitive_absence_is_caught_and_recorded_and_modal_retries_only_run_the_call_again():
    def caught(path: Path) -> list[str]:
        return sorted(ast.unparse(n.type) for n in ast.walk(ast.parse(path.read_text())) if isinstance(n, ast.ExceptHandler))

    package = HERE / "stage_b"
    # every handler of the unit pipeline: three record a definitive absence, one notes that susie was not run
    assert caught(package / "pipeline.py") == ["LDReferenceError", "SourceAbsent", "SourceAbsent", "SourceAbsent"]
    # collect: two record a definitive absence; a 416 restarts the download; a corrupt gzip is reported as corrupt
    assert caught(package / "collect.py") == ["CORRUPT", "RetryableSourceError", "SourceAbsent", "SourceAbsent"]
    wrapper = ast.parse((HERE / "modal_stage_b.py").read_text())
    for name in ("collect_file", "run_unit", "checkpointed_unit"):
        fn = next(n for n in ast.walk(wrapper) if isinstance(n, ast.FunctionDef) and n.name == name)
        assert not [n for n in ast.walk(fn) if isinstance(n, (ast.Try, ast.While, ast.For))], name   # nothing counts or absorbs failures
    marker = next(n for n in ast.walk(ast.parse((package / "status.py").read_text()))
                  if isinstance(n, ast.FunctionDef) and n.name == "marked")
    handler = next(n for n in ast.walk(marker) if isinstance(n, ast.ExceptHandler))
    assert isinstance(handler.body[-1], ast.Raise) and handler.body[-1].exc is None      # the note is written, the error re-raised


# ---- other bytes under the same name, size and ETag ------------------------------------------------------

def test_a_file_served_with_one_other_byte_under_the_same_name_size_and_etag_makes_the_old_checkpoint_refused(world, tmp_path):
    root, path = tmp_path / "vol", f"/decode/{DECODE_KEY}"
    unit = decode_unit(world).model_copy(update={"sentinel": SENTINEL, "outcomes": (
        OutcomeSpec(accession="F_ok", source="finngen", n_case=1000, n_control=9000, risk_coded=True),)})
    assert [(t.source, t.key) for t in collect_tasks([unit])] == [("decode", "0_0"), ("decode_smp", "0_0")]
    for t in collect_tasks([unit]):
        collect_one(t, world.sources(root / "cache"), root, checkpoint_bytes=2048)
    first, old = read_record(root, "decode", "0_0"), volume_collect_digest(root, unit)
    store = DirStore(root / "units" / unit.unit_key)
    done = process_unit(unit, FakeFetcher(), StubBackend(H4), store, fp(unit, collect=old), TOOLS, old)
    assert done["collect_sha256"] == old and done["outcomes"]["F_ok"]["coloc_run"] is True

    data = world.remote.files[path]
    world.remote.etags[path] = world.remote.etag(path)                      # the ETag the listing pins is still served
    world.remote.files[path] = same_size_other_bytes(data)
    assert len(world.remote.files[path]) == len(data) and sum(a != b for a, b in zip(data, world.remote.files[path])) == 1
    record_path(root, "decode", "0_0").unlink()                              # the file is collected again from the source
    (root / first.path).unlink()
    second = collect(world, root, "decode", "0_0")
    assert (second.name, second.bytes, second.etag) == (first.name, first.bytes, first.etag) == (
        unit.pqtl_listing.name, unit.pqtl_listing.size, unit.pqtl_listing.etag)
    assert {k for k, v in first.model_dump().items() if getattr(second, k) != v} - {"utc"} == {"sha256"}
    new = volume_collect_digest(root, unit)
    assert new != old and fp(unit, collect=new) != fp(unit, collect=old)

    before, fetcher = {p.name: p.read_bytes() for p in sorted(store.root.iterdir())}, FakeFetcher()
    with pytest.raises(StaleCheckpointError, match="not resumed"):
        process_unit(unit, fetcher, StubBackend(H4), store, fp(unit, collect=new), TOOLS, new)
    assert fetcher.calls == {} and {p.name: p.read_bytes() for p in sorted(store.root.iterdir())} == before
    with pytest.raises(StaleCheckpointError):
        collect_unit_dir(unit, store.root, fp(unit, collect=new), new)       # assembly recomputes the digest from the records
    fresh = process_unit(unit, FakeFetcher(), StubBackend(H4), DirStore(root / "units_again" / unit.unit_key),
                         fp(unit, collect=new), TOOLS, new)
    assert fresh["collect_sha256"] == new and fresh["fingerprint"] != done["fingerprint"]


# ---- the deCODE download endpoint: a redirect, or a JSON body carrying the address -------------------

def test_the_default_decode_address_puts_the_token_and_the_encoded_key_in_the_query(tmp_path):
    assert DECODE_FILE_URL == "https://download.decode.is/s3/download?token={token}&file={key}"
    listed = CollectTask(source="decode", key="1_1", name="1_1_GENE_Protein name+x.txt.gz")
    resolved = RemoteSources(tmp_path / "cache", "to/ken", lambda: None).resolve(listed)      # the registered addresses
    assert resolved.url() == "https://download.decode.is/s3/download?token=to%2Fken&file=1_1_GENE_Protein%20name%2Bx.txt.gz"
    assert resolved.source_url == "https://download.decode.is/s3/download?token=<token>&file=1_1_GENE_Protein%20name%2Bx.txt.gz"


def test_a_json_reply_carrying_the_address_is_followed_and_resumed_like_a_redirect(world, tmp_path):
    root, path = tmp_path / "vol", f"/decode_smp/{DECODE_KEY}"
    assert "decode_smp" in world.remote.json_links and "decode" not in world.remote.json_links
    world.remote.faults[path] = [("break", len(world.remote.files[path]) // 3)]
    record = collect(world, root, "decode_smp", "0_0")
    assert (root / record.path).read_bytes() == world.remote.files[path]
    fetched = [(p, bool(rng)) for _m, p, rng in world.remote.log if p == path]
    assert fetched == [(path, False), (path, True)]                       # resumed by range through a second link
    text = record_path(root, "decode_smp", "0_0").read_text()
    assert "signed" not in text and TOKENS["decode_smp"] not in text
    assert record.source_url == f"{world.remote.base}/decode_smp/s3/download?token=<token>&file={DECODE_KEY}"
    world.remote.json_routes[("GET", "/nolink")] = lambda q, b, h: (200, {"message": "try later"})
    broken = RemoteFile(name="f.gz", source_url="x", what="decode 0_0", url=lambda: f"{world.remote.base}/nolink", link_reply=True)
    with pytest.raises(RetryableSourceError) as err:
        broken.open_at(0)
    assert err.value.kind == "protocol"


def test_eqtl_catalogue_paths_are_read_over_https_on_the_same_host():
    assert https_path("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000266/QTD000266.all.tsv.gz") == (
        "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000266/QTD000266.all.tsv.gz")
    assert https_path("http://127.0.0.1:1/eqtlcat/x.tsv.gz") == "http://127.0.0.1:1/eqtlcat/x.tsv.gz"


# ---- purge: the whole files are working copies -------------------------------------------------------

def test_purge_deletes_the_whole_files_only_when_every_unit_is_finished_and_keeps_the_records(world, tmp_path):
    root = tmp_path / "vol"
    tasks = collect_tasks(world.units)
    records = {(t.source, t.key): collect_one(t, world.sources(root / "cache"), root) for t in tasks if t.key != "GCST90000002"}
    for u in world.units[:-1]:
        (root / "units" / u.unit_key).mkdir(parents=True)
        (root / "units" / u.unit_key / "result.json").write_text("{}")
        (root / "units" / u.unit_key / "pqtl.tsv.gz").write_bytes(b"extract")
    raw = sorted(p for p in (root / "raw").rglob("*") if p.is_file())
    held = sum(p.stat().st_size for p in raw)
    assert len(raw) == 5
    with pytest.raises(CollectError, match="1 units have no result.json"):
        purge_raw(root, world.units)
    last = world.units[-1]
    (root / "units" / last.unit_key).mkdir()
    (root / "units" / last.unit_key / "result.json").write_text("{}")
    with pytest.raises(CollectError, match="no collect record for gwas_catalog GCST90000002"):
        purge_raw(root, world.units)
    indexed = CollectTask(source="gwas_catalog", key="GCST90000002")
    record_path(root, "gwas_catalog", "GCST90000002").write_text(json.dumps(
        {"status": "remote_indexed", "source": indexed.source, "key": indexed.key, "name": "x.h.tsv.gz"}))
    stray = root / "raw" / "decode" / "not_recorded.txt.gz"
    stray.write_bytes(b"x")
    with pytest.raises(CollectError, match="is the file of no collect record"):
        purge_raw(root, world.units)
    stray.unlink()
    assert all(p.is_file() for p in raw)                                 # every refusal left every file in place
    (root / "raw" / "decode" / "leftover.txt.gz.part").write_bytes(b"partial")
    out = purge_raw(root, world.units)
    assert [p for p in (root / "raw").rglob("*") if p.is_file()] == []
    assert len(out["deleted"]) == 6 and out["bytes_deleted"] == held + len(b"partial")
    assert {d["sha256"] for d in out["deleted"] if d["sha256"]} == {r.sha256 for r in records.values() if r.status == "collected"}
    assert {k: read_record(root, *k) for k in records} == records        # size, sha256, ETag and time are still on record
    assert (root / "units" / world.units[0].unit_key / "pqtl.tsv.gz").read_bytes() == b"extract"
    assert json.loads((root / PURGED_NAME).read_text())["deleted"] == out["deleted"]
    assert purge_raw(root, world.units)["deleted"] == out["deleted"]     # a second purge deletes nothing more
    with pytest.raises(CollectError, match="is missing although"):
        collected_file(root, records[("decode", "0_0")])


# ---- status ----------------------------------------------------------------------------------------

def test_status_tells_finished_failed_pending_and_stranded_apart(world, tmp_path):
    root = tmp_path / "vol"
    tasks = collect_tasks(world.units)
    collect(world, root, "decode", "0_0")
    collect(world, root, "gwas_catalog", NO_DIR)
    world.remote.tokens["decode_smp"] = "smp-two"
    with pytest.raises(RetryableSourceError):
        marked(root, "collect", "decode_smp__0_0", "decode_smp/0_0", lambda: collect(world, root, "decode_smp", "0_0"), lambda: None)
    decode_key, ukb_key = (next(u.unit_key for u in world.units if u.source == s) for s in ("decode", "ukbppp"))
    (root / "units" / ukb_key).mkdir(parents=True)
    (root / "units" / ukb_key / "result.json").write_text("{}")
    (root / "units" / decode_key).mkdir()
    state = volume_state(root)
    assert state["errors"]["collect"]["decode_smp/0_0"]["kind"] == "auth" and state["units_done"] == [ukb_key]
    failed = {"call_id": "fc-1", "state": "failed", **error_of(RetryableSourceError("auth", "decode_smp 0_0: HTTP 403"))}
    calls = {"collect": {"decode/0_0": {"call_id": "fc-0", "state": "done"}, "decode_smp/0_0": failed,
                         f"ukbppp/{OID}": {"call_id": "fc-2", "state": "pending"}},
             "units": {decode_key: {"call_id": "fc-3", "state": "pending"}}}
    live = status_report(tasks, world.units, state, calls, {"deployed": True})
    assert live["files"]["decode/0_0"] == {"state": "done", "record": "collected"}
    assert live["files"]["decode_smp/0_0"] == {"state": "failed", "call_id": "fc-1", "error_class": "RetryableSourceError",
                                               "kind": "auth", "message": "[auth] decode_smp 0_0: HTTP 403"}
    assert live["files"][f"ukbppp/{OID}"]["state"] == "pending" and live["files"][f"gwas_catalog/{DISTINCT}"]["state"] == "not_started"
    assert live["files_by_source"]["gwas_catalog"] == {"absent": 1, "not_started": 3}
    assert live["units_by_source"] == {"decode": {"pending": 1}, "interval": {"not_started": 1}, "ukbppp": {"done": 1}}
    assert live["failures_by_error_class"] == {"RetryableSourceError:auth": 1}
    dead = status_report(tasks, world.units, state, calls, {"deployed": False})   # the same volume, the app gone
    assert dead["files"][f"ukbppp/{OID}"]["state"] == "stranded" and dead["units"][decode_key]["state"] == "stranded"
    assert dead["failures_by_error_class"] == {"RetryableSourceError:auth": 1, "stranded": 2}
    assert json.dumps(dead)                                              # the report is plain JSON
    marked(root, "collect", "decode_smp__0_0", "decode_smp/0_0", lambda: None, lambda: None)   # a later success clears the note
    assert volume_state(root)["errors"]["collect"] == {}
