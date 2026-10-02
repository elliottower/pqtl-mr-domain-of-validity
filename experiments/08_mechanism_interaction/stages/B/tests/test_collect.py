"""The collect phase and the error classification against a local HTTP server (fake_remote.py):
resumable downloads, one record per file, absence against every other fault, and the analyze-side
reads of collected files. The server is the synthetic one of synthetic_unit.py without its
tabix-indexed files, so nothing here needs htslib."""
import gzip
import json
import subprocess
import zlib

import pytest
import requests
from fake_remote import FakeRemote, FakeSynapse
from synthetic_unit import (CHROM, DECODE_KEY, DISTINCT, EXCLUDED_INDEX, GRCH37, N_VARIANTS, NO_DIR, NO_FILE, OID, SENTINEL_INDEX,
                            SHARED, TAR_NAME, TOKENS, build_world, rsid)

from stage_b import collect as collect_module
from stage_b import remote as remote_module
from stage_b.collect import PURGED_NAME, collect_one, collect_tasks, collected_file, download, purge_raw, read_record, record_path
from stage_b.fetch import DECODE_FILE_URL, RemoteFile, RemoteSources, https_path
from stage_b.remote import attempt, check_status, classify, http
from stage_b.schemas import (WINDOW_WIDE, CollectError, CollectTask, OutcomeSpec, RetryableSourceError, SourceAbsent)
from stage_b.status import error_of, marked, status_report, volume_state


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


def test_absence_is_recorded_for_a_404_a_file_the_listing_lacks_and_an_accession_without_a_harmonised_file(world, tmp_path):
    root = tmp_path / "vol"
    no_dir, no_file = collect(world, root, "gwas_catalog", NO_DIR), collect(world, root, "gwas_catalog", NO_FILE)
    assert (no_dir.status, no_dir.detail) == ("absent", f"GWAS Catalog {NO_DIR} listing: HTTP 404")
    assert (no_file.status, no_file.detail) == ("absent", f"GWAS Catalog {NO_FILE}: 0 harmonised files")
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


def test_ensembl_400_is_an_absence_only_when_the_body_says_the_id_is_not_known(world, tmp_path):
    fetcher, remote = world.fetcher(tmp_path / "vol"), world.remote
    sentinel = decode_unit(world).sentinel
    assert fetcher.positions(sentinel) == {"GRCh38": sentinel.pos, "GRCh37": sentinel.pos - 12_345}
    unknown = sentinel.model_copy(update={"rsid": "rs1"})               # the server answers 400 "ID 'rs1' not found"
    assert fetcher.positions(unknown) == {"GRCh38": sentinel.pos, "GRCh37": None}
    remote.json_prefixes.insert(0, ("GET", "/ensembl/GRCh37/variation/human/", lambda rest, q, b, h: (400, {"error": "page size too large"})))
    with pytest.raises(RetryableSourceError) as err:
        fetcher.positions(sentinel)
    assert err.value.kind == "protocol" and str(err.value) == f"[protocol] Ensembl GRCh37 {sentinel.rsid}: HTTP 400"
    vep = ("POST", "/ensembl/GRCh38/vep/human/id")
    assert fetcher.vep([sentinel.rsid], "GRCh38")[0]["id"] == sentinel.rsid
    remote.json_routes[vep] = lambda q, b, h: (400, {"error": "No variant found with ID 'rs1'"})
    with pytest.raises(SourceAbsent, match="the id is not known to Ensembl"):
        fetcher.vep(["rs1"], "GRCh38")
    remote.json_routes[vep] = lambda q, b, h: (400, {"error": "Bad request: POST message exceeds the limit"})
    with pytest.raises(RetryableSourceError) as err:
        fetcher.vep([sentinel.rsid], "GRCh38")
    assert err.value.kind == "protocol"
    remote.json_routes[vep] = lambda q, b, h: (503, {"error": "not found"})          # the words alone are not a 400
    with pytest.raises(RetryableSourceError) as err:
        fetcher.vep([sentinel.rsid], "GRCh38")
    assert err.value.kind == "server"


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
    (root / "units" / "ukbppp__OID00000").mkdir(parents=True)
    (root / "units" / "ukbppp__OID00000" / "result.json").write_text("{}")
    (root / "units" / "decode__0_0").mkdir()
    state = volume_state(root)
    assert state["errors"]["collect"]["decode_smp/0_0"]["kind"] == "auth" and state["units_done"] == ["ukbppp__OID00000"]
    failed = {"call_id": "fc-1", "state": "failed", **error_of(RetryableSourceError("auth", "decode_smp 0_0: HTTP 403"))}
    calls = {"collect": {"decode/0_0": {"call_id": "fc-0", "state": "done"}, "decode_smp/0_0": failed,
                         f"ukbppp/{OID}": {"call_id": "fc-2", "state": "pending"}},
             "units": {"decode__0_0": {"call_id": "fc-3", "state": "pending"}}}
    live = status_report(tasks, world.units, state, calls, {"deployed": True})
    assert live["files"]["decode/0_0"] == {"state": "done", "record": "collected"}
    assert live["files"]["decode_smp/0_0"] == {"state": "failed", "call_id": "fc-1", "error_class": "RetryableSourceError",
                                               "kind": "auth", "message": "[auth] decode_smp 0_0: HTTP 403"}
    assert live["files"][f"ukbppp/{OID}"]["state"] == "pending" and live["files"][f"gwas_catalog/{DISTINCT}"]["state"] == "not_started"
    assert live["files_by_source"]["gwas_catalog"] == {"absent": 1, "not_started": 3}
    assert live["units_by_source"] == {"decode": {"pending": 1}, "interval": {"not_started": 1}, "ukbppp": {"done": 1}}
    assert live["failures_by_error_class"] == {"RetryableSourceError:auth": 1}
    dead = status_report(tasks, world.units, state, calls, {"deployed": False})   # the same volume, the app gone
    assert dead["files"][f"ukbppp/{OID}"]["state"] == "stranded" and dead["units"]["decode__0_0"]["state"] == "stranded"
    assert dead["failures_by_error_class"] == {"RetryableSourceError:auth": 1, "stranded": 2}
    assert json.dumps(dead)                                              # the report is plain JSON
    marked(root, "collect", "decode_smp__0_0", "decode_smp/0_0", lambda: None, lambda: None)   # a later success clears the note
    assert volume_state(root)["errors"]["collect"] == {}
