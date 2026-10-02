"""Error classification for every network access of stage B (collect.py, fetch.py).

The frozen plan gives a consequence only for a source that fails: a regional or outcome file that
is unavailable makes the hypothesis inconclusive. A link that expired, a token that was not
renewed or a dropped connection is not the source failing, so two classes are kept apart:

    SourceAbsent          HTTP 404 or 410. The callers add the other definitive cases: a file the
                          source's listing does not name and an accession with no harmonised
                          file. Ensembl has its own rule (fetch.ensembl_json): the one absence is
                          its HTTP 400 whose JSON error is the unknown-identifier message naming a
                          requested rsID; its 404 is a malformed request and raises as `protocol`.
                          Recorded as unavailable.
    RetryableSourceError  everything else: 401/403 (kind `auth`), 429 (`rate_limit`), 5xx
                          (`server`), any other status (`protocol`), a timeout (`timeout`), a
                          connection error or a failed tabix/bcftools call (`connection`), a
                          truncated or corrupt gzip or tar stream (`corrupt`). Raised, never
                          recorded: the file or unit stays unfinished.

`attempt` retries a retryable fault inside the call with backoff and then raises it; an `auth`
fault is raised at once, because waiting does not renew a credential. A source that keeps failing
is never converted to SourceAbsent, here or by any caller: after any number of retryable faults the
file or unit is unfinished, `status` shows it, and stage B is not assembled. Messages name the source by
a label the caller passes and the exception by its type only: a URL can carry a token, so no URL
and no exception text reaches a message, a log or a file.
"""
import gzip
import http.client
import subprocess
import tarfile
import time
import zlib
from collections.abc import Callable
from typing import TypeVar

import requests
import urllib3

from stage_b.schemas import RetryableSourceError, SourceAbsent

T = TypeVar("T")

ABSENT_STATUS = (404, 410)
ATTEMPTS = 4
WAIT_S = 5.0        # first wait; doubled after each failed attempt
TIMEOUTS = (requests.Timeout, urllib3.exceptions.TimeoutError, TimeoutError, subprocess.TimeoutExpired)
CORRUPT = (EOFError, zlib.error, gzip.BadGzipFile, tarfile.ReadError)
TRANSPORT = (requests.RequestException, urllib3.exceptions.HTTPError, http.client.HTTPException, OSError,
             subprocess.SubprocessError)


def status_kind(code: int) -> str:
    if code in (401, 403):
        return "auth"
    if code == 429:
        return "rate_limit"
    return "server" if 500 <= code < 600 else "protocol"


def check_status(code: int, what: str) -> None:
    """Nothing for 2xx; SourceAbsent for 404/410; RetryableSourceError for any other status."""
    if code in ABSENT_STATUS:
        raise SourceAbsent(f"{what}: HTTP {code}")
    if not 200 <= code < 300:
        raise RetryableSourceError(status_kind(code), f"{what}: HTTP {code}")


def classify(err: BaseException, what: str) -> SourceAbsent | RetryableSourceError:
    """The stage B error for an exception a network or decompression call raised."""
    if isinstance(err, (SourceAbsent, RetryableSourceError)):
        return err
    if isinstance(err, CORRUPT):
        return RetryableSourceError("corrupt", f"{what}: {type(err).__name__}")
    if isinstance(err, TIMEOUTS):
        return RetryableSourceError("timeout", f"{what}: {type(err).__name__}")
    return RetryableSourceError("connection", f"{what}: {type(err).__name__}")


def attempt(call: Callable[[], T], what: str, attempts: int | None = None) -> T:
    """`call()`, with a retryable fault retried up to `attempts` times (default ATTEMPTS), waiting
    WAIT_S, 2 WAIT_S, ... between tries. SourceAbsent and an `auth` fault are raised at once; the
    last retryable fault is raised when the tries are spent. Any other exception passes through."""
    last: RetryableSourceError | None = None
    n = ATTEMPTS if attempts is None else attempts
    for i in range(n):
        try:
            return call()
        except (*CORRUPT, *TRANSPORT, SourceAbsent, RetryableSourceError) as err:
            classified = classify(err, what)
            if isinstance(classified, SourceAbsent):
                raise classified from None
            if classified.kind == "auth":
                raise classified from None
            last = classified
        if i + 1 < n:
            time.sleep(WAIT_S * (2 ** i))
    if last is None:
        raise ValueError("attempt needs at least one try")
    raise last from None


def http(method: str, url: str, what: str, **kwargs) -> requests.Response:
    """One request whose status is classified: a 2xx response, or SourceAbsent / RetryableSourceError.
    Transport errors propagate for `attempt` to classify."""
    r = requests.request(method, url, **kwargs)
    if not 200 <= r.status_code < 300:
        r.close()
        check_status(r.status_code, what)
    return r


def http_json(method: str, url: str, what: str, **kwargs):
    """The JSON body of a classified request, retried by `attempt`."""
    def get():
        with http(method, url, what, **kwargs) as r:
            return r.json()
    return attempt(get, what)
