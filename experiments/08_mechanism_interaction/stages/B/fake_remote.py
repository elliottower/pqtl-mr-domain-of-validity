"""A local HTTP server standing in for stage B's sources, for the tests and the dry run. It serves
only what a test or synthetic_unit.py puts in it; nothing here reaches a real source.

What it reproduces of the real ones:
- files with ETag, Last-Modified and Content-Length, HTTP Range requests (206 with Content-Range,
  416 past the end), so resumable downloads, tabix and bcftools work against it; the ETag is the
  MD5 of the bytes unless `etags` fixes one for the path, which is how a source that serves other
  bytes under an unchanged name, size and ETag is simulated;
- deCODE's download endpoint, `/<source>/s3/download?token=<token>&file=<Key>`: 403 for any token
  but the folder's current one (`tokens`), so a link can expire and be re-issued; 404 for a file
  the folder does not hold; otherwise it hands over a signed address `/signed/<signature>/...`,
  either by a 302 redirect or, for a source listed in `json_links`, in a JSON body;
- pre-signed links, `/synapse/<signature>/<file>` and `/signed/<signature>/<source>/<file>`: 403
  for a signature that is not current;
- JSON endpoints registered by the caller (`json_routes`, `json_prefixes`); a route that returns
  `bytes` has them sent as the body unchanged, so a reply that is not JSON can be served;
- faults, queued per path in `faults` and consumed one per request: ("status", code) answers with
  that status; ("break", n) declares the full length, sends n bytes and drops the connection.

`log` keeps (method, path, Range header) of every request, with tokens and signatures taken out
(`/<source>/s3/download` for the endpoint, `/<source>/<file>` for the file behind a signed address),
so a test can show which files were read from the network and which were not.
"""
import hashlib
import json
import socket
import threading
from collections.abc import Callable, Mapping
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlsplit

JsonRoute = Callable[[dict, bytes, Mapping], tuple[int, object]]
LAST_MODIFIED = "Wed, 30 Sep 2026 00:00:00 GMT"


class FakeRemote:
    def __init__(self):
        self.files: dict[str, bytes] = {}
        self.etags: dict[str, str] = {}        # path -> an ETag served whatever the bytes are
        self.json_routes: dict[tuple[str, str], JsonRoute] = {}
        self.json_prefixes: list[tuple[str, str, Callable[[str, dict, bytes, Mapping], tuple[int, object]]]] = []
        self.tokens: dict[str, str] = {}
        self.json_links: set[str] = set()      # sources whose download endpoint answers with JSON, not a redirect
        self.signatures: set[str] = set()
        self._signed = 0
        self.faults: dict[str, list[tuple[str, int]]] = {}
        self.log: list[tuple[str, str, str]] = []
        self._server: ThreadingHTTPServer | None = None

    # ---- lifecycle --------------------------------------------------------------------------------
    def start(self) -> str:
        remote = self

        class Handler(BaseHTTPRequestHandler):
            protocol_version = "HTTP/1.1"

            def log_message(self, format, *args):  # noqa: A002  (the base class's signature)
                return

            def do_GET(self):
                remote._serve(self, "GET")

            def do_HEAD(self):
                remote._serve(self, "HEAD")

            def do_POST(self):
                remote._serve(self, "POST")

        self._server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self._server.daemon_threads = True
        self._server.handle_error = lambda request, client_address: None   # a client dropping a connection is expected here
        threading.Thread(target=self._server.serve_forever, daemon=True).start()
        return self.base

    def stop(self) -> None:
        if self._server is not None:
            self._server.shutdown()
            self._server.server_close()
            self._server = None

    def __enter__(self) -> "FakeRemote":
        self.start()
        return self

    def __exit__(self, *exc) -> None:
        self.stop()

    @property
    def base(self) -> str:
        if self._server is None:
            raise RuntimeError("the fake remote is not started")
        return f"http://127.0.0.1:{self._server.server_address[1]}"

    # ---- content ----------------------------------------------------------------------------------
    def etag(self, path: str) -> str:
        return self.etags.get(path) or hashlib.md5(self.files[path], usedforsecurity=False).hexdigest()

    def paths_read(self, since: int = 0) -> list[str]:
        """The paths requested from position `since` of the log on."""
        return [path for _method, path, _range in self.log[since:]]

    # ---- requests ---------------------------------------------------------------------------------
    def _logical(self, path: str) -> tuple[int | None, str]:
        """(status to refuse with, path with the credential taken out)."""
        parts = path.split("/")
        if len(parts) >= 5 and parts[1] == "signed":
            return (None if parts[2] in self.signatures else 403), "/" + "/".join(parts[3:])
        if len(parts) >= 4 and parts[1] == "synapse":
            return (None if parts[2] in self.signatures else 403), "/" + "/".join(["synapse", *parts[3:]])
        return None, path

    def _send(self, h: BaseHTTPRequestHandler, status: int, body: bytes, headers: Mapping[str, str] | None = None,
              head: bool = False) -> None:
        h.send_response(status)
        for k, v in {"Content-Length": str(len(body)), **(headers or {})}.items():
            h.send_header(k, v)
        h.end_headers()
        if not head:
            h.wfile.write(body)

    def _download_endpoint(self, h: BaseHTTPRequestHandler, source: str, query: dict) -> None:
        """deCODE's `/s3/download`: check the folder token, then hand over a signed address."""
        token, key = query.get("token", [""])[0], query.get("file", [""])[0]
        if self.tokens.get(source) != token:
            return self._send(h, 403, b"")
        if f"/{source}/{key}" not in self.files:
            return self._send(h, 404, b"")
        self._signed += 1
        signature = f"signed{self._signed:04d}"
        self.signatures.add(signature)
        link = f"{self.base}/signed/{signature}/{source}/{key}"
        if source in self.json_links:
            return self._send(h, 200, json.dumps({"file": key, "url": link}).encode(), {"Content-Type": "application/json"})
        return self._send(h, 302, b"", {"Location": link})

    def _serve(self, h: BaseHTTPRequestHandler, method: str) -> None:
        url = urlsplit(h.path)
        refused, path = self._logical(url.path)
        self.log.append((method, path, h.headers.get("Range", "")))
        length = int(h.headers.get("Content-Length") or 0)
        body = h.rfile.read(length) if length else b""
        if refused is not None:
            return self._send(h, refused, b"")
        if path.endswith("/s3/download"):
            return self._download_endpoint(h, path.split("/")[1], parse_qs(url.query))
        fault = self.faults[path].pop(0) if self.faults.get(path) else None
        if fault is not None and fault[0] == "status":
            return self._send(h, fault[1], b"")
        route = self.json_routes.get((method, path))
        if route is None:
            for m, prefix, fn in self.json_prefixes:
                if m == method and path.startswith(prefix):
                    route = lambda q, b, hd, fn=fn, rest=path[len(prefix):]: fn(rest, q, b, hd)  # noqa: E731
                    break
        if route is not None:
            status, obj = route(parse_qs(url.query), body, h.headers)
            return self._send(h, status, obj if isinstance(obj, bytes) else json.dumps(obj).encode(),
                              {"Content-Type": "application/json"})
        if path not in self.files:
            return self._send(h, 404, b"")
        data = self.files[path]
        headers = {"ETag": f'"{self.etag(path)}"', "Last-Modified": LAST_MODIFIED, "Accept-Ranges": "bytes"}
        start, end, status = 0, len(data) - 1, 200
        rng = h.headers.get("Range", "")
        if rng.startswith("bytes="):
            a, _, b = rng[len("bytes="):].partition("-")
            start = int(a) if a else max(len(data) - int(b), 0)
            end = min(int(b), len(data) - 1) if a and b else len(data) - 1
            if start >= len(data):
                return self._send(h, 416, b"", {"Content-Range": f"bytes */{len(data)}"})
            status, headers["Content-Range"] = 206, f"bytes {start}-{end}/{len(data)}"
        chunk = data[start:end + 1]
        if fault is not None and fault[0] == "break" and method == "GET":
            h.send_response(status)
            for k, v in {"Content-Length": str(len(chunk)), **headers}.items():
                h.send_header(k, v)
            h.end_headers()
            h.wfile.write(chunk[:fault[1]])
            h.wfile.flush()
            h.connection.shutdown(socket.SHUT_RDWR)
            h.close_connection = True
            return None
        return self._send(h, status, chunk, headers, head=method == "HEAD")


class FakeSynapse:
    """fetch.SynapseLike over a FakeRemote: every `file` call issues a new pre-signed link and
    withdraws the earlier ones, as a link that has expired. `refuse` simulates a rejected token."""

    def __init__(self, remote: FakeRemote, folders: Mapping[str, Mapping[str, str]], refuse: BaseException | None = None):
        self.remote, self.folders, self.refuse = remote, folders, refuse
        self.calls = 0

    def children(self, parent: str) -> list[dict]:
        if self.refuse is not None:
            raise self.refuse
        return [{"name": name, "id": entity} for name, entity in sorted(self.folders.get(parent, {}).items())]

    def file(self, entity_id: str) -> dict:
        if self.refuse is not None:
            raise self.refuse
        name = next(n for kids in self.folders.values() for n, e in kids.items() if e == entity_id)
        self.calls += 1
        signature = f"sig{self.calls:04d}"
        self.remote.signatures.clear()
        self.remote.signatures.add(signature)
        data = self.remote.files[f"/synapse/{name}"]
        return {"name": name, "size": len(data), "md5": hashlib.md5(data, usedforsecurity=False).hexdigest(),
                "url": f"{self.remote.base}/synapse/{signature}/{name}"}
