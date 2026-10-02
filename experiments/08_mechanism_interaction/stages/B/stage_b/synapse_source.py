"""Synapse behind the two calls the collect phase needs (fetch.SynapseLike). Imports synapseclient,
which is installed in the Modal image only.

A file is downloaded through a pre-signed link from Synapse's file-handle service
(POST /fileHandle/batch), not through synapseclient's own downloader, so the collect phase resumes
it with HTTP Range requests and checkpoints it like every other file (stage_b/collect.py). The link
is a credential of a few minutes' life: it is requested again at every (re)connection and is never
written anywhere. The same response gives the file's name, size and MD5.
"""
import json
from collections.abc import Callable
from typing import TypeVar

import synapseclient
from synapseclient.core.exceptions import SynapseAuthenticationError, SynapseError, SynapseHTTPError

from stage_b.remote import check_status
from stage_b.schemas import RetryableSourceError

T = TypeVar("T")


def classified(call: Callable[[], T], what: str) -> T:
    """`call()`, with Synapse's errors mapped to stage B's classes (stage_b/remote.py)."""
    try:
        return call()
    except SynapseAuthenticationError:
        raise RetryableSourceError("auth", f"{what}: Synapse refused the credential") from None
    except SynapseHTTPError as err:
        code = getattr(getattr(err, "response", None), "status_code", None)
        if isinstance(code, int):
            check_status(code, what)
        raise RetryableSourceError("protocol", f"{what}: {type(err).__name__}") from None
    except SynapseError as err:
        raise RetryableSourceError("protocol", f"{what}: {type(err).__name__}") from None


class SynapseSource:
    def __init__(self, token: str):
        self.syn = synapseclient.Synapse(silent=True)
        classified(lambda: self.syn.login(authToken=token), "Synapse login")

    def children(self, parent: str) -> list[dict]:
        return classified(lambda: [{"name": k["name"], "id": k["id"]}
                                   for k in self.syn.getChildren(parent, includeTypes=["file"])], f"Synapse {parent} listing")

    def file(self, entity_id: str) -> dict:
        what = f"Synapse {entity_id}"
        entity = classified(lambda: self.syn.restGET(f"/entity/{entity_id}"), what)
        request = {"includeFileHandles": True, "includePreSignedURLs": True, "includePreviewPreSignedURLs": False,
                   "requestedFiles": [{"fileHandleId": entity["dataFileHandleId"], "associateObjectId": entity_id,
                                       "associateObjectType": "FileEntity"}]}
        reply = classified(lambda: self.syn.restPOST("/fileHandle/batch", body=json.dumps(request),
                                                     endpoint=self.syn.fileHandleEndpoint), what)
        got = reply["requestedFiles"][0]
        if "preSignedURL" not in got or "fileHandle" not in got:
            raise RetryableSourceError("auth" if got.get("failureCode") == "UNAUTHORIZED" else "protocol",
                                       f"{what}: no download link ({got.get('failureCode', 'no failure code')})")
        handle = got["fileHandle"]
        return {"name": handle["fileName"], "size": handle.get("contentSize"), "md5": handle.get("contentMd5", ""),
                "url": got["preSignedURL"]}
