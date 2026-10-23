"""Shared helpers: polite cached GET for the public InterPro and UniProt REST APIs (stdlib only).

Raw answers are cached outside the repo (default ~/.cache/adhesionPred-panther, override with PANTHER_CACHE).
Status 204/404 are cached as an empty answer. Run with `python3 -I`.
"""

import hashlib
import json
import os
import time
import urllib.error
import urllib.request
from pathlib import Path

CACHE = Path(os.environ.get("PANTHER_CACHE", Path.home() / ".cache" / "adhesionPred-panther"))
CACHE.mkdir(parents=True, exist_ok=True)
UA = "adhesionPred-scoping/1.0 (jason.stajich@ucr.edu; read-only)"
DELAY = 0.25


def get(url, want_headers=False, tries=7):
    """Return text (and the x-total-results header when want_headers) for url, cached on disk."""
    key = hashlib.sha1(url.encode()).hexdigest()
    path = CACHE / key
    if path.exists():
        d = json.loads(path.read_text())
        return (d["body"], d["total"]) if want_headers else d["body"]
    for k in range(tries):
        try:
            time.sleep(DELAY)
            with urllib.request.urlopen(
                urllib.request.Request(url, headers={"User-Agent": UA}), timeout=90
            ) as r:
                body = r.read().decode()
                total = r.headers.get("x-total-results", "")
            break
        except urllib.error.HTTPError as e:
            if e.code in (204, 404):
                body, total = "", "0"
                break
            time.sleep(3 * 2**k)
        except Exception:
            time.sleep(3 * 2**k)
    else:
        raise SystemExit("failed: " + url)
    path.write_text(json.dumps({"url": url, "body": body, "total": total}))
    return (body, total) if want_headers else body


def get_json(url):
    body = get(url)
    return json.loads(body) if body.strip() else {}
