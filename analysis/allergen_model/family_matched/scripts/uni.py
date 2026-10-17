"""Shared UniProt REST helpers (stdlib only)."""

import re
import time
import urllib.parse
import urllib.request

BASE = "https://rest.uniprot.org/uniprotkb/search"
PF_RE = re.compile(r"PF\d{5}")


def get_all(query, fields, size=500, retries=4):
    """Return list of dict rows (TSV parsed) for a UniProt query, following Link pagination."""
    url = (
        BASE
        + "?"
        + urllib.parse.urlencode({"query": query, "format": "tsv", "fields": fields, "size": size})
    )
    rows, header, release = [], None, None
    while url:
        for attempt in range(retries):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "family_matched/1.0"})
                with urllib.request.urlopen(req, timeout=120) as r:
                    text = r.read().decode("utf-8")
                    release = r.headers.get("X-UniProt-Release", release)
                    link = r.headers.get("Link", "")
                break
            except Exception:
                if attempt == retries - 1:
                    raise
                time.sleep(3 * (attempt + 1))
        lines = text.rstrip("\n").split("\n")
        if header is None:
            header = lines[0].split("\t")
            body = lines[1:]
        else:
            body = lines[1:]
        for ln in body:
            if ln:
                rows.append(dict(zip(header, ln.split("\t"), strict=False)))
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else None
    return rows, release


def pfams(s):
    return PF_RE.findall(s or "")
