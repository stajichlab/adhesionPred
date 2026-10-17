"""Shared helpers for the family-matched negative build (stdlib only)."""

import csv
import time
import urllib.parse
import urllib.request

BASE = "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model"
FAM = BASE + "/allkingdom/family"
UA = "adhesionPred-family-negatives/1.0"
PLANT_TYPES = {"Aero Plant", "Food Plant", "Gliadin"}


def http_get(url, retries=5):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as resp:
                body = resp.read().decode("utf-8")
                return body, resp.headers
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"GET failed after {retries} tries: {url[:200]}: {last}")


def uniprot_search_all(query, fields):
    """Return all rows (list of dicts) for a UniProt search, following pagination."""
    params = {"query": query, "fields": fields, "format": "tsv", "size": "500"}
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)
    rows, header_line = [], None
    total = None
    while url:
        body, headers = http_get(url)
        lines = body.splitlines()
        if header_line is None:
            header_line = lines[0].split("\t")
            total = headers.get("X-Total-Results")
        for line in lines[1:]:
            if line:
                vals = line.split("\t")
                rows.append(dict(zip(header_line, vals, strict=False)))
        url = None
        link = headers.get("Link") or headers.get("link")
        if link and 'rel="next"' in link:
            url = link.split(";")[0].strip().lstrip("<").rstrip(">")
        time.sleep(0.2)
    return rows, total


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def positives():
    """Return list of dicts: accession, organism, length, protein_name, type, group."""
    types = {}
    for r in read_tsv(BASE + "/allergenonline_all/allergenonline_fungi.tsv"):
        types.setdefault(r["uniprot"], set()).add(r["type"])
    out = []
    for r in read_tsv(BASE + "/allergenonline_all/allergenonline_new.tsv"):
        acc = r["accession"]
        tset = types.get(acc, set())
        group = "plant" if (tset & PLANT_TYPES) else "animal"
        out.append(
            {
                "accession": acc,
                "organism": r["organism"],
                "length": r["length"],
                "protein_name": r["protein_name"],
                "type": ";".join(sorted(tset)),
                "group": group,
            }
        )
    return out
