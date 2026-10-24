#!/usr/bin/env python3
"""Fetch UniProt sequence + lineage + signal/GPI/repeat features for accessions named in committed tables.

Read-only public UniProt REST (https://rest.uniprot.org/uniprotkb/accessions), 100 accessions per request,
1 s between requests. Output goes to a cache OUTSIDE the repo (default ~/.cache/adhesionPred/aga1_class_scoping,
override with AGA1_SCOPING_CACHE). Nothing fetched is committed except what later scripts derive (ids, features).
Accessions already in the cache are skipped, so a re-run only fetches what is missing. Lineage is needed for every
accession (clade), so local sequences.faa files are NOT used to skip a fetch.
"""

import argparse
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import CACHE, all_named_accessions, uniprot_cache  # noqa: E402

FIELDS = "accession,reviewed,gene_primary,organism_id,organism_name,length,sequence,lineage,ft_signal,ft_lipid,ft_repeat"


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    pat = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
    named = all_named_accessions()
    skipped = sorted(x for x in named if not pat.match(x))
    if skipped:
        print("not UniProt accessions, not fetched:", skipped)
    acc = sorted({x for x in named if pat.match(x)} - set(uniprot_cache()))
    print(f"{len(acc)} accessions to fetch; cache {CACHE}")
    for i in range(0, len(acc), 100):
        out = CACHE / f"uniprot_batch_{acc[i]}_{len(acc[i:i + 100])}.tsv"
        if out.exists():
            continue
        qs = urllib.parse.urlencode(
            {
                "accessions": ",".join(acc[i : i + 100]),
                "format": "tsv",
                "fields": FIELDS,
                "size": 100,
            }
        )
        req = urllib.request.Request(
            f"https://rest.uniprot.org/uniprotkb/accessions?{qs}",
            headers={"User-Agent": "adhesionPred-scoping (read-only)"},
        )
        with urllib.request.urlopen(req, timeout=60) as r:
            out.write_bytes(r.read())
        print("fetched", out.name)
        time.sleep(1.0)


if __name__ == "__main__":
    main()
