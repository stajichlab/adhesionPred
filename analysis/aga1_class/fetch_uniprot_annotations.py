"""Fetch public UniProt annotations (domains, keywords, location) for a list of accessions (stdlib only; needs network).

Reads one accession per line, asks rest.uniprot.org (accessions endpoint, batches) for: gene names, protein name, length, Pfam and InterPro
cross-references, domain feature names, keywords, subcellular location, review status, and writes one TSV with the UniProt release the
answer came from. The text is public annotation, not sequence. Downloaded data is untrusted: read it with `python3 -I`.

  python3 -I analysis/aga1_class/fetch_uniprot_annotations.py --accessions LIST.txt --out annotations.tsv [--batch 200]

Refuses an existing --out. FATAL: and exit 1 if UniProt does not return every accession (the missing ones are listed), or on a network
or HTTP error after three tries. Accessions that UniProt has merged or deleted are therefore a hard stop, not a silent gap.
"""

import argparse
import csv
import io
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

URL = "https://rest.uniprot.org/uniprotkb/accessions"
FIELDS = "accession,gene_names,protein_name,length,xref_pfam,xref_interpro,ft_domain,keyword,cc_subcellular_location,reviewed"
OUT_COLS = [
    "accession", "gene_names", "protein_name", "length", "pfam", "interpro", "domains", "keywords", "subcellular_location",
    "reviewed", "uniprot_release",
]  # fmt: skip


class FetchError(Exception):
    pass


def domain_names(ft):
    """'DOMAIN 87..249; /note="PA14"; /evidence=...' entries -> 'PA14' (semicolon-joined, order kept)."""
    return ";".join(re.findall(r'/note="([^"]+)"', ft))


def clean_location(text):
    text = re.sub(r"^SUBCELLULAR LOCATION:\s*", "", text)
    return re.sub(r"\s*\{[^}]*\}", "", text).strip()


def parse_tsv(text, release):
    rows = {}
    for r in csv.DictReader(io.StringIO(text), delimiter="\t"):
        rows[r["Entry"]] = {
            "accession": r["Entry"], "gene_names": r["Gene Names"], "protein_name": r["Protein names"], "length": r["Length"],
            "pfam": r["Pfam"].rstrip(";"), "interpro": r["InterPro"].rstrip(";"), "domains": domain_names(r["Domain [FT]"]),
            "keywords": r["Keywords"], "subcellular_location": clean_location(r["Subcellular location [CC]"]),
            "reviewed": r["Reviewed"], "uniprot_release": release,
        }  # fmt: skip
    return rows


def fetch_batch(batch, tries=3):
    q = urllib.parse.urlencode({"accessions": ",".join(batch), "fields": FIELDS, "format": "tsv"})
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(f"{URL}?{q}", timeout=120) as resp:
                return resp.read().decode("utf-8"), resp.headers.get("X-UniProt-Release", "NA")
        except OSError as e:
            last = e
            time.sleep(2 * (i + 1))
    raise FetchError(f"UniProt request failed after {tries} tries: {last}")


def fetch(accessions, batch_size=200):
    rows, releases = {}, set()
    for i in range(0, len(accessions), batch_size):
        text, release = fetch_batch(accessions[i : i + batch_size])
        releases.add(release)
        rows.update(parse_tsv(text, release))
    missing = [a for a in accessions if a not in rows]
    if missing:
        raise FetchError(f"UniProt did not return {len(missing)} accessions, first: {missing[:10]}")
    if len(releases) != 1:
        raise FetchError(f"the batches came from several UniProt releases: {sorted(releases)}")
    return [rows[a] for a in accessions]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--accessions", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch", type=int, default=200)
    a = ap.parse_args(argv)
    if Path(a.out).exists():
        print(f"FATAL: {a.out} exists; nothing is overwritten", file=sys.stderr)
        return 1
    try:
        acc = Path(a.accessions).read_text().split()
        if not acc or len(set(acc)) != len(acc):
            raise FetchError("the accession list is empty or has duplicates")
        rows = fetch(acc, a.batch)
    except (FetchError, OSError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {a.out}: {len(rows)} proteins, UniProt release {rows[0]['uniprot_release']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
