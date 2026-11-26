"""Gene-model check for round 1 targets that lack a signal peptide (for example M7WU28, M7WQR1).

For each accession: the UniProt record (length, first residues, annotated signal peptide, ORF name, proteome) and the other
Rhodotorula proteins in the same UniRef50 cluster (accession, organism, length, signal peptide, first residues). A longer member
from another strain with a signal peptide suggests the NP11 gene model starts too late. A few dozen HTTPS requests; HPCC login node.

  python3 -I analysis/crossclade_domain_search/gene_model_check.py --acc M7WU28 --acc M7WQR1 --out /path/NEW_FILE.tsv
"""

import argparse
import csv
import io
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

UNIPROT = "https://rest.uniprot.org"
UA = {"User-Agent": "crossclade-domain-search/2 (python urllib)"}
RHODOTORULA_TAXID = 5533
NTERM = 40


class CheckError(Exception):
    pass


def http_get(url):
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.read().decode()
    except urllib.error.HTTPError as e:
        raise CheckError(f"{url}: HTTP {e.code}") from e
    except urllib.error.URLError as e:
        raise CheckError(f"{url}: {e}") from e


def _pos(loc, key):
    v = loc.get(key, {}).get("value")
    return "?" if v is None else str(v)


def entry_summary(entry):
    """Fields of a UniProtKB JSON entry used here."""
    if "sequence" not in entry:
        raise CheckError(f"{entry.get('primaryAccession', '?')}: no sequence (obsolete entry?)")
    seq = entry["sequence"]["value"]
    sig = [f for f in entry.get("features", []) if f.get("type") == "Signal"]
    orf = [o["value"] for g in entry.get("genes", []) for o in g.get("orfNames", [])]
    return {
        "length": len(seq),
        "nterm": seq[:NTERM],
        "signal": f"{_pos(sig[0]['location'], 'start')}-{_pos(sig[0]['location'], 'end')}"
        if sig
        else "none",
        "orf": ",".join(orf) or "NA",
        "organism": entry.get("organism", {}).get("scientificName", "NA"),
    }


def uniref50_id(acc, fetch=http_get):
    """The UniRef50 cluster holding acc; tries the `uniprot_id` field, then `uniprotkb` (both appear in UniProt's field lists)."""
    last, answered = None, False
    for field in ("uniprot_id", "uniprotkb"):
        q = urllib.parse.quote(f"{field}:{acc} AND identity:0.5")
        try:
            data = json.loads(fetch(f"{UNIPROT}/uniref/search?query={q}&fields=id&format=json"))
        except CheckError as e:
            last = e
            continue
        answered = True
        ids = [r["id"] for r in data.get("results", [])]
        if ids:
            return ids[0]
    if last is not None and not answered:
        raise last
    return None


def cluster_members(cluster, fetch=http_get):
    q = urllib.parse.quote(f"uniref_cluster_50:{cluster} AND taxonomy_id:{RHODOTORULA_TAXID}")
    tsv = fetch(
        f"{UNIPROT}/uniprotkb/search?query={q}&fields=accession,organism_name,length,ft_signal,sequence&format=tsv&size=500"
    )
    return list(csv.DictReader(io.StringIO(tsv), delimiter="\t"))


def check(acc, fetch=http_get):
    e = entry_summary(json.loads(fetch(f"{UNIPROT}/uniprotkb/{acc}.json")))
    rows = [
        {
            "query": acc,
            "member": acc,
            "organism": e["organism"],
            "length": e["length"],
            "signal": e["signal"],
            "orf": e["orf"],
            "uniref50": "",
            "nterm": e["nterm"],
        }
    ]
    cl = uniref50_id(acc, fetch)
    rows[0]["uniref50"] = cl or "NA"
    if cl:
        for m in cluster_members(cl, fetch):
            if m["Entry"] == acc:
                continue
            sig = m.get("Signal peptide", "")
            rows.append(
                {
                    "query": acc,
                    "member": m["Entry"],
                    "organism": m["Organism"],
                    "length": int(m["Length"]),
                    "signal": sig.split(";")[0].replace("SIGNAL ", "").replace("..", "-")
                    if sig
                    else "none",
                    "orf": "NA",
                    "uniref50": cl,
                    "nterm": m["Sequence"][:NTERM],
                }
            )
    return rows


def main(argv=None, fetch=http_get):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--acc", action="append", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    try:
        out = Path(a.out)
        if out.exists():
            raise CheckError(f"{out} exists; nothing is overwritten")
        rows = [r for acc in a.acc for r in check(acc, fetch)]
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {out}: {len(rows)} rows")
    except CheckError as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
