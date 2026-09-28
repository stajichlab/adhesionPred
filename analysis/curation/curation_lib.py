"""Shared helpers for building curated label tables under data/curated/<category>/.

Network access to www.ebi.ac.uk (QuickGO), rest.uniprot.org and query-api.iedb.org is required.
"""

import csv
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CURATED = REPO / "data" / "curated"
FUNGI = 4751
# Tier-T5 reference genomes (UniProt reference-proteome taxon ids)
REFERENCE_TAXA = {
    559292: "Saccharomyces cerevisiae S288C",
    237561: "Candida albicans SC5314",
    284593: "Nakaseomyces glabratus CBS138",
    498019: "Candidozyma auris B8441",
    284812: "Schizosaccharomyces pombe 972h-",
    330879: "Aspergillus fumigatus Af293",
}
UNIPROT_FIELDS = (
    "accession,reviewed,gene_primary,protein_name,organism_name,organism_id,length,"
    "ft_signal,keyword,xref_pfam"
)
EXPERIMENTAL_ECO = "ECO:0000269"


def log(msg):
    print(msg, file=sys.stderr)


def get(url, accept="application/json", tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": accept})
            return urllib.request.urlopen(req, timeout=120).read().decode()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(5)


def quickgo_experimental(go_terms, taxon=FUNGI):
    """Experimental GO annotations to go_terms (with is_a/part_of descendants) in taxon.

    Returns {uniprot_accession: [(go_id, evidence, reference, assigned_by, qualifier), ...]}.
    """
    base = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
    ann = defaultdict(list)
    for go in go_terms:
        page = 1
        while True:
            q = urllib.parse.urlencode(
                {
                    "goId": go,
                    "goUsage": "descendants",
                    "goUsageRelationships": "is_a,part_of",
                    "taxonId": taxon,
                    "taxonUsage": "descendants",
                    "evidenceCode": EXPERIMENTAL_ECO,
                    "evidenceCodeUsage": "descendants",
                    "limit": 200,
                    "page": page,
                }
            )
            d = json.loads(get(f"{base}?{q}"))
            for r in d["results"]:
                acc = r["geneProductId"].split(":", 1)[1]
                if re.fullmatch(r"[A-Z0-9]{6,10}", acc):
                    ann[acc].append(
                        (
                            r["goId"],
                            r["goEvidence"],
                            r["reference"],
                            r["assignedBy"],
                            r["qualifier"],
                        )
                    )
            if page >= d["pageInfo"]["total"]:
                break
            page += 1
    return ann


def _tsv(text):
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def uniprot_accessions(accs, fields=UNIPROT_FIELDS):
    rows = []
    accs = sorted(set(accs))
    for i in range(0, len(accs), 90):
        q = urllib.parse.urlencode(
            {"accessions": ",".join(accs[i : i + 90]), "fields": fields, "format": "tsv"}
        )
        rows += _tsv(get(f"https://rest.uniprot.org/uniprotkb/accessions?{q}", "text/plain"))
    return {r["Entry"]: r for r in rows}


def uniprot_search(query, size=25, fields=UNIPROT_FIELDS):
    q = urllib.parse.urlencode({"query": query, "fields": fields, "format": "tsv", "size": size})
    return _tsv(get(f"https://rest.uniprot.org/uniprotkb/search?{q}", "text/plain"))


def uniprot_stream(query, fields=UNIPROT_FIELDS):
    """All hits for a query (UniProt stream endpoint, for large result sets)."""
    q = urllib.parse.urlencode({"query": query, "fields": fields, "format": "tsv"})
    return _tsv(get(f"https://rest.uniprot.org/uniprotkb/stream?{q}", "text/plain"))


def read_seeds(path):
    lines = [ln for ln in open(path) if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def is_surface(u):
    return bool(u["Signal peptide"]) or bool(
        re.search(r"GPI-anchor|Cell wall|Secreted", u["Keywords"])
    )


def base_row(u, source):
    return {
        "accession": u["Entry"],
        "reviewed": u["Reviewed"],
        "gene": u["Gene Names (primary)"],
        "protein_name": u["Protein names"][:120],
        "organism": u["Organism"],
        "taxon_id": u["Organism (ID)"],
        "length": u["Length"],
        "signal_peptide": "yes" if u["Signal peptide"] else "no",
        "gpi_anchor": "yes" if "GPI-anchor" in u["Keywords"] else "no",
        "surface": "yes" if is_surface(u) else "no",
        "pfam": u["Pfam"].rstrip(";"),
        "source": source,
    }


def go_summary(hits):
    codes = sorted({h[1] for h in hits})
    pmids = sorted({h[2].split(":", 1)[1] for h in hits if h[2].startswith("PMID:")})
    terms = sorted({h[0] for h in hits})
    quals = sorted({h[4] for h in hits})
    return {
        "evidence_codes": ",".join(codes),
        "pmids": ";".join(pmids),
        "go_terms": ",".join(terms),
        "go_qualifiers": ",".join(quals),
    }


def apply_overrides(out, path):
    """Apply documented manual reclassifications (accession, cls, evidence_level, note)."""
    if not path.exists():
        return
    for ov in read_seeds(path):
        r = out.get(ov["accession"])
        if r is None:
            log(f"  override for absent accession {ov['accession']}")
            continue
        r.update(cls=ov["cls"], evidence_level=ov["evidence_level"])
        r["evidence_summary"] = (
            f"{ov['note']} [manual override; was: {r.get('evidence_summary', '')}]"
        )


def write_table(rows, path, cols, order_key):
    rows = sorted(rows, key=order_key)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {len(rows)} rows to {path.relative_to(REPO)}")
