#!/usr/bin/env python3
"""Pull known fungal antigen proteins from the IEDB Query API for use as the
cross-reactivity confounder panel in 02_score_antigens.py.

This REPLACES the old hand-assembled iedb_antigens.fa, whose provenance was
undocumented (86 UniProt records, no build script, no query recorded).

Data source: https://query-api.iedb.org (PostgREST over IEDB's curated
epitope/assay database; same backing data as the bulk "database_export_v3"
download, just queryable instead of a multi-GB CSV dump).

Query, in order:
  1. bcell_search (B-cell/antibody assays -- what matters for a serodiagnostic
     antigen, as opposed to tcell_search which is cellular/MHC-restricted).
  2. source_organism_iri_search contains the organism taxon (default
     NCBITaxon:4751 = Fungi kingdom; the _iri_search field is pre-expanded to
     the full NCBI taxonomy lineage, so filtering on the kingdom-level ID
     matches every fungal species below it, verified against real records).
  3. host_organism_iri_search contains one of the host taxa (default human
     NCBITaxon:9606 + mouse NCBITaxon:10090 -- the hosts actually relevant to
     a human serodiagnostic; NOT the full Mammalia clade, which would also
     pull in horse/dog/sheep vaccine-response records).
  4. structure_type = 'Linear peptide' (excludes non-peptidic/carbohydrate
     epitopes, which have no protein sequence to fetch).
  5. parent_source_antigen_iri IS NOT NULL (epitope maps to a real source
     protein, not just a free-text antigen description).
  6. qualitative_measure in the Positive family (excludes assays that came
     back Negative -- IEDB records negative results too).

Each surviving row's parent_source_antigen_iri is a "UNIPROT:<accession>"
string (verified against real data). Distinct accessions are deduplicated,
then their full-length sequences are fetched from UniProt (IEDB stores only
the short epitope peptide, not the full source protein).

Known limitation, checked directly against the API with every filter
relaxed: IEDB has NO curated Histoplasma or Blastomyces linear-peptide
protein antigens at all, and only 7 rows total for Paracoccidioides. The
resulting panel is dominated by Aspergillus fumigatus and Candida albicans.
That is a real gap in IEDB's curation, not a bug in this script -- see
docs/reports/ for the literature-curated Coccidioides-specific antigens that
fill part of this gap (not this script's job).

Usage:
    python3 00_download_iedb_antigens.py --out iedb_antigens.fa --meta iedb_antigens_meta.tsv
"""

import argparse
import csv
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

IEDB_API = "https://query-api.iedb.org/bcell_search"
UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
POSITIVE_MEASURES = ("Positive", "Positive-Low", "Positive-Intermediate", "Positive-High")
PAGE_SIZE = 1000
UNIPROT_BATCH = 90


def fetch_iedb_rows(organism_taxon: str, host_taxa: list[str]) -> list[dict]:
    """Page through bcell_search for the fungal-antigen filter described above."""
    host_or = ",".join(f"host_organism_iri_search.cs.{{NCBITaxon:{t}}}" for t in host_taxa)
    params = {
        "select": "parent_source_antigen_iri,parent_source_antigen_name,"
        "parent_source_antigen_source_org_name,host_organism_name,"
        "qualitative_measure,disease_names,reference_iri",
        "source_organism_iri_search": f"cs.{{NCBITaxon:{organism_taxon}}}",
        "structure_type": "eq.Linear peptide",
        "parent_source_antigen_iri": "not.is.null",
        "qualitative_measure": f"in.({','.join(POSITIVE_MEASURES)})",
        "or": f"({host_or})",
    }
    rows = []
    offset = 0
    while True:
        url = f"{IEDB_API}?{urllib.parse.urlencode(params, safe='(){},.')}"
        req = urllib.request.Request(
            url, headers={"Range-Unit": "items", "Range": f"{offset}-{offset + PAGE_SIZE - 1}"}
        )
        with urllib.request.urlopen(req, timeout=60) as resp:
            import json

            page = json.loads(resp.read())
        rows.extend(page)
        if len(page) < PAGE_SIZE:
            break
        offset += PAGE_SIZE
    return rows


def dedup_by_protein(rows: list[dict]) -> dict[str, dict]:
    """One record per parent_source_antigen_iri, keeping the accession only
    when it is a real UniProt-mapped protein (skips ChEBI/other non-protein
    antigen IDs, which have no sequence to fetch)."""
    proteins: dict[str, dict] = {}
    n_non_uniprot = 0
    for r in rows:
        iri = r["parent_source_antigen_iri"]
        if not iri or not iri.startswith("UNIPROT:"):
            n_non_uniprot += 1
            continue
        acc = iri.split(":", 1)[1]
        rec = proteins.setdefault(
            acc,
            {
                "accession": acc,
                "name": r["parent_source_antigen_name"],
                "source_organism": r["parent_source_antigen_source_org_name"],
                "hosts": set(),
                "qualitative_measures": set(),
                "diseases": set(),
                "n_assay_rows": 0,
            },
        )
        rec["hosts"].add(r["host_organism_name"] or "")
        rec["qualitative_measures"].add(r["qualitative_measure"] or "")
        for d in r.get("disease_names") or []:
            rec["diseases"].add(d)
        rec["n_assay_rows"] += 1
    print(
        f"  {n_non_uniprot} assay rows skipped (antigen not UniProt-mapped, e.g. carbohydrate/ChEBI)",
        file=sys.stderr,
    )
    return proteins


def fetch_uniprot_fasta(accessions: list[str]) -> str:
    """Batch-fetch full-length sequences. Chunked to keep query URLs short."""
    fasta_chunks = []
    for i in range(0, len(accessions), UNIPROT_BATCH):
        batch = accessions[i : i + UNIPROT_BATCH]
        query = " OR ".join(f"accession:{a}" for a in batch)
        url = f"{UNIPROT_STREAM}?query=({urllib.parse.quote(query)})&format=fasta"
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=60) as resp:
                    fasta_chunks.append(resp.read().decode())
                break
            except urllib.error.URLError:
                if attempt == 2:
                    raise
                time.sleep(2)
    return "".join(fasta_chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--organism-taxon", default="4751", help="NCBI taxon id, default Fungi kingdom")
    ap.add_argument(
        "--host-taxa",
        default="9606,10090",
        help="comma-separated NCBI taxon ids, default human+mouse",
    )
    ap.add_argument("--out", required=True, help="output FASTA of full-length antigen proteins")
    ap.add_argument("--meta", required=True, help="output TSV of per-protein IEDB provenance")
    args = ap.parse_args()

    host_taxa = args.host_taxa.split(",")
    print(
        f"Querying IEDB: organism=NCBITaxon:{args.organism_taxon}, hosts={host_taxa}",
        file=sys.stderr,
    )
    rows = fetch_iedb_rows(args.organism_taxon, host_taxa)
    print(f"  {len(rows)} positive B-cell linear-peptide assay rows", file=sys.stderr)

    proteins = dedup_by_protein(rows)
    print(f"  {len(proteins)} distinct UniProt-mapped source proteins", file=sys.stderr)

    accessions = sorted(proteins)
    print(f"Fetching {len(accessions)} sequences from UniProt...", file=sys.stderr)
    fasta_text = fetch_uniprot_fasta(accessions)
    Path(args.out).write_text(fasta_text)
    n_written = fasta_text.count(">")
    print(f"  wrote {n_written} sequences to {args.out}", file=sys.stderr)
    if n_written < len(accessions):
        print(
            f"  WARNING: {len(accessions) - n_written} accessions had no UniProt hit (obsolete/merged/demerged entries)",
            file=sys.stderr,
        )

    with open(args.meta, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "accession",
                "name",
                "source_organism",
                "hosts",
                "qualitative_measures",
                "diseases",
                "n_assay_rows",
            ]
        )
        for acc in accessions:
            r = proteins[acc]
            w.writerow(
                [
                    acc,
                    r["name"],
                    r["source_organism"],
                    ";".join(sorted(r["hosts"])),
                    ";".join(sorted(r["qualitative_measures"])),
                    ";".join(sorted(r["diseases"])),
                    r["n_assay_rows"],
                ]
            )
    print(f"  wrote metadata to {args.meta}", file=sys.stderr)


if __name__ == "__main__":
    main()
