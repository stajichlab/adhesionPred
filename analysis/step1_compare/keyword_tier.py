"""D10: the T-c keyword tier for V-kw with test proteins removed. Standard library only.

A T-c row is removed when it is a test protein (spec 2.3, Q4 and Q7). Reasons, first match wins:
1. literature_accession   accession is a row of eurotiomycetes_seeds.tsv
2. spombe_taxon           taxon is S. pombe (284812 or 4896); S. pombe is test-only (Q6)
3. heldout_accession      accession is a gene of a non-training truth source
4. no_sequence            no sequence was found, so the hash rule cannot be checked
5. literature_hash        exact cleaned sequence equals a literature row's sequence
6. heldout_hash           exact cleaned sequence equals a non-training truth protein
"""

import csv
import urllib.parse

import seqhash

SPOMBE_TAXA = frozenset({"284812", "4896"})
TC_SOURCE = "uniprot_surface_kw"
KEYWORD_COLUMNS = ("accession", "gene", "genome", "taxon_id", "length", "seq_sha256", "tier")
REMOVED_COLUMNS = ("accession", "gene", "genome", "taxon_id", "reason", "matched")


def read_seed_accessions(path) -> dict[str, str]:
    """Return accession -> gene for rows of eurotiomycetes_seeds.tsv that have an accession."""
    with open(path, encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    out = {}
    for row in csv.DictReader(lines, delimiter="\t"):
        query = (row.get("uniprot_query") or "").strip()
        if query.startswith("accession:"):
            out[query.split(":", 1)[1].strip()] = row["gene"]
    return out


def accession_fasta_urls(accessions: list[str], batch: int = 100) -> list[str]:
    urls = []
    for start in range(0, len(accessions), batch):
        query = " OR ".join(f"(accession:{a})" for a in accessions[start : start + batch])
        params = {"query": query, "format": "fasta", "size": "500"}
        urls.append(f"https://rest.uniprot.org/uniprotkb/search?{urllib.parse.urlencode(params)}")
    return urls


def build_keyword_tier(
    surface_rows: list[dict[str, str]],
    seq_by_acc: dict[str, str],
    seed_accessions: dict[str, str],
    heldout_accessions: dict[str, str],
    heldout_hashes: dict[str, str],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Return (kept T-c rows, removal log). heldout_* values name the matched test protein."""
    seed_hashes = {seqhash.seq_sha256(seq_by_acc[a]): a for a in seed_accessions if a in seq_by_acc}
    kept, removed = [], []
    for row in surface_rows:
        if row["source"] != TC_SOURCE:
            continue
        acc = row["accession"]
        seq = seq_by_acc.get(acc)
        digest = seqhash.seq_sha256(seq) if seq else ""
        reason, matched = "", ""
        if acc in seed_accessions:
            reason, matched = "literature_accession", seed_accessions[acc]
        elif row["taxon_id"] in SPOMBE_TAXA:
            reason, matched = "spombe_taxon", row["taxon_id"]
        elif acc in heldout_accessions:
            reason, matched = "heldout_accession", heldout_accessions[acc]
        elif not seq:
            reason = "no_sequence"
        elif digest in seed_hashes:
            reason, matched = "literature_hash", seed_hashes[digest]
        elif digest in heldout_hashes:
            reason, matched = "heldout_hash", heldout_hashes[digest]
        base = {k: row[k] for k in ("accession", "gene", "genome", "taxon_id")}
        if reason:
            removed.append({**base, "reason": reason, "matched": matched})
        else:
            length = str(len(seqhash.clean(seq)))
            kept.append({**base, "length": length, "seq_sha256": digest, "tier": "T-c"})
    return kept, removed
