#!/usr/bin/env python3.12
"""T0.3: weak-negative candidates from reviewed Swiss-Prot fungal proteins (spec section 3.2).

Filters (UniProt query, run on 2026-10-10 unless the output says otherwise):
  reviewed:true, taxonomy 4751 (Fungi), NOT keyword KW-0020 (Allergen), NOT cc_allergen:allergic, NOT antigen, NOT Cancer,
  length 10-850, existence:1 (evidence at protein level), annotation_score 4-5 (score above 3).
  Taxa 9606 (human) and 10239 (virus) are excluded by the Fungi filter already; the check is kept in the output as a count of 0.

Output (in --outdir):
  candidates.faa   FASTA, IDs N_<accession>. Input to clustering (T0.4) and pairing search (T0.5).
  candidates.tsv   accession, organism_id, organism_name, genus, length, annotation_score, sha256 of the sequence.
  counts.md        counts per genus and the total, plus the query and the run date.

The UniProt stream is the source; nothing is taken from memory. Re-running gives the same set only if UniProt is unchanged;
the output records the UniProt release from the response headers.
Run: /usr/bin/python3.12 01_candidates.py --outdir <dir>
"""

import argparse
import collections
import datetime
import hashlib
import os
import sys
import urllib.parse
import urllib.request

QUERY = (
    "reviewed:true AND taxonomy_id:4751 AND NOT keyword:KW-0020 AND NOT cc_allergen:allergic "
    "AND NOT antigen AND NOT Cancer AND length:[10 TO 850] AND existence:1 AND (annotation_score:4 OR annotation_score:5)"
)
FIELDS = "accession,organism_id,organism_name,length,annotation_score,sequence"


def fetch():
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        {"query": QUERY, "format": "tsv", "fields": FIELDS}
    )
    req = urllib.request.Request(
        url, headers={"User-Agent": "adhesionPred-allergen-candidates/1.0"}
    )
    with urllib.request.urlopen(req, timeout=600) as resp:
        release = resp.headers.get("X-UniProt-Release", "unknown")
        text = resp.read().decode("utf-8")
    return text, release


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument(
        "--exclude",
        help="positives.tsv: accessions in it are removed from candidates (spec 3.2, no candidate may be a positive)",
    )
    ap.add_argument(
        "--remove",
        help="file with one accession per line: candidates removed for near-identity to a positive (owner decision, 2026-10-10: fident >= 0.9)",
    )
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    excluded_pos = set()
    if args.exclude:
        with open(args.exclude, encoding="utf-8") as fh:
            next(fh)
            excluded_pos = {ln.split("\t")[0] for ln in fh if ln.strip()}
    removed_near = set()
    if args.remove:
        with open(args.remove, encoding="utf-8") as fh:
            removed_near = {ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")}

    text, release = fetch()
    lines = text.rstrip("\n").split("\n")
    if len(lines) < 2:
        sys.exit("FATAL: UniProt returned no rows; the query or the service is wrong")
    header = lines[0].split("\t")
    idx = {h: i for i, h in enumerate(header)}
    for need in ("Entry", "Organism (ID)", "Organism", "Length", "Annotation", "Sequence"):
        if need not in idx:
            sys.exit(f"FATAL: column '{need}' missing; got {header}")

    rows, seen, removed_pos = [], set(), []
    for ln in lines[1:]:
        f = ln.split("\t")
        acc = f[idx["Entry"]]
        if acc in seen:
            sys.exit(f"FATAL: duplicate accession {acc}")
        seen.add(acc)
        if acc in excluded_pos:
            removed_pos.append(acc)
            continue
        if acc in removed_near:
            removed_pos.append(acc + " (near-identical to a positive)")
            continue
        seq = f[idx["Sequence"]].strip()
        if not seq:
            sys.exit(f"FATAL: {acc} has no sequence")
        org_id = f[idx["Organism (ID)"]]
        if org_id in ("9606", "10239"):
            sys.exit(f"FATAL: {acc} is from excluded taxon {org_id}")
        name = f[idx["Organism"]]
        genus = name.split(" ")[0] if name else "unknown"
        rows.append(
            {
                "acc": acc,
                "org_id": org_id,
                "org": name,
                "genus": genus,
                "length": f[idx["Length"]],
                "score": f[idx["Annotation"]],
                "seq": seq,
                "sha": hashlib.sha256(seq.encode()).hexdigest(),
            }
        )

    # Write FASTA with N_ prefix; wrap at 60 columns.
    with open(os.path.join(args.outdir, "candidates.faa"), "w") as fa:
        for r in rows:
            fa.write(f">N_{r['acc']}\n")
            for i in range(0, len(r["seq"]), 60):
                fa.write(r["seq"][i : i + 60] + "\n")
    with open(os.path.join(args.outdir, "candidates.tsv"), "w") as tsv:
        tsv.write(
            "accession\torganism_id\torganism_name\tgenus\tlength\tannotation_score\tsequence_sha256\n"
        )
        for r in rows:
            tsv.write(
                "\t".join(
                    [r["acc"], r["org_id"], r["org"], r["genus"], r["length"], r["score"], r["sha"]]
                )
                + "\n"
            )

    by_genus = collections.Counter(r["genus"] for r in rows)
    today = datetime.date.today().isoformat()
    with open(os.path.join(args.outdir, "counts.md"), "w") as md:
        md.write("# T0.3 weak-negative candidates\n\n")
        md.write(f"- Run date: {today}\n- UniProt release: {release}\n")
        md.write(
            f"- Query: `{QUERY}`\n- Candidates: {len(rows)}\n- Excluded taxa 9606 and 10239 found: 0 (checked)\n"
        )
        md.write(
            f"- Removed because they are positives (--exclude): {len(removed_pos)} {removed_pos}\n\n"
        )
        md.write("| Genus | Candidates |\n|---|---|\n")
        for g, n in by_genus.most_common():
            md.write(f"| {g} | {n} |\n")
    print(f"candidates={len(rows)} genera={len(by_genus)} release={release}")


if __name__ == "__main__":
    main()
