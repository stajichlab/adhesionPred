#!/usr/bin/env python3.12
"""T0.2: fungal allergen positives, union by UniProt accession (spec section 3.1, decisions Q1 and union).

Sources (tables in analysis/allergen_scoping/, built 2026-10-04):
  iuis      iuis_fungal_isoallergens.tsv, column AccUniProt (WHO/IUIS fungal allergens)
  kw0020    uniprot_kw0020_fungi.tsv, column Entry (UniProt keyword KW-0020, Fungi, reviewed)

One row per accession. Sequences come from UniProt for every member, so both sources use the same sequence.
evidence_level:
  iuis                   listed by WHO/IUIS (expert allergen list)
  kw0020_curated         KW-0020 only, name without "homolog"
  similarity_only        KW-0020 only, name with "homolog" (spec: reported separately)
Where a protein is in IUIS and KW-0020, evidence_level is "iuis" (the stronger source).

Output (in --outdir):
  positives.faa   FASTA, IDs P_<accession>
  positives.tsv   accession, source flags, evidence_level, organism, length, sequence_sha256
  counts.md       counts by flag and evidence level, plus the structural checks

Structural checks (fail if violated): no duplicate accession; every accession fetched; no empty sequence.
Run: /usr/bin/python3.12 02_positives.py --tables <analysis/allergen_scoping dir> --outdir <dir>
"""

import argparse
import csv
import datetime
import hashlib
import os
import sys
import urllib.parse
import urllib.request

FIELDS = "accession,organism_name,length,sequence,protein_name"


def read_iuis(path):
    accs = set()
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for a in (r.get("AccUniProt") or "").replace(" ", "").split(","):
                if a and a.lower() != "nan":
                    accs.add(a)
    return accs


def read_kw(path):
    out = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = (r.get("Entry") or "").strip()
            if acc:
                out[acc] = r.get("Protein names") or ""
    return out


def fetch_batch(accs):
    query = " OR ".join(f"accession:{a}" for a in accs)
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        {"query": query, "format": "tsv", "fields": FIELDS}
    )
    req = urllib.request.Request(url, headers={"User-Agent": "adhesionPred-allergen-positives/1.0"})
    with urllib.request.urlopen(req, timeout=600) as resp:
        release = resp.headers.get("X-UniProt-Release", "unknown")
        return resp.read().decode("utf-8").rstrip("\n").split("\n"), release


def fetch_entry(acc):
    """Look up one accession. Returns (ok, record_or_None, primary_accession, entry_type).

    A merged accession is answered with a redirect to the current entry; urllib follows it, and the JSON
    primaryAccession gives the current ID. An inactive entry has no sequence."""
    import json

    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    req = urllib.request.Request(url, headers={"User-Agent": "adhesionPred-allergen-positives/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            d = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        return False, None, acc, f"http {e.code}"
    primary = d.get("primaryAccession", acc)
    etype = d.get("entryType", "")
    seq = (d.get("sequence") or {}).get("value", "")
    org = ((d.get("organism") or {}).get("scientificName")) or ""
    names = (
        d.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value", "")
    )
    if not seq or "Inactive" in etype:
        return True, None, primary, etype
    rec = {"org": org, "length": str(len(seq)), "seq": seq, "name": names}
    return True, rec, primary, etype


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    iuis = read_iuis(os.path.join(args.tables, "iuis_fungal_isoallergens.tsv"))
    kw = read_kw(os.path.join(args.tables, "uniprot_kw0020_fungi.tsv"))
    union = sorted(iuis | set(kw))
    if len(union) != len(set(union)):
        sys.exit("FATAL: duplicate accession in the union")
    print(f"iuis={len(iuis)} kw0020={len(kw)} shared={len(iuis & set(kw))} union={len(union)}")

    # Fetch in batches of 40 accessions to keep the query URL short.
    fetched, release = {}, "unknown"
    for i in range(0, len(union), 40):
        lines, release = fetch_batch(union[i : i + 40])
        if len(lines) < 2:
            continue
        header = lines[0].split("\t")
        idx = {h: j for j, h in enumerate(header)}
        for need in ("Entry", "Organism", "Length", "Sequence", "Protein names"):
            if need not in idx:
                sys.exit(f"FATAL: column '{need}' missing from UniProt response: {header}")
        for ln in lines[1:]:
            f = ln.split("\t")
            fetched[f[idx["Entry"]]] = {
                "org": f[idx["Organism"]],
                "length": f[idx["Length"]],
                "seq": f[idx["Sequence"]].strip(),
                "name": f[idx["Protein names"]],
            }

    # Accessions the batch query did not return: look each up by its own entry. A merged accession redirects to its
    # current primary accession (mapped and recorded). An inactive entry has no sequence and is rejected with the reason.
    missing = [a for a in union if a not in fetched]
    mapped, rejected = {}, []
    for a in missing:
        ok, rec, primary, etype = fetch_entry(a)
        if ok and rec and rec["seq"]:
            if primary != a:
                mapped[a] = primary
                if primary in fetched:
                    continue  # already present under its primary accession
            fetched[primary] = rec
        else:
            rejected.append((a, f"entry type {etype or 'unknown'}, no sequence"))
    empty = [a for a in union if a in fetched and not fetched[a]["seq"]]
    if empty:
        sys.exit(f"FATAL: empty sequence for {empty[:5]}")
    rows, done_keys = [], set()
    for a in union:
        key = mapped.get(a, a)  # a merged accession is written under its current primary accession
        if key in done_keys or key not in fetched:
            continue
        done_keys.add(key)
        in_iuis = a in iuis or any(m == key and s in iuis for s, m in mapped.items())
        in_kw = a in kw or any(m == key and s in kw for s, m in mapped.items())
        name_src = kw.get(a) or kw.get(key) or ""
        if in_iuis:
            level = "iuis"
        elif "homolog" in name_src.lower():
            level = "similarity_only"
        else:
            level = "kw0020_curated"
        rec = fetched[key]
        rows.append(
            {
                "acc": key,
                "iuis": in_iuis,
                "kw0020": in_kw,
                "level": level,
                "org": rec["org"],
                "length": rec["length"],
                "seq": rec["seq"],
                "name": rec["name"],
                "sha": hashlib.sha256(rec["seq"].encode()).hexdigest(),
                "merged_from": ";".join(sorted(s for s, m in mapped.items() if m == key)),
            }
        )

    with open(os.path.join(args.outdir, "positives.faa"), "w") as fa:
        for r in rows:
            fa.write(f">P_{r['acc']}\n")
            for i in range(0, len(r["seq"]), 60):
                fa.write(r["seq"][i : i + 60] + "\n")
    with open(os.path.join(args.outdir, "positives.tsv"), "w") as tsv:
        tsv.write(
            "accession\tin_iuis\tin_kw0020\tevidence_level\tmerged_from\torganism\tlength\tprotein_name\tsequence_sha256\n"
        )
        for r in rows:
            tsv.write(
                "\t".join(
                    [
                        r["acc"],
                        str(r["iuis"]).lower(),
                        str(r["kw0020"]).lower(),
                        r["level"],
                        r["merged_from"],
                        r["org"],
                        r["length"],
                        r["name"],
                        r["sha"],
                    ]
                )
                + "\n"
            )
    with open(os.path.join(args.outdir, "rejected.tsv"), "w") as rj:
        rj.write("accession\treason\n")
        for a, why in rejected:
            rj.write(f"{a}\t{why}\n")

    levels = {}
    for r in rows:
        levels[r["level"]] = levels.get(r["level"], 0) + 1
    today = datetime.date.today().isoformat()
    with open(os.path.join(args.outdir, "counts.md"), "w") as md:
        md.write("# T0.2 fungal allergen positives (union by accession)\n\n")
        md.write(f"- Run date: {today}\n- UniProt release: {release}\n")
        md.write(
            f"- IUIS accessions: {len(iuis)}\n- KW-0020 entries: {len(kw)}\n- Shared: {len(iuis & set(kw))}\n"
        )
        md.write(
            f"- Union: {len(union)}\n- Fetched from UniProt: {len(fetched)}\n- Rows written: {len(rows)}\n"
        )
        md.write(f"- Accessions not returned by UniProt: {len(missing)} {missing[:10]}\n\n")
        md.write("| Evidence level | Rows |\n|---|---|\n")
        for k in ("iuis", "kw0020_curated", "similarity_only"):
            md.write(f"| {k} | {levels.get(k, 0)} |\n")
    print(f"rows={len(rows)} missing={len(missing)} levels={levels}")


if __name__ == "__main__":
    main()
