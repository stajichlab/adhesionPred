#!/usr/bin/env python3.12
"""COMPARE fungal allergen rows (to_import/allergen/compare.tsv, no header) mapped to UniProt.

File columns (no header, 8 fields): species, type, protein name, allergen name (may be blank), accession,
length, year, blank. Placeholder IDs (COMPARE00xxx) have no database accession and are reported, not mapped.
Accessions are GenBank-style protein IDs, mapped with the UniProt ID-mapping service (EMBL-GenBank-DDBJ_CDS).
The reviewed status and the year are kept. A protein already in analysis/allergen_model/positives/positives.tsv
is flagged, not added.

Outputs (in --outdir): compare_fungi.tsv, compare_new.faa (P_ IDs), compare_new.tsv, counts.md, unmapped.tsv.
Run: /usr/bin/python3.12 05_compare_fungi.py --compare <compare.tsv> --positives <positives.tsv> --outdir <dir>
"""

import argparse
import csv
import datetime
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "ao", os.path.join(HERE, "03_allergenonline_fungi.py")
)
ao = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ao)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compare", required=True)
    ap.add_argument("--positives", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    with open(args.compare, encoding="utf-8-sig", errors="replace", newline="") as fh:
        for f in csv.reader(fh, delimiter="\t"):
            if not f or not any(x.strip() for x in f):
                continue
            f = (f + [""] * 8)[:8]
            rows.append(
                {
                    "species": f[0],
                    "type": f[1],
                    "name": f[2],
                    "allergen": f[3].strip(),
                    "acc": f[4].strip(),
                    "length": f[5],
                    "year": f[6],
                }
            )
    print(f"COMPARE rows={len(rows)}")

    pos = {}
    with open(args.positives, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pos[r["accession"]] = r["evidence_level"]

    placeholders = [r for r in rows if r["acc"].startswith("COMPARE")]
    real = [r for r in rows if not r["acc"].startswith("COMPARE")]
    by_src = {}
    for r in real:
        src, qid = ao.classify(r["acc"])
        r["_src"], r["_qid"] = src, qid
        by_src.setdefault(src, set()).add(qid)
    mapping = {src: ao.idmap(src, sorted(ids)) for src, ids in by_src.items()}

    resolved = set()
    for r in real:
        hit = mapping.get(r["_src"], {}).get(r["_qid"])
        r["_uniprot"] = hit[0] if hit else ""
        r["_reviewed"] = ("yes" if hit and hit[1] else "no") if hit else ""
        if hit:
            resolved.add(hit[0])
    seqs = ao.fetch_sequences(sorted(resolved))
    new_ids = sorted(a for a in resolved if a not in pos and a in seqs)

    with open(os.path.join(args.outdir, "compare_fungi.tsv"), "w") as fh:
        fh.write(
            "species\ttype\tprotein_name\tallergen_name\tcompare_accession\tlength\tyear\tuniprot\treviewed\t"
            "in_positives\tnew\n"
        )
        for r in real:
            u = r["_uniprot"]
            fh.write(
                "\t".join(
                    [
                        r["species"],
                        r["type"],
                        r["name"],
                        r["allergen"],
                        r["acc"],
                        r["length"],
                        r["year"],
                        u,
                        r["_reviewed"],
                        "yes" if u in pos else "no",
                        "yes" if u in new_ids else "no",
                    ]
                )
                + "\n"
            )
        for r in placeholders:
            fh.write(
                "\t".join(
                    [
                        r["species"],
                        r["type"],
                        r["name"],
                        r["allergen"],
                        r["acc"],
                        r["length"],
                        r["year"],
                        "",
                        "",
                        "no",
                        "no",
                    ]
                )
                + "\n"
            )
    with open(os.path.join(args.outdir, "unmapped.tsv"), "w") as fh:
        fh.write("compare_accession\tspecies\treason\n")
        for r in real:
            if not r["_uniprot"]:
                fh.write(f"{r['acc']}\t{r['species']}\tno UniProt mapping\n")
        for r in placeholders:
            fh.write(f"{r['acc']}\t{r['species']}\tplaceholder ID, no sequence\n")
    with open(os.path.join(args.outdir, "compare_new.faa"), "w") as fh:
        for a in new_ids:
            fh.write(f">P_{a}\n")
            s = seqs[a]["seq"]
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")
    with open(os.path.join(args.outdir, "compare_new.tsv"), "w") as fh:
        fh.write("accession\torganism\tlength\tprotein_name\n")
        for a in new_ids:
            fh.write(f"{a}\t{seqs[a]['org']}\t{seqs[a]['length']}\t{seqs[a]['name']}\n")

    mapped = sum(1 for r in real if r["_uniprot"])
    in_pos = len({r["_uniprot"] for r in real if r["_uniprot"] in pos})
    today = datetime.date.today().isoformat()
    with open(os.path.join(args.outdir, "counts.md"), "w") as fh:
        fh.write("# COMPARE fungal rows mapped to UniProt\n\n")
        fh.write(
            f"- Run date: {today}\n- Source: to_import/allergen/compare.tsv (no header; {len(rows)} rows)\n"
        )
        fh.write(f"- Rows with a database accession: {len(real)}; mapped to UniProt: {mapped}\n")
        fh.write(f"- Placeholder COMPARE IDs (no sequence): {len(placeholders)}\n")
        fh.write(f"- Distinct UniProt proteins: {len(resolved)}; already in positives: {in_pos}\n")
        fh.write(f"- New proteins written: {len(new_ids)}\n")
        fh.write(
            "- Evidence: the file has no evidence column. COMPARE's own criteria require literature IgE "
            "binding (stated on the COMPARE site, not checked here). Treat as unverified until the source papers are checked.\n"
        )
    print(
        f"real={len(real)} mapped={mapped} unique={len(resolved)} in_positives={in_pos} new={len(new_ids)} placeholders={len(placeholders)}"
    )


if __name__ == "__main__":
    main()
