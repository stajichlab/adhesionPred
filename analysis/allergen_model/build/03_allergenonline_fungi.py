#!/usr/bin/env python3.12
"""Fungal subset of AllergenOnline (version 24, to_import/allergen/AllergenOnline.csv), mapped to UniProt.

Selection:
  - rows with Type "Aero Fungi" or "Food Fungi", and
  - rows of other types whose species name is a fungal genus or species (for example contact allergens of
    Malassezia, Candida, Trichophyton, and the Af293 structure entries).
Accessions are mapped to UniProtKB with the UniProt ID-mapping service, by accession type:
  PDB chains (for example 8GKP_D) -> PDB; RefSeq protein (XP_, NP_) -> RefSeq_Protein;
  UniProt accessions -> UniProtKB_AC-ID; other protein accessions -> EMBL-GenBank-DDBJ_CDS.
Sequences come from UniProt. Proteins already in analysis/allergen_model/positives/positives.tsv are reported,
not added again.

Outputs (in --outdir):
  allergenonline_fungi.tsv      one row per selected entry: source fields, mapped UniProt, reviewed status,
                                already_in_positives, new (yes/no)
  allergenonline_new.faa        FASTA for new UniProt proteins, IDs P_<accession> (same prefix rule as positives)
  allergenonline_new.tsv        metadata for the new proteins
  counts.md                     counts by type, mapping success, reviewed status, and overlap
  unmapped.tsv                  entries that did not map, with the reason
Run: /usr/bin/python3.12 03_allergenonline_fungi.py --csv <AllergenOnline.csv> --positives <positives.tsv> --outdir <dir>
"""

import argparse
import csv
import datetime
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "adhesionPred-allergenonline/1.0"}
FUNGAL = re.compile(
    r"Aspergill|Penicill|Alternaria|Cladosporium|Candida|Malassezia|Trichophyton|Arthroderma|Coccidioides|"
    r"Histoplasma|Fusarium|Saccharomyces|Mucor|Rhizopus|Botrytis|Curvularia|Epicoccum|Stachybotrys|"
    r"Aureobasidium|Paracoccid|Blastomyces|Ustilago|Puccinia|Hormodendrum|Drechslera|Helminthosporium|"
    r"Cryptococcus|Emericella|Neurospora|Pyricularia|Ascochyta|Sporobolomyces|Trichoderma|Verticillium|"
    r"Phoma|Wallemia|Acremonium|Chaetomium|Schizophyllum|Agaricus|Lentinula|Thermomyces|Pichia|Fungi|"
    r"Yeast|Fungus|mold|Cladosporium|Alternaria|Af293|Tritirachium|Mycosphaerella|Ophiostoma"
)
UNIPROT_AC = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9](\.\d+)?$|^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](\.\d+)?$"
)
PDB_CHAIN = re.compile(r"^[0-9][A-Za-z0-9]{3}_[A-Za-z0-9]+$")
REFSEQ = re.compile(r"^(XP|NP|WP|YP)_\d+(\.\d+)?$")


def classify(acc):
    a = acc.strip()
    if PDB_CHAIN.match(a):
        return "PDB", a.split("_")[0].upper()
    if REFSEQ.match(a):
        return "RefSeq_Protein", a
    if UNIPROT_AC.match(a):
        return "UniProtKB_AC-ID", a.split(".")[0]
    return "EMBL-GenBank-DDBJ_CDS", a.split(".")[0]


def http_json(url, data=None, timeout=120):
    body = urllib.parse.urlencode(data).encode() if data is not None else None
    req = urllib.request.Request(url, data=body, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def idmap(src, ids):
    """Run one ID-mapping job for a source database; return {query_id: primary_accession or None, ...} and reviewed flag."""
    out = {}
    for i in range(0, len(ids), 400):
        chunk = ids[i : i + 400]
        try:
            job = http_json(
                "https://rest.uniprot.org/idmapping/run",
                {"from": src, "to": "UniProtKB", "ids": ",".join(chunk)},
            )["jobId"]
        except urllib.error.HTTPError as e:
            body = e.read().decode("utf-8", "replace")[:300]
            sys.exit(
                f"FATAL: idmapping {src} batch starting {chunk[0]} (n={len(chunk)}) HTTP {e.code}: {body}"
            )
        for _ in range(120):
            st = http_json(f"https://rest.uniprot.org/idmapping/status/{job}")
            if st.get("jobStatus") in (None, "FINISHED"):
                break
            time.sleep(2)
        url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job}?size=500"
        while url:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as resp:
                d = json.loads(resp.read().decode("utf-8"))
                link = resp.headers.get("Link", "")
            for r in d.get("results", []):
                to = r.get("to")
                if isinstance(to, dict):
                    out.setdefault(
                        r["from"],
                        (to["primaryAccession"], "reviewed" in to.get("entryType", "").lower()),
                    )
                else:
                    out.setdefault(r["from"], None)
            m = re.search(r'<([^>]+)>;\s*rel="next"', link)
            url = m.group(1) if m else None
    return out


def fetch_sequences(acc_list):
    """Sequences and lengths from UniProt for the mapped accessions (one stream per 40)."""
    seqs = {}
    for i in range(0, len(acc_list), 40):
        batch = acc_list[i : i + 40]
        query = " OR ".join(f"accession:{a}" for a in batch)
        url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
            {
                "query": query,
                "format": "tsv",
                "fields": "accession,organism_name,length,sequence,protein_name",
            }
        )
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=300) as resp:
            lines = resp.read().decode("utf-8").rstrip("\n").split("\n")
        if len(lines) < 2:
            continue
        for ln in lines[1:]:
            f = ln.split("\t")
            seqs[f[0]] = {"org": f[1], "length": f[2], "seq": f[3].strip(), "name": f[4]}
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--positives", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument(
        "--all-kingdoms",
        action="store_true",
        help="select all rows, not only fungal (owner decision Q1)",
    )
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.csv, encoding="utf-8-sig", errors="replace", newline="") as fh:
        rows = list(csv.DictReader(fh))
    rows = [{k.replace("\xa0", " ").strip(): v for k, v in r.items()} for r in rows]
    if args.all_kingdoms:
        sel = rows
    else:
        sel = [
            r
            for r in rows
            if r["Type"] in ("Aero Fungi", "Food Fungi")
            or FUNGAL.search(r["Species"] + " " + r["Common"])
        ]
    print(f"AllergenOnline rows={len(rows)} fungal selected={len(sel)}")

    pos = {}
    with open(args.positives, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pos[r["accession"]] = r["evidence_level"]

    # Map each accession, grouped by source database.
    by_src = {}
    for r in sel:
        src, qid = classify(r["Accession"])
        by_src.setdefault(src, set()).add(qid)
    mapping = {}
    for src, ids in by_src.items():
        print(f"mapping {len(ids)} IDs from {src}")
        mapping[src] = idmap(src, sorted(ids))

    resolved = {}
    for r in sel:
        src, qid = classify(r["Accession"])
        hit = mapping.get(src, {}).get(qid)
        r["_src"], r["_qid"] = src, qid
        r["_uniprot"] = hit[0] if hit else ""
        r["_reviewed"] = ("yes" if hit and hit[1] else "no") if hit else ""
        if hit:
            resolved.setdefault(hit[0], None)
    seqs = fetch_sequences(sorted(resolved))
    print(f"UniProt sequences fetched: {len(seqs)} of {len(resolved)} mapped accessions")

    new_ids = [a for a in sorted(resolved) if a not in pos and a in seqs]
    with open(os.path.join(args.outdir, "allergenonline_fungi.tsv"), "w") as fh:
        fh.write(
            "species\tcommon\tiuis_allergen\ttype\tallergenicity\tsource_accession\tsource_db\tuniprot\treviewed\tin_positives\tnew\n"
        )
        for r in sel:
            u = r["_uniprot"]
            fh.write(
                "\t".join(
                    [
                        r["Species"],
                        r["Common"],
                        r["IUIS\xa0Allergen".replace("\xa0", " ")]
                        if "IUIS\xa0Allergen" in r
                        else r.get("IUIS Allergen", ""),
                        r["Type"],
                        r["Allergenicity"],
                        r["Accession"],
                        r["_src"],
                        u,
                        r["_reviewed"],
                        "yes" if u in pos else "no",
                        "yes" if (u in new_ids) else "no",
                    ]
                )
                + "\n"
            )
    unm = [r for r in sel if not r["_uniprot"]]
    with open(os.path.join(args.outdir, "unmapped.tsv"), "w") as fh:
        fh.write("source_accession\tspecies\tsource_db\treason\n")
        for r in unm:
            fh.write(f"{r['Accession']}\t{r['Species']}\t{r['_src']}\tno UniProt mapping\n")
    with open(os.path.join(args.outdir, "allergenonline_new.faa"), "w") as fh:
        for a in new_ids:
            fh.write(f">P_{a}\n")
            s = seqs[a]["seq"]
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")
    with open(os.path.join(args.outdir, "allergenonline_new.tsv"), "w") as fh:
        fh.write("accession\torganism\tlength\tprotein_name\n")
        for a in new_ids:
            fh.write(f"{a}\t{seqs[a]['org']}\t{seqs[a]['length']}\t{seqs[a]['name']}\n")

    types = {}
    for r in sel:
        types[r["Type"]] = types.get(r["Type"], 0) + 1
    mapped = sum(1 for r in sel if r["_uniprot"])
    overlap = sum(1 for r in sel if r["_uniprot"] and r["_uniprot"] in pos)
    today = datetime.date.today().isoformat()
    with open(os.path.join(args.outdir, "counts.md"), "w") as fh:
        fh.write("# AllergenOnline fungal subset (version 24, mapped to UniProt)\n\n")
        fh.write(
            f"- Run date: {today}\n- Source: to_import/allergen/AllergenOnline.csv ({len(rows)} rows)\n"
        )
        fh.write(
            f"- Fungal rows selected: {len(sel)}\n- Mapped to UniProt: {mapped}\n- Unmapped: {len(sel) - mapped}\n"
        )
        fh.write(
            f"- Distinct UniProt accessions: {len(resolved)}\n- Sequences fetched: {len(seqs)}\n"
        )
        fh.write(
            f"- Already in positives (147 set): {overlap}\n- New proteins written: {len(new_ids)}\n\n"
        )
        fh.write("| Type | Rows |\n|---|---|\n")
        for t, n in sorted(types.items(), key=lambda x: -x[1]):
            fh.write(f"| {t} | {n} |\n")
    print(
        f"selected={len(sel)} mapped={mapped} unique_uniprot={len(resolved)} in_positives={overlap} new={len(new_ids)}"
    )


if __name__ == "__main__":
    main()
