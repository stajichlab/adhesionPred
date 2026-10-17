# ruff: noqa
"""Build controls table draft from candidates.tsv. UniProt JSON (cached), NCBI taxonomy dumps.

Usage: build.py CANDIDATES OUTDIR
Writes OUTDIR/controls.draft.tsv, OUTDIR/sequences.faa, OUTDIR/uniprot_features.tsv
"""

import csv
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "..", "dl", "json")
os.makedirs(CACHE, exist_ok=True)
TAX = "/srv/projects/db/taxonomy"
FAMILY_PFAM = {}  # filled from family_table
REPO = "/bigdata/stajichlab/jstajich/projects/adhesionPred"
for r in csv.DictReader(open(f"{REPO}/data/sorting_hat/family_table.tsv"), delimiter="\t"):
    FAMILY_PFAM[r["pfam_acc"]] = r["name"]

ADH = re.compile(r"adhe|attach|adher|agglutin|flocc", re.I)
CLADES = [
    (147537, "Saccharomycotina"),
    (451866, "Taphrinomycotina"),
    (33183, "Onygenales"),
    (5042, "Eurotiales"),
    (147538, "Pezizomycotina_other"),
    (5204, "Basidiomycota"),
]


def load_tax():
    parent, rank, name = {}, {}, {}
    with open(f"{TAX}/nodes.dmp") as f:
        for line in f:
            p = line.split("\t|\t")
            parent[int(p[0])] = int(p[1])
            rank[int(p[0])] = p[2]
    with open(f"{TAX}/names.dmp") as f:
        for line in f:
            p = line.rstrip("\t|\n").split("\t|\t")
            if p[3] == "scientific name":
                name[int(p[0])] = p[1]
    return parent, rank, name


def species_of(t, parent, rank):
    seen = t
    while t in parent and t != 1:
        if rank[t] == "species":
            return t
        t = parent[t]
    raise ValueError(f"no species above {seen}")


def lineage(t, parent):
    out = []
    while t in parent and t != 1:
        out.append(t)
        t = parent[t]
    return out


def clade_of(t, parent):
    lin = set(lineage(t, parent))
    for tid, nm in CLADES:
        if tid in lin:
            return nm
    return "other_fungi" if 4751 in lin else "non_fungal"


def get_json(acc):
    fn = f"{CACHE}/{acc}.json"
    if not os.path.exists(fn):
        x = subprocess.run(
            ["curl", "-s", "-m", "60", f"https://rest.uniprot.org/uniprotkb/{acc}.json"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        open(fn, "w").write(x)
        time.sleep(0.2)
    return json.load(open(fn))


def ncbi_seq(acc):
    fn = f"{CACHE}/{acc}.ncbi.fa"
    if not os.path.exists(fn):
        x = subprocess.run(
            [
                "curl",
                "-s",
                "-m",
                "60",
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={acc}&rettype=fasta&retmode=text",
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        open(fn, "w").write(x)
    lines = open(fn).read().splitlines()
    return lines[0], "".join(lines[1:])


def ev_str(evs):
    out = []
    for e in evs or []:
        s = e.get("evidenceCode", "")
        if e.get("source"):
            s += f"|{e['source']}:{e.get('id', '')}"
        out.append(s)
    return out


def main():
    cand, outdir = sys.argv[1], sys.argv[2]
    parent, rank, name = load_tax()
    rows, feats = [], []
    seqs = {}
    for c in csv.DictReader(open(cand), delimiter="\t"):
        acc = c["accession"]
        r = {
            "accession": acc,
            "stratum": c["stratum"],
            "origin": c["origin"],
            "manual_pmid": c.get("manual_pmid") or "",
            "manual_quote": c.get("manual_quote") or "",
            "cand_notes": c.get("notes") or "",
        }
        if re.match(r"^[A-Z]{2,3}_?\d+(\.\d+)?$", acc) and not re.match(
            r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", acc
        ):
            hdr, s = ncbi_seq(acc)
            seqs[acc] = s
            m = re.search(r"\[(.+)\]", hdr)
            org = m.group(1) if m else ""
            tid = {"Blastomyces dermatitidis ER-3": 559297}.get(org)
            r.update(
                gene=hdr.split()[1] if len(hdr.split()) > 1 else "",
                reviewed="ncbi",
                org_taxid=tid,
                func="",
                func_exp_pmids="",
                func_ev="",
                pfam="",
                gpi="",
                signal="",
                tm="",
                st_feature="",
                repeat_feature="",
                subcell="",
                go_adh="",
                length=len(s),
                protein=hdr[:120],
                cat="",
                cat_exp_pmids="",
            )
        else:
            j = get_json(acc)
            if "primaryAccession" not in j:
                print("MISSING", acc, file=sys.stderr)
                continue
            if j["primaryAccession"] != acc:
                print("SECONDARY", acc, "->", j["primaryAccession"], file=sys.stderr)
            seqs[acc] = j["sequence"]["value"]
            genes = j.get("genes", [{}])
            g = genes[0] if genes else {}
            gene = (
                g.get("geneName", {}).get("value")
                or (g.get("orderedLocusNames") or [{}])[0].get("value")
                or (g.get("orfNames") or [{}])[0].get("value")
                or ""
            )
            pname = j.get("proteinDescription", {}).get("recommendedName", {}).get(
                "fullName", {}
            ).get("value") or (j.get("proteinDescription", {}).get("submissionNames") or [{}])[
                0
            ].get("fullName", {}).get("value", "")
            funcs, fev = [], []
            for cm in j.get("comments", []):
                if cm.get("commentType") == "FUNCTION":
                    for t in cm.get("texts", []):
                        funcs.append(t["value"])
                        fev += ev_str(t.get("evidences"))
            subc = []
            for cm in j.get("comments", []):
                if cm.get("commentType") == "SUBCELLULAR LOCATION":
                    for sl in cm.get("subcellularLocations", []):
                        loc = sl.get("location", {})
                        subc.append(
                            loc.get("value", "")
                            + "{"
                            + ",".join(e.split("|")[0] for e in ev_str(loc.get("evidences")))
                            + "}"
                        )
            exp_pm = sorted(
                {e.split("PubMed:")[1] for e in fev if e.startswith("ECO:0000269|PubMed")}
            )
            cat, cev = [], []
            for cm in j.get("comments", []):
                if cm.get("commentType") == "CATALYTIC ACTIVITY":
                    rx = cm.get("reaction", {})
                    cat.append(rx.get("name", ""))
                    cev += ev_str(rx.get("evidences"))
            cat_pm = sorted(
                {e.split("PubMed:")[1] for e in cev if e.startswith("ECO:0000269|PubMed")}
            )
            pf = [x["id"] for x in j.get("uniProtKBCrossReferences", []) if x["database"] == "Pfam"]
            go_adh = [
                x["id"]
                + ":"
                + next((p["value"] for p in x.get("properties", []) if p["key"] == "GoTerm"), "")
                for x in j.get("uniProtKBCrossReferences", [])
                if x["database"] == "GO"
                and ADH.search(
                    next((p["value"] for p in x.get("properties", []) if p["key"] == "GoTerm"), "")
                )
            ]
            ft = j.get("features", [])

            def fts(typ, pat=None):
                out = []
                for f in ft:
                    if f["type"] == typ and (
                        pat is None or re.search(pat, f.get("description", ""), re.I)
                    ):
                        out.append(
                            f"{f['location']['start']['value']}-{f['location']['end']['value']}:{f.get('description','')}"
                        )
                return out

            gpi = fts("Lipidation", "GPI")
            st = fts("Compositional bias", r"ser|thr") + fts("Region", r"linker")
            rep = fts("Repeat")
            r.update(
                gene=gene,
                protein=pname,
                reviewed=j.get("entryType", "")[:20],
                org_taxid=j["organism"]["taxonId"],
                func=" ".join(funcs),
                func_exp_pmids=";".join(exp_pm),
                func_ev=";".join(sorted(set(e.split("|")[0] for e in fev))),
                pfam=";".join(pf),
                gpi=";".join(gpi),
                signal=";".join(fts("Signal")),
                tm=str(len(fts("Transmembrane"))),
                st_feature=";".join(st),
                repeat_feature=str(len(rep)),
                subcell=";".join(subc),
                go_adh=";".join(go_adh),
                length=len(seqs[acc]),
                cat=" | ".join(cat)[:300],
                cat_exp_pmids=";".join(cat_pm),
            )
        sq = seqs[acc]
        best = 0.0
        for i in range(0, max(1, len(sq) - 29)):
            w = sq[i : i + 30]
            best = max(best, sum(w.count(x) for x in "ST") / 30)
        r["max_ST30"] = f"{best:.2f}"
        sp = species_of(r["org_taxid"], parent, rank)
        r["taxon_id"] = sp
        r["species"] = name[sp]
        r["strain_taxid"] = r["org_taxid"]
        r["clade"] = clade_of(sp, parent)
        r["family_pfam"] = ";".join(
            f"{p}({FAMILY_PFAM[p]})" for p in r["pfam"].split(";") if p in FAMILY_PFAM
        )
        r["adh_text_flag"] = "yes" if ADH.search(r["func"]) or r["go_adh"] else ""
        rows.append(r)
    keys = list(rows[0].keys())
    with open(f"{outdir}/candidates_annotated.tsv", "w") as f:
        w = csv.DictWriter(f, keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with open(f"{outdir}/sequences.faa", "w") as f:
        for r in rows:
            f.write(f">{r['accession']}\n{seqs[r['accession']]}\n")
    print(len(rows), "rows")


main()
