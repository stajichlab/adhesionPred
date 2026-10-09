# ruff: noqa
"""Assemble curation_evidence.tsv for the 61 unlabelled hydrophobin calls (evidence only, no decisions).
Usage: assemble.py <hydrophobin_truth dir> <work dir> <upjson dir> <literature_gene_map.tsv or '-'> <out tsv>
Reads: evidence_sheets.tsv, run_list.tsv, positives_t2t3.faa (174), work/lineA.tsv, work/q61_vs_pos.tsv,
work/pos_vs_<proteome>.tsv, work/sp61.tsv, work/pfam61.domtbl, work/up/hits_<proteome>.tsv, work/fdb/hits_<proteome>.tsv.
"""

import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

D, W, UJ, LIT, OUT = (Path(x) if x != "-" else None for x in sys.argv[1:6])
HYD_PFAM = {"PF01185", "PF06766", "PF28987", "PF22354", "PF29785", "PF29802", "PF29465"}
FMT = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle".split()


def read_fasta(p):
    seqs, name = {}, None
    for line in open(p):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    return seqs


def read_blast(p):
    hits = defaultdict(lambda: defaultdict(list))  # q -> s -> [hsp dicts]
    order = defaultdict(list)
    if not Path(p).exists():
        return hits, order
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        h = dict(zip(FMT, f))
        for k in ("pident", "bitscore"):
            h[k] = float(h[k])
        for k in ("length", "qlen", "slen", "qstart", "qend", "sstart", "send"):
            h[k] = int(h[k])
        if h["sseqid"] not in hits[h["qseqid"]]:
            order[h["qseqid"]].append(h["sseqid"])
        hits[h["qseqid"]][h["sseqid"]].append(h)
    return hits, order


def merged_cov(ivs, L):
    ivs = sorted((min(a, b), max(a, b)) for a, b in ivs)
    tot, cs, ce = 0, None, None
    for a, b in ivs:
        if cs is None or a > ce + 1:
            if cs is not None:
                tot += ce - cs + 1
            cs, ce = a, b
        else:
            ce = max(ce, b)
    if cs is not None:
        tot += ce - cs + 1
    return 100.0 * tot / L


def pair_stats(hsps):
    """identity of the top HSP; merged coverage of query and subject over all HSPs of the pair; total bitscore."""
    top = max(hsps, key=lambda h: h["bitscore"])
    qc = merged_cov([(h["qstart"], h["qend"]) for h in hsps], top["qlen"])
    sc = merged_cov([(h["sstart"], h["send"]) for h in hsps], top["slen"])
    return top["pident"], qc, sc, top["bitscore"], top["evalue"], top


def best(hitd):
    if not hitd:
        return None, None
    s = max(hitd, key=lambda s: max(h["bitscore"] for h in hitd[s]))
    return s, hitd[s]


sheet = list(csv.DictReader(open(D / "evidence_sheets.tsv"), delimiter="\t"))
lineA = {r["key"]: r for r in csv.DictReader(open(W / "lineA.tsv"), delimiter="\t")}
pos_names = set(read_fasta(D / "positives_t2t3.faa"))
truth = {r["accession"]: r for r in csv.DictReader(open(D / "truth_all.tsv"), delimiter="\t")}
q2pos, _ = read_blast(W / "q61_vs_pos.tsv")
sp, _ = read_blast(W / "sp61.tsv")

# Pfam domtbl
pfam = defaultdict(list)
for line in open(W / "pfam61.domtbl"):
    if line.startswith("#"):
        continue
    f = line.split()
    pfam[f[3]].append(
        dict(
            name=f[0],
            acc=f[1].split(".")[0],
            qlen=int(f[5]),
            ievalue=f[12],
            score=f[13],
            env=(int(f[19]), int(f[20])),
        )
    )

# literature map (optional)
lit = defaultdict(list)
if LIT and LIT.exists():
    for r in csv.DictReader(open(LIT), delimiter="\t"):
        if r.get("proteome_id"):
            lit[(r["proteome"], r["proteome_id"])].append(r)


# UniProt JSON
def upinfo(acc):
    p = UJ / f"{acc}.json"
    if not p.exists():
        return None
    d = json.load(open(p))
    pd = d.get("proteinDescription", {})
    nm = pd.get("recommendedName") or (pd.get("submissionNames") or [{}])[0]
    full = nm.get("fullName", {})
    ev = [
        f"{e.get('evidenceCode', '')}"
        + (f" {e.get('source', '')}:{e.get('id', '')}" if e.get("source") else "")
        for e in full.get("evidences", [])
    ]
    genes = []
    for g in d.get("genes", []):
        genes += [g.get("geneName", {}).get("value", "")] + [
            o["value"] for o in g.get("orderedLocusNames", []) + g.get("orfNames", [])
        ]
    pmids = []
    for c in d.get("references", []):
        pmids += [
            x["id"]
            for x in c["citation"].get("citationCrossReferences", [])
            if x["database"] == "PubMed"
        ]
    exp = any(
        e.get("evidenceCode") == "ECO:0000269"
        for c in d.get("comments", [])
        for t in c.get("texts", [])
        for e in t.get("evidences", [])
    )
    return dict(
        acc=d["primaryAccession"],
        reviewed="Swiss-Prot" in d["entryType"] and "unreviewed" not in d["entryType"],
        entrytype=d["entryType"],
        name=full.get("value", ""),
        name_ev=";".join(ev) or "none given",
        genes=",".join(x for x in genes if x),
        pmids=",".join(pmids),
        pe=d.get("proteinExistence", ""),
        exp_comment=exp,
    )


rows = []
for r in sheet:
    p, i = r["proteome"], r["id"]
    key = f"{p}__{i}"
    o = {"proteome": p, "id": i, "kind": r["kind"], "length": r["length"]}
    # line A
    o["line_A"] = lineA[key]["line_A"]
    o["A_detail"] = lineA[key]["A_detail"]
    # line B
    o["line_B"] = "yes" if r["strict_families"].strip() else "no"
    o["B_detail"] = f"strict_families={r['strict_families'] or 'none'} (from evidence_sheets.tsv)"
    # line D: reciprocal best hit vs 174 known hydrophobins
    s, hs = best(q2pos.get(key, {}))
    if s is None:
        o["line_D"] = "no"
        o["D_detail"] = "no BLASTP hit to the 174 known hydrophobins at E<=1e-3"
    else:
        pid, qc, sc, bits, ev, _ = pair_stats(hs)
        rev, _o = read_blast(W / f"pos_vs_{p}.tsv")
        rs, rhs = best(rev.get(s, {}))
        recip = rs == i
        tier = truth.get(s, {}).get("tier", "?")
        sp_of = truth.get(s, {}).get("species", "?")
        ok = recip and pid >= 40 and qc >= 70 and sc >= 70
        o["line_D"] = "yes" if ok else "no"
        o["D_detail"] = (
            f"best hit {s} ({truth.get(s, {}).get('entry_name', '?')}, tier {tier}, {sp_of}): identity {pid:.1f}% (top HSP), "
            f"coverage query {qc:.0f}% / subject {sc:.0f}% (merged HSPs), bits {bits:.0f}, E {ev}; "
            f"reverse best hit of {s} in {p}: {rs or 'none'} -> {'reciprocal' if recip else 'not reciprocal'}; "
            f"one-to-one ortholog clause not checked"
        )
    # line E: Swiss-Prot + Pfam-A
    es = []
    E_yes = False
    s, hs = best(sp.get(key, {}))
    if s is None:
        es.append("Swiss-Prot 2023_03: no hit at E<=1e-3")
    else:
        pid, qc, sc, bits, ev, top = pair_stats(hs)
        acc = s.split("|")[1] if "|" in s else s
        title = top["stitle"]
        hydro_annot = acc in pos_names or re.search(r"hydrophob|rodlet", title, re.I) is not None
        uncharacterized = re.search(r"uncharacteri[sz]ed", title, re.I) is not None
        strong = pid >= 40 and qc >= 70
        es.append(
            f"Swiss-Prot best: {title[:140]} identity {pid:.1f}% query cov {qc:.0f}% subject cov {sc:.0f}% E {ev}; "
            f"annotation={'hydrophobin' if hydro_annot else ('uncharacterized' if uncharacterized else 'non-hydrophobin')}; "
            f"strong(>=40%id,>=70%qcov)={'yes' if strong else 'no'}"
        )
        if strong and not hydro_annot and not uncharacterized:
            E_yes = True
    doms = pfam.get(key, [])
    L = int(r["length"])
    if not doms:
        es.append("Pfam-A 38.2 --cut_ga: no domain")
    else:
        parts = []
        for dd in sorted(doms, key=lambda x: x["env"][0]):
            c = 100 * (dd["env"][1] - dd["env"][0] + 1) / L
            parts.append(f"{dd['name']}({dd['acc']}) env {dd['env'][0]}-{dd['env'][1]} {c:.0f}%")
        nonh = [dd for dd in doms if dd["acc"] not in HYD_PFAM]
        nonh_cov = merged_cov([dd["env"] for dd in nonh], L) if nonh else 0.0
        es.append(
            f"Pfam-A 38.2 --cut_ga: {'; '.join(parts)}; non-hydrophobin domain coverage {nonh_cov:.0f}% of length"
        )
        if nonh_cov > 50:
            E_yes = True
    o["line_E"] = "yes" if E_yes else "no"
    o["E_detail"] = " || ".join(es)
    o["other_pfam"] = (
        ",".join(sorted({dd["name"] for dd in doms if dd["acc"] not in HYD_PFAM})) or ""
    )
    # UniProt mapping
    if p in ("Afum_Af293", "Afum_A1163") and i.startswith(("tr|", "sp|")):
        acc = i.split("|")[1]
        how = "proteome FASTA is UniProt; accession taken from the ID"
        ident = qc = sc = 100.0
    else:
        uh, _ = read_blast(W / f"up/hits_{p}.tsv")
        s, hs = best(uh.get(key, {}))
        acc = None
        how = "BLASTP vs UniProtKB 2026_03 entries of the taxon; no hit"
        if s:
            ident, qc, sc, _b, _e, top = pair_stats(hs)
            ox = re.search(r"OX=(\d+)", top["stitle"])
            how = f"BLASTP vs UniProtKB 2026_03 (taxon set for {p}); best {s.split('|')[1]} OX={ox.group(1) if ox else '?'} identity {ident:.1f}% query cov {qc:.0f}% subject cov {sc:.0f}%"
            if ident >= 95 and qc >= 90 and sc >= 90:
                acc = s.split("|")[1]
                how += " -> mapped (>=95% id, >=90% coverage of both)"
            elif ident >= 95 and (qc >= 90 or sc >= 90):
                how += f" -> PARTIAL (>=95% id but <90% coverage of one sequence; not mapped; candidate entry {s.split('|')[1]})"
            else:
                how += " -> not mapped (below >=95% id / >=90% coverage of both)"
    info = upinfo(acc) if acc else None
    o["uniprot_accession"] = acc or ""
    o["uniprot_name"] = info["name"] if info else ""
    o["uniprot_entry"] = (
        f"{info['entrytype']}; name evidence {info['name_ev']}; genes {info['genes']}; PubMed {info['pmids'] or 'none'}; "
        f"PE {info['pe']}; experimental (ECO:0000269) comment={'yes' if info['exp_comment'] else 'no'}"
        if info
        else ""
    )
    o["uniprot_mapping"] = how
    # FungiDB
    fh, _ = read_blast(W / f"fdb/hits_{p}.tsv")
    s, hs = best(fh.get(key, {}))
    if not Path(W / f"fdb/hits_{p}.tsv").exists():
        o["fungidb_product"] = ""
        o["fungidb_mapping"] = "no FungiDB release-31 proteome for this species"
    elif s is None:
        o["fungidb_product"] = ""
        o["fungidb_mapping"] = "FungiDB release-31: no hit"
    else:
        ident, qc, sc, _b, _e, top = pair_stats(hs)
        prod = re.search(r"gene_product=([^|]+)", top["stitle"])
        org = re.search(r"organism=([^|]+)", top["stitle"])
        gene = re.search(r"gene=([^|\s]+)", top["stitle"])
        mapped = ident >= 95 and qc >= 90 and sc >= 90
        partial = (not mapped) and ident >= 95 and (qc >= 90 or sc >= 90)
        o["fungidb_product"] = prod.group(1).strip() if (prod and mapped) else ""
        o["fungidb_mapping"] = (
            f"FungiDB-31 {org.group(1).strip() if org else '?'} {gene.group(1) if gene else s}: identity {ident:.1f}% query cov {qc:.0f}% subject cov {sc:.0f}%"
            + (
                " -> mapped (>=95% id, >=90% coverage of both)"
                if mapped
                else (
                    " -> PARTIAL (>=95% id, <90% coverage of one; not mapped)"
                    if partial
                    else " -> not mapped"
                )
            )
            + (
                ""
                if mapped
                else f" (product of best hit: {prod.group(1).strip() if prod else '?'})"
            )
        )
    # line C
    csrc = []
    C = "no"
    if info:
        hyd = re.search(r"hydrophob|rodlet", info["name"], re.I)
        if hyd and (info["reviewed"] or "PubMed" in info["name_ev"] or info["exp_comment"]):
            C = "yes"
            csrc.append(
                f"UniProt {info['acc']} '{info['name']}' ({info['entrytype']}; name evidence {info['name_ev']})"
            )
        elif hyd:
            csrc.append(
                f"UniProt {info['acc']} '{info['name']}' (unreviewed; name evidence {info['name_ev']}; no PubMed name evidence, no ECO:0000269 comment): does not count"
            )
    if o["fungidb_product"] and re.search(r"hydrophob|rodlet", o["fungidb_product"], re.I):
        csrc.append(
            f"FungiDB-31 gene_product '{o['fungidb_product']}' (basis of the product name not checked)"
        )
        C = "yes"
    for lr in lit.get((p, i), []):
        if lr["mapped"] == "yes":
            csrc.append(
                f"paper {lr['paper']} gene {lr['paper_gene_id']} (identity {lr['identity']}%, {lr['coverage']}; {lr['how_obtained']})"
            )
            C = "yes"
        elif lr["mapped"] == "partial":
            csrc.append(
                f"PARTIAL match, does not meet the mapping rule: paper {lr['paper']} gene {lr['paper_gene_id']} identity {lr['identity']}%, {lr['coverage']}"
            )
        else:
            csrc.append(
                f"best paper gene hit below rule: {lr['paper']} {lr['paper_gene_id']} identity {lr['identity']}%, {lr['coverage']}"
            )
    if C == "no":
        csrc.insert(
            0,
            "no mapped UniProt, FungiDB or paper record names it a hydrophobin"
            + ("" if LIT else "; literature not yet checked"),
        )
    o["line_C"] = C if LIT else ("yes" if C == "yes" else "not checked")
    o["C_source"] = " | ".join(csrc)
    rows.append(o)

cols = [
    "proteome",
    "id",
    "kind",
    "length",
    "line_A",
    "A_detail",
    "line_B",
    "B_detail",
    "line_C",
    "C_source",
    "line_D",
    "D_detail",
    "line_E",
    "E_detail",
    "other_pfam",
    "uniprot_accession",
    "uniprot_name",
    "uniprot_entry",
    "uniprot_mapping",
    "fungidb_product",
    "fungidb_mapping",
]
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
from collections import Counter

for c in ("line_A", "line_B", "line_C", "line_D", "line_E"):
    print(c, dict(Counter(r[c] for r in rows)))
print(
    "mapped to UniProt:",
    sum(1 for r in rows if r["uniprot_accession"]),
    "; other Pfam domain:",
    sum(1 for r in rows if r["other_pfam"]),
)
