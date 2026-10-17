# ruff: noqa
"""Assemble controls.tsv, clusters.tsv and counts from the annotated candidates and MMseqs2 output.

Usage: assemble.py WORKDIR OUTDIR REJECT_ACCS(comma) LEVEL_OVERRIDES(acc=N2,...)
"""

import collections
import csv
import sys

W, OUT = sys.argv[1], sys.argv[2]
DROP = set(sys.argv[3].split(",")) if len(sys.argv) > 3 and sys.argv[3] else set()
LEVEL = (
    dict(x.split("=") for x in sys.argv[4].split(",")) if len(sys.argv) > 4 and sys.argv[4] else {}
)
REPO = "/bigdata/stajichlab/jstajich/projects/adhesionPred"
DATE = "2026-10-08"
CURATOR = "claude-opus-5-5 (task 07 agent)"

fam = {
    r["pfam_acc"]: r
    for r in csv.DictReader(open(f"{REPO}/data/sorting_hat/family_table.tsv"), delimiter="\t")
}
ADH_FAM = {k for k, v in fam.items() if v["module"] == "pfam_adhesion"}

rows = [
    r
    for r in csv.DictReader(open(f"{W}/candidates_annotated.tsv"), delimiter="\t")
    if r["accession"] not in DROP
]

# clusters (combined with positives)
members = collections.defaultdict(list)
rep_of = {}
for line in open(f"{W}/clu_cluster.tsv"):
    rep, m = line.rstrip("\n").split("\t")
    members[rep].append(m)
    rep_of[m] = rep
cid = {rep: f"HN{i:04d}" for i, rep in enumerate(sorted(members), 1)}
mixed = {
    rep
    for rep, ms in members.items()
    if any(m.startswith("pos|") for m in ms)
    and any(m.startswith("neg|") and m[4:] not in DROP for m in ms)
}

# seed membership and seed homology
seed_members = collections.defaultdict(set)
for line in open(f"{W}/../dl/seed/seed_members.tsv"):
    pf, acc, _ = line.rstrip("\n").split("\t")
    seed_members[acc].add(pf)
seed_hits = collections.defaultdict(set)
for line in open(f"{W}/seed_hits.m8"):
    q, t = line.split("\t")[:2]
    pid = float(line.split("\t")[2])
    seed_hits[q].add((t.split("|")[0], pid))
sow = {line.split("\t")[0] for line in open(f"{W}/sow_hits.m8")}

ORDER = ["none", "unknown", "partial", "tuned_on_truth", "in_reference"]
out = []
for r in rows:
    a = r["accession"]
    # evidence and quote
    if r["manual_pmid"]:
        src = f"PMID:{r['manual_pmid']}"
        if r["reviewed"] != "ncbi":
            src += f"; UniProtKB:{a}"
        else:
            src += f"; NCBI Protein:{a}"
        quote = r["manual_quote"]
        level = "N1"
    elif r["func_exp_pmids"]:
        src = (
            f"UniProtKB:{a} FUNCTION (ECO:0000269; PMID:{','.join(r['func_exp_pmids'].split(';'))})"
        )
        quote = r["func"]
        level = "N1"
    elif r["cat_exp_pmids"]:
        src = f"UniProtKB:{a} CATALYTIC ACTIVITY (ECO:0000269; PMID:{','.join(r['cat_exp_pmids'].split(';'))})"
        quote = "CATALYTIC ACTIVITY: " + r["cat"]
        level = "N1"
    else:
        ev = r["func_ev"] or (
            "FUNCTION without evidence code" if r["func"] else "no FUNCTION comment"
        )
        src = f"UniProtKB:{a} (FUNCTION evidence: {ev}; SUBCELLULAR LOCATION)"
        quote = (r["func"] + " " if r["func"] else "") + (
            "SUBCELLULAR LOCATION: " + r["subcell"] if r["subcell"] else ""
        )
        level = "N2"
    level = LEVEL.get(a, level)
    quote = " ".join(quote.split())
    if len(quote) > 400:
        quote = quote[:397] + "..."
    # serves
    pf = [p for p in r["pfam"].split(";") if p]
    adh_dom = [p for p in pf if p in ADH_FAM]
    other_dom = [p for p in pf if p in fam and p not in ADH_FAM]
    serves = []
    if r["stratum"] in ("secretory_non_surface", "mucin_sensor"):
        serves.append("step1")
    nrep = int(r["repeat_feature"] or 0)
    if (
        nrep >= 2
        or float(r["max_ST30"]) >= 0.5
        or r["stratum"] in ("repeat_non_adhesin", "st_linker_enzyme", "mucin_sensor")
    ):
        serves.append("repeat")
    if adh_dom:
        serves.append("family_domain")
    if r["stratum"] not in ("secretory_non_surface",):
        serves.append("adhesion_level")
    # leakage per module
    lk_step1 = "none"
    lk_repeat = "partial" if a in sow else "none"
    exact = seed_members.get(a, set()) & ADH_FAM
    hom = {p for p, _ in seed_hits.get(a, set()) if p in ADH_FAM}
    lk_fd = "in_reference" if exact else ("partial" if hom else "none")
    served_lk = [
        lk_step1 if "step1" in serves else "none",
        lk_repeat if "repeat" in serves else "none",
        lk_fd if "family_domain" in serves else "none",
    ]
    leak = max(served_lk, key=ORDER.index)
    notes = []
    if r["cand_notes"]:
        notes.append(r["cand_notes"])
    if r["reviewed"] != "ncbi" and not r["reviewed"].startswith("UniProtKB reviewed"):
        notes.append("UniProt entry is unreviewed (TrEMBL).")
    rep = rep_of.get("neg|" + a)
    if rep in mixed:
        pos = [m[4:] for m in members[rep] if m.startswith("pos|")]
        notes.append(
            f"FLAG: same 30% cluster as adhesin positive(s) {','.join(pos)}; drop or handle at cluster level."
        )
    other_hits = sorted({p for p, _ in seed_hits.get(a, set()) if p not in ADH_FAM})
    if other_hits:
        notes.append(
            "Homolog of Pfam seed of other families: "
            + ",".join(f"{p}({fam[p]['name']},{fam[p]['module']})" for p in other_hits)
            + "."
        )
    if other_dom:
        notes.append(
            "Has Pfam "
            + ",".join(f"{p}({fam[p]['name']},{fam[p]['module']})" for p in other_dom)
            + "."
        )
    if r["species"] != r.get("org_name", r["species"]):
        pass
    out.append(
        {
            "accession": a,
            "gene": r["gene"],
            "species": r["species"],
            "taxon_id": r["taxon_id"],
            "label": "negative",
            "stratum": r["stratum"],
            "evidence_level": level,
            "source": src,
            "quote": quote,
            "tuned": "no",
            "leakage": leak,
            "curator": CURATOR,
            "date": DATE,
            "notes": " ".join(notes),
            "serves": ",".join(serves),
            "clade": r["clade"],
            "strain_taxid": r["strain_taxid"],
            "protein_name": r["protein"],
            "origin": r["origin"],
            "cluster_id": cid[rep],
            "mixed_cluster_with_positive": "yes" if rep in mixed else "no",
            "leakage_step1": lk_step1,
            "leakage_repeat": lk_repeat,
            "leakage_family_domain": lk_fd,
            "family_domain_pfam": ",".join(adh_dom),
            "gpi_feature": "yes" if r["gpi"] else "no",
            "signal_peptide": "yes" if r["signal"] else "no",
            "max_ST30": r["max_ST30"],
            "uniprot_repeat_features": r["repeat_feature"],
            "uniprot_status": r["reviewed"],
        }
    )
cols = list(out[0].keys())
with open(f"{OUT}/controls.tsv", "w") as f:
    w = csv.DictWriter(f, cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(out)
keep = {o["accession"] for o in out}
with open(f"{OUT}/clusters.tsv", "w") as f:
    f.write("accession\tcluster_id\tset\n")
    for rep, ms in sorted(members.items()):
        for m in ms:
            kind, acc = m.split("|", 1)
            if kind == "neg" and acc not in keep:
                continue
            f.write(
                f"{acc}\t{cid[rep]}\t{'hard_negative' if kind == 'neg' else 'adhesin_positive'}\n"
            )
seqs = {}
cur = None
for line in open(f"{W}/sequences.faa"):
    if line.startswith(">"):
        cur = line[1:].strip()
        seqs[cur] = ""
    else:
        seqs[cur] += line.strip()
with open(f"{OUT}/sequences.faa", "w") as f:
    for o in out:
        f.write(f">{o['accession']}\n{seqs[o['accession']]}\n")
print(
    len(out),
    "rows;",
    len({o["cluster_id"] for o in out}),
    "clusters; mixed:",
    sum(o["mixed_cluster_with_positive"] == "yes" for o in out),
)
