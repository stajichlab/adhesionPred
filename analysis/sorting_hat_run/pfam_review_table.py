#!/usr/bin/env python3
"""One review table per proteome: every Pfam family hit with the facts a reviewer needs.

Plan 2, Task 13 step 3. Columns: family, protein, annotation (FASTA header), length, n_domains,
best_domain_score, best_i_evalue, domain_cys (Cys in the envelope of the best-scoring domain), n_tm (TMHMM helices),
signal_peptide (rule R0 call, when the module table
exists), draft_member (from make_draft_members.py), review_class (empty: a reviewer fills it).
The classes to use: false_domain_hit, uncharacterised_true_member, receptor_like, known_member.

Usage: pfam_review_table.py --domtbl D --fasta F --tm-table T [--sp-module step1_rule@R0.tsv.gz]
         [--members members.tsv] --out OUT.tsv
"""

import argparse
import collections
import csv
import gzip
import sys


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_fasta(p):
    seqs, heads, name, buf = {}, {}, None, []
    with opener(p) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).rstrip("*")
                name = line[1:].split()[0]
                heads[name] = line[1:]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).rstrip("*")
    return seqs, heads


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--tm-table", required=True)
    ap.add_argument("--sp-module")
    ap.add_argument("--members")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    seqs, heads = read_fasta(a.fasta)
    tm = {}
    with open(a.tm_table) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            tm[r["protein_id"]] = r
    sp = {}
    if a.sp_module:
        with opener(a.sp_module) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                sp[r["id"]] = r
    members = collections.defaultdict(set)
    if a.members:
        with open(a.members) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                members[r["pfam_acc"]].add(r["protein_id"])
    hits = collections.defaultdict(list)
    with open(a.domtbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split(None, 22)
            acc = c[4].split(".")[0]
            hits[(acc, c[3], c[0])].append(
                {"score": float(c[13]), "ie": float(c[12]), "env": (int(c[19]), int(c[20]))}
            )
    rows = []
    for (acc, name, prot), ds in sorted(hits.items()):
        best = max(ds, key=lambda d: d["score"])
        s = seqs.get(prot, "")
        e0, e1 = best["env"]
        t = tm.get(prot, {})
        row = {
            "family": f"{acc} {name}",
            "protein": prot,
            "annotation": heads.get(prot, "")[:120],
            "length": len(s),
            "n_domains": len(ds),
            "best_domain_score": f"{best['score']:.1f}",
            "best_i_evalue": f"{best['ie']:.2g}",
            "domain_cys": s[e0 - 1 : e1].count("C"),
            "n_tm": t.get("pred_hel", "NA"),
            "signal_peptide": sp.get(prot, {}).get("call", "NA"),
            "draft_member": "yes" if prot in members.get(acc, set()) else "no",
            "review_class": "",
        }
        rows.append(row)
    cols = list(rows[0]) if rows else []
    with open(a.out, "w") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{a.out}: {len(rows)} family-protein hits", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
