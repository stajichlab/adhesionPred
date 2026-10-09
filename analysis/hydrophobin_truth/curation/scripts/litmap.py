# ruff: noqa
"""Map paper gene IDs (via UniProt sequences) to proteome IDs by BLASTP.
Usage: litmap.py <lit_genes.tsv> <lit_seqs.faa> <blastdb dir> <q61 ids tsv: evidence_sheets.tsv> <out tsv>
Needs blastp on PATH. Best hit per (paper gene, proteome) by bitscore. Coverage merged over HSPs.
mapped = 'yes' if identity >= 95 and coverage >= 90 of BOTH sequences; 'partial' if >= 95 identity and >= 90 of one; else 'no'.
"""

import csv
import subprocess
import sys
import tempfile
from collections import defaultdict

genes, seqf, dbdir, sheet, out = sys.argv[1:6]
seqs, name = {}, None
for line in open(seqf):
    if line.startswith(">"):
        name = line[1:].split("|")[1]
        seqs[name] = ""
    else:
        seqs[name] += line.strip()
in61 = {(r["proteome"], r["id"]) for r in csv.DictReader(open(sheet), delimiter="\t")}


def cov(ivs, L):
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


rows = []
for g in csv.DictReader(open(genes), delimiter="\t"):
    for p in g["proteomes"].split(","):
        with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
            fh.write(f">{g['uniprot']}\n{seqs[g['uniprot']]}\n")
        res = subprocess.run(
            [
                "blastp",
                "-query",
                fh.name,
                "-db",
                f"{dbdir}/{p}",
                "-evalue",
                "1e-5",
                "-max_target_seqs",
                "5",
                "-outfmt",
                "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore",
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        hits = defaultdict(list)
        for line in res.splitlines():
            f = line.split("\t")
            hits[f[1]].append(f)
        o = dict(
            paper=f"{g['paper']} (PMID {g['pmid']})",
            paper_gene_id=g["paper_gene_id"],
            proteome=p,
            proteome_id="",
            identity="",
            coverage="",
            mapped="no",
            how_obtained=f"{g['id_source']}; sequence UniProt {g['uniprot']}; BLASTP 2.14.0+ vs {p} proteome",
            in_61="",
        )
        if hits:
            s = max(hits, key=lambda k: max(float(x[11]) for x in hits[k]))
            hs = hits[s]
            top = max(hs, key=lambda x: float(x[11]))
            pid = float(top[2])
            qc = cov([(int(x[6]), int(x[7])) for x in hs], int(top[4]))
            sc = cov([(int(x[8]), int(x[9])) for x in hs], int(top[5]))
            m = (
                "yes"
                if pid >= 95 and qc >= 90 and sc >= 90
                else ("partial" if pid >= 95 and (qc >= 90 or sc >= 90) else "no")
            )
            o.update(
                proteome_id=s,
                identity=f"{pid:.1f}",
                coverage=f"paper protein {qc:.0f}% / proteome protein {sc:.0f}%",
                mapped=m,
                in_61="yes" if (p, s) in in61 else "no",
            )
        else:
            o["how_obtained"] += "; no hit at E<=1e-5"
        rows.append(o)
cols = [
    "paper",
    "paper_gene_id",
    "proteome",
    "proteome_id",
    "identity",
    "coverage",
    "mapped",
    "in_61",
    "how_obtained",
]
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(
    len(rows),
    "rows;",
    sum(r["mapped"] == "yes" for r in rows),
    "mapped yes;",
    sum(r["mapped"] == "partial" for r in rows),
    "partial;",
    sum(r["in_61"] == "yes" and r["mapped"] != "no" for r in rows),
    "mapped/partial to one of the 61",
)
