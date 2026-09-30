#!/usr/bin/env python3
"""Map Pro/Cys-rich (SOWgp/BAD1-type) class-2a candidates to C. immitis RS and join RS1 TPM.

The spherule/mycelium RNA-seq (RS1_kallisto.TPM.csv) is quantified on C. immitis RS
transcripts only. Candidates from the long-read strains and Silveira are mapped to RS by
DIAMOND blastp best hit. Low-complexity masking is off, because the repeat arrays are
Pro-rich and masking removes most of the signal.

Run on a compute node, not the login node (diamond makedb is killed there):
    srun -p short -c 4 --mem=8G --time=0:30:00 \
        bash -lc 'module load diamond/2.1.24 && python3 04_procys_rs_expression.py'
The RS-genome tblastn check for the unmapped Pro-rich family (report section 4.5) was run
by hand and is not part of this script.
Outputs:
    procys_rs_map.tsv   one row per candidate: best RS hit, identity, coverage
    procys_rs_tpm.tsv   one row per RS locus: TPM per replicate, condition means, log2FC
"""

import csv
import math
import statistics as st
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis._common.paths import COCCI_LONGREAD, COCCI_PANGENOME, SPHERULE_RNASEQ  # noqa: E402

HERE = Path(__file__).resolve().parent
CLASS = "ProCys_rich(SOWgp_BAD1_type)"
RS_FA = COCCI_PANGENOME / "input_run2" / "CimmitisRS_FungiDB.fasta"
TPM = SPHERULE_RNASEQ / "RS1_kallisto.TPM.csv"


def read_fasta(path):
    name, seq, hdr = None, [], ""
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, hdr, "".join(seq)
            hdr = line[1:].rstrip()
            name, seq = hdr.split()[0], []
        else:
            seq.append(line.strip())
    if name:
        yield name, hdr, "".join(seq)


def proteome(strain):
    if strain.endswith("_FungiDB"):
        return COCCI_PANGENOME / "input_run2" / f"{strain}.fasta"
    return next(COCCI_LONGREAD.glob(f"*/{strain}.proteins.fa"))


cands = [
    r
    for r in csv.DictReader(open(HERE / "class2a_candidates.tsv"), delimiter="\t")
    if r["comp_class"] == CLASS
]
want = {}
for r in cands:
    want.setdefault(r["strain"], set()).add(r["protein"])
seqs = {}
for strain, ids in want.items():
    for name, _, s in read_fasta(proteome(strain)):
        if name in ids:
            seqs[(strain, name)] = s
missing = [(r["strain"], r["protein"]) for r in cands if (r["strain"], r["protein"]) not in seqs]
if missing:
    sys.exit(f"sequences not found: {missing}")

product = {}
for name, hdr, _ in read_fasta(RS_FA):
    f = dict(x.split("=", 1) for x in hdr.split(" | ")[1:] if "=" in x)
    product[name] = f.get("gene_product", "")

with tempfile.TemporaryDirectory() as td:
    q = Path(td) / "q.faa"
    with open(q, "w") as fh:
        for (strain, name), s in seqs.items():
            fh.write(f">{strain}|{name}\n{s}\n")
    subprocess.run(
        ["diamond", "makedb", "--in", str(RS_FA), "-d", f"{td}/rs", "--quiet"], check=True
    )
    out = Path(td) / "hits.m8"
    subprocess.run(
        [
            "diamond",
            "blastp",
            "-q",
            str(q),
            "-d",
            f"{td}/rs",
            "-o",
            str(out),
            "--ultra-sensitive",
            "--masking",
            "0",
            "--max-target-seqs",
            "5",
            "--evalue",
            "1e-5",
            "--quiet",
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qlen",
            "slen",
            "evalue",
            "bitscore",
            "qcovhsp",
            "scovhsp",
        ],
        check=True,
    )
    hits = {}
    for line in open(out):
        f = line.rstrip("\n").split("\t")
        hits.setdefault(f[0], []).append(f)

tpm = {}
with open(TPM) as fh:
    rd = csv.reader(fh)
    header = next(rd)[1:]
    for row in rd:
        tpm[row[0].strip('"')] = [float(x) for x in row[1:]]

maprows = []
for r in cands:
    key = f"{r['strain']}|{r['protein']}"
    h = hits.get(key, [])
    best = h[0] if h else None
    maprows.append(
        {
            "strain": r["strain"],
            "protein": r["protein"],
            "length": r["length"],
            "rep_period": r["rep_period"],
            "rep_n_copies": r["rep_n_copies"],
            "pct_pro": r["pct_pro"],
            "pct_cys": r["pct_cys"],
            "rs_best_hit": best[1] if best else "",
            "pident": best[2] if best else "",
            "qcov": best[8] if best else "",
            "scov": best[9] if best else "",
            "evalue": best[6] if best else "",
            "second_hit": h[1][1] if len(h) > 1 else "",
            "second_bits": h[1][7] if len(h) > 1 else "",
            "best_bits": best[7] if best else "",
        }
    )

with open(HERE / "procys_rs_map.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(maprows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(maprows)

loci = {}
for m in maprows:
    if m["rs_best_hit"]:
        loci.setdefault(m["rs_best_hit"], []).append(m)

rows = []
for rs, ms in loci.items():
    tx = rs.rsplit("-p", 1)[0]
    v = tpm.get(tx)
    rec = {
        "rs_protein": rs,
        "rs_product": product.get(rs, ""),
        "n_candidates": len(ms),
        "strains": ",".join(sorted({m["strain"] for m in ms})),
        "min_pident": min(float(m["pident"]) for m in ms),
    }
    if v:
        rec.update({h: round(x, 2) for h, x in zip(header, v)})
        myc, s48, s8d = st.mean(v[0:2]), st.mean(v[2:4]), st.mean(v[4:6])
        rec.update(
            mycelia_mean=round(myc, 2),
            spherule48h_mean=round(s48, 2),
            spherule8d_mean=round(s8d, 2),
            log2FC_48h_vs_myc=round(math.log2((s48 + 1) / (myc + 1)), 2),
            log2FC_8d_vs_myc=round(math.log2((s8d + 1) / (myc + 1)), 2),
        )
    rows.append(rec)
rows.sort(key=lambda r: -r.get("spherule48h_mean", -1))

fields = [
    "rs_protein",
    "rs_product",
    "n_candidates",
    "strains",
    "min_pident",
    *header,
    "mycelia_mean",
    "spherule48h_mean",
    "spherule8d_mean",
    "log2FC_48h_vs_myc",
    "log2FC_8d_vs_myc",
]
with open(HERE / "procys_rs_tpm.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", restval="NA")
    w.writeheader()
    w.writerows(rows)

print(
    f"{len(cands)} candidates -> {len(loci)} RS loci; "
    f"{sum(1 for m in maprows if not m['rs_best_hit'])} with no RS hit"
)
