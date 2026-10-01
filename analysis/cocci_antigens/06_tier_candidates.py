import csv
import math
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis._common.paths import FUNGI5K_DUCKDB  # noqa: E402

DB = str(FUNGI5K_DUCKDB)
RS = "FA2214EC"
con = duckdb.connect(DB, read_only=True)
con.execute("PRAGMA memory_limit='3GB'")
best = {}
for line in open("cocci_antigens/idmap.m8"):
    f = line.rstrip("\n").split("\t")
    q, t, pid = f[0], f[1], float(f[2])
    if q not in best or pid > best[q][1]:
        best[q] = (t, pid)
sp = set(
    con.execute(
        f"SELECT DISTINCT replace(protein_id,'.protein','') p FROM signalp WHERE protein_id LIKE '{RS}%'"
    )
    .fetchdf()
    .p
)
tm = set(
    con.execute(
        f"SELECT DISTINCT replace(protein_id,'.protein','') p FROM tmhmm WHERE protein_id LIKE '{RS}%'"
    )
    .fetchdf()
    .p
)
pf = con.execute(
    f"SELECT replace(protein_id,'.protein','') p, pfam_id FROM pfam WHERE protein_id LIKE '{RS}%'"
).fetchdf()

doms = defaultdict(set)
for p, d in zip(pf.p, pf.pfam_id, strict=False):
    doms[p].add(d)
ln = con.execute(
    f"SELECT replace(protein_id,'.protein','') p, length FROM gene_proteins WHERE protein_id LIKE '{RS}%'"
).fetchdf()
length = dict(zip(ln.p, ln.length, strict=False))
tpm = {}
with open("RS1_kallisto.TPM.csv") as f:
    rd = csv.reader(f)
    next(rd)
    for row in rd:
        v = [float(x) for x in row[1:]]
        tpm[row[0].strip('"')] = {
            "myc": st.mean(v[0:2]),
            "s48": st.mean(v[2:4]),
            "s8d": st.mean(v[4:6]),
        }


def lfc(d):
    return math.log2((d["s48"] + 1) / (d["myc"] + 1))


rank = [
    x
    for x in csv.DictReader(open("cocci_antigens/cocci_antigen_ranking.tsv"), delimiter="\t")
    if x["is_orthogroup_representative"] == "yes"
]


def g(x, k):
    return float(x[k])


out = []
for x in rank:
    fid = best.get(x["protein"], ("", 0))[0]
    t = tpm.get(x["protein"].replace("-p1", ""))
    x.update(
        fungi5k=fid,
        annotated="yes" if fid else "no",
        signal_peptide="yes" if fid and fid in sp else ("no" if fid else "unknown"),
        tm_helix="yes" if fid and fid in tm else ("no" if fid else "unknown"),
        pfam=";".join(sorted(doms.get(fid, []))),
        length=length.get(fid, ""),
        mycelia_TPM=round(t["myc"], 1) if t else "",
        spherule48h_TPM=round(t["s48"], 1) if t else "",
        log2FC_spherule48h=round(lfc(t), 2) if t else "",
    )
    out.append(x)
sec = [
    x
    for x in out
    if x["signal_peptide"] == "yes"
    and g(x, "prevalence") >= 0.95
    and g(x, "max_fungal_crossreact_pid") == 0
    and g(x, "human_homolog_pid") == 0
]
tier1 = [x for x in sec if x["log2FC_spherule48h"] != "" and float(x["log2FC_spherule48h"]) > 1]
print(f"TIER 1 (secreted + universal + specific + spherule-induced): {len(tier1)}")
print(f"TIER 2 (secreted + universal + specific, any expression):     {len(sec)}")
tier1.sort(key=lambda z: -float(z["spherule48h_TPM"] or 0))
print(f"\n{'protein':<21}{'len':>5}{'TM':>4}{'myc':>8}{'sph48':>8}{'lfc':>7}{'prev':>6}  pfam")
for x in tier1:
    print(
        f"{x['protein'][:20]:<21}{str(x['length']):>5}{x['tm_helix'][:1].upper():>4}"
        f"{float(x['mycelia_TPM']):>8.1f}{float(x['spherule48h_TPM']):>8.1f}"
        f"{float(x['log2FC_spherule48h']):>7.2f}{g(x, 'prevalence'):>6.2f}  {x['pfam'][:40] or '(none)'}"
    )
cols = [
    "protein",
    "orthogroup",
    "anchor",
    "fungi5k",
    "length",
    "signal_peptide",
    "tm_helix",
    "pfam",
    "prevalence",
    "copy_number_mean",
    "copy_number_cv",
    "max_fungal_crossreact_pid",
    "human_homolog_pid",
    "mycelia_TPM",
    "spherule48h_TPM",
    "log2FC_spherule48h",
    "score",
]
for name, rows in [
    ("cocci_antigens/TIER1_candidates.tsv", tier1),
    ("cocci_antigens/TIER2_candidates.tsv", sec),
]:
    with open(name, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {name}")
