"""VWD (von Willebrand factor type D, PF00094) family profile across chytrids:
copy number, length, signal peptide/TM status, partner domains, call status.
Writes vwd_chytrids.tsv.

Run on HPCC (needs function.duckdb): python 04_vwd_profile.py
"""

import csv

import duckdb
import pandas as pd

DB = "/bigdata/stajichlab/shared/projects/Fungi_5k/functionalDB/function.duckdb"
RES = "/bigdata/stajichlab/jstajich/projects/adhesionPred/results"
SP = {
    "Bd_JAM81": ("F2ADAE73", "Batrachochytrium_dendrobatidis_JAM81"),
    "Bsal_AMFP13": ("F61BA062", "Batrachochytrium_salamandrivorans_AMFP13"),
    "Hpolyrhiza": ("F37F92FF", "Homolaphlyctis_polyrhiza_JEL142"),
    "Spunctatus": ("F931F230", "Spizellomyces_punctatus_DAOM_BR117"),
}
con = duckdb.connect(DB, read_only=True)
con.execute("PRAGMA memory_limit='3GB'")
allr = []
for tag, (lt, fn) in SP.items():
    calls = {
        r["id"]: float(r["probability_adhesion"])
        for r in csv.DictReader(open(f"{RES}/{fn}.adhesion_predict.csv"))
    }
    pf = con.execute(f"""SELECT replace(protein_id,'.protein','') protein_id,
                              string_agg(DISTINCT pfam_id, chr(59)) dom, count(*) nd
                       FROM pfam WHERE protein_id LIKE '{lt}%' GROUP BY protein_id""").fetchdf()
    ln = con.execute(
        f"SELECT replace(protein_id,'.protein','') protein_id, length FROM gene_proteins WHERE protein_id LIKE '{lt}%'"
    ).fetchdf()
    sp = set(
        con.execute(
            f"SELECT DISTINCT replace(protein_id,'.protein','') protein_id FROM signalp WHERE protein_id LIKE '{lt}%'"
        )
        .fetchdf()
        .protein_id
    )
    tm = set(
        con.execute(
            f"SELECT DISTINCT replace(protein_id,'.protein','') protein_id FROM tmhmm WHERE protein_id LIKE '{lt}%'"
        )
        .fetchdf()
        .protein_id
    )
    v = pf[pf.dom.str.contains("VWD", na=False)].merge(ln, on="protein_id", how="left")
    tot = len(ln)
    print(f"=== {tag}: proteome {tot}, VWD proteins {len(v)} ({len(v) / tot * 1000:.2f} per 1000)")
    if not len(v):
        continue
    v["signal_peptide"] = v.protein_id.isin(sp)
    v["tm"] = v.protein_id.isin(tm)
    v["called"] = v.protein_id.isin(calls)
    v["p"] = v.protein_id.map(calls)
    v["species"] = tag
    print(
        f"    called {v.called.sum()}/{len(v)}; signal peptide {v.signal_peptide.sum()}; TM {v.tm.sum()}"
    )
    print(
        f"    length median {v.length.median():.0f} range {v.length.min():.0f}-{v.length.max():.0f}"
    )
    print("    partner domains (with VWD):")
    from collections import Counter

    c = Counter(d for s in v.dom for d in s.split(";") if d != "VWD")
    for d, n in c.most_common(10):
        print(f"       {n:>3}  {d}")
    allr.append(v)
if allr:
    pd.concat(allr)[
        ["species", "protein_id", "length", "nd", "signal_peptide", "tm", "called", "p", "dom"]
    ].sort_values(["species", "length"], ascending=[True, False]).to_csv(
        "vwd_chytrids.tsv", sep="\t", index=False
    )
    print("\nwrote vwd_chytrids.tsv")
