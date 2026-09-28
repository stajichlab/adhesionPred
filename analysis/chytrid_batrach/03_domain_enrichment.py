"""Systematic Pfam-domain enrichment among adhesion calls, per chytrid genome
(Fisher exact, BH-corrected). Writes domain_enrichment_chytrids.tsv.

Run on HPCC (needs function.duckdb): python 03_domain_enrichment.py
"""

import csv
from collections import defaultdict

import duckdb
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

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
out = []
for tag, (lt, fn) in SP.items():
    calls = {r["id"] for r in csv.DictReader(open(f"{RES}/{fn}.adhesion_predict.csv"))}
    tot = con.execute(
        f"SELECT count(*) FROM gene_proteins WHERE protein_id LIKE '{lt}%'"
    ).fetchone()[0]
    pf = con.execute(f"""SELECT DISTINCT replace(protein_id,'.protein','') protein_id, pfam_id
                       FROM pfam WHERE protein_id LIKE '{lt}%'""").fetchdf()
    bydom = defaultdict(set)
    for pid, dom in zip(pf.protein_id, pf.pfam_id):
        bydom[dom].add(pid)
    rows = []
    for dom, prots in bydom.items():
        if len(prots) < 3:
            continue
        k = len(prots & calls)
        if k == 0:
            continue
        odds, p = fisher_exact(
            [[k, len(prots) - k], [len(calls) - k, tot - len(prots) - (len(calls) - k)]],
            alternative="greater",
        )
        rows.append(
            {
                "species": tag,
                "domain": dom,
                "n_proteins": len(prots),
                "n_called": k,
                "frac_called": k / len(prots),
                "genome_call_rate": len(calls) / tot,
                "odds_ratio": odds,
                "p": p,
            }
        )
    if rows:
        d = pd.DataFrame(rows)
        d["q"] = multipletests(d.p, method="fdr_bh")[1]
        out.append(d)
res = pd.concat(out).sort_values(["species", "q", "odds_ratio"], ascending=[True, True, False])
res.to_csv("domain_enrichment_chytrids.tsv", sep="\t", index=False, float_format="%.4g")
for sp in res.species.unique():
    d = res[(res.species == sp) & (res.q < 0.05)]
    print(
        f"\n=== {sp}: {len(d)} domains enriched at q<0.05 (call rate {d.genome_call_rate.iloc[0]:.2%})"
        if len(d)
        else f"\n=== {sp}: none at q<0.05"
    )
    for _, r in d.head(15).iterrows():
        print(
            f"   {r.domain:<22} {r.n_called:>3}/{r.n_proteins:<4} ({r.frac_called:>5.0%})  OR={r.odds_ratio:>7.1f}  q={r.q:.2g}"
        )
