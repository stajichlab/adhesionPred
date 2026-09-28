"""CBM18 (Pfam Chitin_bind_1, PF00187) vs adhesion calls in Batrachochytrium and
two non-pathogenic chytrid outgroups. Reports Fisher enrichment, copy number and length.

Run on HPCC (needs function.duckdb): python 01_cbm18_vs_calls.py
"""

import csv
import sys
from pathlib import Path

import duckdb
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis._common.paths import FUNGI5K_DUCKDB, REPO_ROOT  # noqa: E402

DB = str(FUNGI5K_DUCKDB)
RES = str(REPO_ROOT / "results")
SP = {
    "Batrachochytrium_dendrobatidis_JAM81": ("F2ADAE73", 7021),
    "Batrachochytrium_salamandrivorans_AMFP13": ("F61BA062", 16260),
    "Homolaphlyctis_polyrhiza_JEL142": ("F37F92FF", None),
    "Spizellomyces_punctatus_DAOM_BR117": ("F931F230", None),
}
con = duckdb.connect(DB, read_only=True)
con.execute("PRAGMA memory_limit='3GB'")
for name, (lt, _) in SP.items():
    try:
        calls = {
            r["id"]: float(r["probability_adhesion"])
            for r in csv.DictReader(open(f"{RES}/{name}.adhesion_predict.csv"))
        }
    except FileNotFoundError:
        print(f"=== {name}: no result file")
        continue
    tot = con.execute(
        f"SELECT count(*) FROM gene_proteins WHERE protein_id LIKE '{lt}%'"
    ).fetchone()[0]
    pf = con.execute(f"""SELECT protein_id, string_agg(DISTINCT pfam_id, chr(59)) AS dom, count(*) n
                       FROM pfam WHERE protein_id LIKE '{lt}%' GROUP BY protein_id""").fetchdf()
    pf["protein_id"] = pf.protein_id.str.replace(".protein", "", regex=False)
    cb = pf[pf.dom.str.contains("Chitin_bind_1", na=False)].copy()
    ln = con.execute(f"""SELECT replace(protein_id,'.protein','') protein_id, length FROM gene_proteins
                       WHERE protein_id LIKE '{lt}%'""").fetchdf()
    cb = cb.merge(ln, on="protein_id", how="left")
    cb["called"] = cb.protein_id.isin(calls)
    cb["p"] = cb.protein_id.map(calls)
    k = int(cb.called.sum())
    N = len(cb)
    if N:
        odds, pv = fisher_exact(
            [[k, N - k], [len(calls) - k, tot - N - (len(calls) - k)]], alternative="greater"
        )
        print(f"=== {name}  proteome={tot} calls={len(calls)} ({len(calls) / tot:.2%})")
        print(
            f"    CBM18 proteins={N}, called={k} ({k / N:.0%})  enrichment OR={odds:.1f}, Fisher p={pv:.2g}"
        )
        print(
            f"    length: median {cb.length.median():.0f}, range {cb.length.min():.0f}-{cb.length.max():.0f}  "
            f"(#CBM18 copies/protein: {cb.n.min()}-{cb.n.max()})"
        )
        for _, r in cb.sort_values("length", ascending=False).iterrows():
            print(
                f"      {'CALL' if r.called else '    '} {str(r.p)[:5]:>5} len={int(r.length):>5} copies={r.n:>2}  {r.dom[:70]}"
            )
    else:
        print(f"=== {name}  proteome={tot} calls={len(calls)}: no CBM18")
