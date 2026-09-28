"""CBM18 domain architecture (alone / +chitin deacetylase / +tyrosinase / +GH18),
signal-peptide status and call status; writes cbm18_batrachochytrium.tsv.

Run on HPCC (needs function.duckdb): python 02_cbm18_architecture.py
"""

import csv
import sys
from pathlib import Path

import duckdb
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis._common.paths import FUNGI5K_DUCKDB, REPO_ROOT  # noqa: E402

DB = str(FUNGI5K_DUCKDB)
RES = str(REPO_ROOT / "results")
SP = {
    "Bd_JAM81": ("F2ADAE73", "Batrachochytrium_dendrobatidis_JAM81"),
    "Bsal_AMFP13": ("F61BA062", "Batrachochytrium_salamandrivorans_AMFP13"),
}
con = duckdb.connect(DB, read_only=True)
con.execute("PRAGMA memory_limit='3GB'")
print("signalp cols:", [d[0] for d in con.execute("select * from signalp limit 1").description])
rows = []
for tag, (lt, fn) in SP.items():
    calls = {
        r["id"]: r["probability_adhesion"]
        for r in csv.DictReader(open(f"{RES}/{fn}.adhesion_predict.csv"))
    }
    pf = con.execute(f"""SELECT protein_id, string_agg(DISTINCT pfam_id, chr(59)) AS dom, count(*) n
                       FROM pfam WHERE protein_id LIKE '{lt}%' GROUP BY protein_id""").fetchdf()
    pf["protein_id"] = pf.protein_id.str.replace(".protein", "", regex=False)
    cb = pf[pf.dom.str.contains("Chitin_bind_1", na=False)].copy()
    sp = con.execute(
        f"SELECT DISTINCT replace(protein_id,'.protein','') protein_id FROM signalp WHERE protein_id LIKE '{lt}%'"
    ).fetchdf()
    ln = con.execute(
        f"SELECT replace(protein_id,'.protein','') protein_id, length FROM gene_proteins WHERE protein_id LIKE '{lt}%'"
    ).fetchdf()
    cb = cb.merge(ln, on="protein_id", how="left")
    cb["signal_peptide"] = cb.protein_id.isin(set(sp.protein_id))
    cb["called_adhesion"] = cb.protein_id.isin(calls)
    cb["probability"] = cb.protein_id.map(calls)
    cb["species"] = tag

    def arch(d):
        return (
            "CBM18_only"
            if d.replace("Chitin_bind_1", "").strip(";") == ""
            else "CBM18+tyrosinase"
            if "Tyrosinase" in d
            else "CBM18+chitin_deacetylase"
            if "Polysacc_deac" in d
            else "CBM18+GH18"
            if "Glyco_hydro_18" in d
            else "CBM18+other"
        )

    cb["architecture"] = cb.dom.map(arch)
    rows.append(cb)
df = pd.concat(rows)
print("\nsignal peptide by architecture and call status:")
print(
    df.groupby(["species", "architecture"])
    .agg(
        n=("protein_id", "size"),
        signal_peptide=("signal_peptide", "sum"),
        called=("called_adhesion", "sum"),
        median_len=("length", "median"),
        copies=("n", "max"),
    )
    .to_string()
)
df.rename(columns={"n": "cbm18_copies", "dom": "pfam_domains"})[
    [
        "species",
        "protein_id",
        "architecture",
        "cbm18_copies",
        "length",
        "signal_peptide",
        "called_adhesion",
        "probability",
        "pfam_domains",
    ]
].sort_values(["species", "architecture", "length"], ascending=[True, True, False]).to_csv(
    "cbm18_batrachochytrium.tsv", sep="\t", index=False
)
print("\nwrote cbm18_batrachochytrium.tsv")
