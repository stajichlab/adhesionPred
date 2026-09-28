"""Known Batrachochytrium virulence/pathogenicity families (M36 fungalysin, S41,
CBM18, tyrosinase, VWD, CAP, AA1 laccase) vs adhesion calls and signal peptides.
Shows the classifier ignores the secreted protease expansions entirely.

Run on HPCC (needs function.duckdb): python 05_known_virulence_families.py
"""

import csv

import duckdb

DB = "/bigdata/stajichlab/shared/projects/Fungi_5k/functionalDB/function.duckdb"
RES = "/bigdata/stajichlab/jstajich/projects/adhesionPred/results"
SP = {
    "Bd_JAM81": ("F2ADAE73", "Batrachochytrium_dendrobatidis_JAM81"),
    "Bsal_AMFP13": ("F61BA062", "Batrachochytrium_salamandrivorans_AMFP13"),
    "Hpolyrhiza": ("F37F92FF", "Homolaphlyctis_polyrhiza_JEL142"),
    "Spunctatus": ("F931F230", "Spizellomyces_punctatus_DAOM_BR117"),
}
FAM = {
    "M36 fungalysin": "Peptidase_M36",
    "CBM18": "Chitin_bind_1",
    "tyrosinase": "Tyrosinase",
    "VWD": "VWD",
    "CAP/PR-1": "^CAP$",
    "Cu-oxidase/AA1": "Cu-oxidase",
    "crinkler": "CRN|Crinkler",
    "S41 protease": "Peptidase_S41",
    "aspartyl protease": "^Asp$",
    "chitin deacetylase": "Polysacc_deac",
}
con = duckdb.connect(DB, read_only=True)
con.execute("PRAGMA memory_limit='3GB'")
print(f"{'family':<22}" + "".join(f"{t:>20}" for t in SP))
rows = {}
for tag, (lt, fn) in SP.items():
    calls = {r["id"] for r in csv.DictReader(open(f"{RES}/{fn}.adhesion_predict.csv"))}
    pf = con.execute(f"""SELECT replace(protein_id,'.protein','') protein_id,
                       string_agg(DISTINCT pfam_id, chr(59)) dom FROM pfam
                       WHERE protein_id LIKE '{lt}%' GROUP BY protein_id""").fetchdf()
    sp = set(
        con.execute(
            f"SELECT DISTINCT replace(protein_id,'.protein','') protein_id FROM signalp WHERE protein_id LIKE '{lt}%'"
        )
        .fetchdf()
        .protein_id
    )
    rows[tag] = (pf, calls, sp)
for label, pat in FAM.items():
    cells = []
    for _tag, (pf, calls, sp) in rows.items():
        h = pf[pf.dom.str.contains(pat, na=False, regex=True, case=False)]
        k = len(set(h.protein_id) & calls)
        s = len(set(h.protein_id) & sp)
        cells.append(f"{len(h):>4} ({k:>3} call,{s:>3} SP)")
    print(f"{label:<22}" + "".join(f"{c:>20}" for c in cells))
