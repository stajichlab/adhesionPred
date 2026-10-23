"""Step 3: size, fungal clade spread, reviewed count, cross-kingdom spread and overlaps per PANTHER family.

Usage: python3 -I s3_family_summary.py REPO_ROOT ACC_TSV PFAM_TSV OUT_DIR [MIN_FRAC]
Families = every family named in accession_to_panther.tsv plus every family that holds at least MIN_FRAC (default 0.10)
of the sample of some panel Pfam family in pfam_to_panther.tsv. Counts are UniProt REST (UniProtKB cross-references of
PANTHER; the UniProt release lags InterPro, see the README), the family name and the all-taxa count are InterPro REST.
The clade rule is that of analysis/wall_fusion_truth/map_kiss.py::clades_from_lineage (fine level), applied to the UniProt
lineage ids with ranks. Fungal streams above MAX_STREAM proteins are not streamed (counts only).
Writes family_summary.tsv, family_clade_counts.tsv (long format) and family_member_overlap.tsv.
"""

import collections
import csv
import re
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import get, get_json  # noqa: E402

root, acc_tsv, pfam_tsv, outdir = (Path(x) for x in sys.argv[1:5])
MIN_FRAC = float(sys.argv[5]) if len(sys.argv) > 5 else 0.10
MAX_STREAM = 60000
U = "https://rest.uniprot.org/"

fams = collections.defaultdict(set)  # family -> why
with acc_tsv.open(newline="") as fh:
    acc_rows = list(csv.DictReader(fh, delimiter="\t"))
for r in acc_rows:
    for f in r["panther_family"].split(";"):
        if f:
            fams[f].add("accession")
with pfam_tsv.open(newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["panther_family"] and float(r["frac_of_sample"] or 0) >= MIN_FRAC:
            fams[r["panther_family"]].add("pfam:" + r["pfam_acc"])
ours = {r["accession"]: r for r in acc_rows}

panel = set()  # Pfam accessions of the panel
with pfam_tsv.open(newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        panel.add(r["pfam_acc"])
panel_family_table = set()
with (root / "data/sorting_hat/family_table.tsv").open(newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["active"] == "yes":
            panel_family_table.add(r["pfam_acc"])
panel_models = set()
with (root / "analysis/wall_fusion_report/wall_fusion_models.tsv").open(newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        panel_models.add(r["id"])
shipped = panel_family_table | panel_models  # the Pfam panel actually used by the project

RANK = re.compile(r"(\d+) \((\w[\w ]*)\)")


def count(query):
    _, t = get(
        f"{U}uniprotkb/search?query={urllib.parse.quote(query)}&format=tsv&fields=accession&size=1",
        want_headers=True,
    )
    return int(t or 0)


names_cache = {}


def taxname(ids):
    need = [i for i in ids if i not in names_cache]
    for k in range(0, len(need), 40):
        ch = need[k : k + 40]
        q = " OR ".join("tax_id:" + x for x in ch)
        d = get_json(
            f"{U}taxonomy/search?query={urllib.parse.quote(q)}&fields=id,scientific_name&format=json&size=100"
        )
        for r in d.get("results", []):
            names_cache[str(r["taxonId"])] = r["scientificName"]
    return {i: names_cache.get(i, i) for i in ids}


def fine_clade(lineage):
    lin = {}
    for i, rk in RANK.findall(lineage):
        lin.setdefault(rk, i)
    phylum, sub, cls = lin.get("phylum", ""), lin.get("subphylum", ""), lin.get("class", "")
    if phylum == "4890":  # Ascomycota
        return cls if sub == "147538" else sub  # Pezizomycotina -> class, else subphylum
    if phylum == "5204":  # Basidiomycota
        return cls
    return phylum


KINGDOMS = {"Metazoa": "33208", "Viridiplantae": "33090", "Bacteria": "2", "Archaea": "2157"}
sumrows, clade_rows, ov_rows = [], [], []
for fam in sorted(fams):
    meta = get_json(f"https://www.ebi.ac.uk/interpro/api/entry/panther/{fam}")
    m = meta.get("metadata", {})
    name = (m.get("name") or {}).get("name", "")
    ipr = m.get("integrated") or ""
    n_ip_all = (m.get("counters") or {}).get("proteins", "")
    base = f"xref:panther-{fam}"
    n_all = count(base)
    n_fung = count(f"({base}) AND (taxonomy_id:4751)")
    n_rev = count(f"({base}) AND (taxonomy_id:4751) AND (reviewed:true)")
    n_nonf = count(f"({base}) AND (NOT taxonomy_id:4751)")
    kg = {k: count(f"({base}) AND (taxonomy_id:{v})") for k, v in KINGDOMS.items()}
    clade_n, pf_any, pf_shipped, ctrl, miss = collections.Counter(), 0, 0, [], []
    streamed = n_fung <= MAX_STREAM and n_fung > 0
    pf_by = collections.Counter()
    if streamed:
        q = urllib.parse.quote(f"({base}) AND (taxonomy_id:4751)")
        tsv = get(
            f"{U}uniprotkb/stream?query={q}&format=tsv&fields=accession,reviewed,organism_id,lineage_ids,xref_pfam"
        )
        rows = [ln.split("\t") for ln in tsv.strip().split("\n")[1:]]
        lin_cache = {}
        raw = []
        for c in rows:
            c += [""] * (5 - len(c))
            fc = fine_clade(c[3])
            raw.append((c[0], c[1], fc, {x for x in c[4].split(";") if x}))
        nm = taxname(sorted({fc for _, _, fc, _ in raw if fc}))
        for acc, rev, fc, pfs in raw:
            clade_n[nm.get(fc, fc) if fc else "(no fine clade)"] += 1
            if pfs & panel:
                pf_any += 1
            if pfs & shipped:
                pf_shipped += 1
            for p in pfs & panel:
                pf_by[p] += 1
            if acc in ours:
                ov_rows.append(
                    [
                        fam,
                        acc,
                        ours[acc]["source_set"],
                        ours[acc]["gene"],
                        ours[acc]["slice"],
                        nm.get(fc, fc),
                        rev,
                    ]
                )
        n_fung = len(
            raw
        )  # the stream is the reference count (the search count differs by under 1 percent)
        n_rev = sum(1 for _, rev, _, _ in raw if rev == "reviewed")
        for cl, n in clade_n.most_common():
            clade_rows.append([fam, cl, n])
    n_cl5 = sum(1 for k, v in clade_n.items() if v >= 5 and not k.startswith("("))
    top = ";".join(f"{k}:{v}" for k, v in clade_n.most_common(5))
    nonf = sum(kg.values())
    kd_present = sum(1 for v in kg.values() if v >= 100)
    n_all = n_fung + n_nonf
    if n_all >= 20000 or (kd_present >= 3 and n_all >= 5000):
        flag = "broad"
    elif kd_present >= 2 or n_fung / max(n_all, 1) < 0.5:
        flag = "intermediate"
    else:
        flag = "specific"
    in_ours = sum(1 for r in ov_rows if r[0] == fam)
    sumrows.append(
        [
            fam,
            name,
            ipr,
            "+".join(sorted(fams[fam])),
            n_ip_all,
            n_all,
            n_fung,
            n_rev,
            n_fung - n_rev,
            n_cl5,
            top,
            kg["Metazoa"],
            kg["Viridiplantae"],
            kg["Bacteria"],
            kg["Archaea"],
            f"{n_fung / n_all:.3f}" if n_all else "",
            flag,
            pf_any if streamed else "",
            pf_shipped if streamed else "",
            ";".join(f"{k}:{v}" for k, v in pf_by.most_common(4)),
            in_ours,
            "yes" if streamed else "no_counts_only",
        ]
    )
    print(fam, name, n_all, n_fung, n_rev, flag, flush=True)

outdir.mkdir(parents=True, exist_ok=True)
with (outdir / "family_summary.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(
        [
            "panther_family",
            "name",
            "interpro_entry",
            "why_listed",
            "n_all_taxa_interpro",
            "n_all_taxa_uniprot",
            "n_fungal",
            "n_reviewed_fungal",
            "n_unreviewed_fungal",
            "n_fine_clades_ge5",
            "top_clades",
            "n_metazoa",
            "n_viridiplantae",
            "n_bacteria",
            "n_archaea",
            "frac_fungal",
            "broad_or_specific",
            "n_fungal_with_panel_pfam",
            "n_fungal_with_shipped_pfam",
            "top_panel_pfam_counts",
            "n_our_accessions_in_family",
            "fungal_list_streamed",
        ]
    )
    w.writerows(sumrows)
with (outdir / "family_clade_counts.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["panther_family", "clade_fine", "n_fungal_proteins"])
    w.writerows(clade_rows)
with (outdir / "family_member_overlap.tsv").open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(
        [
            "panther_family",
            "accession",
            "source_set",
            "gene",
            "slice",
            "clade_fine",
            "uniprot_status",
        ]
    )
    w.writerows(ov_rows)
