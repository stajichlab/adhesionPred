"""Hand-check of PTHR46430 and PTHR35192: subfamilies, kingdom counts, domains, reviewed members, Kiss mapping.

Usage: python3 -I s1_family_tables.py REPO_ROOT OUT_DIR
Calls (public, read-only, cached outside the repo; see the README of the data folder):
  InterPro  entry/panther/{FAM}, entry/panther/{FAM}/?subfamilies, entry/{pfam|interpro}/{id}
  UniProt   uniprotkb/stream?query=xref:panther-{FAM}&format=tsv&fields=... (all taxa), plus a reviewed-only stream with
            function and location comments.
Writes (ids and metadata only, no sequences): family_totals.tsv, subfamily_sizes.tsv, domain_counts.tsv,
reviewed_members.tsv, kiss_mapping.tsv, kiss_in_family.tsv, fungal_clades_by_subfamily.tsv.
"""

import collections
import csv
import re
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import get, get_json  # noqa: E402

root, out = Path(sys.argv[1]), Path(sys.argv[2])
FAMS = ["PTHR46430", "PTHR35192"]
U = "https://rest.uniprot.org/"
IP = "https://www.ebi.ac.uk/interpro/api/"
RANK = re.compile(r"(\d+) \((\w[\w ]*)\)")
FIELDS = "accession,reviewed,organism_name,organism_id,lineage,lineage_ids,xref_pfam,xref_interpro,protein_name,gene_names,length,xref_panther"


def stream(query, fields):
    q = urllib.parse.quote(query)
    tsv = get(f"{U}uniprotkb/stream?query={q}&format=tsv&fields={fields}")
    lines = tsv.strip("\n").split("\n")
    head = lines[0].split("\t")
    return [
        {h: (ln.split("\t") + [""] * len(head))[i] for i, h in enumerate(head)} for ln in lines[1:]
    ]


def fine_clade(lineage_ids):
    lin = {}
    for i, rk in RANK.findall(lineage_ids):
        lin.setdefault(rk, i)
    phylum, sub, cls = lin.get("phylum", ""), lin.get("subphylum", ""), lin.get("class", "")
    if phylum == "4890":
        return cls if sub == "147538" else sub
    if phylum == "5204":
        return cls
    return phylum


def top_groups(lineage):
    """(domain, kingdom) from the UniProt lineage text."""
    dom = king = ""
    for name, rk in re.findall(r"([^,]+?) \((\w[\w ]*)\)", lineage):
        name = name.strip()
        if rk == "domain" or rk == "superkingdom":
            dom = name
        if rk == "kingdom":
            king = name
    return dom, king


names = {}


def taxname(i):
    if i not in names:
        d = get_json(
            f"{U}taxonomy/search?query=tax_id:{i}&fields=id,scientific_name&format=json&size=1"
        )
        r = d.get("results", [])
        names[i] = r[0]["scientificName"] if r else i
    return names[i]


def entry_name(db, acc):
    d = get_json(f"{IP}entry/{db}/{acc}")
    return (d.get("metadata", {}).get("name") or {}).get("name", "") if d else ""


def split_ids(x):
    return [t for t in x.split(";") if t]


def pct(a, b):
    return f"{a / b:.3f}" if b else ""


kiss = {}  # accession -> (cluster, gene)
with (root / "data/controls/wall-fusion-kiss2019/clusters.tsv").open(newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        kiss[r["accession"]] = [r["cluster_id"], ""]
with (root / "data/controls/wall-fusion-kiss2019/controls.tsv").open(newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        kiss.setdefault(r["accession"], [r["cluster_id"], ""])[1] = r["gene"]
with (root / "docs/reports/data/kiss_baseline_2026-10-22/clusters_missed.tsv").open(
    newline=""
) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        gs = r["genes"].split(",")
        accs = r["accessions"].split(",")
        for i, a in enumerate(accs):
            if not a.startswith("jgi:"):
                kiss.setdefault(a, [r["cluster_id"], gs[i] if len(gs) == len(accs) else r["genes"]])

tot_rows, sub_rows, dom_rows, rev_rows, kiss_rows, kissin_rows, clade_rows = (
    [],
    [],
    [],
    [],
    [],
    [],
    [],
)
for fam in FAMS:
    meta = get_json(f"{IP}entry/panther/{fam}")["metadata"]
    subs = get_json(f"{IP}entry/panther/{fam}/?subfamilies").get("results", [])
    subname = {s["metadata"]["accession"]: s["metadata"]["name"] for s in subs}
    rows = stream(f"xref:panther-{fam}", FIELDS)
    for r in rows:
        sf = [t for t in split_ids(r["PANTHER"]) if t.startswith(fam + ":")]
        r["sub"] = sf[0] if sf else "(family only)"
        r["fungal"] = "Fungi (kingdom)" in r["Taxonomic lineage"]
        r["domain"], r["kingdom"] = top_groups(r["Taxonomic lineage"])
        r["clade"] = fine_clade(r["Taxonomic lineage (Ids)"]) if r["fungal"] else ""
    fung = [r for r in rows if r["fungal"]]
    # family totals by domain/kingdom
    dk = collections.Counter((r["domain"], r["kingdom"] or "(none)") for r in rows)
    ndom = collections.Counter(r["domain"] for r in rows)
    for (d, k), n in sorted(dk.items(), key=lambda x: -x[1]):
        tot_rows.append(
            [
                fam,
                meta["name"]["name"],
                len(rows),
                meta["counters"]["proteins"],
                d,
                k,
                n,
                pct(n, len(rows)),
                ndom[d],
            ]
        )
    # subfamilies
    for sf in sorted({r["sub"] for r in rows}):
        m = [r for r in rows if r["sub"] == sf]
        mf = [r for r in m if r["fungal"]]
        cl = collections.Counter(r["clade"] for r in mf)
        dm = collections.Counter(r["domain"] for r in m)
        sub_rows.append(
            [
                fam,
                sf,
                subname.get(sf, ""),
                len(m),
                len(mf),
                pct(len(mf), len(m)),
                sum(r["Reviewed"] == "reviewed" for r in m),
                sum(r["Reviewed"] == "reviewed" for r in mf),
                sum(1 for c, n in cl.items() if n >= 5 and c),
                dm.get("Bacteria", 0),
                dm.get("Archaea", 0),
                dm.get("Eukaryota", 0) - len(mf),
            ]
        )
        for c, n in cl.most_common():
            clade_rows.append([fam, sf, taxname(c) if c else "(no fine clade)", n])
    # domains: fungal, per subfamily and for all of the family
    scopes = [("fungal_all_subfamilies", True, None)]
    scopes += [(f"fungal_{x}", True, x) for x in sorted({r["sub"] for r in fung})]
    scopes += [("non_fungal_all_subfamilies", False, None)]
    for scope, want_fungal, want_sub in scopes:
        grp = [
            r
            for r in rows
            if (r["fungal"] == want_fungal)
            and (want_sub is None or not want_fungal or r["sub"] == want_sub)
        ]
        if not grp:
            continue
        for db, col in (("pfam", "Pfam"), ("interpro", "InterPro")):
            c = collections.Counter(i for r in grp for i in set(split_ids(r[col])))
            none = sum(1 for r in grp if not split_ids(r[col]))
            for i, n in c.most_common():
                if n / len(grp) >= 0.02 or n >= 20:
                    dom_rows.append(
                        [fam, scope, db, i, entry_name(db, i), n, len(grp), pct(n, len(grp))]
                    )
            dom_rows.append(
                [
                    fam,
                    scope,
                    db,
                    "(none)",
                    "no entry of this database",
                    none,
                    len(grp),
                    pct(none, len(grp)),
                ]
            )
    # reviewed members, all taxa, with function and location comments
    rv = stream(
        f"(xref:panther-{fam}) AND (reviewed:true)", FIELDS + ",cc_function,cc_subcellular_location"
    )
    for r in rv:
        sf = [t for t in split_ids(r["PANTHER"]) if t.startswith(fam + ":")]
        dom, king = top_groups(r["Taxonomic lineage"])
        fn = re.sub(r"\s+", " ", r.get("Function [CC]", ""))[:300]
        loc = re.sub(r"\s+", " ", r.get("Subcellular location [CC]", ""))[:160]
        rev_rows.append(
            [
                fam,
                sf[0] if sf else "(family only)",
                r["Entry"],
                r["Gene Names"].split(" ")[0] if r["Gene Names"] else "",
                r["Protein names"][:120],
                re.sub(r" \(.*", "", r["Organism"]),
                r["Organism (ID)"],
                king or dom,
                r["Length"],
                r["Pfam"],
                fn,
                loc,
            ]
        )
    # Kiss mapping and any other Kiss accession in the family
    idx = {r["Entry"]: r for r in rows}
    for a, (cl, g) in sorted(kiss.items(), key=lambda x: x[1][0]):
        if a in idx:
            r = idx[a]
            kissin_rows.append(
                [
                    fam,
                    r["sub"],
                    subname.get(r["sub"], ""),
                    a,
                    cl,
                    g,
                    re.sub(r" \(.*", "", r["Organism"]),
                    r["clade"] and taxname(r["clade"]),
                    r["Reviewed"],
                    r["Pfam"],
                ]
            )


def write(name, head, rows):
    with (out / name).open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(head)
        w.writerows(rows)


out.mkdir(parents=True, exist_ok=True)
write(
    "family_totals.tsv",
    [
        "panther_family",
        "family_name",
        "n_uniprot",
        "n_interpro",
        "domain_of_life",
        "kingdom",
        "n",
        "frac_of_family",
        "n_in_domain",
    ],
    tot_rows,
)
write(
    "subfamily_sizes.tsv",
    [
        "panther_family",
        "subfamily",
        "subfamily_name",
        "n_all",
        "n_fungal",
        "fungal_share",
        "n_reviewed",
        "n_reviewed_fungal",
        "fine_clades_ge5_fungal",
        "n_bacteria",
        "n_archaea",
        "n_other_eukaryota",
    ],
    sub_rows,
)
write(
    "fungal_clades_by_subfamily.tsv",
    ["panther_family", "subfamily", "fine_clade", "n_fungal"],
    clade_rows,
)
write(
    "domain_counts.tsv",
    [
        "panther_family",
        "scope",
        "database",
        "entry",
        "entry_name",
        "n_with_entry",
        "n_in_scope",
        "fraction",
    ],
    dom_rows,
)
write(
    "reviewed_members.tsv",
    [
        "panther_family",
        "subfamily",
        "accession",
        "gene",
        "protein_name",
        "organism",
        "organism_taxid",
        "kingdom_or_domain",
        "length",
        "pfam",
        "function_comment_300ch",
        "location_comment_160ch",
    ],
    rev_rows,
)
write(
    "kiss_in_family.tsv",
    [
        "panther_family",
        "subfamily",
        "subfamily_name",
        "accession",
        "kiss_cluster",
        "gene",
        "organism",
        "fine_clade",
        "uniprot_status",
        "pfam",
    ],
    kissin_rows,
)
print("done", len(rev_rows), "reviewed;", len(kissin_rows), "kiss rows")
