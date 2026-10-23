"""Step 1: per-accession PANTHER family and subfamily, InterPro entries and Pfam entries (InterPro REST).

Usage: python3 -I s1_accession_to_panther.py REPO_ROOT OUT_TSV
Accessions: (a) clusters_missed.tsv of the Kiss baseline, (c) controls.tsv of the wall-fusion control set.
Keys that are not UniProt accessions (jgi:...) are skipped and counted. No sequences are fetched.
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import get_json  # noqa: E402

root, out = Path(sys.argv[1]), Path(sys.argv[2])
rows = {}  # accession -> dict
skipped = []
with (root / "docs/reports/data/kiss_baseline_2026-10-22/clusters_missed.tsv").open(
    newline=""
) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        accs = r["accessions"].split(",")
        genes = r["genes"].split(",")
        for i, acc in enumerate(accs):
            if acc.startswith("jgi:"):
                skipped.append(("missed44", acc))
                continue
            g = genes[i] if len(genes) == len(accs) else r["genes"]
            d = rows.setdefault(
                acc, {"sets": set(), "gene": g, "slice": r["slice"], "cluster": r["cluster_id"]}
            )
            d["sets"].add("missed44")
with (root / "data/controls/wall-fusion-kiss2019/controls.tsv").open(newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        acc = r["accession"]
        if acc.startswith("jgi:"):
            skipped.append(("controls", acc))
            continue
        d = rows.setdefault(
            acc,
            {"sets": set(), "gene": r["gene"], "slice": r["stratum"], "cluster": r["cluster_id"]},
        )
        d["sets"].add("controls112")
        d["gene"] = r["gene"]

B = "https://www.ebi.ac.uk/interpro/api/entry/{db}/protein/uniprot/{acc}"


def entries(db, acc):
    d = get_json(B.format(db=db, acc=acc))
    return d.get("results", [])


cols = [
    "accession",
    "source_set",
    "slice",
    "kiss_cluster",
    "gene",
    "panther_family",
    "panther_subfamily",
    "panther_family_name",
    "panther_score",
    "interpro_entries",
    "pfam_entries",
]
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for acc in sorted(rows):
        d = rows[acc]
        pan = entries("panther", acc)
        fam, sub, name, score = [], [], [], []
        for e in pan:
            m = e["metadata"]
            fam.append(m["accession"])
            name.append(m["name"])
            for p in e["proteins"]:
                for loc in p["entry_protein_locations"]:
                    if loc.get("subfamily"):
                        sub.append(loc["subfamily"]["accession"])
                    score.append(str(loc.get("score", "")))
        ipr = [
            f'{e["metadata"]["accession"]}={e["metadata"]["name"]}'
            for e in entries("interpro", acc)
        ]
        pf = [f'{e["metadata"]["accession"]}={e["metadata"]["name"]}' for e in entries("pfam", acc)]
        w.writerow(
            [
                acc,
                "+".join(sorted(d["sets"])),
                d["slice"],
                d["cluster"],
                d["gene"],
                ";".join(fam),
                ";".join(sub),
                ";".join(name),
                ";".join(score),
                ";".join(ipr),
                ";".join(pf),
            ]
        )
print(len(rows), "accessions; skipped jgi keys:", len(skipped), file=sys.stderr)
for s in skipped:
    print("skipped", *s, file=sys.stderr)
