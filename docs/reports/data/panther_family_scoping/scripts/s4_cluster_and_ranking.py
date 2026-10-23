"""Step 4: join the tables into (i) the 44 missed clusters with their PANTHER family, (ii) coverage of the 112 controls
by PANTHER, (iii) a mechanical screen of the families as independent test sets.

Usage: python3 -I s4_cluster_and_ranking.py REPO_ROOT DATA_DIR
Reads accession_to_panther.tsv and family_summary.tsv in DATA_DIR and the PR 6 step 6.0 table scan_per_cluster.tsv.
Screen (all thresholds are stated in the report, section 6; they are a way to sort, not a measured quality):
  size >= 500 fungal proteins; fine clades with >= 5 proteins >= 8; fungal share of all UniProt members >= 0.90;
  at least 5 reviewed fungal proteins; at least 1 accession of ours in the family.
"""

import collections
import csv
import sys
from pathlib import Path

root, data = Path(sys.argv[1]), Path(sys.argv[2])


def rd(p):
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


acc = {r["accession"]: r for r in rd(data / "accession_to_panther.tsv")}
fam = {r["panther_family"]: r for r in rd(data / "family_summary.tsv")}
scan = rd(root / "docs/reports/data/pr6_pfam_scan_2026-10-22/scan_per_cluster.tsv")
missed = {
    r["cluster_id"]: r
    for r in rd(root / "docs/reports/data/kiss_baseline_2026-10-22/clusters_missed.tsv")
}


def w(path, header, rows):
    with path.open("w", newline="") as fh:
        t = csv.writer(fh, delimiter="\t", lineterminator="\n")
        t.writerow(header)
        t.writerows(rows)


rows = []
for s in scan:
    cid = s["cluster_id"]
    accs = missed[cid]["accessions"].split(",")
    real = [a for a in accs if not a.startswith("jgi:")]
    fams = collections.Counter()
    for a in real:
        for f in acc[a]["panther_family"].split(";"):
            if f:
                fams[f] += 1
    top = fams.most_common(1)[0][0] if fams else ""
    f = fam.get(top, {})
    rows.append(
        [
            cid,
            s["slice"],
            s["genes"],
            s["class"],
            s["other_pfam_families"].split(":")[0]
            if s["other_pfam_families"]
            else s["panel_models"],
            len(accs),
            len(real),
            len(accs) - len(real),
            ";".join(f"{k}x{v}" for k, v in fams.items()),
            f.get("name", ""),
            f.get("n_fungal", ""),
            f.get("broad_or_specific", ""),
            f.get("n_fine_clades_ge5", ""),
        ]
    )
w(
    data / "missed44_cluster_panther.tsv",
    [
        "cluster_id",
        "slice",
        "genes",
        "pfam_scan_class",
        "pfam_hit_first",
        "n_accessions",
        "n_uniprot_accessions",
        "n_jgi_keys_skipped",
        "panther_families",
        "panther_name",
        "n_fungal_in_family",
        "broad_or_specific",
        "n_fine_clades_ge5",
    ],
    rows,
)

# controls by slice
cov = collections.defaultdict(lambda: [0, 0])
for a in acc.values():
    if "controls112" in a["source_set"]:
        cov[a["slice"]][0] += 1
        cov[a["slice"]][1] += 1 if a["panther_family"] else 0
w(
    data / "controls_panther_coverage.tsv",
    ["slice", "n_controls_with_uniprot_accession", "n_with_panther_family"],
    [[k, v[0], v[1]] for k, v in sorted(cov.items())],
)

out = []
for f, r in fam.items():
    if r["fungal_list_streamed"] != "yes":
        continue
    n = int(r["n_fungal"])
    fr = float(r["frac_fungal"] or 0)
    flags = [
        n >= 500,
        int(r["n_fine_clades_ge5"]) >= 8,
        fr >= 0.90,
        int(r["n_reviewed_fungal"]) >= 5,
        int(r["n_our_accessions_in_family"]) >= 1,
    ]
    shipped = int(r["n_fungal_with_shipped_pfam"] or 0)
    out.append(
        [
            f,
            r["name"],
            n,
            r["n_reviewed_fungal"],
            r["n_fine_clades_ge5"],
            r["frac_fungal"],
            r["broad_or_specific"],
            r["n_our_accessions_in_family"],
            sum(flags),
            n - shipped,
            r["why_listed"],
        ]
    )
out.sort(key=lambda x: (-x[8], -x[2]))
w(
    data / "test_set_screen.tsv",
    [
        "panther_family",
        "name",
        "n_fungal",
        "n_reviewed_fungal",
        "n_fine_clades_ge5",
        "frac_fungal_of_all_taxa",
        "broad_or_specific",
        "n_our_accessions_in_family",
        "n_screen_criteria_met_of_5",
        "n_fungal_without_shipped_pfam",
        "why_listed",
    ],
    out,
)
