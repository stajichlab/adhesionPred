"""Join (no API calls): the three Kiss accessions of the two families with the committed Pfam scan of step 6.0 and the
PANTHER subfamily found by s1. Usage: python3 -I s3_kiss_mapping.py REPO_ROOT OUT_DIR"""

import csv
import sys
from pathlib import Path

root, out = Path(sys.argv[1]), Path(sys.argv[2])


def rows(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


scan = {
    r["cluster_id"]: r
    for r in rows(root / "docs/reports/data/pr6_pfam_scan_2026-10-22/scan_per_cluster.tsv")
}
res = []
for r in rows(out / "kiss_in_family.tsv"):
    s = scan[r["kiss_cluster"]]
    res.append(
        [
            r["kiss_cluster"],
            r["gene"],
            r["accession"],
            r["organism"],
            r["panther_family"],
            r["subfamily"],
            r["subfamily_name"],
            s["class"],
            s["proposal"],
            s["weak_families"] if "weak_families" in s else "",
            s["other_pfam_families"],
        ]
    )
with open(out / "kiss_mapping.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(
        [
            "kiss_cluster",
            "gene",
            "accession",
            "organism",
            "panther_family",
            "subfamily",
            "subfamily_name",
            "pfam_scan_class",
            "pfam_scan_proposal",
            "weak_families",
            "other_pfam_hits",
        ]
    )
    w.writerows(res)
print(res)
