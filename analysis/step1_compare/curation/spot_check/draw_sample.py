"""Draw the owner's spot-check sample from the two curation passes (seeded, reproducible).

Pool: all draft rows except the 3 H99 pass-2 rows the independent reviewer rejected.
The sheet hides the reviewer verdicts so the owner judges the rows blind.
"""

import csv
import random

SEED = 20261003
N = 15
REJECTED_H99 = {22, 24, 33}  # file line numbers in h99_pass2/draft_rows.tsv

pool = []
for tag, path in (("H99", "h99_pass2/draft_rows.tsv"), ("UM", "umay_pilot/draft_rows.tsv")):
    with open(path) as fh:
        for line_no, row in enumerate(csv.DictReader(fh, delimiter="\t"), 2):
            if tag == "H99" and line_no in REJECTED_H99:
                continue
            row["row_id"] = f"{tag}:{line_no}"
            pool.append(row)

random.Random(SEED).shuffle(pool)
sample = sorted(pool[:N], key=lambda r: r["row_id"])
cols = [
    "row_id", "species", "symbol", "uniprot_accession", "go_term", "evidence_code",
    "expected_label", "selected_by_predictor", "pmid", "evidence_note",
]  # fmt: skip
with open("spot_check/spot_check_sheet.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(cols + ["owner_verdict", "owner_comment"])
    for r in sample:
        w.writerow([r[c] for c in cols] + ["", ""])
print(len(pool), len(sample))
