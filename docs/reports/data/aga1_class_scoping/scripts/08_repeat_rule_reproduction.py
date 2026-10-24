#!/usr/bin/env python3
"""Recompute the existing repeat call on the committed truth_v3 tables from the locally re-run detectors, so that the
numbers in the report are reproduced rather than copied. Also the exploratory z_seq >= Z extra rule (zseq_rule_scan.py).
Writes repeat_rule_reproduction.tsv. Stdlib only. Checks: S288C 0.435 / 0.987 and C. albicans 0.462 / 0.976 are the
statuses in docs/class_descriptions.yaml; agreement with the committed per-protein calls of curated_rows_mapped.tsv.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import OUT, RCT, read_tsv, truth_tables, write_tsv  # noqa: E402


def main():
    F = {r["accession"]: r for r in read_tsv(OUT / "features.tsv")}
    rows = []
    for prot, truth in truth_tables().items():
        for z in (None, 5.0, 4.5, 4.0, 3.5):
            tp = fn = fp = tn = 0
            gained, added = [], []
            for t in truth:
                f = F[t["accession"]]
                base = int(f["call_repeat"])
                call = base or (
                    z is not None and f["r14_z_seq"] != "" and float(f["r14_z_seq"]) >= z
                )
                lab = int(t["label"])
                tp += lab and call
                fn += lab and not call
                fp += (not lab) and call
                tn += (not lab) and not call
                if z is not None and call and not base:
                    (gained if lab else added).append(t["gene"] or t["accession"])
            rows.append({"truth": prot, "z_seq_min": "none" if z is None else z, "positives": tp + fn, "sensitivity": round(tp / (tp + fn), 3), "negatives": fp + tn, "false_positives": fp, "specificity": round(tn / (tn + fp), 3), "gained_positives": ";".join(gained), "n_added_false_positives": len(added)})  # fmt: skip
    write_tsv(OUT / "repeat_rule_reproduction.tsv", rows, list(rows[0]))
    cm = read_tsv(RCT / "curated_rows_mapped.tsv")
    agree = sum(1 for x in cm if x["accession"] in F and (x["tandem_repeat_protein"] == "called") == bool(int(F[x["accession"]]["call_repeat"])))  # fmt: skip
    n = sum(1 for x in cm if x["accession"] in F)
    print(
        f"curated_rows_mapped agreement {agree}/{n} (rows without a local feature row: {len(cm) - n})"
    )


if __name__ == "__main__":
    main()
