#!/usr/bin/python3.12
"""Apply curation_rubric.md to curation/curation_evidence.tsv. Tier T4: curator judgement from evidence, not owner-reviewed.

Rules (in order): hydrophobin if C, or A with B or D, unless E also holds and C does not (then unresolved: conflict);
not_hydrophobin if E and not C; otherwise unresolved. Writes curator_decisions.tsv. Usage: apply_rubric.py --dir analysis/hydrophobin_truth
"""

import argparse
import csv
import sys
from pathlib import Path


def decide(r):
    a, b, c, dd, e = (r[f"line_{x}"] == "yes" for x in "ABCDE")
    if c:
        return "hydrophobin", "C: named a hydrophobin in a paper or curated record"
    if a and (b or dd):
        if e:
            return (
                "unresolved",
                "conflict: hydrophobin evidence (A with B or D) and another function (E)",
            )
        return "hydrophobin", "A with " + " and ".join(x for x, v in (("B", b), ("D", dd)) if v)
    if e:
        return "not_hydrophobin", "E: evidence of another function"
    lines = [x for x, v in zip("ABCDE", (a, b, c, dd, e), strict=True) if v]
    return "unresolved", "lines holding: " + (",".join(lines) or "none")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    rows = list(csv.DictReader(open(d / "curation/curation_evidence.tsv"), delimiter="\t"))
    out = []
    for r in rows:
        dec, why = decide(r)
        out.append(
            {
                "proteome": r["proteome"],
                "id": r["id"],
                "kind": r["kind"],
                "tier": "T4",
                "decision": dec,
                "why": why,
                "A": r["line_A"],
                "B": r["line_B"],
                "C": r["line_C"],
                "D": r["line_D"],
                "E": r["line_E"],
                "E_detail": r["E_detail"][:200],
                "uniprot_name": r["uniprot_name"],
                "note": "not owner-reviewed",
            }
        )
    with open(d / "curator_decisions.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    from collections import Counter

    print(Counter((o["kind"], o["decision"]) for o in out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
