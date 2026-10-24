"""Join the converter table with the readable columns of docs/reports/data/aga1_class_scoping/proteins.tsv (stdlib).

Usage: python3 -I build_record_table.py --predictions CONVERTED.tsv --out RECORD.tsv   (refuses an existing --out)
`role_group` is `role`, except that class_positive is split by tier (class_positive_A, class_positive_B).
Fails loudly (FATAL:, exit 1) if the two tables do not hold exactly the same accessions.
"""

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROTEINS = ROOT / "docs/reports/data/aga1_class_scoping/proteins.tsv"
KEEP = [
    "gene", "organism", "clade_coarse", "role", "pos_tier", "curated_class",
    "uniprot_gpi_feature", "uniprot_signal_feature",
]  # fmt: skip


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build(pred_rows, protein_rows):
    prot = {r["accession"]: r for r in protein_rows}
    if len(prot) != len(protein_rows):
        raise SystemExit("FATAL: duplicate accession in proteins.tsv")
    ids = [r["accession"] for r in pred_rows]
    if len(set(ids)) != len(ids) or set(ids) != set(prot):
        raise SystemExit("FATAL: prediction table and proteins.tsv do not hold the same accessions")
    out = []
    for r in pred_rows:
        p = prot[r["accession"]]
        if p["role"] == "class_positive" and p["pos_tier"] not in ("A", "B"):
            raise SystemExit(f"FATAL: class_positive without tier A/B: {r['accession']}")
        group = f"class_positive_{p['pos_tier']}" if p["role"] == "class_positive" else p["role"]
        head = {"accession": r["accession"], **{k: p[k] for k in KEEP[:3]}, "role_group": group}
        head.update({k: p[k] for k in KEEP[3:]})
        out.append({**head, **{k: v for k, v in r.items() if k != "accession"}})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if Path(a.out).exists():
        print(f"FATAL: {a.out} exists; nothing is overwritten", file=sys.stderr)
        return 1
    try:
        rows = build(read(a.predictions), read(PROTEINS))
    except SystemExit as e:
        print(e, file=sys.stderr)
        return 1
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {a.out}: {len(rows)} proteins")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
