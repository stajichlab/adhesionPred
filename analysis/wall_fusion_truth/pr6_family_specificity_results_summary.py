"""Combine the two PR 6 family-specificity runs into one table (stdlib, offline, deterministic).

Reads `<record>/run1` and `<record>/run2` (tables copied back from HPCC) and writes `<record>/summary.tsv`:
one row per family with the counts of `family_specificity.tsv` plus three derived columns. A proteome hit whose
protein id is one of the family's own Kiss members is counted in `own_in_proteomes` (the member is itself
inside a kept proteome); `kept_hits_excl_own` leaves those out. Nothing here is a sensitivity or a specificity.

Usage: python3 -I analysis/wall_fusion_truth/pr6_family_specificity_results_summary.py [record_dir]
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RECORD = ROOT / "docs/reports/data/pr6_family_specificity_results_2026-10-24"
PER = 10_000
COLUMNS = [
    "run",
    "pfam_acc",
    "name",
    "in_family_table",
    "n_own_members",
    "own_hits",
    "neg_hits_n1",
    "neg_hits_n2",
    "other_positive_hits_adhesion",
    "other_positive_hits_wall",
    "other_positive_hits_fusion",
    "kept_hits_total",
    "kept_per10k",
    "own_in_proteomes",
    "kept_hits_excl_own",
    "kept_per10k_excl_own",
]


def _read(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _protein_id(header_id: str) -> str:
    parts = header_id.split("|")
    return parts[1] if len(parts) >= 3 else header_id


def summarize(record: Path = RECORD) -> list[dict]:
    out = []
    for run in ("run1", "run2"):
        own = {
            (r["pfam_acc"], r["accession"])
            for r in _read(record / run / "positives_hits.tsv")
            if r["own_or_other"] == "own"
        }
        in_prot: dict[str, int] = {}
        for r in _read(record / run / "proteome_hits.tsv"):
            if (r["pfam_acc"], _protein_id(r["protein"])) in own:
                in_prot[r["pfam_acc"]] = in_prot.get(r["pfam_acc"], 0) + 1
        for r in _read(record / run / "family_specificity.tsv"):
            total = int(r["kept_hits_total"])
            n_prot = int(r["kept_n_proteins"])
            excl = total - in_prot.get(r["pfam_acc"], 0)
            row = {k: r[k] for k in COLUMNS if k in r}
            row.update(
                run=run,
                own_in_proteomes=in_prot.get(r["pfam_acc"], 0),
                kept_hits_excl_own=excl,
                kept_per10k_excl_own=f"{PER * excl / n_prot:.2f}",
            )
            out.append(row)
    out.sort(key=lambda x: (x["run"], -float(x["kept_per10k"]), x["pfam_acc"]))
    return out


def main(argv: list[str]) -> None:
    record = Path(argv[1]) if len(argv) > 1 else RECORD
    rows = summarize(record)
    with (record / "summary.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in COLUMNS])
    print(f"wrote {record / 'summary.tsv'}: {len(rows)} families")


if __name__ == "__main__":
    main(sys.argv)
