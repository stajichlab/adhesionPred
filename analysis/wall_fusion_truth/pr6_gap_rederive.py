#!/usr/bin/env python3
"""Re-derive the verdicts of a PR 6 gap-search round from its COMMITTED tables under the fixed clade rule.

Round 2 (UniRef50, job 30262455) counted the '(no clade)' bucket (fungal hits without an assigned fine clade) as a
clade. The rule now counts named clades only (pr6_gap_summary.training_rule). This script reads the committed run record
(`gap_homologs.tsv.gz`, `gap_per_query.tsv`), applies the fixed rule and writes a small table next to the old verdicts.
It never edits the run record. Standard library only; run with `python3 -I`.

    python3 -I analysis/wall_fusion_truth/pr6_gap_rederive.py RECORD_DIR OUT_DIR
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
HEADER = [
    "query",
    "genes",
    "method",
    "n_outside_seq",
    "n_unassigned_seq",
    "n_clades_committed",
    "n_named_clades_fixed",
    "rule_met_committed",
    "rule_met_fixed",
    "strict_committed",
    "strict_fixed",
    "status_committed",
    "status_fixed",
    "changed",
]


def load_gap():
    spec = importlib.util.spec_from_file_location("pr6_gap_summary", HERE / "pr6_gap_summary.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pr6_gap_summary"] = mod
    spec.loader.exec_module(mod)
    return mod


def read_tsv(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def rederive(record: Path) -> list[list]:
    gap = load_gap()
    per_query = read_tsv(record / "gap_per_query.tsv")
    hits = read_tsv(record / "gap_homologs.tsv.gz")
    out = []
    for q in per_query:
        rules = {}
        for method in ("phmmer", "jackhmmer"):
            rows = [
                {
                    "where": h["where"],
                    "clade_fine": h["clade_fine"],
                    # the species is not kept in the record; the taxid is a stand-in, only the status is used here
                    "species_taxid": h["taxid"],
                }
                for h in hits
                if h["query"] == q["query"] and h["method"] == method
            ]
            rules[method] = gap.training_rule(rows)
        status_fixed, _ = gap.target_status(rules)
        for method in ("phmmer", "jackhmmer"):
            r = rules[method]
            unassigned = r["per_clade"].get(gap.NO_CLADE, 0)
            out.append(
                [
                    q["query"],
                    q["genes"],
                    method,
                    r["n_outside_seq"],
                    unassigned,
                    int(q[f"{method}_outside_clades"]),
                    r["n_outside_clades"],
                    q[f"{method}_rule_met"],
                    "yes" if r["rule_met"] else "no",
                    q[f"{method}_strict_5_in_2_clades"],
                    "yes" if r["strict_per_clade_met"] else "no",
                    q["status"],
                    status_fixed,
                    "yes" if status_fixed != q["status"] else "no",
                ]
            )
    return out


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        print(__doc__)
        return 2
    record, out = Path(argv[0]), Path(argv[1])
    if out.exists() and any(out.iterdir()):
        print(f"FATAL: {out} exists and is not empty; nothing is overwritten", file=sys.stderr)
        return 1
    rows = rederive(record)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "rederived_verdicts.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(rows)
    print(f"wrote {out / 'rederived_verdicts.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
