"""Margin tables for the custom Ecm33 `all` model, from COMMITTED tables only (no search is run).

Usage: python3 -I analysis/wall_fusion_truth/ecm33_margin_tables.py OUT_DIR
Writes small TSVs to OUT_DIR and prints the numbers used by docs/reports/2026-10-24-ecm33-margin-report.md.
The cutoff is read, never changed.
"""

import csv
import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "docs/reports/data/pr6_ecm33_build_2026-10-23"
EVAL = ROOT / "docs/reports/data/pr6_ecm33_independent_eval_2026-10-24"
CUT = ROOT / "data/custom_models/ecm33_all_2026-10-23/cutoff.json"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def heldout_margins():
    """Per unit and mode: held-out custom scores against that unit's own cutoff."""
    rows = []
    for unit_dir in sorted((BUILD / "ecm33_build").iterdir()):
        if unit_dir.name == "all":
            continue
        for mode in ("cluster30", "none"):
            d = unit_dir / mode
            cutoff = json.loads((d / "cutoff.json").read_text())["cutoff"]
            scores = [float(r["score"]) for r in read_tsv(d / "hits_custom.tsv")]
            rows.append(
                {
                    "unit": unit_dir.name,
                    "mode": mode,
                    "n_hits": len(scores),
                    "cutoff": cutoff,
                    "min_score": min(scores),
                    "median_score": round(statistics.median(scores), 1),
                    "min_margin": round(min(scores) - cutoff, 1),
                }
            )
    return rows


def proteome_hits():
    cutoff = json.loads(CUT.read_text())["cutoff"]
    rows = read_tsv(BUILD / "ecm33_build/all/none/proteome_hits.tsv")
    return cutoff, [(r["proteome"], r["id"], float(r["score"])) for r in rows]


def eval_scores():
    cutoff = json.loads(CUT.read_text())["cutoff"]
    rows = read_tsv(EVAL / "search/hits_custom.tsv")
    n_clusters = len(read_tsv(EVAL / "extract/test_clusters.tsv"))
    return cutoff, sorted(float(r["score"]) for r in rows), n_clusters


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    hm = heldout_margins()
    with open(out / "heldout_margins.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(hm[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(hm)
    cutoff, ph = proteome_hits()
    above = [(p, i, s) for p, i, s in ph if s >= cutoff]
    below = [(p, i, s) for p, i, s in ph if s < cutoff]
    with open(out / "proteome_margins.tsv", "w", newline="") as fh:
        fh.write("proteome\tid\tscore\tmargin_to_cutoff\tside\n")
        for p, i, s in ph:
            fh.write(
                f"{p}\t{i}\t{s}\t{round(s - cutoff, 1)}\t{'at_or_above' if s >= cutoff else 'below'}\n"
            )
    c, sc, n = eval_scores()
    with open(out / "eval_histogram.txt", "w") as fh:
        for lo in range(25, 500, 25):
            k = sum(lo <= s < lo + 25 for s in sc)
            fh.write(f"{lo:3d}-{lo + 24:3d} bits | {'#' * k} {k}\n")
    print("cutoff", c, "proteome hits at/above", len(above), "below", len(below))
    print("eval hits", len(sc), "clusters", n, "lowest", sc[:8])
    print("within10", sum(s < c + 10 for s in sc), "within20", sum(s < c + 20 for s in sc))
    print("median", statistics.median(sc), "max", sc[-1])


if __name__ == "__main__":
    main(sys.argv[1])
