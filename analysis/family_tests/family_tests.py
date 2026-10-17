#!/usr/bin/python3.12
"""Test every family of family_table.tsv against curated positives, hard negatives and the scan proteomes.

Inputs (made by hmmsearch --cut_ga with the family HMMs, see README.md): domtbl files for the adhesin positives,
the shared hard negatives and the scan proteomes. Output: family_tests.tsv (one row per family) in --out.
This is a descriptive test of each family: counts, no estimated sensitivity or specificity (see README.md).
Usage: family_tests.py --work _workdir/family_tests --out docs/reports/data/sorting_hat/family_tests
"""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def hits(path):
    """{family accession without version: set of target ids} from a hmmsearch --domtblout file."""
    out = defaultdict(set)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        out[f[4].split(".")[0]].add(f[0])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    work, out = Path(a.work), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    fam = list(
        csv.DictReader(open(ROOT / "data/sorting_hat/family_table.tsv", newline=""), delimiter="\t")
    )
    ev = {}
    for r in csv.DictReader(
        open(ROOT / "data/curated/adhesins/adhesins.tsv", newline=""), delimiter="\t"
    ):
        if r["cls"] == "adhesin":
            ev[r["accession"]] = r["evidence_level"]
    neg = {}
    for r in csv.DictReader(
        open(ROOT / "data/controls/hard-negatives-shared/controls.tsv", newline=""), delimiter="\t"
    ):
        neg[r["accession"]] = r["evidence_level"]
    cl = {}
    for r in csv.DictReader(
        open(ROOT / "data/controls/hard-negatives-shared/clusters.tsv", newline=""), delimiter="\t"
    ):
        if r["set"] == "hard_negative":
            cl[r["accession"]] = r["cluster_id"]
    hp, hn = hits(work / "adhesin_positives.domtbl"), hits(work / "hard_negatives.domtbl")
    scan = defaultdict(set)  # family -> {proteome|protein}
    sp, rep = (
        set(),
        set(),
    )  # "proteome|protein" with signal_peptide_protein[R0] / tandem_repeat_protein called
    for lt in sorted((ROOT / "_workdir/fungi_scan").glob("*/out/calls.long.tsv.gz")):
        name = lt.parts[-3]
        with gzip.open(lt, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["value"] != "called":
                    continue
                if r["call"] == "signal_peptide_protein":
                    sp.add(f"{name}|{r['protein']}")
                elif r["call"] == "tandem_repeat_protein":
                    rep.add(f"{name}|{r['protein']}")
    for dom in sorted((ROOT / "_workdir/fungi_scan").glob("*/raw/pfam/domtbl.txt")):
        name = dom.parts[-4]
        for acc, ids in hits(dom).items():
            scan[acc] |= {f"{name}|{i}" for i in ids}
    rows = []
    for f in fam:
        acc = f["pfam_acc"]
        pos = {e: len([p for p in hp.get(acc, ()) if ev.get(p) == e]) for e in ("E1", "E2", "E3")}
        npos = {e: sum(1 for v in ev.values() if v == e) for e in ("E1", "E2", "E3")}
        n1 = {p for p in hn.get(acc, ()) if neg.get(p) == "N1"}
        n2 = {p for p in hn.get(acc, ()) if neg.get(p) == "N2"}
        s = scan.get(acc, set())
        rows.append(
            {
                "pfam_acc": acc,
                "name": f["name"],
                "active": f["active"],
                "pos_E1_hit": f"{pos['E1']}/{npos['E1']}",
                "pos_E2_hit": f"{pos['E2']}/{npos['E2']}",
                "pos_E3_hit": f"{pos['E3']}/{npos['E3']}",
                "neg_N1_hit": len(n1),
                "neg_N1_clusters_hit": len({cl[p] for p in n1 if p in cl}),
                "neg_N2_hit": len(n2),
                "scan_proteins": len(s),
                "scan_proteomes": len({x.split("|")[0] for x in s}),
                "scan_with_signal_peptide": len(s & sp),
                "scan_with_repeat_call": len(s & rep),
            }
        )
    with open(out / "family_tests.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} families")


if __name__ == "__main__":
    main()
