#!/usr/bin/python3.12
"""Measure the hydrophobin calls against the truth tables (task H6 of the hydrophobin plan).

Variants: pfam_only (hydrophobin_domain), pfam_or_cys8 (hydrophobin_domain or cys8_pattern hit), cys8_only.
Per proteome: confusion counts, sensitivity and specificity with cluster-bootstrap 95% intervals, plus lists of
misses, called proteins with no label (candidate false positives) and rescue-only calls (cys8 hit, no Pfam call).
Negatives are ASSUMED (absent annotation). Reads plain or .gz inputs.

Usage: measure.py --dir analysis/hydrophobin_truth --out OUTDIR
"""

import argparse
import csv
import gzip
import math
import random
import sys
from pathlib import Path

Z = 1.959964


def wilson(k, n):
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def keep_rescue(outcomes, floor=5, min_lower=0.5):
    """Rule 6.3: outcomes are one per cluster of rescue-only calls: hydrophobin, not_hydrophobin, unresolved."""
    resolved = [o for o in outcomes if o in ("hydrophobin", "not_hydrophobin")]
    k = sum(1 for o in resolved if o == "hydrophobin")
    n = len(resolved)
    if n < floor:
        return "no evidence", k, n
    return ("keep" if wilson(k, n)[0] >= min_lower else "no evidence"), k, n


def confusion(labels, called):
    tp = sum(1 for p, y in labels.items() if y == 1 and p in called)
    fn = sum(1 for p, y in labels.items() if y == 1 and p not in called)
    fp = sum(1 for p, y in labels.items() if y == 0 and p in called)
    tn = sum(1 for p, y in labels.items() if y == 0 and p not in called)
    return {
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
        "sens": tp / (tp + fn) if tp + fn else float("nan"),
        "spec": tn / (tn + fp) if tn + fp else float("nan"),
    }


def boot_ci(labels, clusters, called, n_boot=2000, seed=1):
    """Cluster bootstrap 95% interval for sensitivity and specificity."""
    by = {}
    for p in labels:
        by.setdefault(clusters[p], []).append(p)
    ids = sorted(by)
    rng = random.Random(seed)
    sens, spec = [], []
    for _ in range(n_boot):
        tp = fn = fp = tn = 0
        for c in (rng.choice(ids) for _ in ids):
            for p in by[c]:
                y, hit = labels[p], p in called
                tp += y == 1 and hit
                fn += y == 1 and not hit
                fp += y == 0 and hit
                tn += y == 0 and not hit
        if tp + fn:
            sens.append(tp / (tp + fn))
        if tn + fp:
            spec.append(tn / (tn + fp))

    def q(v):
        v = sorted(v)
        return (v[int(0.025 * len(v))], v[int(0.975 * len(v)) - 1]) if v else (float("nan"),) * 2

    return q(sens), q(spec)


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_tsv(p):
    with _open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def headers(fasta):
    out = {}
    with _open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                k = line[1:].split()[0]
                out[k] = line[1:].strip()[:100]
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d, out = Path(a.dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    runs = read_tsv(d / "run_list.tsv")
    summary, lists = [], []
    for r in runs:
        name, w = r["name"], Path(r["workdir"])
        cal = d / "calibration" / name / "truth_all.tsv"
        truth = read_tsv(cal) if cal.exists() else []
        labels = {t["id"]: int(t["label"]) for t in truth}
        clusters = {t["id"]: t["cluster"] for t in truth}
        calls = {
            x["protein"]
            for x in read_tsv(w / "out_hyd/calls.long.tsv.gz")
            if x["call"] == "hydrophobin_domain" and x["value"] == "called"
        }
        cys = {x["id"]: x for x in read_tsv(w / "modules/cys8_pattern.tsv.gz")}
        cys_hit = {p for p, x in cys.items() if x["hit"] == "1"}
        cys_match = {p for p, x in cys.items() if x["pattern_match"] == "1"}
        hdr = headers(r["fasta"])
        variants = {"pfam_only": calls, "pfam_or_cys8": calls | cys_hit, "cys8_only": cys_hit}
        for v, called in variants.items():
            row = {"proteome": name, "variant": v, "called": len(called)}
            if labels:
                c = confusion(labels, called)
                (sl, sh), (pl, ph) = boot_ci(labels, clusters, called)
                row.update(
                    c,
                    sens_lo=sl,
                    sens_hi=sh,
                    spec_lo=pl,
                    spec_hi=ph,
                    n_clusters=len(set(clusters.values())),
                )
            summary.append(row)
        pos = {p for p, y in labels.items() if y == 1}
        for p in sorted(pos - calls):
            lists.append(
                (name, "miss_pfam", p, p in cys_hit, p in cys_match, cys.get(p, {}), hdr.get(p, ""))
            )
        for p in sorted(calls - pos):
            lists.append(
                (
                    name,
                    "pfam_call_not_positive",
                    p,
                    p in cys_hit,
                    p in cys_match,
                    cys.get(p, {}),
                    hdr.get(p, ""),
                )
            )
        for p in sorted(cys_hit - calls):
            lists.append((name, "rescue_only", p, True, True, cys.get(p, {}), hdr.get(p, "")))
        for p in sorted(cys_match - cys_hit - calls):
            lists.append(
                (name, "pattern_no_signal_peptide", p, False, True, cys.get(p, {}), hdr.get(p, ""))
            )
    keys = [
        "proteome",
        "variant",
        "called",
        "tp",
        "fn",
        "fp",
        "tn",
        "sens",
        "sens_lo",
        "sens_hi",
        "spec",
        "spec_lo",
        "spec_hi",
        "n_clusters",
    ]
    with open(out / "summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(summary)
    with open(out / "review_lists.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "proteome",
                "list",
                "id",
                "cys8_hit",
                "pattern_match",
                "length",
                "n_cys",
                "sets",
                "header",
            ]
        )
        for name, kind, p, h, pm, c, hd in lists:
            w.writerow(
                [
                    name,
                    kind,
                    p,
                    int(h),
                    int(pm),
                    c.get("length", ""),
                    c.get("n_cys", ""),
                    c.get("sets", ""),
                    hd,
                ]
            )
    print(f"wrote {out}/summary.tsv and review_lists.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
