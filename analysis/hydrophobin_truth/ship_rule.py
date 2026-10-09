#!/usr/bin/python3.12
"""L8b: the ship rule of the hydrophobin extended-level spec (6.2 to 6.5) from stored results. Writes ship_decision.json.

Conditions: recall (L5), cost in every test proteome (L7a), hard-negative rate per group (L5), precision of the extra (relaxed-only) calls
by cluster with T4 outcomes (not owner-reviewed). Needs mmseqs on PATH. Usage: ship_rule.py --dir analysis/hydrophobin_truth
"""

import argparse
import csv
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
Z = 1.959964


def wilson_lower(k, n):
    if n == 0:
        return 0.0
    p = k / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return max(0.0, c - h)


def cluster_outcomes(clusters, decisions):
    by = {}
    for p, c in clusters.items():
        by.setdefault(c, set()).add(decisions[p])
    return {c: (next(iter(s)) if len(s) == 1 else "unresolved") for c, s in by.items()}


def precision_rule(outcomes, floor=5, min_lower=0.5):
    resolved = [o for o in outcomes.values() if o in ("hydrophobin", "not_hydrophobin")]
    k = sum(1 for o in resolved if o == "hydrophobin")
    n = len(resolved)
    lower = wilson_lower(k, n)
    if n < floor:
        return "no evidence", k, n, lower
    return ("pass" if lower >= min_lower else "fail"), k, n, lower


def ship(c):
    return bool(c["recall"] and c["cost"] and c["hard_negatives"] and c["precision"] == "pass")


def cluster_seqs(seqs):
    names = list(seqs)
    with tempfile.TemporaryDirectory() as tmp:
        faa = Path(tmp) / "x.faa"
        with open(faa, "w") as f:
            for i, k in enumerate(names):
                f.write(f">s{i}\n{seqs[k]}\n")
        pre = Path(tmp) / "c"
        subprocess.run(
            [
                "mmseqs",
                "easy-cluster",
                str(faa),
                str(pre),
                str(Path(tmp) / "w"),
                "--min-seq-id",
                "0.3",
                "-c",
                "0.5",
            ],
            check=True,
            capture_output=True,
        )
        cl = {}
        for line in open(f"{pre}_cluster.tsv"):
            rep, mem = line.split()
            cl[names[int(mem[1:])]] = names[int(rep[1:])]
    return cl


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--reports", default="docs/reports/data/sorting_hat/hydrophobin_ext")
    a = ap.parse_args()
    d, rep = Path(a.dir), Path(a.reports)
    l5 = json.load(open(rep / "l5_results.json"))
    cost = list(csv.DictReader(open(rep / "relaxed_cost.tsv"), delimiter="\t"))
    test_cost = [r for r in cost if r["part"] == "test"]
    cost_ok = all(r["within_limit"] == "True" for r in test_cost)
    hn = l5["hard_negatives_test"]
    hn_ok = all(v["relaxed_passes"] is not False for v in hn.values())
    extras = [
        r
        for r in csv.DictReader(open(d / "extra_calls.tsv"), delimiter="\t")
        if r["kind"] == "unlabelled_extra"
    ]
    dec = {
        (r["proteome"], r["id"]): r["decision"]
        for r in csv.DictReader(open(d / "curator_decisions.tsv"), delimiter="\t")
    }
    seqs = {f"{r['proteome']}|{r['id']}": r["sequence"] for r in extras}
    cl = cluster_seqs(seqs)
    decisions = {k: dec[tuple(k.split("|", 1))] for k in seqs}
    outcomes = cluster_outcomes(cl, decisions)
    prec, k, n, lower = precision_rule(outcomes)
    cond = {
        "recall": bool(l5["relaxed"]["recall_rule_pass"]),
        "cost": cost_ok,
        "hard_negatives": hn_ok,
        "precision": prec,
    }
    result = {
        "ship": ship(cond),
        "conditions": cond,
        "cost_test_proteomes": {
            r["proteome"]: [r["unlabelled_per_10000"], r["within_limit"]] for r in test_cost
        },
        "precision": {
            "basis": "relaxed-only unlabelled extra calls, clustered across proteomes (MMseqs2 30%, 0.5); outcomes are T4 (not owner-reviewed)",
            "extra_calls": len(extras),
            "clusters": len(outcomes),
            "hydrophobin": sum(1 for o in outcomes.values() if o == "hydrophobin"),
            "not_hydrophobin": sum(1 for o in outcomes.values() if o == "not_hydrophobin"),
            "unresolved": sum(1 for o in outcomes.values() if o == "unresolved"),
            "resolved": n,
            "wilson_lower_bound": round(lower, 3),
        },
        "note": "hydrophobin_extended is added to categories.yaml only if ship is true",
    }
    json.dump(result, open(d / "ship_decision.json", "w"), indent=1)
    print(json.dumps(result, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
