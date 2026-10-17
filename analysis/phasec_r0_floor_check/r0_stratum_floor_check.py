#!/usr/bin/env python3
"""S1 gate of the sorting-hat review plan (D16, Q4): do N-int plus N-sec alone meet the status floors
for the R0 call in each Phase C test species?

Read-only. Standard library only. Reads three files of the Phase C work directory and writes one TSV
(the path given by --out; default: next to this script). It changes nothing in the work directory.

Floors (docs/agent-tasks/COMMON-RULES.md rule 7, status.py): at least 20 negatives, at least 20
independent clusters of negatives, and a 95% interval half-width of at most 0.10 for specificity.

How the numbers are made (same rules as calibration/phasec.py, written out here):
- Proteins and clusters: rows of eval_table.tsv.gz with origin go, homology_only no, class neg, in the
  species' source; the cluster comes from clusters.tsv.gz (join on seq_sha256).
- A protein counts as "curated negative" when its stratum list holds N-int or N-sec and does not hold
  PM-TM. Proteins with both are counted apart (column n_both) and are NOT in the curated numbers.
- FPR per stratum comes from metrics.json (test_sets[key].truth[truth].metrics[stratum][variant]
  [candidate].fpr). The combined false-positive count is the sum of FPR * n over N-int and N-sec.
- Interval of the combined specificity: Wilson on the protein count, then widened with the Wilson
  interval on the cluster count (union), as widest_with_clusters does. The bootstrap interval of
  metrics.json cannot be combined across strata, so this interval is Wilson only. It is shown beside
  the pooled specificity that the current R0 status files carry (all negatives, as in phasec.py).
"""

import argparse
import csv
import gzip
import json
import math
import sys
from pathlib import Path

SETS = {
    "S1:Scer_SGD": "Scer_SGD",
    "S1:Calb_CGD": "Calb_CGD",
    "S3-Eurotiomycetes:Afum_ASPFU": "Afum_ASPFU",
    "S3-Eurotiomycetes:Anid_EMENI": "Anid_EMENI",
    "S3-Basidiomycota:Cneo_H99_GOA": "Cneo_H99_GOA",
    "S3-Basidiomycota:Umay_MYCMD": "Umay_MYCMD",
}
CURATED = ("N-int", "N-sec")
PM_TM = "PM-TM"
MIN_N, MIN_CLUSTERS, MAX_HALF_WIDTH = 20, 20, 0.10
Z = 1.959964


def wilson(k, n):
    p = k / n
    denom = 1 + Z * Z / n
    centre = (p + Z * Z / (2 * n)) / denom
    half = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def wilson_p(p, n):
    denom = 1 + Z * Z / n
    centre = (p + Z * Z / (2 * n)) / denom
    half = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def open_text(path):
    return (
        gzip.open(path, "rt", encoding="utf-8")
        if str(path).endswith(".gz")
        else open(path, encoding="utf-8")
    )


def split_list(cell):
    return {
        x.strip()
        for x in cell.replace("[", "").replace("]", "").replace("'", "").split(",")
        if x.strip()
    }


def load_clusters(path):
    cluster_of = {}
    with open_text(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if cluster_of.setdefault(r["seq_sha256"], r["cluster_id"]) != r["cluster_id"]:
                sys.exit(f"STOP: {r['seq_sha256']} is in two clusters")
    return cluster_of


def count_by_set(eval_table, cluster_of):
    """{source: {group: (proteins, clusters set)}}; groups: all, curated, pmtm, both."""
    wanted = set(SETS.values())
    acc = {s: {g: [0, set()] for g in ("all", "curated", "pmtm", "both")} for s in wanted}
    with open_text(eval_table) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        need = {"seq_sha256", "origin", "class", "homology_only", "source_ids", "stratum"}
        missing = need - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"STOP: eval table lacks columns {sorted(missing)}")
        for r in reader:
            if r["origin"] != "go" or r["homology_only"] != "no" or r["class"] != "neg":
                continue
            strata = split_list(r["stratum"])
            has_cur = bool(strata & set(CURATED))
            has_pm = PM_TM in strata
            group = (
                "both"
                if has_cur and has_pm
                else "curated"
                if has_cur
                else "pmtm"
                if has_pm
                else None
            )
            cluster = cluster_of.get(r["seq_sha256"])
            if cluster is None:
                sys.exit(f"STOP: {r['seq_sha256']} is not in the cluster file")
            for source in split_list(r["source_ids"]) & wanted:
                for g in ("all", group):
                    if g:
                        acc[source][g][0] += 1
                        acc[source][g][1].add(cluster)
    return acc


def rate(metrics, stratum, variant, candidate, key):
    cell = metrics.get(stratum, {}).get(variant, {}).get(candidate)
    if not cell or cell.get(key) is None or cell[key].get("value") is None:
        return None
    return cell[key]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--workdir", required=True, help="the phasec directory of the step1 work directory"
    )
    ap.add_argument("--variant", default="V-go")
    ap.add_argument("--candidate", default="R0")
    ap.add_argument("--truth", default="direct")
    ap.add_argument("--out", default=str(Path(__file__).with_name("r0_stratum_floor_check.tsv")))
    a = ap.parse_args()
    wd = Path(a.workdir)
    metrics_all = json.loads((wd / "metrics.json").read_text())["test_sets"]
    cluster_of = load_clusters(wd / "clusters.tsv.gz")
    counts = count_by_set(wd / "eval_table.tsv.gz", cluster_of)

    cols = [
        "set",
        "n_neg_all",
        "clusters_all",
        "spec_all",
        "lo_all",
        "hi_all",
        "halfw_all",
        "n_curated",
        "clusters_curated",
        "n_pmtm_only",
        "n_both",
        "fp_curated",
        "spec_curated",
        "lo_curated",
        "hi_curated",
        "halfw_curated",
        "n_ok",
        "clusters_ok",
        "halfw_ok",
        "FLOORS_MET",
        "metrics_n_N-int",
        "metrics_n_N-sec",
        "n_mismatch_note",
    ]
    rows = []
    for key, source in SETS.items():
        ts = metrics_all[key]
        t = ts["truth"][a.truth]
        m = t["metrics"]
        n_json = t["n"]
        c = counts[source]
        n_all, cl_all = c["all"][0], len(c["all"][1])
        fpr_all = rate(m, "all", a.variant, a.candidate, "fpr")
        spec_all = lo_all = hi_all = hw_all = ""
        if fpr_all:
            lo, hi = 1 - fpr_all["hi"], 1 - fpr_all["lo"]
            v = 1 - fpr_all["value"]
            if cl_all >= 1:
                wl, wh = wilson_p(v, cl_all)
                lo, hi = min(lo, wl, v), max(hi, wh, v)
            spec_all, lo_all, hi_all, hw_all = v, lo, hi, (hi - lo) / 2
        fp, n_cur, notes = 0.0, 0, []
        for s in CURATED:
            f = rate(m, s, a.variant, a.candidate, "fpr")
            ns = n_json.get(s, {}).get("neg", 0)
            if f and ns:
                fp += f["value"] * ns
                n_cur += ns
        mine = c["curated"][0]
        if n_cur != mine:
            notes.append(
                f"strata n from metrics.json ({n_cur}) differs from eval-table count ({mine})"
            )
        cl_cur = len(c["curated"][1])
        if n_cur >= 1:
            k = int(round(fp))
            lo, hi = wilson(k, n_cur)
            v = 1 - fp / n_cur
            lo, hi = 1 - hi, 1 - lo
            if cl_cur >= 1:
                wl, wh = wilson_p(v, cl_cur)
                lo, hi = min(lo, wl, v), max(hi, wh, v)
            hw = (hi - lo) / 2
            n_ok, cl_ok, hw_ok = n_cur >= MIN_N, cl_cur >= MIN_CLUSTERS, hw <= MAX_HALF_WIDTH + 1e-9
            met = n_ok and cl_ok and hw_ok
            row = [
                key,
                n_all,
                cl_all,
                spec_all,
                lo_all,
                hi_all,
                hw_all,
                n_cur,
                cl_cur,
                c["pmtm"][0],
                c["both"][0],
                round(fp, 3),
                v,
                lo,
                hi,
                hw,
                n_ok,
                cl_ok,
                hw_ok,
                met,
                n_json.get("N-int", {}).get("neg", ""),
                n_json.get("N-sec", {}).get("neg", ""),
                "; ".join(notes),
            ]
        else:
            row = [
                key,
                n_all,
                cl_all,
                spec_all,
                lo_all,
                hi_all,
                hw_all,
                0,
                0,
                c["pmtm"][0],
                c["both"][0],
                "",
                "",
                "",
                "",
                "",
                False,
                False,
                False,
                False,
                n_json.get("N-int", {}).get("neg", ""),
                n_json.get("N-sec", {}).get("neg", ""),
                "no N-int or N-sec negatives",
            ]
        rows.append(row)

    def fmt(x):
        return f"{x:.4f}" if isinstance(x, float) else str(x)

    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(x) for x in r) + "\n")
    print("\t".join(cols))
    for r in rows:
        print("\t".join(fmt(x) for x in r))
    print(f"\nwritten: {a.out}")
    print("FLOORS_MET for every set:", all(r[cols.index("FLOORS_MET")] is True for r in rows))


if __name__ == "__main__":
    main()
