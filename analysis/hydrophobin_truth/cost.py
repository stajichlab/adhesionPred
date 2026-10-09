#!/usr/bin/python3.12
"""L7a: cost of the relaxed level from the module tables (no core run needed).

Extra calls = hydrophobin_relaxed hit not called by pfam_hydrophobin and not by pfam_hsba and not labelled (T1/T2/T3 mapped or LP).
Reported per 10,000 proteins for the test proteomes (out of sample for the cutoff) and the tuning proteomes (in sample); HsbA overlaps
(E11) and labelled extra calls in their own columns. Usage: cost.py --dir analysis/hydrophobin_truth --out OUTDIR
"""

import argparse
import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def classify(relaxed, strict, hsba, labelled):
    strict_also = sorted(relaxed & strict)
    rest = relaxed - strict
    hs = sorted(rest & hsba)
    rest2 = rest - hsba
    lab = sorted(rest2 & labelled)
    unl = sorted(rest2 - labelled)
    return {
        "strict_also": strict_also,
        "hsba_overlap": hs,
        "labelled_extra": lab,
        "unlabelled_extra": unl,
    }


def per_10000(n, size):
    return 0.0 if size == 0 else n * 10000.0 / size


def read_module(path):
    with gzip.open(path, "rt") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d, out = Path(a.dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pd = load("prefreeze_data")
    runs = list(csv.DictReader(open(d / "run_list.tsv"), delimiter="\t"))
    split = {
        r["proteome"]: r["part"]
        for r in csv.DictReader(open(d / "proteome_split.tsv"), delimiter="\t")
    }
    pmap = list(csv.DictReader(open(d / "proteome_map.tsv"), delimiter="\t"))
    lp = {
        k.split("|")[1]: v for k, v in pd.read_fasta(d / "literature/jensen_sequences.faa").items()
    }
    thr = json.load(open(d / "freeze.json"))["thresholds"][
        "cost_limit_per_10000_in_each_test_proteome"
    ]
    rows, extras = [], []
    for r in runs:
        name, w = r["name"], Path(r["workdir"])
        seqs = pd.read_fasta(r["fasta"])
        mod = read_module(w / "modules/hydrophobin_relaxed.tsv.gz")
        relaxed = {i for i, x in mod.items() if x["hit"] == "1"}
        strict = pd.read_hits(w / "modules/pfam_hydrophobin.tsv.gz")
        hsba = pd.read_hits(w / "modules/pfam_hsba.tsv.gz")
        labelled = {x["proteome_id"] for x in pmap if x["proteome"] == name} | pd.exact_matches(
            lp, seqs
        )
        c = classify(relaxed, strict, hsba, labelled)
        n = len(seqs)
        rate = per_10000(len(c["unlabelled_extra"]), n)
        rows.append(
            {
                "proteome": name,
                "part": split[name],
                "proteins": n,
                "relaxed_hits": len(relaxed),
                "strict_also": len(c["strict_also"]),
                "hsba_overlap": len(c["hsba_overlap"]),
                "labelled_extra": len(c["labelled_extra"]),
                "unlabelled_extra": len(c["unlabelled_extra"]),
                "unlabelled_per_10000": round(rate, 2),
                "limit": thr,
                "within_limit": rate <= thr,
            }
        )
        for kind in ("hsba_overlap", "labelled_extra", "unlabelled_extra"):
            for i in c[kind]:
                x = mod[i]
                extras.append(
                    [name, i, kind, x["score"], x["n_cys"], x["query"], len(seqs[i]), seqs[i]]
                )
    with open(out / "relaxed_cost.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(d / "extra_calls.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["proteome", "id", "kind", "relaxed_score", "n_cys", "best_model", "length", "sequence"]
        )
        w.writerows(extras)
    for r in rows:
        print(
            r["proteome"],
            r["part"],
            "n",
            r["proteins"],
            "relaxed",
            r["relaxed_hits"],
            "strict_also",
            r["strict_also"],
            "hsba",
            r["hsba_overlap"],
            "labelled_extra",
            r["labelled_extra"],
            "unlabelled",
            r["unlabelled_extra"],
            "per10k",
            r["unlabelled_per_10000"],
            "OK" if r["within_limit"] else "OVER",
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
