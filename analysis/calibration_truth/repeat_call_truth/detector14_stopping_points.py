#!/usr/bin/env python3
"""Where does detector 14 stop for the truth positives that it misses? (diagnosis only)

For each truth positive that the combined repeat call misses, run the stages of
analysis/cocci_repeats/14_repeat_detect_general.py directly: the best period, its sequence-level
z-score (gate Z_MIN) and composition-corrected score (gate MIN_SCORE), and the region that
periodic_region() extracts at the default cut. Compare the region with the span of the UniProt
`Repeat` features. Nothing is changed or written except the output table.

Usage: detector14_stopping_points.py --detector FILE --work WORKDIR_ROOT --truth DIR --prefix truth_v2 --out FILE
"""

import argparse
import csv
import importlib.util
import json
import sys
import urllib.request

SPECIES = ("Scer_S288C", "Calb_SC5314")


def load(path):
    spec = importlib.util.spec_from_file_location("detector14", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).rstrip("*")
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf).rstrip("*")
    return seqs


def uniprot_span(acc):
    d = json.load(
        urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60)
    )
    reps = [
        (f["location"]["start"]["value"], f["location"]["end"]["value"])
        for f in d["features"]
        if f["type"] == "Repeat"
    ]
    return len(reps), (max(b for _, b in reps) - min(a for a, _ in reps) + 1) / d["sequence"][
        "length"
    ]


def called(row):
    return (
        bool(row)
        and int(row["rep_period"]) > 0
        and float(row["rep_coverage"]) >= 0.25
        and float(row["rep_n_copies"]) >= 2.5
    )


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--detector", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--prefix", default="truth_v2")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = load(a.detector)
    out = []
    for prot in SPECIES:
        seqs = fasta(f"{a.work}/{prot}.faa")
        r02 = {
            r["protein"]: r
            for r in csv.DictReader(
                open(f"{a.work}/{prot}/raw/repeats/repeat02.tsv"), delimiter="\t"
            )
        }
        r14 = {
            r["protein"]: r
            for r in csv.DictReader(
                open(f"{a.work}/{prot}/raw/repeats/repeat14.tsv"), delimiter="\t"
            )
        }
        for t in csv.DictReader(open(f"{a.truth}/{a.prefix}.{prot}.annotated.tsv"), delimiter="\t"):
            if t["label"] != "1" or called(r02.get(t["id"])) or called(r14.get(t["id"])):
                continue
            s = seqs[t["id"]]
            if len(s) < 80:
                continue
            idx = d.encode(s)
            exp = d.expected_rate(d.background(idx), d.SIM)
            scores = d.period_scan(idx, d.SIM, exp)
            if not scores:
                continue
            raw = {q: v[1] for q, v in scores.items()}
            zs = d.z_scan(raw)
            p, _, _ = d.pick_period({q: (zs[q], raw[q]) for q in scores})
            sc = scores[p][0]
            d.REGION_THRESHOLD = 0.5
            span = d.periodic_region(idx, p, d.SIM, exp)
            n_rep, uni_frac = uniprot_span(t["accession"])
            gate = (
                "z_seq below Z_MIN"
                if zs[p] < d.Z_MIN
                else "score below MIN_SCORE"
                if sc < d.MIN_SCORE
                else "region test or coverage"
            )
            out.append(
                {
                    "proteome": prot,
                    "gene": t["gene"] or t["accession"],
                    "length": len(s),
                    "best_period": p,
                    "z_seq": f"{zs[p]:.2f}",
                    "score": f"{sc:.3f}",
                    "stops_at": gate,
                    "region_fraction_at_cut_0.5": f"{(span[1] - span[0]) / len(s):.2f}",
                    "uniprot_repeats": n_rep,
                    "uniprot_span_fraction": f"{uni_frac:.2f}",
                }
            )
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} missed positives -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
