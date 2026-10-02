#!/usr/bin/env python3
"""Choose the proteins that go into the PF28404 family tree.

Tips: every Onygenales hit of length <= --max-len, plus up to --per-order proteomes (fixed seed)
from every other order that has hits, taking each proteome's best-scoring hit of length <=
--max-len. Hits longer than --max-len are multi-domain proteins and are left out of the tree.
Writes tips.tsv and tips.faa (header `<label>|<protein_id>`). Reads hits.tsv and hits.faa of the
work directory.
"""

import argparse
import csv
import random
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", required=True)
    ap.add_argument("--max-len", type=int, default=400)
    ap.add_argument("--per-order", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20261002)
    a = ap.parse_args()
    work = Path(a.work)
    with open(work / "hits.tsv") as fh:
        hits = [h for h in csv.DictReader(fh, delimiter="\t") if int(h["length"]) <= a.max_len]
    seqs, name = {}, None
    for ln in open(work / "hits.faa"):
        if ln.startswith(">"):
            name = ln[1:].strip()
            seqs[name] = ""
        else:
            seqs[name] += ln.strip()
    keep = [h for h in hits if h["order"] == "Onygenales"]
    rng = random.Random(a.seed)
    by_order = {}
    for h in hits:
        if h["order"] != "Onygenales":
            by_order.setdefault((h["class"], h["order"]), {}).setdefault(h["label"], []).append(h)
    for key in sorted(by_order):
        labels = sorted(by_order[key])
        for lab in rng.sample(labels, min(a.per_order, len(labels))):
            keep.append(max(by_order[key][lab], key=lambda h: float(h["score"])))
    with open(work / "tips.tsv", "w") as fh:
        fh.write("tip\tlabel\tprotein_id\tclass\torder\tlength\tscore\n")
        for h in keep:
            fh.write(
                f"{h['label']}|{h['protein_id']}\t{h['label']}\t{h['protein_id']}\t{h['class']}\t{h['order']}\t{h['length']}\t{h['score']}\n"
            )
    with open(work / "tips.faa", "w") as fh:
        for h in keep:
            n = f"{h['label']}|{h['protein_id']}"
            fh.write(f">{n}\n{seqs[n]}\n")
    print(f"{len(keep)} tips ({sum(1 for h in keep if h['order'] == 'Onygenales')} Onygenales)")


if __name__ == "__main__":
    main()
