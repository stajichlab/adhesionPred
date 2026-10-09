#!/usr/bin/python3.12
"""Task L3: hard-negative sets (rules in hard_negative_queries.md, written before the fetch).

Needs mmseqs and hmmsearch on PATH (module load mmseqs2/17-b804f hmmer/3.4).
Usage: make_hard_negatives.py --dir analysis/hydrophobin_truth
"""

import argparse
import csv
import gzip
import json
import random
import re
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

SEED = 20261008
HYD_NAME = re.compile(r"hydrophobin|rodlet", re.I)
HYD_MODELS = {
    "Hydrophobin",
    "Hydrophobin_2",
    "Eas",
    "DewD",
    "Hyd1F",
    "Hydrophobin_D",
    "Hydrophobin_like",
}
PFAM_GROUPS = {
    "CFEM": "xref:pfam-PF05730",
    "cerato_platanin": "xref:pfam-PF07249",
    "HsbA": "xref:pfam-PF12296",
}
CELLWALL_Q = "(xref:pfam-PF00399 OR gene:CCW12 OR gene:CCW14)"
SECRETED_Q = "(keyword:KW-0964 OR keyword:KW-0800) AND length:[40 TO 250]"
MODELS_HMM = "_workdir/sorting_hat/calibration/hydrophobin_discovery/models.hmm"


def small_secreted_ok(r):
    s = r["sequence"].upper()
    return 40 <= len(s) <= 250 and s.count("C") >= 6 and not HYD_NAME.search(r["name"] or "")


def remove_hydrophobin_like(groups, hits, names):
    drop = set(hits) | {a for a, n in names.items() if HYD_NAME.search(n or "")}
    kept = {g: [r for r in rows if r["accession"] not in drop] for g, rows in groups.items()}
    removed = sorted(
        a for rows in groups.values() for r in rows for a in [r["accession"]] if a in drop
    )
    return kept, sorted(set(removed))


def sample_unreviewed(pool, n, seed):
    rng = random.Random(seed)
    pool = sorted(pool, key=lambda r: r["accession"])
    return pool if len(pool) <= n else rng.sample(pool, n)


def split_parts(clusters, seed=SEED):
    ids = sorted(set(clusters.values()))
    rng = random.Random(seed)
    rng.shuffle(ids)
    test = set(ids[: len(ids) // 2])
    return {p: ("test" if c in test else "tuning") for p, c in clusters.items()}


def fetch(query, reviewed):
    q = f"({query}) AND taxonomy_id:4751 AND reviewed:{'true' if reviewed else 'false'}"
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        {
            "query": q,
            "format": "tsv",
            "fields": "accession,id,protein_name,organism_name,length,sequence",
        }
    )
    with urllib.request.urlopen(url, timeout=300) as r:
        text = r.read().decode()
    rows = []
    for x in csv.DictReader(text.splitlines(), delimiter="\t"):
        rows.append(
            {
                "accession": x["Entry"],
                "name": x["Protein names"],
                "species": x["Organism"],
                "sequence": x["Sequence"],
                "length": int(x["Length"]),
                "source": "reviewed" if reviewed else "unreviewed_sample",
            }
        )
    return rows


def hmm_hits(rows, tmp):
    faa = Path(tmp) / "all.faa"
    with open(faa, "w") as f:
        for r in rows:
            f.write(f">{r['accession']}\n{r['sequence']}\n")
    dom = Path(tmp) / "dom.tbl"
    subprocess.run(
        [
            "hmmsearch",
            "--cut_ga",
            "--noali",
            "-o",
            "/dev/null",
            "--domtblout",
            str(dom),
            MODELS_HMM,
            str(faa),
        ],
        check=True,
        capture_output=True,
    )
    out = set()
    for line in open(dom):
        if line[0] != "#":
            c = line.split(None, 22)
            if c[3] in HYD_MODELS:
                out.add(c[0])
    return out


def cluster(seqs, tmp):
    names = list(seqs)
    faa = Path(tmp) / "c.faa"
    with open(faa, "w") as f:
        for i, k in enumerate(names):
            f.write(f">s{i}\n{seqs[k]}\n")
    pre = Path(tmp) / "clu"
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
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    groups = {}
    for g, q in PFAM_GROUPS.items():
        rev = fetch(q, True)
        unrev = sample_unreviewed(fetch(q, False), 300, SEED)
        groups[g] = rev + unrev
    cw_rev = fetch(CELLWALL_Q, True)
    cw_unrev = sample_unreviewed(fetch("xref:pfam-PF00399", False), 300, SEED)
    groups["cell_wall_cys"] = cw_rev + cw_unrev
    groups["small_secreted_cys"] = [r for r in fetch(SECRETED_Q, True) if small_secreted_ok(r)]
    allrows = {r["accession"]: r for rows in groups.values() for r in rows}
    with tempfile.TemporaryDirectory() as tmp:
        hits = hmm_hits(list(allrows.values()), tmp)
    names = {k: r["name"] for k, r in allrows.items()}
    kept, removed = remove_hydrophobin_like(groups, hits, names)
    with open(d / "removed_hydrophobin_like.tsv", "w") as fh:
        fh.write("accession\treason\n")
        for acc in removed:
            fh.write(f"{acc}\t{'pfam_hydrophobin_hit' if acc in hits else 'name'}\n")
    out_rows = []
    counts = {}
    for g, rows in kept.items():
        with tempfile.TemporaryDirectory() as tmp:
            cl = cluster({r["accession"]: r["sequence"] for r in rows}, tmp)
        part = split_parts(cl)
        for r in rows:
            out_rows.append(
                (
                    g,
                    r["accession"],
                    r["source"],
                    r["species"],
                    cl[r["accession"]],
                    part[r["accession"]],
                    r["sequence"],
                )
            )
        counts[g] = {
            "reviewed": sum(1 for r in rows if r["source"] == "reviewed"),
            "unreviewed_sample": sum(1 for r in rows if r["source"] != "reviewed"),
            "clusters": len(set(cl.values())),
            "tuning": sum(1 for p in part.values() if p == "tuning"),
            "test": sum(1 for p in part.values() if p == "test"),
        }
    with gzip.open(d / "hard_negatives.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["group", "accession", "source", "species", "cluster", "part", "sequence"])
        w.writerows(out_rows)
    json.dump(counts, open(d / "hard_negative_counts.json", "w"), indent=1)
    print(json.dumps(counts, indent=1))
    print("removed as hydrophobin-like:", len(removed))
    return 0


if __name__ == "__main__":
    sys.exit(main())
