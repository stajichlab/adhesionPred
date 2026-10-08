#!/usr/bin/python3.12
"""Inputs for `cellsurface_sorting_hat_calibrate truth` per proteome (task H2c).

For each proteome with at least --min-pos mapped T1/T2 positives:
  label 1  mapped T1 or T2 truth protein
  excluded mapped T3 protein (a rule-only label is neither positive nor negative)
  label 0  every other protein of the proteome (ASSUMED negative: absent annotation, not measured absence)
All rows are clustered with MMseqs2 (30% identity, coverage 0.5), so every negative has a cluster. The
development/test split is made once, by cluster, with a fixed seed, before any measurement. The status
tables are written from the test part.

Usage: make_calibration_inputs.py --dir DIR --proteome-dir-config taxa.tsv [--min-pos 3]
Needs mmseqs on PATH (module load mmseqs2/17-b804f). Reads plain or .gz FASTA.
"""

import argparse
import csv
import gzip
import random
import subprocess
import sys
import tempfile
from pathlib import Path

FUNGI5K = "/bigdata/stajichlab/shared/projects/Fungi_5k/input"
SORTING = "_workdir/sorting_hat"
SEED = 20261008


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with _open(path) as fh:
        for line in fh:
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


def build_rows(proteome_ids, positives, excluded):
    missing = sorted(set(positives) - set(proteome_ids))
    if missing:
        raise ValueError(f"positives not in the proteome: {missing}")
    rows = []
    for p in proteome_ids:
        if p in excluded and p not in positives:
            continue
        rows.append((p, 1 if p in positives else 0))
    return rows


def attach_clusters(rows, clusters):
    miss = [p for p, _ in rows if p not in clusters]
    if miss:
        raise ValueError(f"{len(miss)} row(s) have no cluster, for example {miss[:3]}")
    return [(p, label, clusters[p]) for p, label in rows]


def split_parts(clusters, seed=SEED, dev_fraction=0.3):
    ids = sorted(set(clusters.values()))
    rng = random.Random(seed)
    rng.shuffle(ids)
    dev = set(ids[: max(1, round(len(ids) * dev_fraction))])
    return {p: ("dev" if c in dev else "test") for p, c in clusters.items()}


def cluster(seqs, tmp):
    faa = Path(tmp) / "rows.faa"
    names = list(seqs)  # MMseqs2 rewrites IDs such as sp|ACC|NAME, so cluster under neutral IDs
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
            str(Path(tmp) / "work"),
            "--min-seq-id",
            "0.3",
            "-c",
            "0.5",
            "--threads",
            "4",
        ],
        check=True,
        capture_output=True,
    )
    cl = {}
    for line in open(f"{pre}_cluster.tsv"):
        rep, member = line.split()
        cl[names[int(member[1:])]] = names[int(rep[1:])]
    return cl


def proteome_path(row):
    f = row["file"]
    if f.startswith("Fungi_5k "):
        return f"{FUNGI5K}/{f.split(' ', 1)[1]}"
    return f"{SORTING}/{f}"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--dir", required=True)
    ap.add_argument("--min-pos", type=int, default=3)
    a = ap.parse_args()
    d = Path(a.dir)
    truth = {r["accession"]: r for r in csv.DictReader(open(d / "truth_all.tsv"), delimiter="\t")}
    taxa = list(csv.DictReader(open(d / "taxa.tsv"), delimiter="\t"))
    mp = list(csv.DictReader(open(d / "proteome_map.tsv"), delimiter="\t"))
    split_path = d / "split.tsv"  # the committed copy is split.tsv.gz
    if split_path.exists():
        raise SystemExit(f"{split_path} exists: the split is made once")
    split_rows = []
    for t in taxa:
        name = t["proteome"]
        mapped = [r for r in mp if r["proteome"] == name]
        pos = {
            r["proteome_id"] for r in mapped if truth[r["truth_accession"]]["tier"] in ("T1", "T2")
        }
        t3 = {r["proteome_id"] for r in mapped if truth[r["truth_accession"]]["tier"] == "T3"}
        if len(pos) < a.min_pos or not t["taxon_id"]:
            print(f"{name}: {len(pos)} positives, skipped", flush=True)
            continue
        seqs = read_fasta(proteome_path(t))
        rows = build_rows(list(seqs), pos, t3)
        with tempfile.TemporaryDirectory() as tmp:
            cl = cluster({p: seqs[p] for p, _ in rows}, tmp)
        full = attach_clusters(rows, cl)
        part = split_parts({p: c for p, _, c in full})
        out = d / "calibration" / name
        out.mkdir(parents=True, exist_ok=True)
        for tag, keep in (("all", None), ("test", "test"), ("dev", "dev")):
            with open(out / f"truth_{tag}.tsv", "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["id", "label", "cluster"])
                for p, label, c in full:
                    if keep is None or part[p] == keep:
                        w.writerow([p, label, c])
        for p, _label, c in full:
            split_rows.append([name, p, c, part[p], SEED])
        npos = {
            k: sum(1 for p, lab, _ in full if lab == 1 and part[p] == k) for k in ("dev", "test")
        }
        print(
            f"{name}: {len(pos)} positives, {sum(1 for _, lab in rows if lab == 0)} assumed negatives, "
            f"{len({c for _, _, c in full})} clusters; positives dev {npos['dev']} test {npos['test']}",
            flush=True,
        )
    with open(split_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["proteome", "id", "cluster", "part", "seed"])
        w.writerows(split_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
