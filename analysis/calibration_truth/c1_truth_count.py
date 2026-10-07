"""Count the curated truth for the repeat call (decision C1), by clusters.

Reads data/curated/adhesins/adhesins.tsv, fetches the sequences from UniProt (live), clusters
them with MMseqs2 (30% identity, 50% coverage, default coverage mode) and prints counts of rows,
sequences and independent clusters per class.

This is a count. It is not a calibration. Re-running can change the numbers, because UniProt
changes. Usage: python3.12 c1_truth_count.py --out DIR   (MMseqs2 must be on PATH)
"""

import argparse
import csv
import subprocess
import time
import urllib.parse
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def fetch(accs):
    seqs = {}
    for i in range(0, len(accs), 40):
        q = " OR ".join(f"accession:{a}" for a in accs[i : i + 40])
        url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=" + urllib.parse.quote(
            q
        )
        text = subprocess.run(
            ["curl", "-s", "-m", "60", url], capture_output=True, text=True, check=True
        ).stdout
        cur = None
        for line in text.splitlines():
            if line.startswith(">"):
                cur = line.split("|")[1]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
        time.sleep(0.3)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--table", default=str(REPO / "data/curated/adhesins/adhesins.tsv"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(args.table), delimiter="\t"))
    wanted = ("adhesin", "hard_negative", "surface_other_adhesion_phenotype")
    accs = sorted({r["accession"] for r in rows if r["accession"] and r["cls"] in wanted})
    seqs = fetch(accs)
    missing = [a for a in accs if a not in seqs]
    with open(out / "all.faa", "w") as f:
        for a, s in seqs.items():
            f.write(f">{a}\n{s}\n")
    subprocess.run(
        [
            "mmseqs",
            "easy-cluster",
            str(out / "all.faa"),
            str(out / "clu"),
            str(out / "tmp"),
            "--min-seq-id",
            "0.3",
            "-c",
            "0.5",
            "-v",
            "1",
            "--threads",
            "4",
        ],
        check=True,
        capture_output=True,
    )
    rep = {}
    for line in open(out / "clu_cluster.tsv"):
        a, b = line.rstrip("\n").split("\t")
        rep[b] = a

    def count(sel, label):
        got = {r["accession"] for r in sel if r["accession"] in rep}
        clusters = {rep[a] for a in got}
        print(f"{label}\trows={len(sel)}\tsequences={len(got)}\tclusters={len(clusters)}")
        return clusters

    adh = [r for r in rows if r["cls"] == "adhesin"]
    e1 = [r for r in adh if r["evidence_level"] == "E1"]
    e12 = [r for r in adh if r["evidence_level"] in ("E1", "E2")]
    print(f"accessions without a sequence: {len(missing)} {missing[:10]}")
    count(e1, "adhesin E1")
    pos = count(e12, "adhesin E1+E2")
    count([r for r in e12 if r["reviewed"] == "reviewed"], "adhesin E1+E2, reviewed")
    count([r for r in e1 if r["needs_review"] == "no"], "adhesin E1, needs_review=no")
    hn = [r for r in rows if r["cls"] == "hard_negative"]
    neg = count(hn, "hard negative N1+N2")
    count([r for r in hn if r["evidence_level"] == "N1"], "hard negative N1")
    count(
        [r for r in rows if r["cls"] == "surface_other_adhesion_phenotype"],
        "surface_other_adhesion_phenotype",
    )
    print(f"clusters that hold both an E1+E2 adhesin and a hard negative: {len(pos & neg)}")

    def text(r):
        return " ".join((r["family"], r["evidence_summary"], r["protein_name"])).lower()

    count(
        [r for r in e12 if "repeat" in text(r) or "tandem" in text(r)],
        "adhesin E1+E2 with 'repeat' or 'tandem' in family, summary or name",
    )
    print("rows with an empty family column among E1+E2:", sum(1 for r in e12 if not r["family"]))


if __name__ == "__main__":
    main()
