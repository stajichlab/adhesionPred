"""Cluster the adhesion benchmark's proteins and recompute the fold list on clusters, not proteins.

Why this exists: `build_benchmark.py` counts proteins, and protein counts flatter the benchmark.
The 32 Saccharomycetes positives are largely FLO/ALS paralogues, which collapse to a handful of
independent observations once sequences are clustered. A fold list built on protein counts can
therefore claim a clade is scorable when it rests on one or two families. Clusters are the unit the
splits must use, so a cluster never straddles a train/test boundary.

Two steps, because only the first needs the network:

1. ``fetch`` - download one sequence per accession in ``proteins.tsv`` (UniProt REST, or NCBI eutils
   for a RefSeq/GenBank accession) and write ``sequences.faa``. Resumable: an accession already in
   the output file is not fetched again. Login node; a few hundred small requests.
2. ``cluster`` - run MMseqs2 ``easy-cluster`` at the given identity and coverage, then rewrite
   ``folds.tsv`` with cluster counts and add ``clusters.tsv`` (accession -> cluster representative).
   Seconds for a few hundred sequences.

Usage::

    python3 -I analysis/curation/cluster_benchmark.py fetch   --benchmark data/benchmark/adhesion_v1
    python3 -I analysis/curation/cluster_benchmark.py cluster --benchmark data/benchmark/adhesion_v1 \\
        --mmseqs mmseqs --tmpdir "$SCRATCH/mmseqs_benchmark"

``cluster`` writes ``folds.tsv``, ``clusters.tsv`` and refreshes ``SHA256SUMS``; it leaves
``proteins.tsv`` untouched, so the membership of the set never changes here.
"""

import argparse
import csv
import importlib.util
import json
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

_spec = importlib.util.spec_from_file_location(
    "build_benchmark", Path(__file__).parent / "build_benchmark.py"
)
_bb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_bb)

BuildError = _bb.BuildError
UNIPROT = "https://rest.uniprot.org/uniprotkb"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UA = {"User-Agent": "adhesionpred-benchmark/1 (python urllib)"}
SCORED_TIERS = ("positive", "hard_negative")


class FetchError(Exception):
    pass


def http_get(url, retries=4, backoff=2.0):
    """GET with retries on timeouts, 429 and 5xx. A 404 raises at once: it will not get better."""
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode()
        except urllib.error.HTTPError as e:
            if e.code in (404, 400):
                raise FetchError(f"{url}: HTTP {e.code}") from e
            last = e
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last = e
        time.sleep(backoff * (attempt + 1))
    raise FetchError(f"{url}: {last}")


def parse_fasta(text):
    """{id: sequence} from FASTA text; the id is the first whitespace-delimited token."""
    out, name, chunks = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        elif name:
            chunks.append(line.strip())
    if name:
        out[name] = "".join(chunks)
    return out


def fetch_sequence(accession, db, fetch=http_get):
    """One sequence, from UniProt or NCBI. Returns '' when the accession has no record."""
    acc = accession.strip()
    try:
        if db.lower().startswith("uniprot") or "_" not in acc.split(".")[0]:
            text = fetch(f"{UNIPROT}/{acc}.fasta")
            if text.strip():
                seqs = parse_fasta(text)
                if seqs:
                    return next(iter(seqs.values()))
    except FetchError:
        pass
    try:
        text = fetch(f"{EUTILS}?db=protein&id={acc}&rettype=fasta&retmode=text")
        seqs = parse_fasta(text)
        return next(iter(seqs.values()), "")
    except FetchError:
        return ""


def read_existing(path):
    if not path.exists():
        return {}
    return parse_fasta(path.read_text())


def cmd_fetch(a):
    bench = Path(a.benchmark)
    rows = _bb.read_tsv(bench / "proteins.tsv")
    wanted = [r for r in rows if r["tier"] in SCORED_TIERS]
    out = bench / "sequences.faa"
    have = read_existing(out)

    missing, failed = [], []
    for r in wanted:
        if r["accession"] not in have:
            missing.append(r)
    print(f"{len(wanted)} scored proteins; {len(have)} already present; fetching {len(missing)}")

    # Append as we go, so an interrupted run loses nothing.
    with open(out, "a") as fh:
        for i, r in enumerate(missing, 1):
            seq = fetch_sequence(r["accession"], r["accession_db"] or "UniProtKB")
            if not seq:
                failed.append(r["accession"])
                print(f"  [{i}/{len(missing)}] {r['accession']}: no sequence", file=sys.stderr)
                continue
            fh.write(f">{r['accession']}\n")
            for j in range(0, len(seq), 60):
                fh.write(seq[j : j + 60] + "\n")
            fh.flush()
            if i % 25 == 0:
                print(f"  [{i}/{len(missing)}]")
            time.sleep(a.pause)

    have = read_existing(out)
    print(f"\nwrote {out}: {len(have)} sequences")
    if failed:
        (bench / "sequences_missing.txt").write_text("\n".join(failed) + "\n")
        print(f"{len(failed)} without a sequence; listed in sequences_missing.txt", file=sys.stderr)
    return 0


def run_mmseqs(mmseqs, faa, tmpdir, min_seq_id, coverage):
    """{member: representative} from `mmseqs easy-cluster`."""
    tmpdir = Path(tmpdir)
    tmpdir.mkdir(parents=True, exist_ok=True)
    prefix = tmpdir / "clu"
    cmd = [
        mmseqs,
        "easy-cluster",
        str(faa),
        str(prefix),
        str(tmpdir / "tmp"),
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(coverage),
        "--cov-mode",
        "0",
        "-v",
        "1",
    ]
    print(" ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise BuildError(f"mmseqs failed ({proc.returncode}): {proc.stderr[-2000:]}")
    tsv = Path(f"{prefix}_cluster.tsv")
    if not tsv.exists():
        raise BuildError(f"{tsv} was not written")
    members = {}
    for line in tsv.read_text().splitlines():
        rep, member = line.split("\t")[:2]
        members[member] = rep
    return members


def fold_rows_from_clusters(rows, clusters, groups, min_pos, min_neg):
    """Like build_benchmark.plan_folds, but a clade's count is its number of distinct clusters."""
    by_class = {}
    unclustered = []
    for r in rows:
        if r["tier"] not in SCORED_TIERS:
            continue
        rep = clusters.get(r["accession"])
        if rep is None:
            unclustered.append(r["accession"])
            continue
        by_class.setdefault(r["clade_class"], {"positive": set(), "hard_negative": set()})
        by_class[r["clade_class"]][r["tier"]].add(rep)

    counts = {c: {t: len(s) for t, s in tiers.items()} for c, tiers in by_class.items()}
    missing_group = sorted(c for c in counts if c not in groups)
    if missing_group:
        raise BuildError(f"classes with no group: {missing_group}")

    folds, pooled = [], {}
    for clade, n in sorted(counts.items()):
        if n["positive"] >= min_pos and n["hard_negative"] >= min_neg:
            folds.append(
                {
                    "fold": clade,
                    "level": "class",
                    "members": clade,
                    "n_positive": n["positive"],
                    "n_hard_negative": n["hard_negative"],
                    "scored": "yes",
                }
            )
        else:
            g = pooled.setdefault(groups[clade], {"positive": 0, "hard_negative": 0, "members": []})
            g["positive"] += n["positive"]
            g["hard_negative"] += n["hard_negative"]
            g["members"].append(clade)

    for group, n in sorted(pooled.items()):
        ok = n["positive"] >= min_pos and n["hard_negative"] >= min_neg
        folds.append(
            {
                "fold": group,
                "level": "group",
                "members": ",".join(sorted(n["members"])),
                "n_positive": n["positive"],
                "n_hard_negative": n["hard_negative"],
                "scored": "yes" if ok else "no",
            }
        )
    folds.sort(key=lambda f: (-f["n_positive"], f["fold"]))
    return folds, unclustered


def cmd_cluster(a):
    bench = Path(a.benchmark)
    faa = bench / "sequences.faa"
    if not faa.exists():
        raise BuildError(f"{faa} is missing; run the fetch step first")
    rows = _bb.read_tsv(bench / "proteins.tsv")
    groups = _bb.read_clade_groups(a.clade_groups)

    clusters = run_mmseqs(a.mmseqs, faa, a.tmpdir, a.min_seq_id, a.coverage)
    folds, unclustered = fold_rows_from_clusters(
        rows, clusters, groups, a.min_fold_positives, a.min_fold_negatives
    )

    with open(bench / "clusters.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(("accession", "cluster_representative"))
        for acc in sorted(clusters):
            w.writerow((acc, clusters[acc]))

    cols = ("fold", "level", "members", "n_positive", "n_hard_negative", "scored")
    with open(bench / "folds.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cols), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(folds)

    (bench / "clustering.json").write_text(
        json.dumps(
            {
                "min_seq_id": a.min_seq_id,
                "coverage": a.coverage,
                "cov_mode": 0,
                "n_sequences": len(parse_fasta(faa.read_text())),
                "n_clusters": len(set(clusters.values())),
                "unclustered": unclustered,
                "mmseqs_version": subprocess.run(
                    [a.mmseqs, "version"], capture_output=True, text=True
                ).stdout.strip(),
            },
            indent=2,
        )
        + "\n"
    )

    sums = []
    for path in sorted(bench.iterdir()):
        if path.is_file() and path.name != "SHA256SUMS":
            sums.append(f"{_bb.sha256_file(path)}  {path.name}")
    (bench / "SHA256SUMS").write_text("\n".join(sums) + "\n")

    print("fold\tlevel\tpositive\thard_negative\tscored")
    for f in folds:
        print(
            f"{f['fold']}\t{f['level']}\t{f['n_positive']}\t{f['n_hard_negative']}\t{f['scored']}"
        )
    print(f"\n{len(set(clusters.values()))} clusters from {len(clusters)} sequences")
    if unclustered:
        print(f"WARNING: {len(unclustered)} scored proteins have no sequence", file=sys.stderr)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="download one sequence per scored accession")
    f.add_argument("--benchmark", required=True)
    f.add_argument("--pause", type=float, default=0.15, help="seconds between requests")

    c = sub.add_parser("cluster", help="cluster the sequences and rewrite folds.tsv")
    c.add_argument("--benchmark", required=True)
    c.add_argument("--mmseqs", default="mmseqs")
    c.add_argument("--tmpdir", required=True, help="scratch directory for MMseqs2")
    c.add_argument("--min-seq-id", type=float, default=0.3)
    c.add_argument("--coverage", type=float, default=0.5)
    c.add_argument("--clade-groups", default=str(ROOT / "data/benchmark/clade_groups.tsv"))
    c.add_argument("--min-fold-positives", type=int, default=_bb.MIN_FOLD_POSITIVES)
    c.add_argument("--min-fold-negatives", type=int, default=_bb.MIN_FOLD_NEGATIVES)

    a = ap.parse_args(argv)
    try:
        return cmd_fetch(a) if a.cmd == "fetch" else cmd_cluster(a)
    except (BuildError, FetchError, OSError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
