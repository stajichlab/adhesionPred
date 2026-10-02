#!/usr/bin/env python3
"""Phase C step 2: MMseqs2 clusters, split tables, T-c removals, maximum identity.

Reads $STEP1_WORKDIR/phasec/eval_table.tsv.gz, eval_literature.tsv and eval_sequences.fasta.gz
(08; their SHA-256 must equal build_run.json) and species.tsv (its SHA-256 must equal the one
08 used). Runs MMseqs2 (spec 3.4 items 3 and 6):

  easy-cluster  --min-seq-id 0.3 -c 0.5 --cov-mode 0           over all table and literature
                                                               sequences
  easy-search   -s 7.5 -c 0.5 --cov-mode 0                     per S2 and S3 split: query =
                --format-output query,target,fident            test proteins, target = GO
                                                               training proteins (same for
                                                               V-go and V-kw)

Writes to $STEP1_WORKDIR/phasec/ (mmseqs.log goes to logs/mmseqs.log there, also after a STOP):

  clusters.tsv.gz       seq_sha256 -> cluster_id (the representative's hash)
  split_members.tsv.gz  one row per split, fold and sequence: part, origin, class, cluster
  tc_removed.tsv        one row per T-c row removed from a split (rule a, b or c)
  max_identity.tsv.gz   one row per S2/S3 split and test sequence; no hit = below 0.3
  splits_run.json       MMseqs2 version and commands, seed, counts, input hashes

STOP (exit 2, no output): stale 08 outputs; another species.tsv than 08 used; `mmseqs version`
fails (exit 132 means an AVX2 build on a CPU without AVX2: run on partition epyc or set
STEP1_MMSEQS to a non-AVX2 binary); an MMseqs2 command fails; a sequence without a cluster; a
split that breaks a leakage rule (splits.check_*).
"""

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import evalio
import manifest
import paths
import runinfo
import splits
import truth_table

OUTPUT_NAMES = (
    "clusters.tsv.gz",
    "split_members.tsv.gz",
    "tc_removed.tsv",
    "max_identity.tsv.gz",
    "splits_run.json",
)
NON_AVX2_MMSEQS = "/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs"
CLUSTER_ARGS = ("--min-seq-id", "0.3", "-c", "0.5", "--cov-mode", "0")
SEARCH_ARGS = (
    "-s",
    "7.5",
    "-c",
    "0.5",
    "--cov-mode",
    "0",
    "--format-output",
    "query,target,fident",
)


def mmseqs_version(mmseqs: str) -> str:
    try:
        done = subprocess.run([mmseqs, "version"], capture_output=True, text=True)
    except OSError as exc:
        raise evalio.StopError(
            f"cannot run {mmseqs} ({exc}); module load MMseqs2/17-b804f"
        ) from exc
    if done.returncode != 0:
        raise evalio.StopError(
            f"`{mmseqs} version` exited with {done.returncode} (132 or -4 = illegal instruction: an "
            "AVX2 build on a CPU without AVX2; run on partition epyc or set STEP1_MMSEQS to "
            f"the non-AVX2 binary {NON_AVX2_MMSEQS})"
        )
    return done.stdout.strip()


def run_mmseqs(argv: list[str], log_path: Path) -> None:
    with open(log_path, "a") as log:
        done = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)
    if done.returncode != 0:
        raise evalio.StopError(f"mmseqs {argv[1]} exited with {done.returncode}; see {log_path}")


def write_fasta(path: Path, hashes, seqs: dict) -> None:
    with open(path, "w") as handle:
        for h in sorted(set(hashes)):
            handle.write(f">{h}\n{seqs[h]}\n")


def read_fasta_gz(path: Path) -> dict[str, str]:
    seqs, head = {}, None
    for line in gzip.decompress(Path(path).read_bytes()).decode().splitlines():
        if line.startswith(">"):
            head = line[1:].strip()
            seqs[head] = ""
        elif head is not None:
            seqs[head] += line.strip()
    return seqs


def run(work: Path, species_path: Path, mmseqs: str, tmp: Path, threads: int, arguments=()):
    """Run in a fresh directory made inside `tmp`; only that directory is removed (never `tmp`)."""
    tmp = Path(tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    fresh = Path(tempfile.mkdtemp(prefix="splits_", dir=tmp))
    try:
        return _run(work, species_path, mmseqs, fresh, threads, arguments)
    finally:
        shutil.rmtree(fresh, ignore_errors=True)


def _run(work: Path, species_path: Path, mmseqs: str, tmp: Path, threads: int, arguments=()):
    out = evalio.out_dir(work)
    build = evalio.require_current(
        out,
        "build_run.json",
        ("eval_table.tsv.gz", "eval_literature.tsv", "eval_sequences.fasta.gz"),
        "08_build_eval_tables.py",
    )
    if manifest.sha256_file(species_path) != build["input_sha256"]["species.tsv"]:
        raise evalio.StopError(f"{species_path} differs from the species.tsv that 08 used")
    version = mmseqs_version(mmseqs)
    # the log lives on /bigdata, not in the temp dir that run() removes, so a STOP can name it
    log_path = out / "logs" / "mmseqs.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("")
    seqs = read_fasta_gz(out / "eval_sequences.fasta.gz")
    fasta = tmp / "all.fasta"
    write_fasta(fasta, seqs, seqs)
    commands = []
    argv = [mmseqs, "easy-cluster", str(fasta), str(tmp / "clu"), str(tmp / "w_clu"), *CLUSTER_ARGS]
    argv += ["--threads", str(threads)]
    commands.append(["easy-cluster", *CLUSTER_ARGS])
    run_mmseqs(argv, log_path)
    raw = tmp / "clu_cluster.tsv"
    cluster_of = splits.read_cluster_tsv(raw.read_text().splitlines())
    splits.check_clusters(cluster_of, seqs)
    table = truth_table.read_tsv(out / "eval_table.tsv.gz")
    lit = truth_table.read_tsv(out / "eval_literature.tsv")
    defs = splits.split_definitions(truth_table.read_tsv(species_path))
    members, removed = splits.build(table, lit, cluster_of, defs, evalio.SEED)
    splits.check_test_truth(members)
    splits.check_no_shared_hash(members)
    splits.check_no_cluster_spans(members)
    identity = []
    for d in defs:
        if not d.test_sources:
            continue
        rows = [m for m in members if m["split_id"] == d.split_id]
        query = [m["seq_sha256"] for m in rows if m["part"] in ("test", "test_lit")]
        target = [m["seq_sha256"] for m in rows if m["part"] == "train"]
        q, t, hits = (
            tmp / f"{d.split_id}.q.fasta",
            tmp / f"{d.split_id}.t.fasta",
            tmp / f"{d.split_id}.m8",
        )
        write_fasta(q, query, seqs)
        write_fasta(t, target, seqs)
        argv = [mmseqs, "easy-search", str(q), str(t), str(hits), str(tmp / f"w_{d.split_id}")]
        argv += [*SEARCH_ARGS, "--threads", str(threads)]
        run_mmseqs(argv, log_path)
        best = splits.parse_hits(hits.read_text().splitlines())
        identity += splits.identity_rows(d.split_id, query, best)
    commands.append(["easy-search", *SEARCH_ARGS])
    clusters = [{"seq_sha256": h, "cluster_id": cluster_of[h]} for h in sorted(cluster_of)]
    log = {
        "all_sources": build["all_sources"],
        "truth_set_sha256": build["truth_set_sha256"],
        "input_sha256": {
            n: build["outputs_sha256"][n]
            for n in ("eval_table.tsv.gz", "eval_literature.tsv", "eval_sequences.fasta.gz")
        }
        | {"species.tsv": manifest.sha256_file(species_path)},
        "mmseqs_version": version,
        "mmseqs_commands": commands,
        "cluster_tsv_sha256": manifest.sha256_file(raw),
        "seed": evalio.SEED,
        "sequences": len(seqs),
        "clusters": len(set(cluster_of.values())),
        "splits": [d.split_id for d in defs],
        "members_by_split_fold_part": {
            f"{s}|{f}|{p}": n
            for (s, f, p), n in sorted(
                Counter((m["split_id"], m["fold"], m["part"]) for m in members).items()
            )
        },
        "tc_removed_by_split_fold_rule": {
            f"{s}|{f}|{r}": n
            for (s, f, r), n in sorted(
                Counter((x["split_id"], x["fold"], x["rule"]) for x in removed).items()
            )
        },
        "identity_below_0.3_by_split": dict(
            sorted(Counter(r["split_id"] for r in identity if r["below_0.3"] == "yes").items())
        ),
        "git_commit": runinfo.git_commit(),
        "library_versions": evalio.library_versions(),
        "arguments": list(arguments),
    }
    evalio.write_outputs(
        out,
        {
            "clusters.tsv.gz": lambda p: truth_table.write_tsv(p, splits.CLUSTER_COLUMNS, clusters),
            "split_members.tsv.gz": lambda p: truth_table.write_tsv(
                p, splits.MEMBER_COLUMNS, members
            ),
            "tc_removed.tsv": lambda p: truth_table.write_tsv(p, splits.REMOVED_COLUMNS, removed),
            "max_identity.tsv.gz": lambda p: truth_table.write_tsv(
                p, splits.IDENTITY_COLUMNS, identity
            ),
        },
        "splits_run.json",
        log,
    )
    return log


def main(argv=None) -> int:
    import os

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--mmseqs", default=os.environ.get("STEP1_MMSEQS", "mmseqs"))
    parser.add_argument("--tmp-dir", required=True, help="node-local work dir, e.g. $SCRATCH/x")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        log = run(
            work,
            Path(args.species),
            args.mmseqs,
            Path(args.tmp_dir),
            args.threads,
            list(argv) if argv is not None else sys.argv[1:],
        )
    except (evalio.StopError, splits.SplitError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    print(f"sequences={log['sequences']} clusters={log['clusters']} mmseqs={log['mmseqs_version']}")
    for key, n in log["tc_removed_by_split_fold_rule"].items():
        print(f"  tc_removed {key}={n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
