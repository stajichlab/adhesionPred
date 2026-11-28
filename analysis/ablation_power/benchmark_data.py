"""Load the frozen adhesion benchmark and build leave-one-clade-out splits for the power study.

This module decides nothing about components. It exists so that the power calibration
(`run_power.py`) and any later ablation harness read the benchmark the same way.

Two rules it enforces, both from the 2026-11-27 pre-flight checks:

  STRICT CLUSTER RULE (decision D-A). When fold F is held out, every protein whose MMseqs cluster
  contains any member of F is dropped from TRAINING. Five clusters span folds in adhesion_v1, so
  without this a held-out clade is still represented in training by homologues at >= 30% identity --
  exactly the leakage leave-clade-out exists to prevent.

  TRAINING POOL (decision D-D) must be stated, never implied. "all" trains on every protein outside
  the held-out fold, including the six folds marked scored=no; "scored" trains only on the three
  scorable folds. The choice changes the leakage drop counts (Pezizomycotina 5 vs 2) and brings in
  ~26 positives from clades that are never tested, where positivity is confounded with clade. The
  caller passes it explicitly and it is recorded in every output.

Scoring is at the PROTEIN level and resampling is by CLUSTER (decision D-B): the FLO/ALS/EPA
paralogue families would otherwise count as independent evidence in a confidence interval.
"""

import csv
import hashlib
from pathlib import Path

import numpy as np

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
SCORED_TIERS = ("positive", "hard_negative")
TRAIN_POOLS = ("all", "scored")

# Feature-set definitions. The baseline "bar" is the best of this declared family, not
# composition+length by default: measured on adhesion_v1, log length ALONE scores a pooled
# AP of 0.690 against 0.551 for the 21-feature model, so the 21-feature model is not a bar.
FEATURE_SETS = {
    "length": (20,),
    "composition": tuple(range(20)),
    "composition_length": tuple(range(21)),
}


class BenchmarkError(Exception):
    pass


class Benchmark:
    """The scored part of a frozen benchmark, with folds, clusters and baseline features."""

    def __init__(self, accessions, y, fold, cluster, features, version, digest, train_pool):
        self.accessions = accessions
        self.y = y
        self.fold = fold
        self.cluster = cluster
        self.features = features
        self.version = version
        self.digest = digest
        self.train_pool = train_pool

    def __len__(self):
        return len(self.accessions)

    @property
    def scored_folds(self):
        """Folds that may be held out, in a stable order."""
        return self._scored_folds

    def subset(self, idx):
        idx = np.asarray(idx, dtype=int)
        return Benchmark(
            [self.accessions[i] for i in idx],
            self.y[idx],
            [self.fold[i] for i in idx],
            [self.cluster[i] for i in idx],
            self.features[idx],
            self.version,
            self.digest,
            self.train_pool,
        )


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def verify_checksums(directory):
    """Refuse to run on a benchmark that does not match its own SHA256SUMS.

    A result that cannot name the exact benchmark it came from is not a result.
    """
    directory = Path(directory)
    sums = directory / "SHA256SUMS"
    if not sums.exists():
        raise BenchmarkError(f"{sums} is missing; refusing to run on an unverified benchmark")
    bad = []
    for line in sums.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        want, _, name = line.partition("  ")
        name = name.strip()
        target = directory / name
        if not target.exists():
            bad.append(f"{name}: missing")
        elif sha256_file(target) != want:
            bad.append(f"{name}: digest does not match")
    if bad:
        raise BenchmarkError("benchmark does not match SHA256SUMS: " + "; ".join(bad))
    return sha256_file(sums)


def parse_fasta(path):
    seqs, acc, parts = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if acc is not None:
                seqs[acc] = "".join(parts)
            acc, parts = line[1:].split()[0], []
        elif acc is not None:
            parts.append(line.strip())
    if acc is not None:
        seqs[acc] = "".join(parts)
    return seqs


def composition_features(sequence):
    """20 amino-acid frequencies followed by log length. Unknown residues are ignored."""
    n = len(sequence)
    if n == 0:
        raise BenchmarkError("empty sequence")
    counts = [sequence.count(a) / n for a in AA_ORDER]
    return np.array(counts + [np.log(n)], dtype=float)


def clade_to_fold(folds_rows):
    """Map each clade_class to its fold, and record which folds are scorable."""
    mapping, scored = {}, []
    for row in folds_rows:
        fold = row["fold"]
        if row["scored"] == "yes":
            scored.append(fold)
        for member in row["members"].split(","):
            member = member.strip()
            if member:
                mapping[member] = fold
    return mapping, scored


def load(directory, train_pool="all", verify=True):
    """Read a frozen benchmark directory into a Benchmark."""
    if train_pool not in TRAIN_POOLS:
        raise BenchmarkError(f"train_pool must be one of {TRAIN_POOLS}, got {train_pool!r}")
    directory = Path(directory)
    digest = verify_checksums(directory) if verify else "unverified"

    proteins = {
        r["accession"]: r for r in read_tsv(directory / "proteins.tsv") if r["tier"] in SCORED_TIERS
    }
    clusters = {
        r["accession"]: r["cluster_representative"] for r in read_tsv(directory / "clusters.tsv")
    }
    mapping, scored = clade_to_fold(read_tsv(directory / "folds.tsv"))
    sequences = parse_fasta(directory / "sequences.faa")
    version = (directory / "VERSION").read_text().splitlines()[0].strip()

    accessions, y, fold, cluster, feats = [], [], [], [], []
    for acc in sorted(sequences):
        row = proteins.get(acc)
        if row is None:
            continue
        f = mapping.get(row["clade_class"])
        if f is None:
            # An unmapped clade excludes the row rather than guessing its fold.
            continue
        accessions.append(acc)
        y.append(1 if row["tier"] == "positive" else 0)
        fold.append(f)
        cluster.append(clusters.get(acc, acc))
        feats.append(composition_features(sequences[acc]))

    if not accessions:
        raise BenchmarkError(f"no scored protein with a sequence and a fold in {directory}")

    bench = Benchmark(
        accessions,
        np.array(y, dtype=int),
        fold,
        cluster,
        np.vstack(feats),
        version,
        digest,
        train_pool,
    )
    bench._scored_folds = [f for f in scored if f in set(fold)]
    return bench


def splits(bench):
    """Leave-one-clade-out splits honouring the strict cluster rule and the declared training pool.

    Yields (fold_name, train_index, test_index, n_dropped_for_leakage).
    """
    fold = np.array(bench.fold)
    cluster = np.array(bench.cluster)
    scored = set(bench.scored_folds)

    cluster_folds = {}
    for c, f in zip(cluster, fold, strict=False):
        cluster_folds.setdefault(c, set()).add(f)

    for held in bench.scored_folds:
        test = np.flatnonzero(fold == held)
        eligible = fold != held
        if bench.train_pool == "scored":
            eligible &= np.isin(fold, list(scored))
        leaky = np.array([held in cluster_folds[c] for c in cluster])
        train = np.flatnonzero(eligible & ~leaky)
        dropped = int(np.sum(eligible & leaky))
        yield held, train, test, dropped


def assert_no_leakage(bench):
    """Every split must share no cluster between train and test. Raises if one does."""
    cluster = np.array(bench.cluster)
    for held, train, test, _ in splits(bench):
        shared = set(cluster[train]) & set(cluster[test])
        if shared:
            raise BenchmarkError(
                f"fold {held}: {len(shared)} cluster(s) in both train and test, "
                f"first: {sorted(shared)[:3]}"
            )
    return True
