"""Tests for analysis/ablation_power/benchmark_data.py.

Fixtures are synthetic, so a change to the real curated tables cannot make these pass or fail for
the wrong reason. The exception is the last group, which asserts facts about the REAL frozen
adhesion_v1 that the plan relies on -- those are meant to break if the benchmark is ever edited.
"""

import hashlib
from pathlib import Path

import benchmark_data as bd
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
REAL = ROOT / "data/benchmark/adhesion_v1"


def _write(path, rows, columns):
    path.write_text(
        "\t".join(columns)
        + "\n"
        + "".join("\t".join(str(r[c]) for c in columns) + "\n" for r in rows)
    )


def make_benchmark(tmp_path, proteins, clusters, folds, sequences, with_sums=True):
    """A minimal frozen-benchmark directory."""
    d = tmp_path / "bench"
    d.mkdir()
    _write(
        d / "proteins.tsv",
        proteins,
        ["accession", "tier", "clade_class", "length"],
    )
    _write(d / "clusters.tsv", clusters, ["accession", "cluster_representative"])
    _write(
        d / "folds.tsv",
        folds,
        ["fold", "level", "members", "n_positive", "n_hard_negative", "scored"],
    )
    (d / "sequences.faa").write_text("".join(f">{a}\n{s}\n" for a, s in sequences.items()))
    (d / "VERSION").write_text("test_v1\n")
    if with_sums:
        lines = []
        for name in ("proteins.tsv", "clusters.tsv", "folds.tsv", "sequences.faa", "VERSION"):
            digest = hashlib.sha256((d / name).read_bytes()).hexdigest()
            lines.append(f"{digest}  {name}")
        (d / "SHA256SUMS").write_text("\n".join(lines) + "\n")
    return d


@pytest.fixture
def simple(tmp_path):
    """Two scorable folds; cluster 'shared' deliberately spans both."""
    proteins, clusters, sequences = [], [], {}
    for i in range(4):
        for clade, tier in (
            ("Alpha", "positive"),
            ("Alpha", "hard_negative"),
            ("Beta", "positive"),
            ("Beta", "hard_negative"),
        ):
            acc = f"{clade[0]}{tier[0]}{i}"
            proteins.append({"accession": acc, "tier": tier, "clade_class": clade, "length": 100})
            clusters.append({"accession": acc, "cluster_representative": acc})
            sequences[acc] = "ACDEFGHIKL" * 10
    # one cluster with a member in each fold
    for acc, clade in (("SH_A", "Alpha"), ("SH_B", "Beta")):
        proteins.append({"accession": acc, "tier": "positive", "clade_class": clade, "length": 100})
        clusters.append({"accession": acc, "cluster_representative": "shared"})
        sequences[acc] = "MKVLAAWTYG" * 10
    folds = [
        {
            "fold": "Alpha",
            "level": "class",
            "members": "Alpha",
            "n_positive": 5,
            "n_hard_negative": 4,
            "scored": "yes",
        },
        {
            "fold": "Beta",
            "level": "class",
            "members": "Beta",
            "n_positive": 5,
            "n_hard_negative": 4,
            "scored": "yes",
        },
    ]
    return make_benchmark(tmp_path, proteins, clusters, folds, sequences)


# --------------------------------------------------------------- loading
def test_load_reads_only_scored_tiers(tmp_path):
    proteins = [
        {"accession": "P1", "tier": "positive", "clade_class": "Alpha", "length": 10},
        {"accession": "P2", "tier": "hard_negative", "clade_class": "Alpha", "length": 10},
        {"accession": "P3", "tier": "excluded", "clade_class": "Alpha", "length": 10},
        {"accession": "P4", "tier": "positive_moonlighting", "clade_class": "Alpha", "length": 10},
    ]
    clusters = [
        {"accession": p["accession"], "cluster_representative": p["accession"]} for p in proteins
    ]
    folds = [
        {
            "fold": "Alpha",
            "level": "class",
            "members": "Alpha",
            "n_positive": 1,
            "n_hard_negative": 1,
            "scored": "yes",
        }
    ]
    seqs = {p["accession"]: "ACDEFGHIKL" for p in proteins}
    d = make_benchmark(tmp_path, proteins, clusters, folds, seqs)
    bench = bd.load(d)
    assert set(bench.accessions) == {"P1", "P2"}
    assert list(bench.y) == [1, 0]


def test_load_drops_a_protein_whose_clade_has_no_fold(tmp_path):
    proteins = [
        {"accession": "P1", "tier": "positive", "clade_class": "Alpha", "length": 10},
        {"accession": "P2", "tier": "positive", "clade_class": "Nowhere", "length": 10},
    ]
    clusters = [
        {"accession": p["accession"], "cluster_representative": p["accession"]} for p in proteins
    ]
    folds = [
        {
            "fold": "Alpha",
            "level": "class",
            "members": "Alpha",
            "n_positive": 1,
            "n_hard_negative": 0,
            "scored": "yes",
        }
    ]
    seqs = {"P1": "ACDE", "P2": "ACDE"}
    bench = bd.load(make_benchmark(tmp_path, proteins, clusters, folds, seqs))
    assert bench.accessions == ["P1"]


def test_load_refuses_an_unknown_training_pool(simple):
    with pytest.raises(bd.BenchmarkError):
        bd.load(simple, train_pool="whatever")


def test_load_refuses_a_benchmark_that_fails_its_checksums(simple):
    (simple / "proteins.tsv").write_text(
        "accession\ttier\tclade_class\tlength\nX\tpositive\tAlpha\t1\n"
    )
    with pytest.raises(bd.BenchmarkError, match="SHA256SUMS"):
        bd.load(simple)


def test_load_refuses_a_benchmark_with_no_checksums(tmp_path):
    d = make_benchmark(tmp_path, [], [], [], {}, with_sums=False)
    with pytest.raises(bd.BenchmarkError, match="SHA256SUMS"):
        bd.load(d)


def test_no_verify_allows_an_unchecksummed_benchmark(tmp_path):
    proteins = [
        {"accession": "P1", "tier": "positive", "clade_class": "Alpha", "length": 10},
        {"accession": "P2", "tier": "hard_negative", "clade_class": "Alpha", "length": 10},
    ]
    clusters = [
        {"accession": p["accession"], "cluster_representative": p["accession"]} for p in proteins
    ]
    folds = [
        {
            "fold": "Alpha",
            "level": "class",
            "members": "Alpha",
            "n_positive": 1,
            "n_hard_negative": 1,
            "scored": "yes",
        }
    ]
    d = make_benchmark(
        tmp_path, proteins, clusters, folds, {"P1": "ACDE", "P2": "KLMN"}, with_sums=False
    )
    assert len(bd.load(d, verify=False)) == 2


# --------------------------------------------------------------- features
def test_composition_features_are_frequencies_then_log_length():
    got = bd.composition_features("AAAACCCC")
    assert got[bd.AA_ORDER.index("A")] == pytest.approx(0.5)
    assert got[bd.AA_ORDER.index("C")] == pytest.approx(0.5)
    assert got[20] == pytest.approx(np.log(8))
    assert got[:20].sum() == pytest.approx(1.0)


def test_composition_ignores_unknown_residues_rather_than_failing():
    got = bd.composition_features("AAXX")
    assert got[bd.AA_ORDER.index("A")] == pytest.approx(0.5)
    assert got[20] == pytest.approx(np.log(4))


def test_composition_refuses_an_empty_sequence():
    with pytest.raises(bd.BenchmarkError):
        bd.composition_features("")


# --------------------------------------------------------------- the strict cluster rule (D-A)
def test_no_cluster_is_ever_in_both_train_and_test(simple):
    bench = bd.load(simple)
    assert bd.assert_no_leakage(bench) is True


def test_the_spanning_cluster_is_dropped_from_training(simple):
    """'shared' has a member in each fold, so holding out either must drop the other member."""
    bench = bd.load(simple)
    for held, train, _test, dropped in bd.splits(bench):
        assert dropped == 1, f"fold {held} dropped {dropped}, expected the one shared member"
        assert "shared" not in set(np.array(bench.cluster)[train])


def test_without_the_strict_rule_the_spanning_cluster_would_leak(simple):
    """Shows the rule is doing work: the naive split puts 'shared' on both sides."""
    bench = bd.load(simple)
    fold = np.array(bench.fold)
    cluster = np.array(bench.cluster)
    for held in bench.scored_folds:
        naive_train = np.flatnonzero(fold != held)
        test = np.flatnonzero(fold == held)
        assert set(cluster[naive_train]) & set(cluster[test]) == {"shared"}


def test_assert_no_leakage_raises_when_the_rule_is_broken(simple, monkeypatch):
    bench = bd.load(simple)

    def naive(b):
        f = np.array(b.fold)
        for held in b.scored_folds:
            yield held, np.flatnonzero(f != held), np.flatnonzero(f == held), 0

    monkeypatch.setattr(bd, "splits", naive)
    with pytest.raises(bd.BenchmarkError, match="both train and test"):
        bd.assert_no_leakage(bench)


# --------------------------------------------------------------- the training pool (D-D)
def test_the_training_pool_changes_what_is_trained_on(tmp_path):
    """A protein in an unscored fold is in the 'all' pool and out of the 'scored' pool."""
    proteins, clusters, seqs = [], [], {}
    for clade in ("Alpha", "Beta", "Tiny"):
        for i, tier in enumerate(("positive", "hard_negative")):
            acc = f"{clade}{i}"
            proteins.append({"accession": acc, "tier": tier, "clade_class": clade, "length": 10})
            clusters.append({"accession": acc, "cluster_representative": acc})
            seqs[acc] = "ACDEFGHIKL"
    folds = [
        {
            "fold": c,
            "level": "class",
            "members": c,
            "n_positive": 1,
            "n_hard_negative": 1,
            "scored": "yes" if c != "Tiny" else "no",
        }
        for c in ("Alpha", "Beta", "Tiny")
    ]
    d = make_benchmark(tmp_path, proteins, clusters, folds, seqs)

    sizes = {}
    for pool in ("all", "scored"):
        bench = bd.load(d, train_pool=pool)
        assert bench.scored_folds == ["Alpha", "Beta"]
        sizes[pool] = {held: len(train) for held, train, _, _ in bd.splits(bench)}
    # 'all' carries the two Tiny proteins that 'scored' leaves out
    assert sizes["all"]["Alpha"] == sizes["scored"]["Alpha"] + 2


# --------------------------------------------------------------- facts about the REAL benchmark
@pytest.mark.skipif(not REAL.exists(), reason="adhesion_v1 is not present")
class TestRealBenchmark:
    """These pin the numbers the power study and the plan quote. They SHOULD break on a v2."""

    def test_the_frozen_benchmark_verifies_against_its_own_checksums(self):
        assert bd.verify_checksums(REAL)

    def test_three_folds_are_scorable(self):
        bench = bd.load(REAL)
        assert bench.scored_folds == ["Saccharomycetes", "Pezizomycotina", "Eurotiomycetes"]

    def test_the_strict_rule_holds_on_the_real_benchmark(self):
        assert bd.assert_no_leakage(bd.load(REAL)) is True

    def test_pooled_counts_are_50_positives_and_75_hard_negatives(self):
        bench = bd.load(REAL)
        mask = np.array([f in set(bench.scored_folds) for f in bench.fold])
        assert int(bench.y[mask].sum()) == 50
        assert int((bench.y[mask] == 0).sum()) == 75

    def test_the_leakage_drop_counts_depend_on_the_declared_training_pool(self):
        """The plan's D-A table is only reproducible with train_pool='all'; this pins both."""
        got = {}
        for pool in ("all", "scored"):
            bench = bd.load(REAL, train_pool=pool)
            got[pool] = {held: dropped for held, _tr, _te, dropped in bd.splits(bench)}
        assert got["all"] == {"Saccharomycetes": 4, "Pezizomycotina": 5, "Eurotiomycetes": 4}
        assert got["scored"] == {"Saccharomycetes": 4, "Pezizomycotina": 2, "Eurotiomycetes": 4}

    def test_no_mixed_label_cluster_inside_the_scorable_folds(self):
        """Corrects the plan's P2: the showcase hydrophobin pairs are outside C1 entirely."""
        bench = bd.load(REAL)
        scored = set(bench.scored_folds)
        labels = {}
        for y, f, c in zip(bench.y, bench.fold, bench.cluster, strict=False):
            if f in scored:
                labels.setdefault(c, set()).add(int(y))
        assert [c for c, ys in labels.items() if len(ys) > 1] == []
