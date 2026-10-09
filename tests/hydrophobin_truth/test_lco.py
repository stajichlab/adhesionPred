import importlib.util
import shutil
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/lco.py"
spec = importlib.util.spec_from_file_location("lco", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def aln(n_cys_cols, n_seq=10, extra=3, drop_frac=0.0):
    cols = ["C"] * n_cys_cols + ["A"] * extra
    seqs = []
    for i in range(n_seq):
        s = list(cols)
        if i < int(drop_frac * n_seq):
            s[0] = "-"
        seqs.append("".join(s))
    return seqs


def test_conserved_cysteine_columns_counts_columns_at_or_above_the_fraction():
    assert m.conserved_cys_columns(aln(8), 0.8) == 8
    assert (
        m.conserved_cys_columns(aln(8, drop_frac=0.3), 0.8) == 7
    )  # first column has C in only 70% of sequences


def test_check_refuses_an_alignment_with_fewer_than_eight_conserved_columns():
    assert m.cys_check(aln(8), 0.8)
    assert not m.cys_check(aln(7), 0.8)


def test_training_ids_exclude_the_held_out_cluster():
    clusters = {"a": "c1", "b": "c1", "c": "c2", "d": "c3"}
    assert m.training_ids(clusters, "c1") == ["c", "d"]
    assert m.training_ids(clusters, None) == ["a", "b", "c", "d"]


@pytest.mark.skipif(
    shutil.which("mafft") is None or shutil.which("hmmbuild") is None,
    reason="mafft or hmmer not on PATH",
)
def test_build_runs_mafft_and_hmmbuild_on_a_small_synthetic_family(tmp_path):
    base = "MKLLVAAGCAADDEEFFGGCCHHIIKKLLCMMNNPPQQCRRSSTTVVCWWYYACDDEEFFCGGHHIICKKLLMMNNPPQQRRSSTTVVWWYYC"
    seqs = {f"s{i}": base[:i] + "A" + base[i + 1 :] if base[i] != "C" else base for i in range(6)}
    out = m.build_hmm(seqs, tmp_path, "toy", threads=1)
    assert out["hmm"].exists() and out["hmm"].read_text().startswith("HMMER3")
    assert out["n_cys_columns"] >= 8


def test_the_aligner_command_is_l_ins_i_for_every_fold():
    assert m.MAFFT_ARGS[:2] == ["--localpair", "--maxiterate"]


def test_ordinal_check_accepts_a_split_eighth_cysteine_column():
    # 6 class I-like sequences with C8 at column 20, 4 class II-like with C8 at column 14: same ordinal cysteine, nearby columns
    base = ["-"] * 30

    def seq(cols):
        s = list(base)
        for j in cols:
            s[j] = "C"
        return "".join(s)

    c1 = [2, 4, 5, 9, 12, 15, 16]
    aln = [seq(c1 + [20]) for _ in range(6)] + [seq(c1 + [14]) for _ in range(4)]
    assert m.conserved_cys_columns(aln, 0.8) == 7  # the column rule fails on this alignment
    assert m.ordinal_cys_check(aln, window=10, min_frac=0.8)


def test_ordinal_check_refuses_sequences_with_scattered_cysteines():
    base = ["-"] * 60

    def seq(cols):
        s = list(base)
        for j in cols:
            s[j] = "C"
        return "".join(s)

    good = [2, 4, 5, 9, 12, 15, 16, 20]
    bad = [2, 30, 35, 40, 45, 50, 52, 58]
    aln = [seq(good)] * 6 + [seq(bad)] * 4
    assert not m.ordinal_cys_check(aln, window=10, min_frac=0.8)


def test_ordinal_check_counts_sequences_with_fewer_than_eight_cysteines_as_inconsistent():
    base = ["-"] * 30

    def seq(cols):
        s = list(base)
        for j in cols:
            s[j] = "C"
        return "".join(s)

    full = seq([2, 4, 5, 9, 12, 15, 16, 20])
    short = seq([2, 4, 5, 9, 12, 15])
    assert m.ordinal_cys_check([full] * 9 + [short], window=10, min_frac=0.8)
    assert not m.ordinal_cys_check([full] * 7 + [short] * 3, window=10, min_frac=0.8)


def test_slots_check_accepts_a_cysteine_split_over_distant_columns_and_refuses_a_garbled_alignment():
    def seq(cols, width=80):
        s = ["-"] * width
        for j in cols:
            s[j] = "C"
        return "".join(s)

    c = [2, 4, 5, 9, 12, 15, 16]
    split = [seq(c + [20])] * 6 + [
        seq(c + [60])
    ] * 4  # C8 split between two distant columns, 60% and 40%
    assert m.cys_slots_check(split)
    garbled = [
        seq([10 * k + i for k in range(8)]) for i in range(10)
    ]  # every column holds a cysteine in 1 of 10 sequences
    assert not m.cys_slots_check(garbled)


def test_aligner_commands():
    cmd, out = m.aligner_command("mafft", "x.faa", 4)
    assert cmd[0] == "mafft" and "--localpair" in cmd and out is None
    cmd, out = m.aligner_command("famsa", "x.faa", 4)
    assert cmd[:3] == ["famsa", "-t", "4"] and cmd[-1] == "-" and out is None
    cmd, out = m.aligner_command("muscle5", "x.faa", 4)
    assert cmd[:2] == ["muscle", "-align"] and out == "x.faa.muscle.afa"
    try:
        m.aligner_command("clustalo", "x.faa", 1)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")
