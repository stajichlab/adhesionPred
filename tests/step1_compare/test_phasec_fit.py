"""10_fit_and_score.py on the fixture chain (08 -> 09 with stub MMseqs2 -> 10)."""

import json

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import phasec_fixture as pf  # noqa: E402
import splits  # noqa: E402
import truth_table  # noqa: E402
from conftest import load_phasec  # noqa: E402

SUBSET = "R2,M8,H"  # re-runs fit only these candidates (time); the rest is fitted in the chain


def _run10(work, *extra):
    return load_phasec("10_fit_and_score").main(["--work-dir", str(work), *extra])


def _units(work):
    return json.loads((work / "phasec" / "scores_run.json").read_text())["units"]


def _score_rows(work):
    rows = truth_table.read_tsv(work / "phasec" / "scores.tsv.gz")
    return [r for r in rows if r["candidate"] in SUBSET.split(",")]


def _subset(units):
    return {k: {c: v[c] for c in SUBSET.split(",")} for k, v in units.items()}


def test_scores_cover_every_unit_candidate_and_scored_row(phasec_chain):
    out = phasec_chain["work"] / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    units = {(m["split_id"], m["fold"]) for m in members}
    assert set(_units(phasec_chain["work"])) == {
        f"{s}|{f}|{v}" for s, f in units for v in ("V-go", "V-kw")
    }
    rows = truth_table.read_tsv(out / "scores.tsv.gz")
    n_unique = len(truth_table.read_tsv(phasec_chain["work"] / "phaseb" / "features_unique.tsv.gz"))
    full = [
        r
        for r in rows
        if r["split_id"] == "FULL" and r["variant"] == "V-go" and r["candidate"] == "M8"
    ]
    assert len(full) == n_unique and {r["part"] for r in full} == {"all"}
    s1 = {
        (r["seq_sha256"], r["part"])
        for r in rows
        if r["split_id"] == "S1" and r["candidate"] == "R2"
    }
    want = {(m["seq_sha256"], m["part"]) for m in members
            if m["split_id"] == "S1" and m["part"] in ("test", "test_tc")}  # fmt: skip
    assert s1 == want
    assert all(r["score"] == "" and r["prob"] == "" and r["call"] in ("0", "1")
               for r in rows if r["candidate"] == "R2")  # fmt: skip
    assert all(r["score"] != "" and 0.0 <= float(r["prob"]) <= 1.0
               for r in rows if r["candidate"] == "H")  # fmt: skip


def test_fitted_settings_are_recorded(phasec_chain):
    units = _units(phasec_chain["work"])
    unit = units["S1|0|V-go"]
    assert unit["M8"]["C"] in (0.001, 0.003, 0.01, 0.1, 1.0, 10.0)  # ruling C-12
    assert unit["H"]["h_variant"] in ("M8", "M35", "M8-C", "M35-C")
    assert unit["R2"]["g"] in ("highly_probable", "probable", "weakly")
    assert unit["R2"]["t"] in (0.20, 0.25, 0.30, 0.35, 0.40)  # owner decision 2026-10-01
    assert set(unit["M35"]["inner_pr_auc"]) == {"0.001", "0.003", "0.01", "0.1", "1.0", "10.0"}
    assert len(unit["H"]["inner_pr_auc"]) == 4 * 6  # (ESM variant, C) pairs
    assert units["S1|0|V-kw"]["M8"]["n_train_pos"] > unit["M8"]["n_train_pos"]  # T-c positives


def test_training_rows_keep_homology_only_rows(phasec_chain):
    # spec 3.2: the direct-evidence filter applies to test rows only (11 applies it)
    out = phasec_chain["work"] / "phasec"
    table = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "eval_table.tsv.gz")}
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    train = [m for m in members if m["split_id"] == "FULL" and m["part"] == "train"]
    assert any(table[m["seq_sha256"]]["homology_only"] == "yes" for m in train)
    assert _units(phasec_chain["work"])["FULL|0|V-go"]["B0"]["n_train"] == len(train)


def test_threshold_fit_uses_train_only(phasec_chain, tmp_path):
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    test = [
        m for m in members if m["part"] in ("test", "test_lit") and m["class"] in ("pos", "neg")
    ]
    before = [m["class"] for m in test]
    for m, c in zip(test, np.random.default_rng(0).permutation(before), strict=True):
        m["class"] = str(c)
    assert [m["class"] for m in test] != before  # the outer test labels really changed
    truth_table.write_tsv(out / "split_members.tsv.gz", splits.MEMBER_COLUMNS, members)
    pf.fix_recorded_hash(
        out / "splits_run.json", "split_members.tsv.gz", out / "split_members.tsv.gz"
    )
    assert _run10(fx["work"], "--candidates", SUBSET) == 0
    # t, g, C, the ML threshold, Platt a and b and the H variant: all unchanged
    assert _units(fx["work"]) == _subset(_units(phasec_chain["work"]))
    # the scaler is not stored; identical scores show that it did not change either
    assert _score_rows(fx["work"]) == _score_rows(phasec_chain["work"])


def test_workers_give_the_same_scores(phasec_chain, tmp_path):
    one = pf.copy_work(phasec_chain, tmp_path / "one")
    two = pf.copy_work(phasec_chain, tmp_path / "two")
    assert _run10(one["work"], "--candidates", SUBSET) == 0
    assert _run10(two["work"], "--candidates", SUBSET, "--workers", "2") == 0
    a = (one["work"] / "phasec" / "scores.tsv.gz").read_bytes()
    assert (two["work"] / "phasec" / "scores.tsv.gz").read_bytes() == a
    assert _units(one["work"]) == _subset(_units(phasec_chain["work"]))


def test_unknown_candidate_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    assert _run10(fx["work"], "--candidates", "M8,XGB") == 2
    assert "unknown candidates ['XGB']" in capsys.readouterr().err


def test_stale_input_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    path = fx["work"] / "phasec" / "clusters.tsv.gz"
    path.write_bytes(path.read_bytes() + b"\n")
    before = (fx["work"] / "phasec" / "scores.tsv.gz").read_bytes()
    assert _run10(fx["work"]) == 2
    assert "clusters.tsv.gz differs from the SHA-256 in splits_run.json" in capsys.readouterr().err
    assert (fx["work"] / "phasec" / "scores.tsv.gz").read_bytes() == before  # not replaced


def test_changed_embedding_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    emb = fx["work"] / "phaseb" / "emb" / "esm2_t6_8M_UR50D.nterm.npy"
    arr = np.load(emb)
    arr[0, 0] += 1.0
    np.save(emb, arr)
    assert _run10(fx["work"]) == 2
    assert "differs from embedding_run.json" in capsys.readouterr().err


def test_too_few_positive_clusters_stops(phasec_chain, tmp_path, capsys):
    # Review Focus 5: an S2 training set whose positives fall into 2 clusters
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    for m in members:
        if m["split_id"] == "S2-Calb_CGD" and m["part"] == "train" and m["class"] == "pos":
            m["cluster_id"] = "one_cluster"
    truth_table.write_tsv(out / "split_members.tsv.gz", splits.MEMBER_COLUMNS, members)
    pf.fix_recorded_hash(
        out / "splits_run.json", "split_members.tsv.gz", out / "split_members.tsv.gz"
    )
    assert _run10(fx["work"], "--candidates", "R2,M8") == 2
    err = capsys.readouterr().err
    assert "S2-Calb_CGD|0|V-go: 1 positive and" in err and "clusters" in err
