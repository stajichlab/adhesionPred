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


def _one_positive_cluster(out):
    """S2-Calb_CGD training positives fall into one cluster (members and clusters agree)."""
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    hashes = set()
    for m in members:
        if m["split_id"] == "S2-Calb_CGD" and m["part"] == "train" and m["class"] == "pos":
            hashes.add(m["seq_sha256"])
    for m in members:  # clusters.tsv.gz has one cluster per hash, so change it in every split
        if m["seq_sha256"] in hashes:
            m["cluster_id"] = "one_cluster"
    truth_table.write_tsv(out / "split_members.tsv.gz", splits.MEMBER_COLUMNS, members)
    clusters = truth_table.read_tsv(out / "clusters.tsv.gz")
    for r in clusters:
        if r["seq_sha256"] in hashes:
            r["cluster_id"] = "one_cluster"
    truth_table.write_tsv(out / "clusters.tsv.gz", splits.CLUSTER_COLUMNS, clusters)
    for name in ("split_members.tsv.gz", "clusters.tsv.gz"):
        pf.fix_recorded_hash(out / "splits_run.json", name, out / name)


def test_threshold_fit_uses_train_only(phasec_chain, tmp_path):
    # permute the class of ALL rows that are not training rows (pos, neg, excluded, test_tc,
    # test_lit): every fitted setting and every score of all ten candidates stays the same
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    other = [m for m in members if m["part"] not in ("train", "train_tc")]
    before = [m["class"] for m in other]
    assert {"pos", "neg", "excluded"} <= set(before)
    assert {"test_tc", "test_lit"} <= {m["part"] for m in other}
    for m, c in zip(other, np.random.default_rng(0).permutation(before), strict=True):
        m["class"] = str(c)
    assert [m["class"] for m in other] != before  # the outer test labels really changed
    truth_table.write_tsv(out / "split_members.tsv.gz", splits.MEMBER_COLUMNS, members)
    pf.fix_recorded_hash(
        out / "splits_run.json", "split_members.tsv.gz", out / "split_members.tsv.gz"
    )
    assert _run10(fx["work"]) == 0
    # t, g, C, the ML threshold, Platt a and b and the H variant: all unchanged
    assert _units(fx["work"]) == _units(phasec_chain["work"])
    # the scaler is not stored; identical scores show that it did not change either
    a = (fx["work"] / "phasec" / "scores.tsv.gz").read_bytes()
    assert a == (phasec_chain["work"] / "phasec" / "scores.tsv.gz").read_bytes()


def _tasks(work, candidates=("M8",), clusters=True):
    out = work / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    order = json.loads((out / "splits_run.json").read_text())["splits"]
    hashes = [r["seq_sha256"] for r in truth_table.read_tsv(work / "phaseb" / "features_unique.tsv.gz")]  # fmt: skip
    cl = {r["seq_sha256"]: r["cluster_id"] for r in truth_table.read_tsv(out / "clusters.tsv.gz")}
    mod = load_phasec("10_fit_and_score")
    return mod, members, order, hashes, cl, candidates


def test_tasks_keep_training_and_test_rows_apart(phasec_chain):
    mod, members, order, hashes, cl, cands = _tasks(phasec_chain["work"])
    task_list = mod.tasks(members, order, hashes, cands, cl)
    assert len(task_list) == 2 * len(mod.unit_order(members, order))
    rows = {}
    for m in members:
        rows.setdefault((m["split_id"], m["fold"]), []).append(m)
    by_key = {t["key"]: t for t in task_list}
    for t in task_list:
        unit = rows[(t["split_id"], t["fold"])]
        train = set(t["train"])
        assert len(train) == len(t["train"])
        test = {m["seq_sha256"] for m in unit if m["part"] in ("test", "test_tc", "test_lit")}
        assert not train & test
        tc = {m["seq_sha256"] for m in unit if m["origin"] == "tc"}
        train_tc = {m["seq_sha256"] for m in unit if m["part"] == "train_tc"}
        if t["variant"] == "V-go":
            assert not train & tc  # V-go has no T-c rows
        else:
            go = set(by_key[f"{t['split_id']}|{t['fold']}|V-go"]["train"])
            assert train - go == train_tc and go <= train
    assert any(t["variant"] == "V-kw" and set(t["train"]) - set(by_key[
        t["key"].replace("V-kw", "V-go")]["train"]) for t in task_list)  # fmt: skip


def test_tasks_stop_on_overlap_and_on_wrong_clusters(phasec_chain):
    mod, members, order, hashes, cl, cands = _tasks(phasec_chain["work"])
    bad = {h: "x" + c for h, c in cl.items()}
    with pytest.raises(mod.evalio.StopError, match="differs from clusters.tsv.gz"):
        mod.tasks(members, order, hashes, cands, bad)
    s1 = [m for m in members if m["split_id"] == "S1" and m["fold"] == "0"]
    train = next(m for m in s1 if m["part"] == "train")
    leaked = [*members, {**train, "part": "test"}]
    with pytest.raises(mod.evalio.StopError, match="training and scored parts"):
        mod.tasks(leaked, order, hashes, cands, cl)


def test_a_hash_in_two_scored_parts_is_scored_once(phasec_chain):
    lit_hash = phasec_chain["named"]["LIT_AND_GO"]
    out = phasec_chain["work"] / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    parts = {m["part"] for m in members if m["split_id"] == "S3-Eurotiomycetes"
             and m["seq_sha256"] == lit_hash}  # fmt: skip
    assert parts == {"test", "test_lit"}  # the fixture really has the double row
    seen = {}
    for r in truth_table.read_tsv(out / "scores.tsv.gz"):
        key = (r["split_id"], r["fold"], r["variant"], r["candidate"], r["seq_sha256"])
        assert key not in seen, key
        seen[key] = r["part"]
    assert seen[("S3-Eurotiomycetes", "0", "V-go", "M8", lit_hash)] == "test,test_lit"


def test_result_for_another_unit_stops(phasec_chain, tmp_path):
    mod, members, order, hashes, cl, cands = _tasks(phasec_chain["work"])
    task_list = mod.tasks(members, order, hashes, cands, cl)[:2]
    results = [(t["key"], {}, {}) for t in task_list][::-1]
    with pytest.raises(mod.evalio.StopError, match="where S1|0|V-go was expected"):
        mod.write_scores(tmp_path / "s.tsv.gz", task_list, results)


def test_first_failure_cancels_the_other_units(monkeypatch):
    import concurrent.futures as cf

    mod = load_phasec("10_fit_and_score")
    made = []

    class FakePool:
        def __init__(self, **kw):
            self.cancel = None

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def submit(self, fn, task):
            f = cf.Future()
            if task["key"] == "bad":
                f.set_exception(mod.models.ModelError("bad: broken"))
            made.append(f)
            return f

        def shutdown(self, wait=True, cancel_futures=False):
            self.cancel = cancel_futures
            for f in made:
                if cancel_futures:
                    f.cancel()

    pool = FakePool()
    monkeypatch.setattr(mod, "ProcessPoolExecutor", lambda **kw: pool)
    with pytest.raises(mod.models.ModelError, match="bad: broken"):
        mod.run_tasks([{"key": "bad"}, {"key": "b"}, {"key": "c"}], None, 2, None)
    assert pool.cancel is True and [f.cancelled() for f in made] == [False, True, True]


def test_failing_unit_with_two_workers_names_the_unit(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    _one_positive_cluster(out)
    for name in ("scores.tsv.gz", "scores_run.json"):
        (out / name).unlink()
    assert _run10(fx["work"], "--candidates", "R2,M8", "--workers", "2") == 2
    assert "S2-Calb_CGD|0|V-go: 1 positive and" in capsys.readouterr().err
    assert not (out / "scores.tsv.gz").exists() and not (out / "scores_run.json").exists()


def test_split_run_from_another_build_stops(phasec_chain, tmp_path, capsys):
    # 08 output changed and build_run.json re-pinned, but 09 was not re-run
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    table = out / "eval_table.tsv.gz"
    table.write_bytes(table.read_bytes() + b"\n")
    pf.fix_recorded_hash(out / "build_run.json", "eval_table.tsv.gz", table)
    for name in ("scores.tsv.gz", "scores_run.json"):
        (out / name).unlink()
    assert _run10(fx["work"]) == 2
    err = capsys.readouterr().err
    assert "another eval_table.tsv.gz than build_run.json records" in err
    assert not (out / "scores.tsv.gz").exists()


def test_missing_key_in_a_run_file_is_named(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    path = fx["work"] / "phasec" / "build_run.json"
    obj = json.loads(path.read_text())
    del obj["truth_set_sha256"]
    path.write_text(json.dumps(obj))
    assert _run10(fx["work"], "--candidates", "R0") == 2
    assert "STOP: missing key 'truth_set_sha256' in an input or run file" in capsys.readouterr().err


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
    _one_positive_cluster(out)
    assert _run10(fx["work"], "--candidates", "R2,M8") == 2
    err = capsys.readouterr().err
    assert "S2-Calb_CGD|0|V-go: 1 positive and" in err and "clusters" in err
