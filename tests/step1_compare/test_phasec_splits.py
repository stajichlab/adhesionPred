"""splits.py and 09_make_splits.py: clusters, S1-S3, T-c removals, max identity (spec 3.4)."""

import json
import os
import shutil
import subprocess

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import paths  # noqa: E402
import phasec_fixture as pf  # noqa: E402
import splits  # noqa: E402
import truth_table  # noqa: E402
from conftest import load_phasec  # noqa: E402

STUB = paths.STEP1_DIR.parents[1] / "tests/step1_compare/fixtures/phasec/stub_mmseqs/mmseqs"
SPECIES = truth_table.read_tsv(paths.STEP1_DIR / "species.tsv")


def go(h, cls, source, gene=None, **kw):
    return {"seq_sha256": h, "origin": "go", "class": cls, "source_ids": source,
            "gene_ids": gene or f"g_{h}", "taxon_ids": "", "clades": "", **kw}  # fmt: skip


def tc(h, acc, taxon="1", clade="Saccharomycotina"):
    return {"seq_sha256": h, "origin": "tc", "class": "pos", "source_ids": "T-c",
            "gene_ids": acc, "taxon_ids": taxon, "clades": clade}  # fmt: skip


def test_read_cluster_tsv():
    got = splits.read_cluster_tsv(["a\ta\n", "a\tb\n", "c\tc\n", "\n"])
    assert got == {"a": "a", "b": "a", "c": "c"}
    with pytest.raises(splits.SplitError, match="in two clusters"):
        splits.read_cluster_tsv(["a\ta", "c\ta"])
    with pytest.raises(splits.SplitError, match="expected 2 fields"):
        splits.read_cluster_tsv(["a b"])
    with pytest.raises(splits.SplitError, match="1 sequences without a cluster"):
        splits.check_clusters({"a": "a"}, ["a", "z"])


def test_split_definitions_from_species_tsv():
    defs = {d.split_id: d for d in splits.split_definitions(SPECIES)}
    assert list(defs) == ["S1", "S2-Calb_CGD", "S2-Scer_SGD", "S2-Spom_PomBase",
                          "S3-Basidiomycota", "S3-Eurotiomycetes", "FULL"]  # fmt: skip
    assert defs["S2-Calb_CGD"].train_sources == ("Scer_SGD",)
    assert defs["S2-Calb_CGD"].test_taxa == frozenset({"237561"})
    assert defs["S3-Eurotiomycetes"].test_sources == ("Afum_ASPFU", "Anid_EMENI")
    assert defs["S3-Eurotiomycetes"].literature and not defs["S3-Basidiomycota"].literature
    assert defs["S3-Basidiomycota"].test_sources == ("Cneo_H99_GOA", "Umay_MYCMD")
    tested = {s for d in defs.values() for s in d.test_sources}
    assert not tested & {"Spom_SCHPO-mod", "Cneo_JEC21_GOA", "Cneo_CRYD1"}  # alternate files


def _toy(n_per_class=12):
    table, cluster_of = [], {}
    for src in ("Scer_SGD", "Calb_CGD", "Spom_PomBase", "Afum_ASPFU"):
        for i in range(n_per_class):
            for cls in ("pos", "neg"):
                h = f"{src[:4]}{cls}{i:02d}"
                table.append(go(h, cls, src))
                cluster_of[h] = f"cl_{h}" if i % 4 else f"fam{i % 3}_{cls}"
        h = f"{src[:4]}excl"
        table.append(go(h, "excluded", src))
        cluster_of[h] = f"cl_{h}"
    for i in range(6):
        h = f"tc{i:02d}"
        table.append(tc(h, f"Q{i}"))
        cluster_of[h] = f"cl_{h}"
    return table, cluster_of


def _species():
    return [
        {"source_id": "Scer_SGD", "taxon_id": "559292", "in_clade": "Saccharomycotina", "role": "train"},
        {"source_id": "Calb_CGD", "taxon_id": "237561", "in_clade": "Saccharomycotina", "role": "train"},
        {"source_id": "Spom_PomBase", "taxon_id": "284812", "in_clade": "Taphrinomycotina", "role": "test_species"},
        {"source_id": "Afum_ASPFU", "taxon_id": "330879", "in_clade": "Eurotiomycetes", "role": "test_clade"},
    ]  # fmt: skip


def _build(table, cluster_of, lit=()):
    return splits.build(table, list(lit), cluster_of, splits.split_definitions(_species()), 1)


def test_s1_folds_keep_clusters_whole_and_cover_tc_rows():
    table, cluster_of = _toy()
    members, _ = _build(table, cluster_of)
    s1 = [m for m in members if m["split_id"] == "S1"]
    fold_of = {}
    for m in s1:
        if m["part"] in ("test", "test_tc"):
            assert m["seq_sha256"] not in fold_of
            fold_of[m["seq_sha256"]] = m["fold"]
    by_cluster = {}
    for h, f in fold_of.items():
        by_cluster.setdefault(cluster_of[h], set()).add(f)
    assert all(len(v) == 1 for v in by_cluster.values())
    assert {f"tc{i:02d}" for i in range(6)} <= set(fold_of)  # every T-c row has an oof fold
    assert "Scerexcl" in fold_of and "Spompos00" not in fold_of
    splits.check_no_cluster_spans(members)
    splits.check_no_shared_hash(members)


def test_no_cluster_spans_train_and_test():
    table, cluster_of = _toy()
    members, _ = _build(table, cluster_of)
    s1 = [dict(m) for m in members if m["split_id"] == "S1" and m["fold"] == "0"]
    test = next(m for m in s1 if m["part"] == "test")
    train = next(m for m in s1 if m["part"] == "train")
    train["cluster_id"] = test["cluster_id"]  # plant one spanning cluster
    with pytest.raises(splits.SplitError, match="spans train and test"):
        splits.check_no_cluster_spans(s1)


def test_s2_s3_allow_go_homologs_but_not_tc_cluster_mates():
    table, cluster_of = _toy()
    members, _ = _build(table, cluster_of)
    s3 = [m for m in members if m["split_id"] == "S3-Eurotiomycetes"]
    # fam clusters hold Scer, Calb and Afum rows: GO train and test share them by design
    train_c = {m["cluster_id"] for m in s3 if m["part"] == "train"}
    test_c = {m["cluster_id"] for m in s3 if m["part"] == "test"}
    assert train_c & test_c
    splits.check_no_cluster_spans(members)  # does not raise
    bad = [dict(m) for m in s3]
    next(m for m in bad if m["part"] == "train_tc")["cluster_id"] = next(iter(test_c))
    with pytest.raises(splits.SplitError, match="spans train and test"):
        splits.check_no_cluster_spans(bad)


def test_test_truth_sources():
    table, cluster_of = _toy()
    members, _ = _build(table, cluster_of)
    splits.check_test_truth(members)
    assert not [m for m in members if m["origin"] == "tc" and m["part"] in ("test", "test_lit")]
    routed = [dict(m) for m in members]
    next(m for m in routed if m["origin"] == "tc")["part"] = "test"
    with pytest.raises(splits.SplitError, match="T-c row .* is in a test set"):
        splits.check_test_truth(routed)


def test_tc_excludes_test_proteins():
    table, cluster_of = _toy()
    table.append(tc("Afumpos00", "Q99"))  # same hash as an A. fumigatus test protein
    table.append(tc("tcacc", "g_Afumneg03"))  # accession of a test protein
    cluster_of["tcacc"] = "cl_tcacc"
    # the hash row would be dropped by 08 (GO label wins); here it tests rule a directly
    kept, removed = splits.tc_removals(
        [r for r in table if r["origin"] == "tc"],
        [r for r in table if r["origin"] == "go" and r["source_ids"] == "Afum_ASPFU"],
        cluster_of,
    )
    rules = {r["gene_ids"]: rule for r, rule in removed}
    assert rules == {"Q99": "a_test_protein", "g_Afumneg03": "a_test_protein"}
    assert len(kept) == 6


def test_tc_excludes_test_cluster_mates():
    table, cluster_of = _toy()
    cluster_of["tc00"] = cluster_of["Spompos05"]  # T-c row in the cluster of a S. pombe protein
    members, removed = _build(table, cluster_of)
    got = [(r["split_id"], r["rule"]) for r in removed if r["seq_sha256"] == "tc00"]
    assert got == [("S2-Spom_PomBase", "b_cluster_mate")]
    assert any(m["seq_sha256"] == "tc00" and m["part"] == "train_tc" for m in members
               if m["split_id"] == "S2-Calb_CGD")  # fmt: skip


def test_tc_excludes_test_taxa():
    table, cluster_of = _toy()
    for h, taxon, clade in (
        ("tcalb", "237561", "Saccharomycotina"),
        ("tcfum", "330879", "Eurotiomycetes"),
    ):
        table.append(tc(h, f"A_{h}", taxon, clade))
        cluster_of[h] = f"cl_{h}"
    members, removed = _build(table, cluster_of)
    got = {(r["split_id"], r["seq_sha256"]): r["rule"] for r in removed}
    assert got[("S2-Calb_CGD", "tcalb")] == "c_test_taxon"  # S2: taxon of the test species
    assert got[("S3-Eurotiomycetes", "tcfum")] == "c_test_taxon"  # S3: taxon in the test clade
    assert ("S2-Scer_SGD", "tcalb") not in got and ("S2-Calb_CGD", "tcfum") not in got
    s1_tc = {m["seq_sha256"] for m in members if m["split_id"] == "S1" and m["part"] == "train_tc"}
    assert {"tcalb", "tcfum"} <= s1_tc  # S1 keeps them (spec 4, context note)


def test_hash_in_a_training_and_a_test_source_stops():
    # Review Focus 4: one sequence labelled in S. cerevisiae and in S. pombe
    table, cluster_of = _toy()
    table.append(go("both", "neg", "Scer_SGD,Spom_PomBase"))
    cluster_of["both"] = "cl_both"
    with pytest.raises(splits.SplitError, match="belongs to a training source and a test source"):
        _build(table, cluster_of)


def test_literature_rows_are_tested_in_s3_eurotiomycetes_only():
    table, cluster_of = _toy()
    lit = [{"accession": "P1", "seq_sha256": "lit1", "literature_positive": "yes"},
           {"accession": "P2", "seq_sha256": "lit2", "literature_positive": "no"},
           {"accession": "P3", "seq_sha256": "", "literature_positive": "no"}]  # fmt: skip
    cluster_of.update({"lit1": "cl_lit1", "lit2": "cl_lit2"})
    members, _ = _build(table, cluster_of, lit)
    got = {(m["split_id"], m["seq_sha256"], m["class"]) for m in members if m["part"] == "test_lit"}
    assert got == {("S3-Eurotiomycetes", "lit1", "pos"), ("S3-Eurotiomycetes", "lit2", "excluded")}


def test_max_identity_excludes_self_hits():
    lines = ["q1\tq1\t1.000", "q1\tt1\t0.420", "q1\tt2\t0.610", "q2\tt1\t0.250", "q3\tq3\t1.0"]
    best = splits.parse_hits(lines)
    assert best == {"q1": 0.61, "q2": 0.25}
    rows = {r["seq_sha256"]: r for r in splits.identity_rows("S2-x", ["q1", "q2", "q3"], best)}
    assert rows["q1"]["max_identity"] == "0.6100" and rows["q1"]["below_0.3"] == "no"
    assert rows["q2"]["below_0.3"] == "yes"
    assert rows["q3"]["max_identity"] == "" and rows["q3"]["below_0.3"] == "yes"  # self hit only


# --- 09_make_splits.py with the stub MMseqs2 ---


def _built(tmp_path):
    fx = pf.make_work(tmp_path)
    assert load_phasec("08_build_eval_tables").main(pf.build_argv(fx)) == 0
    return fx


def _run09(fx, tmp_path, mmseqs=STUB):
    argv = ["--work-dir", str(fx["work"]), "--species", str(fx["species"]),
            "--mmseqs", str(mmseqs), "--tmp-dir", str(tmp_path / "scratch")]  # fmt: skip
    return load_phasec("09_make_splits").main(argv)


def test_make_splits_with_the_stub(tmp_path, monkeypatch):
    fx = _built(tmp_path)
    log_file = tmp_path / "stub.log"
    monkeypatch.setenv("STUB_MMSEQS_LOG", str(log_file))
    assert _run09(fx, tmp_path) == 0
    out = fx["work"] / "phasec"
    run = json.loads((out / "splits_run.json").read_text())
    assert run["mmseqs_version"] == "stub-17" and run["seed"] == 20261001
    calls = [line.split("\t") for line in log_file.read_text().splitlines()]
    cluster = next(c for c in calls if c[0] == "easy-cluster")
    assert cluster[4:10] == ["--min-seq-id", "0.3", "-c", "0.5", "--cov-mode", "0"]
    search = [c for c in calls if c[0] == "easy-search"]
    assert len(search) == 5  # S2-Calb_CGD, S2-Scer_SGD, S2-Spom_PomBase, S3-Basidio, S3-Euro
    assert search[0][5:13] == ["-s", "7.5", "-c", "0.5", "--cov-mode", "0", "--format-output",
                               "query,target,fident"]  # fmt: skip
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    assert {m["split_id"] for m in members} == set(run["splits"])
    s1_test = {m["seq_sha256"] for m in members if m["split_id"] == "S1" and m["part"] == "test"}
    table = truth_table.read_tsv(out / "eval_table.tsv.gz")
    train_go = {r["seq_sha256"] for r in table if r["origin"] == "go" and r["roles"] == "train"}
    assert s1_test == train_go  # every training-source GO row is scored out of fold once
    ident = truth_table.read_tsv(out / "max_identity.tsv.gz")
    assert {r["split_id"] for r in ident} == {
        "S2-Calb_CGD",
        "S2-Scer_SGD",
        "S2-Spom_PomBase",
        "S3-Basidiomycota",
        "S3-Eurotiomycetes",
    }


def test_make_splits_is_deterministic(tmp_path):
    fx = _built(tmp_path)
    assert _run09(fx, tmp_path) == 0
    out = fx["work"] / "phasec"
    first = {n: (out / n).read_bytes() for n in ("clusters.tsv.gz", "split_members.tsv.gz")}
    assert _run09(fx, tmp_path) == 0
    assert first == {n: (out / n).read_bytes() for n in first}


def test_stale_input_stops(tmp_path, capsys):
    fx = _built(tmp_path)
    table = fx["work"] / "phasec" / "eval_table.tsv.gz"
    table.write_bytes(table.read_bytes() + b"\n")
    assert _run09(fx, tmp_path) == 2
    assert "eval_table.tsv.gz differs from the SHA-256 in build_run.json" in capsys.readouterr().err
    assert not (fx["work"] / "phasec" / "splits_run.json").exists()


def test_mmseqs_failure_stops_without_output(tmp_path, monkeypatch, capsys):
    fx = _built(tmp_path)
    monkeypatch.setenv("STUB_MMSEQS_FAIL", "easy-search")
    assert _run09(fx, tmp_path) == 2
    assert "easy-search exited with 1" in capsys.readouterr().err
    assert not (fx["work"] / "phasec" / "split_members.tsv.gz").exists()


def test_mmseqs_version_failure_names_avx2(tmp_path, capsys):
    fx = _built(tmp_path)
    bad = tmp_path / "mmseqs"
    bad.write_text("#!/bin/bash\nexit 132\n")
    bad.chmod(0o755)
    assert _run09(fx, tmp_path, bad) == 2
    assert "AVX2" in capsys.readouterr().err


def _real_mmseqs():
    exe = os.environ.get("STEP1_MMSEQS") or shutil.which("mmseqs")
    if not exe:
        return None
    try:
        ok = subprocess.run([exe, "version"], capture_output=True).returncode == 0
    except OSError:
        return None
    return exe if ok else None


def _plant_homolog(fx, out):
    """Make one S. pombe test sequence a ~80% identical copy of an S. cerevisiae training one."""
    import gzip

    table = truth_table.read_tsv(out / "eval_table.tsv.gz")
    train = next(r["seq_sha256"] for r in table if r["source_ids"] == "Scer_SGD"
                 and r["class"] == "neg")  # fmt: skip
    test = next(r["seq_sha256"] for r in table if r["source_ids"] == "Spom_PomBase")
    path = out / "eval_sequences.fasta.gz"
    seqs, head = {}, None
    for line in gzip.decompress(path.read_bytes()).decode().splitlines():
        if line.startswith(">"):
            head = line[1:].strip()
            seqs[head] = ""
        else:
            seqs[head] += line.strip()
    src = seqs[train]
    seqs[test] = "".join("A" if i % 5 == 4 and c != "A" else c for i, c in enumerate(src))
    text = "".join(f">{h}\n{q}\n" for h, q in seqs.items())
    path.write_bytes(gzip.compress(text.encode(), mtime=0))
    pf.fix_recorded_hash(out / "build_run.json", "eval_sequences.fasta.gz", path)


@pytest.mark.skipif(
    _real_mmseqs() is None, reason="needs a working MMseqs2 (module or STEP1_MMSEQS)"
)
def test_make_splits_with_real_mmseqs(tmp_path):
    fx = _built(tmp_path)
    out = fx["work"] / "phasec"
    _plant_homolog(fx, out)
    assert _run09(fx, tmp_path, _real_mmseqs()) == 0
    clusters = truth_table.read_tsv(out / "clusters.tsv.gz")
    seqs = pf.gz_lines(out / "eval_sequences.fasta.gz")
    assert len(clusters) == sum(x.startswith(">") for x in seqs)
    run = json.loads((out / "splits_run.json").read_text())
    assert run["mmseqs_version"] and run["clusters"] >= 1
    ident = [r for r in truth_table.read_tsv(out / "max_identity.tsv.gz") if r["max_identity"]]
    assert any(0.3 < float(r["max_identity"]) <= 1.0 for r in ident)


# --- fix round 1 ---


def _wrapper(tmp_path, keep):
    """A mmseqs wrapper that copies the easy-search FASTA files to `keep`, then runs the stub."""
    keep.mkdir()
    wrap = tmp_path / "wrap_mmseqs"
    wrap.write_text(
        f'#!/bin/bash\nif [ "$1" = easy-search ]; then cp "$2" "{keep}/$(basename "$2")"; '
        f'cp "$3" "{keep}/$(basename "$3")"; fi\nexec "{STUB}" "$@"\n'
    )
    wrap.chmod(0o755)
    return wrap


def _fasta_ids(path):
    return {x[1:].strip() for x in path.read_text().splitlines() if x.startswith(">")}


def test_search_targets_are_train_part_and_queries_include_literature(tmp_path):
    fx = _built(tmp_path)
    keep = tmp_path / "keep"
    assert _run09(fx, tmp_path, _wrapper(tmp_path, keep)) == 0
    members = truth_table.read_tsv(fx["work"] / "phasec" / "split_members.tsv.gz")
    for split_id in ("S2-Calb_CGD", "S3-Eurotiomycetes", "S3-Basidiomycota"):
        rows = [m for m in members if m["split_id"] == split_id]
        assert any(m["part"] == "train_tc" for m in rows)  # the fixture has T-c training rows
        train = {m["seq_sha256"] for m in rows if m["part"] == "train"}
        query = {m["seq_sha256"] for m in rows if m["part"] in ("test", "test_lit")}
        assert _fasta_ids(keep / f"{split_id}.t.fasta") == train
        assert _fasta_ids(keep / f"{split_id}.q.fasta") == query
    lit = {m["seq_sha256"] for m in members
           if m["split_id"] == "S3-Eurotiomycetes" and m["part"] == "test_lit"}  # fmt: skip
    assert lit and lit <= _fasta_ids(keep / "S3-Eurotiomycetes.q.fasta")


def test_tc_rule_order_a_before_c_and_literature_accession():
    table, cluster_of = _toy()
    row = tc("tcboth", "Q77", "237561", "Saccharomycotina")  # taxon rule c would also match
    cluster_of.update({"tcboth": "cl_tcboth", "litx": "cl_litx"})
    lit = [{"seq_sha256": "litx", "origin": "lit", "class": "pos", "gene_ids": "Q77"}]
    kept, removed = splits.tc_removals([row], lit, cluster_of, frozenset({"237561"}))
    assert not kept and [rule for _, rule in removed] == ["a_test_protein"]
    same_taxon = tc("tcboth", "Q78", "237561", "Saccharomycotina")
    test = [r for r in table if r["source_ids"] == "Calb_CGD"][:1] + [
        {**go("t1", "pos", "Calb_CGD"), "gene_ids": "Q78"}
    ]
    cluster_of["t1"] = "cl_t1"
    _, removed = splits.tc_removals([same_taxon], test, cluster_of, frozenset({"237561"}))
    assert [rule for _, rule in removed] == ["a_test_protein"]


def test_literature_accession_removes_tc_row_in_build():
    table, cluster_of = _toy()
    table.append(tc("tclit", "P1", "330879", "Eurotiomycetes"))
    cluster_of["tclit"] = "cl_tclit"
    lit = [{"accession": "P1", "seq_sha256": "lit1", "literature_positive": "yes"}]
    cluster_of["lit1"] = "cl_lit1"
    _, removed = _build(table, cluster_of, lit)
    got = {(r["split_id"], r["rule"]) for r in removed if r["seq_sha256"] == "tclit"}
    assert ("S3-Eurotiomycetes", "a_test_protein") in got


def test_user_files_in_tmp_dir_survive_and_fresh_dir_is_removed(tmp_path):
    fx = _built(tmp_path)
    scratch = tmp_path / "scratch"
    (scratch / "mine").mkdir(parents=True)
    (scratch / "mine" / "keep.txt").write_text("x")
    (scratch / "file.txt").write_text("y")
    assert _run09(fx, tmp_path) == 0
    assert (scratch / "mine" / "keep.txt").read_text() == "x"
    assert (scratch / "file.txt").read_text() == "y"
    assert sorted(p.name for p in scratch.iterdir()) == ["file.txt", "mine"]


def test_fresh_dir_is_removed_after_a_failed_command(tmp_path, monkeypatch):
    fx = _built(tmp_path)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    (scratch / "file.txt").write_text("y")
    monkeypatch.setenv("STUB_MMSEQS_FAIL", "easy-search")
    assert _run09(fx, tmp_path) == 2
    assert [p.name for p in scratch.iterdir()] == ["file.txt"]


@pytest.mark.parametrize("how", ["exit 132", "kill -ILL $$"])
def test_avx2_message_for_exit_132_and_sigill(tmp_path, capsys, how):
    fx = _built(tmp_path)
    bad = tmp_path / "mmseqs_bad"
    bad.write_text(f"#!/bin/bash\n{how}\n")
    bad.chmod(0o755)
    assert _run09(fx, tmp_path, bad) == 2
    err = capsys.readouterr().err
    assert "AVX2" in err and "132 or -4" in err
    code = "132" if how == "exit 132" else "-4"
    assert f"exited with {code}" in err
