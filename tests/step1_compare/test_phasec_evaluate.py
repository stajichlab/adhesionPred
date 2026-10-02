"""11_evaluate.py: metrics.json, labels, prevalence, proteome calls, findings, golden file."""

import gzip
import json
import math
import os

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import findings  # noqa: E402
import metrics as mt  # noqa: E402
import phasec_fixture as pf  # noqa: E402
import splits  # noqa: E402
import truth_table  # noqa: E402
from conftest import TESTS_DIR, load_phasec  # noqa: E402

GOLDEN = TESTS_DIR / "fixtures" / "phasec" / "golden_metrics.json.gz"


def _metrics(work):
    return json.loads((work / "phasec" / "metrics.json").read_text())


def _run11(fx, n=50):
    return load_phasec("11_evaluate").main(pf.eval_argv(fx, n))


def test_metrics_json_has_the_test_sets_and_no_nan(phasec_chain):
    text = (phasec_chain["work"] / "phasec" / "metrics.json").read_text()

    def refuse(name):
        raise AssertionError(f"metrics.json holds {name}")

    m = json.loads(text, parse_constant=refuse)
    assert set(m["test_sets"]) == {
        "S1:all", "S1:Calb_CGD", "S1:Scer_SGD", "S2-Calb_CGD:Calb_CGD", "S2-Scer_SGD:Scer_SGD",
        "S2-Spom_PomBase:Spom_PomBase", "S3-Basidiomycota:clade",
        "S3-Basidiomycota:Cneo_H99_GOA", "S3-Basidiomycota:Umay_MYCMD",
        "S3-Eurotiomycetes:clade", "S3-Eurotiomycetes:Afum_ASPFU",
        "S3-Eurotiomycetes:literature",
    }  # fmt: skip
    s1 = m["test_sets"]["S1:all"]["truth"]["direct"]["n"]
    assert "identity_below_0.3" not in s1  # S1 has no maximum-identity stratum
    s2 = m["test_sets"]["S2-Spom_PomBase:Spom_PomBase"]["truth"]["direct"]["n"]
    assert "identity_below_0.3" in s2
    assert all(ts["label"] in ("estimate", "smoke test") for ts in m["test_sets"].values())


def test_undefined_metrics_are_null(phasec_chain):
    # Review Focus 1: a stratum without positives, a test set without negatives
    m = _metrics(phasec_chain["work"])
    nsec = m["test_sets"]["S1:all"]["truth"]["direct"]["metrics"]["N-sec"]["V-go"]["M8"]
    assert nsec["recall"]["value"] is None and nsec["fpr"]["value"] is not None
    lit = m["test_sets"]["S3-Eurotiomycetes:literature"]["truth"]["direct"]["metrics"]["all"][
        "V-go"
    ]
    assert lit["M8"]["roc_auc"]["value"] is None and lit["M8"]["fpr"]["value"] is None
    assert lit["M8"]["recall"]["value"] is not None
    # review I-1: without negatives precision is 1 by construction; spec 3.2 says recall only
    for c in ("M8", "R2", "B1"):
        assert lit[c]["precision"]["value"] is None, c
    assert lit["M8"]["pr_auc"]["value"] is None
    assert lit["M8"]["precision_at_recall_0.8"]["value"] is None
    assert lit["M8"]["precision_at_recall_0.9"]["value"] is None
    lit_truth = m["test_sets"]["S3-Eurotiomycetes:literature"]["truth"]["direct"]
    vs_rule = lit_truth["vs_rule"]["V-go"]["M8"]["precision_at_rule_recall"]
    assert vs_rule["value"] is None
    # review fix 1: strata with one class only (wall and extracellular-only: positives;
    # N-int, N-sec, PM-TM: negatives) have no precision, PR-AUC or precision at recall
    mets = m["test_sets"]["S1:all"]["truth"]["direct"]["metrics"]
    for stratum in ("wall", "extracellular-only"):
        c = mets[stratum]["V-go"]["M8"]
        assert c["precision"]["value"] is None and c["pr_auc"]["value"] is None, stratum
        assert c["precision_at_recall_0.8"]["value"] is None, stratum
        assert c["recall"]["value"] is not None, stratum
    assert mets["wall"]["V-go"]["M8"]["fpr"]["value"] is None
    nsec_m8 = mets["N-sec"]["V-go"]["M8"]
    assert nsec_m8["precision"]["value"] is None and nsec_m8["pr_auc"]["value"] is None
    assert nsec_m8["fpr"]["value"] is not None
    # a test set with negatives keeps its precision
    s1 = m["test_sets"]["S1:all"]["truth"]["direct"]["metrics"]["all"]["V-go"]["M8"]
    assert s1["precision"]["value"] is not None and s1["pr_auc"]["value"] is not None


def test_direct_filter_applies_to_test_rows_only(phasec_chain):
    out = phasec_chain["work"] / "phasec"
    table = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "eval_table.tsv.gz")}
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    test = [
        m
        for m in members
        if m["split_id"] == "S1" and m["part"] == "test" and m["class"] in ("pos", "neg")
    ]
    direct = [m for m in test if table[m["seq_sha256"]]["homology_only"] == "no"]
    assert len(direct) < len(test)
    n = _metrics(phasec_chain["work"])["test_sets"]["S1:all"]["truth"]
    assert n["all"]["n"]["all"]["pos"] + n["all"]["n"]["all"]["neg"] == len(test)
    assert n["direct"]["n"]["all"]["pos"] + n["direct"]["n"]["all"]["neg"] == len(direct)
    train = [
        m for m in members if m["split_id"] == "S1" and m["fold"] == "0" and m["part"] == "train"
    ]
    assert any(table[m["seq_sha256"]]["homology_only"] == "yes" for m in train)


def test_estimate_label_rule():
    ok = {c: 0.10 for c in findings.LABEL_CANDIDATES}
    assert findings.estimate_label(ok, 20) == "estimate"
    assert findings.estimate_label({**ok, "H": 0.1001}, 20) == "smoke test"
    assert findings.estimate_label({**ok, "R2": 0.0999}, 20) == "estimate"
    assert findings.estimate_label({**ok, "M35": None}, 20) == "smoke test"  # no positives
    assert findings.estimate_label({c: 0.05 for c in ("R2", "M8")}, 50) == "smoke test"  # missing
    assert "B1" not in findings.LABEL_CANDIDATES and "R0" not in findings.LABEL_CANDIDATES
    # count floor (owner decision 2026-10-01): a zero-width interval of 7 positives is a smoke test
    zero = {c: 0.0 for c in findings.LABEL_CANDIDATES}
    assert findings.estimate_label(zero, 7) == "smoke test"
    assert findings.estimate_label(zero, 19) == "smoke test"
    small = {c: 0.04 for c in findings.LABEL_CANDIDATES}
    assert findings.estimate_label(small, 20) == "estimate"
    assert findings.floor_met(20) and not findings.floor_met(19)


def test_label_in_metrics_follows_the_half_widths(phasec_chain):
    for name, ts in _metrics(phasec_chain["work"])["test_sets"].items():
        n_pos = ts["truth"]["direct"]["n"]["all"]["pos"]
        assert ts["n_direct_positives"] == n_pos, name
        assert ts["floor_met"] is (n_pos >= findings.MIN_DIRECT_POSITIVES), name
        assert ts["label"] == findings.estimate_label(ts["recall_half_width"], n_pos), name
        direct = ts["truth"]["direct"]["metrics"]["all"]["V-go"]
        r2 = direct["R2"]["recall"]
        if r2["lo"] is not None:
            assert ts["recall_half_width"]["R2"] == pytest.approx((r2["hi"] - r2["lo"]) / 2)


def test_prevalence_table_formula(phasec_chain):
    # hand-computed: recall 0.8, FPR 0.1, prevalence 0.05 -> 0.04 / (0.04 + 0.095)
    assert mt.precision_at_prevalence(0.8, 0.1, 0.05) == pytest.approx(0.2962962962962963)
    block = _metrics(phasec_chain["work"])["test_sets"]["S1:all"]["truth"]["direct"]
    r2 = block["metrics"]["all"]["V-go"]["R2"]
    for key, pi in (("0.01", 0.01), ("0.05", 0.05), ("0.1", 0.10)):
        want = mt.precision_at_prevalence(r2["recall"]["value"], r2["fpr"]["value"], pi)
        assert block["prevalence"]["V-go"]["R2"][key]["value"] == pytest.approx(float(want))


def test_proteome_calls_use_oof_for_training_hashes(phasec_chain):
    out = phasec_chain["work"] / "phasec"
    calls = truth_table.read_tsv(out / "proteome_calls.tsv.gz")
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    s1 = {
        m["seq_sha256"]: m["fold"]
        for m in members
        if m["split_id"] == "S1" and m["part"] in ("test", "test_tc")
    }
    train = {m["seq_sha256"] for m in members if m["part"] in ("train", "train_tc")}
    scores = {(r["split_id"], r["fold"], r["variant"], r["candidate"], r["seq_sha256"]): r
              for r in truth_table.read_tsv(out / "scores.tsv.gz")}  # fmt: skip
    seen_oof = 0
    for r in calls:
        h = r["seq_sha256"]
        if h in train:
            assert r["score_source"] != "final", h  # a training hash never gets `final`
        if h in s1:
            seen_oof += 1
            assert r["score_source"] == "oof"
            want = scores[("S1", s1[h], "V-go", "M8", h)]
            assert r["score_M8_V-go"] == want["score"] and r["call_M8_V-go"] == want["call"]
        else:
            assert r["score_source"] == "final"
            assert r["score_M8_V-go"] == scores[("FULL", "0", "V-go", "M8", h)]["score"]
    assert seen_oof >= 2  # the fixture puts truth and T-c sequences into the proteome sets
    assert {r["set_id"] for r in calls} == {"Scer_proteome", "Cimm_RS_proteome"}


def test_findings_booleans():
    assert findings.finding_a({"S1:all": 0.995, "S2-x:x": 0.95})["holds"] is False  # saturated
    assert findings.finding_a({"S1:all": 0.97, "S2-x:x": 0.95})["holds"] is True
    assert findings.finding_a({"S1:all": None})["holds"] is False
    assert findings.finding_a({})["holds"] is False
    up = {"value": 0.2, "lo": 0.05, "hi": 0.3}
    zero = {"value": 0.1, "lo": -0.01, "hi": 0.3}
    b = findings.finding_b({"M8": {"B1": up, "R2": up}, "M35": {"B1": up, "R2": zero}})
    assert b["holds_for"] == ["M8"] and b["candidates"]["M35"]["beats_R2"] is False
    c = findings.finding_c(
        {"S2-a:a": b, "S2-b:b": findings.finding_b({"M8": {"B1": zero, "R2": up}})}
    )
    assert c["holds_for"] == []
    c = findings.finding_c({"S2-a:a": b, "S2-b:b": b, "S2-c:c": b})
    assert c["holds_for"] == ["M8"] and c["n_s2"] == 3
    # review M-1: fewer than three S2 sets never holds
    c = findings.finding_c({"S2-a:a": b, "S2-b:b": b})
    assert c["holds_for"] == [] and c["n_s2"] == 2 and c["expected_s2"] == 3


def test_findings_json_reads_metrics_json(phasec_chain):
    m = _metrics(phasec_chain["work"])
    f = json.loads((phasec_chain["work"] / "phasec" / "findings.json").read_text())
    b1 = m["test_sets"]["S1:all"]["truth"]["direct"]["metrics"]["all"]["V-go"]["B1"]["roc_auc"]
    assert f["a_b1_not_saturated"]["b1_roc_auc"]["S1:all"] == b1["value"]
    diff = m["test_sets"]["S1:all"]["truth"]["direct"]["nsec_fpr_at_rule_recall"]["V-go"]["diff"]
    assert f["b_ml_beats_b1_and_r2_on_nsec_s1"]["candidates"]["M8"]["B1"] == diff["M8"]["B1"]
    assert set(f["c_same_under_s2"]["test_sets"]) == {
        "S2-Calb_CGD:Calb_CGD",
        "S2-Scer_SGD:Scer_SGD",
        "S2-Spom_PomBase:Spom_PomBase",
    }


def test_named_panel_and_context(phasec_chain):
    m = _metrics(phasec_chain["work"])
    names = [p["name"] for p in m["named_panel"]]
    assert names[:2] == ["FLO1", "SAG1"] and "LIT2" in names  # hard_negative row appended
    lit2 = next(p for p in m["named_panel"] if p["name"] == "LIT2")
    assert lit2["found"] and lit2["literature_class"] == "hard_negative"
    assert not next(p for p in m["named_panel"] if p["name"] == "FLO1")["found"]
    ctx = m["context"]["onygenales_tc_rows_in_vkw_training"]
    assert ctx.get("S1|0", 0) + ctx.get("S1|1", 0) + ctx.get("S1|2", 0) >= 1
    assert "S3-Eurotiomycetes|0" not in ctx  # removed by rule c (ruling C-7)


def test_stale_input_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    path = fx["work"] / "phasec" / "scores.tsv.gz"
    path.write_bytes(path.read_bytes() + b"\n")
    assert _run11(fx) == 2
    assert "scores.tsv.gz differs from the SHA-256 in scores_run.json" in capsys.readouterr().err


def test_tc_row_in_a_test_set_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    next(m for m in members if m["origin"] == "tc" and m["split_id"] == "S2-Calb_CGD")["part"] = (
        "test"
    )
    truth_table.write_tsv(out / "split_members.tsv.gz", splits.MEMBER_COLUMNS, members)
    pf.fix_recorded_hash(
        out / "splits_run.json", "split_members.tsv.gz", out / "split_members.tsv.gz"
    )
    pf.fix_recorded_hash(out / "scores_run.json", "split_members.tsv.gz",
                         out / "split_members.tsv.gz", key="input_sha256")  # fmt: skip
    assert _run11(fx) == 2
    assert "is in a test set" in capsys.readouterr().err


def _close(a, b, path="$"):
    if isinstance(a, dict):
        assert isinstance(b, dict) and set(a) == set(b), path
        for k in a:
            _close(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            _close(x, y, f"{path}[{i}]")
    elif isinstance(a, float) or isinstance(b, float):
        assert a is not None and b is not None and math.isclose(a, b, rel_tol=0, abs_tol=1e-9), path
    else:
        assert a == b, path


def test_golden_metrics(tmp_path):
    fx = pf.run_chain(tmp_path, 10, load_phasec, candidates="B1,R2,M8")
    assert _run11(fx, 50) == 0
    got = _metrics(fx["work"])
    # hand check: R2 recall on S1:all (direct truth, V-go) from the scores and the tables
    out = fx["work"] / "phasec"
    table = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "eval_table.tsv.gz")}
    calls = {r["seq_sha256"]: r["call"] for r in truth_table.read_tsv(out / "scores.tsv.gz")
             if r["split_id"] == "S1" and r["variant"] == "V-go" and r["candidate"] == "R2"}  # fmt: skip
    pos = [h for h, r in table.items() if r["origin"] == "go" and r["roles"] == "train"
           and r["class"] == "pos" and r["homology_only"] == "no"]  # fmt: skip
    r2 = got["test_sets"]["S1:all"]["truth"]["direct"]["metrics"]["all"]["V-go"]["R2"]["recall"]
    assert r2["value"] == pytest.approx(sum(calls[h] == "1" for h in pos) / len(pos), abs=1e-12)
    if os.environ.get("PHASEC_WRITE_GOLDEN") == "1":
        GOLDEN.write_bytes(
            gzip.compress(json.dumps(got, indent=1, sort_keys=True).encode(), mtime=0)
        )
        pytest.skip("golden file written")
    _close(got, json.loads(gzip.decompress(GOLDEN.read_bytes())))
    run = json.loads((out / "evaluate_run.json").read_text())
    assert run["library_versions"]["numpy"] and run["library_versions"]["sklearn"]


# ---- addendum tests (controller addendum to Task 9) ----


def test_scores_with_a_hash_in_test_and_test_lit(phasec_chain):
    # addendum 1: `part` can be `test,test_lit`; the reader indexes by hash only
    ev = load_phasec("11_evaluate")
    out = phasec_chain["work"] / "phasec"
    h = phasec_chain["named"]["LIT_AND_GO"]
    raw = [r for r in truth_table.read_tsv(out / "scores.tsv.gz")
           if r["seq_sha256"] == h and r["candidate"] == "R2" and r["variant"] == "V-go"]  # fmt: skip
    parts = {(r["split_id"], r["part"]) for r in raw}
    assert ("S3-Eurotiomycetes", "test,test_lit") in parts
    scores = ev.read_scores(out / "scores.tsv.gz")
    d = scores[("S3-Eurotiomycetes", "0", "V-go", "R2")]
    assert len(d["index"]) == len(set(d["index"]))
    row = next(r for r in raw if r["split_id"] == "S3-Eurotiomycetes")
    assert d["call"][d["index"][h]] == (row["call"] == "1")
    # the hash is a row of the GO test set and of the literature test set; both read its score
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    table = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "eval_table.tsv.gz")}
    lit = truth_table.read_tsv(out / "eval_literature.tsv")
    ident = truth_table.read_tsv(out / "max_identity.tsv.gz")
    species = truth_table.read_tsv(phasec_chain["species"])
    sets = {t.name: t for t in ev.test_sets(members, table, lit, ident, species)}
    for name in ("S3-Eurotiomycetes:Afum_ASPFU", "S3-Eurotiomycetes:literature"):
        rows = [r for r in sets[name].rows if r["seq_sha256"] == h]
        assert len(rows) == 1, name
        ev.gather(scores, rows, "V-go", "R2")  # no StopError


def test_link_08_to_09_stops_on_a_stale_split(phasec_chain, tmp_path, capsys):
    # addendum 2: eval_table changes and build_run.json is re-pinned; 09 stays stale
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    rows = truth_table.read_tsv(out / "eval_table.tsv.gz")
    cols = list(rows[0])
    rows[0]["gene_ids"] += "x"
    truth_table.write_tsv(out / "eval_table.tsv.gz", cols, rows)
    pf.fix_recorded_hash(out / "build_run.json", "eval_table.tsv.gz", out / "eval_table.tsv.gz")
    before = {p.name for p in out.iterdir()}
    (out / "metrics.json").unlink()
    assert _run11(fx) == 2
    err = capsys.readouterr().err
    assert "STOP" in err and "eval_table.tsv.gz" in err and "re-run 09" in err
    assert {p.name for p in out.iterdir()} == before - {"metrics.json"}  # no output written


def test_link_09_to_10_stops_on_other_clusters(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    log = json.loads((out / "scores_run.json").read_text())
    log["input_sha256"]["clusters.tsv.gz"] = "0" * 64
    (out / "scores_run.json").write_text(json.dumps(log))
    (out / "metrics.json").unlink()
    assert _run11(fx) == 2
    assert "clusters.tsv.gz" in capsys.readouterr().err
    assert not (out / "metrics.json").exists()


def test_point_value_is_the_unweighted_metric(phasec_chain):
    # addendum 3: row 0 of the evaluated weight matrix is all ones. If the point value came
    # from the first random resample, these would differ.
    out = phasec_chain["work"] / "phasec"
    table = {r["seq_sha256"]: r for r in truth_table.read_tsv(out / "eval_table.tsv.gz")}
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    test = [
        m["seq_sha256"]
        for m in members
        if m["split_id"] == "S1" and m["part"] == "test" and m["class"] in ("pos", "neg")
        and table[m["seq_sha256"]]["homology_only"] == "no"
    ]  # fmt: skip
    cls = {m["seq_sha256"]: m["class"] for m in members
           if m["split_id"] == "S1" and m["part"] == "test"}  # fmt: skip
    scored = {(r["fold"], r["seq_sha256"]): r for r in truth_table.read_tsv(out / "scores.tsv.gz")
              if r["split_id"] == "S1" and r["variant"] == "V-go" and r["candidate"] == "M8"}  # fmt: skip
    fold = {m["seq_sha256"]: m["fold"] for m in members
            if m["split_id"] == "S1" and m["part"] == "test"}  # fmt: skip
    pos = [h for h in test if cls[h] == "pos"]
    neg = [h for h in test if cls[h] == "neg"]
    call = {h: scored[(fold[h], h)]["call"] == "1" for h in test}
    blk = _metrics(phasec_chain["work"])["test_sets"]["S1:all"]["truth"]["direct"]
    m8 = blk["metrics"]["all"]["V-go"]["M8"]
    assert m8["recall"]["value"] == pytest.approx(sum(call[h] for h in pos) / len(pos), abs=1e-12)
    assert m8["fpr"]["value"] == pytest.approx(sum(call[h] for h in neg) / len(neg), abs=1e-12)


def test_proteome_hash_that_is_tc_truth_and_proteome_at_once():
    # addendum 6: one hash = a T-c row, an unlabelled truth gene and two proteome genes
    ev = load_phasec("11_evaluate")

    def d(hs, score):
        return {"index": {h: i for i, h in enumerate(hs)}, "score": np.array([score] * len(hs)),
                "call": np.array([score > 1] * len(hs))}  # fmt: skip

    h1, h2, h3 = "h1", "h2", "h3"
    scores = {}
    for v in ("V-go", "V-kw"):
        scores[("S1", "1", v, "M8")] = d([h1], 0.5)  # h1 is in an S1 test part
        scores[("FULL", "0", v, "M8")] = d([h1, h2, h3], 9.0)
    members = [
        {"set_id": "truth", "source_id": "Scer_SGD", "gene_id": "U1", "seq_sha256": h1},
        {"set_id": "uniprot_kw", "source_id": "uniprot_kw", "gene_id": "Q1", "seq_sha256": h1},
        {"set_id": "Scer_proteome", "source_id": "Scer_proteome", "gene_id": "A", "seq_sha256": h1},
        {"set_id": "Cimm_RS_proteome", "source_id": "Cimm_RS_proteome", "gene_id": "B",
         "seq_sha256": h1},
        {"set_id": "Scer_proteome", "source_id": "Scer_proteome", "gene_id": "C", "seq_sha256": h2},
        {"set_id": "Scer_proteome", "source_id": "Scer_proteome", "gene_id": "D", "seq_sha256": h3},
    ]  # fmt: skip
    table = {h1: {"class": "pos", "label": "P-ext", "origin": "tc"}}
    rows = ev.proteome_rows(members, {"Scer_proteome", "Cimm_RS_proteome"}, table, scores,
                            {h1: "1"}, {h1, h2}, ("M8",))  # fmt: skip
    assert [(r["set_id"], r["gene_id"]) for r in rows] == [
        ("Scer_proteome", "A"), ("Cimm_RS_proteome", "B"), ("Scer_proteome", "C"),
        ("Scer_proteome", "D"),
    ]  # fmt: skip
    src = {r["gene_id"]: r["score_source"] for r in rows}
    assert src == {"A": "oof", "B": "oof", "C": "in_sample", "D": "final"}
    by = {r["gene_id"]: r for r in rows}
    assert by["A"]["score_M8_V-go"] == "0.5" and by["B"]["score_M8_V-go"] == "0.5"
    assert by["C"]["score_M8_V-go"] == "9" and by["D"]["call_M8_V-kw"] == "1"
    assert by["A"]["origin"] == "tc" and by["C"]["origin"] == ""


# ---- fix round 1 tests ----


def _hand_rows():
    """42 rows: 12 positives (6 wall, 6 extracellular-only), 20 N-sec and 10 N-int negatives;
    clusters of 2 rows."""
    rows = []

    def add(i, cls, subset, stratum, cluster):
        rows.append({"label": "", "subset": subset, "stratum": stratum, "d8_class": "",
                     "homology_only": "no", "internal_evidence_htp_only": "no",
                     "source_ids": "x", "gene_ids": "g", "length": 100, "seq_sha256": f"h{i}",
                     "split": "S1", "fold": "0", "class": cls, "cluster": cluster,
                     "part": "test", "below": ""})  # fmt: skip

    i = 0
    for k in range(12):
        add(i, "pos", "wall" if k < 6 else "extracellular-only", "", f"p{k // 2}")
        i += 1
    for k in range(20):
        add(i, "neg", "", "N-sec", f"n{k // 2}")
        i += 1
    for k in range(10):
        add(i, "neg", "", "N-int", f"m{k // 2}")
        i += 1
    return rows


def _hand_scores(rows, spec):
    """spec[cand] = (score_fn(row_index, row) -> float, call_fn(...) -> bool); both variants
    get the same values."""
    out = {}
    for v in ("V-go", "V-kw"):
        for cand, (sf, cf) in spec.items():
            s = np.array([sf(i, r) for i, r in enumerate(rows)])
            out[("S1", "0", v, cand)] = {
                "index": {r["seq_sha256"]: i for i, r in enumerate(rows)},
                "score": s, "prob": s.copy(),
                "call": np.array([cf(i, r) for i, r in enumerate(rows)], dtype=bool),
            }  # fmt: skip
    return out


def _ml_beats_case():
    rows = _hand_rows()
    pos_i = lambda r: int(r["seq_sha256"][1:])  # noqa: E731
    spec = {
        # R2 finds 10 of 12 positives and calls every negative
        "R2": (lambda i, r: 1.0 if r["class"] == "neg" or pos_i(r) < 10 else 0.0,
               lambda i, r: r["class"] == "neg" or pos_i(r) < 10),
        # M8 separates the classes; it calls every positive
        "M8": (lambda i, r: 1.0 if r["class"] == "pos" else 0.0, lambda i, r: r["class"] == "pos"),
        "B1": (lambda i, r: 0.5, lambda i, r: False),
    }  # fmt: skip
    return rows, _hand_scores(rows, spec)


def _ml_loses_case():
    rows = _hand_rows()
    spec = {
        "R2": (lambda i, r: 1.0 if r["class"] == "pos" else 0.0, lambda i, r: r["class"] == "pos"),
        "M8": (lambda i, r: 0.5, lambda i, r: False),
        "B1": (lambda i, r: 1.0 if r["class"] == "pos" else 0.0, lambda i, r: r["class"] == "pos"),
    }  # fmt: skip
    return rows, _hand_scores(rows, spec)


def _hand_block(case, truth="direct", n=200):
    ev = load_phasec("11_evaluate")
    rows, scores = case
    ts = ev.TestSet("S1:all", "S1", "pooled", rows)
    return ev.evaluate_truth(ts, truth, scores, ("B1", "R2", "M8"), n, 1)


def test_hand_built_by_construction_values_are_null():
    block, *_ = _hand_block(_ml_beats_case())
    m = block["metrics"]
    for stratum in ("wall", "extracellular-only"):  # positives only
        c = m[stratum]["V-go"]["M8"]
        assert c["precision"]["value"] is None and c["pr_auc"]["value"] is None
        assert c["precision_at_recall_0.8"]["value"] is None
        assert c["precision_at_recall_0.9"]["value"] is None
        assert c["recall"]["value"] == 1.0 and c["fpr"]["value"] is None
    for stratum in ("N-sec", "N-int"):  # negatives only
        c = m[stratum]["V-go"]["R2"]
        assert c["precision"]["value"] is None and c["pr_auc"]["value"] is None
        assert c["recall"]["value"] is None and c["fpr"]["value"] == 1.0
        assert c["roc_auc"]["value"] is None
    assert m["all"]["V-go"]["M8"]["precision"]["value"] == 1.0  # defined with both classes


def test_variant_effect_of_identical_scores_is_exactly_zero():
    # kills a bootstrap weight matrix that is redrawn per candidate or variant
    block, *_ = _hand_block(_ml_beats_case())
    for stratum in ("all", "wall", "N-sec"):
        for cand in ("R2", "M8", "B1"):
            for name, e in block["variant_effect"][stratum][cand].items():
                if e["value"] is None:
                    continue
                assert e["value"] == 0.0 and e["lo"] == 0.0 and e["hi"] == 0.0, (
                    stratum,
                    cand,
                    name,
                )
    e = block["variant_effect"]["all"]["M8"]["recall"]
    assert e["value"] == 0.0 and e["lo"] == 0.0 and e["hi"] == 0.0


def test_finding_b_hand_built_beats_and_mirror():
    ev = load_phasec("11_evaluate")
    block, *_ = _hand_block(_ml_beats_case())
    diff = block["nsec_fpr_at_rule_recall"]["V-go"]["diff"]
    # comparator N-sec FPR minus ML: R2 calls all N-sec (1.0), B1 puts all at one score (1.0), M8 0.0
    assert diff["M8"]["R2"]["value"] == 1.0 and diff["M8"]["B1"]["value"] == 1.0
    assert diff["M8"]["R2"]["lo"] > 0 and diff["M8"]["B1"]["lo"] > 0
    f = findings.finding_b(diff)
    assert f["holds_for"] == ["M8"] and f["candidates"]["M8"]["beats_B1"] is True
    assert f["candidates"]["M8"]["beats_R2"] is True
    mirror, *_ = _hand_block(_ml_loses_case())
    diff = mirror["nsec_fpr_at_rule_recall"]["V-go"]["diff"]
    assert diff["M8"]["R2"]["value"] == -1.0 and diff["M8"]["B1"]["value"] == -1.0
    f = findings.finding_b(diff)
    assert f["holds_for"] == [] and f["candidates"]["M8"]["beats_R2"] is False
    assert f["candidates"]["M8"]["beats_B1"] is False
    assert ev.findings is findings


def test_comparison_level_is_the_rule_recall_per_resample(monkeypatch):
    # kills "FPR at the ML threshold": the level passed to fpr_at_recall is R2's recall in each
    # resample, not the ML candidate's (M8 recall is 1.0 here, R2's is about 0.83)
    import bootstrap

    ev = load_phasec("11_evaluate")
    rows, scores = _ml_beats_case()
    seen = []
    real = ev.metrics.fpr_at_recall

    def spy(W, y, score, level, mask):
        seen.append(np.array(level, dtype=float))
        return real(W, y, score, level, mask)

    monkeypatch.setattr(ev.metrics, "fpr_at_recall", spy)
    ts = ev.TestSet("S1:all", "S1", "pooled", rows)
    ev.evaluate_truth(ts, "direct", scores, ("B1", "R2", "M8"), 200, 1)
    y = np.array([r["class"] == "pos" for r in rows])
    W = bootstrap.cluster_weights([r["cluster"] for r in rows], 200,
                                  bootstrap.seed_for(1, "S1:all|direct"))  # fmt: skip
    Wx = np.vstack([np.ones((1, len(rows))), W])
    r2 = scores[("S1", "0", "V-go", "R2")]["call"]
    want = np.asarray(mt.recall(Wx, y, r2), dtype=float)
    assert len(seen) == 4  # B1 and M8, two variants
    for level in seen:
        assert np.allclose(level, want, equal_nan=True)
    assert not np.allclose(want, 1.0)


def test_read_scores_stops_on_a_duplicate_hash(tmp_path):
    ev = load_phasec("11_evaluate")
    path = tmp_path / "scores.tsv.gz"
    row = ["S1", "0", "V-go", "R2", "hh", "test", "0.5", "0.5", "1"]
    truth_table.write_tsv(
        path, ev.SCORE_COLUMNS, [dict(zip(ev.SCORE_COLUMNS, row, strict=True))] * 2
    )
    with pytest.raises(ev.evalio.StopError, match="twice"):
        ev.read_scores(path)


def test_evaluate_run_records_every_input_hash(phasec_chain):
    import manifest

    out = phasec_chain["work"] / "phasec"
    got = json.loads((out / "evaluate_run.json").read_text())["input_sha256"]
    for name in ("eval_literature.tsv", "clusters.tsv.gz", "max_identity.tsv.gz"):
        assert got[name] == manifest.sha256_file(out / name), name
    assert got["species.tsv"] == manifest.sha256_file(phasec_chain["species"])


def test_agreement_records_score_source_counts(phasec_chain):
    m = _metrics(phasec_chain["work"])
    counts = m["agreement"]["score_source_counts"]
    assert counts["proteomes"] == m["score_sources"]
    out = phasec_chain["work"] / "phasec"
    table = truth_table.read_tsv(out / "eval_table.tsv.gz")
    for cls in ("pos", "neg"):
        n = sum(1 for r in table if r["origin"] == "go" and r["class"] == cls)
        assert sum(counts["truth"][cls].values()) == n
        assert counts["truth"][cls].get("oof", 0) > 0 and "in_sample" not in counts["truth"][cls]
