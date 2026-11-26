"""Fixture tests for analysis/crossclade_domain_search (offline: no network, no HPCC)."""

import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DIR = ROOT / "analysis/crossclade_domain_search"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


bi = _load("build_inputs")
sm = _load("summarize")


def _write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


# --- queries.tsv -----------------------------------------------------------------------------------------------------------


def test_queries_table_panels_and_pfams():
    rows = bi.read_tsv(DIR / "queries.tsv")
    assert {r["panel"] for r in rows} == {"asco", "basidio"}
    assert len({r["accession"] for r in rows}) == len(rows) == 14
    assert {r["pfam"] for r in rows if r["panel"] == "basidio"} == {"PF21671"}
    # the asco panel holds no basidiomycete protein and the basidio panel only C. neoformans H99 (held out from the NP11 target)
    assert all(r["species"] != "Rtor_NP11" for r in rows)
    assert {r["species"] for r in rows if r["panel"] == "basidio"} == {"Cneo_H99"}


# --- build_inputs pure functions -------------------------------------------------------------------------------------------


def test_pick_location_best_score_then_n_terminal():
    locs = [
        {"fragments": [{"start": 300, "end": 400}], "score": 1e-5},
        {"fragments": [{"start": 50, "end": 120}, {"start": 130, "end": 200}], "score": 1e-20},
        {"fragments": [{"start": 10, "end": 40}], "score": 1e-20},
    ]
    assert bi.pick_location(locs) == (10, 40, 3)
    with pytest.raises(bi.BuildError):
        bi.pick_location([])


def test_crop_pads_and_clips():
    seq = "ABCDEFGHIJ"
    assert bi.crop(seq, 3, 5, 1) == ("BCDEF", 2, 6)
    assert bi.crop(seq, 1, 10, 5) == (seq, 1, 10)
    with pytest.raises(bi.BuildError):
        bi.crop(seq, 5, 11, 0)


def test_query_id_round_trips_to_accession():
    qid = bi.query_id("A0A2H0ZI42", "PF11765", 1, 350)
    assert qid == "A0A2H0ZI42_PF11765_1_350"
    assert sm.query_accession(qid) == "A0A2H0ZI42"
    assert sm.query_accession("dec07_" + qid) == "A0A2H0ZI42"
    assert sm.decoy_query("dec00_" + qid) == qid
    assert not sm.is_decoy(qid) and sm.is_decoy("dec42_" + qid)


def test_decoys_reverse_then_deterministic_window_shuffles():
    sub = "STSTSTAAAAGGGGPPPPNNNQ"
    dec = bi.decoys_for("Q1_PF1_1_22", sub, 3, 5)
    assert [d for d, _ in dec] == [
        "dec00_Q1_PF1_1_22",
        "dec01_Q1_PF1_1_22",
        "dec02_Q1_PF1_1_22",
        "dec03_Q1_PF1_1_22",
    ]
    assert dec[0][1] == sub[::-1]
    for _, seq in dec[1:]:
        assert len(seq) == len(sub)
        # each 5-residue window keeps its composition
        assert all(sorted(seq[i : i + 5]) == sorted(sub[i : i + 5]) for i in range(0, len(sub), 5))
    assert dec == bi.decoys_for("Q1_PF1_1_22", sub, 3, 5)  # deterministic
    assert len({seq for _, seq in dec[1:]}) > 1


def test_pick_location_missing_score_ranks_last():
    locs = [
        {"fragments": [{"start": 5, "end": 9}]},
        {"fragments": [{"start": 50, "end": 90}], "score": 1e-3},
    ]
    assert bi.pick_location(locs) == (50, 90, 2)


def test_parse_fasta_reduces_uniprot_headers():
    text = (
        ">sp|P08640|FLO11_YEAST Flocculin\nMKV\nLL\n>tr|M7WRY1|M7WRY1_RHOT1 x\nAA\n>plain desc\nC\n"
    )
    assert bi.parse_fasta_text(text) == [("P08640", "MKVLL"), ("M7WRY1", "AA"), ("plain", "C")]


def _staged(tmp_path):
    st = tmp_path / "staged"
    (st / "chunks").mkdir(parents=True)
    _write_tsv(
        st / "proteins.tsv",
        ["accession", "species", "chunk", "length"],
        [
            ["X1", "Rtor_NP11", "Rtor_NP11_001", "3"],
            ["X2", "Cneo_H99", "Cneo_H99_001", "2"],
            ["X3", "Afum", "Afum_001", "1"],
        ],
    )
    (st / "chunks/Rtor_NP11_001.faa").write_text(">X1\nMKL\n")
    (st / "chunks/Cneo_H99_001.faa").write_text(">X2\nMA\n")
    (st / "chunks/Afum_001.faa").write_text(">X3\nM\n")
    return st


def test_staged_targets_filters_species_and_checks_lengths(tmp_path):
    st = _staged(tmp_path)
    assert bi.staged_targets(st, {"Rtor_NP11", "Cneo_H99"}) == {
        "X1": ("Rtor_NP11", "MKL"),
        "X2": ("Cneo_H99", "MA"),
    }
    with pytest.raises(bi.BuildError, match="not in"):
        bi.staged_targets(st, {"Rtor_NP11", "Missing_sp"})
    (st / "chunks/Cneo_H99_001.faa").write_text(">X2\nMAA\n")
    with pytest.raises(bi.BuildError, match="length"):
        bi.staged_targets(st, {"Cneo_H99"})


def test_main_refuses_output_inside_repo_or_existing(tmp_path, capsys):
    st = _staged(tmp_path)
    args = [
        "--queries",
        str(DIR / "queries.tsv"),
        "--staged",
        str(st),
        "--species",
        "Rtor_NP11=1130832",
    ]
    assert bi.main(args + ["--out", str(ROOT / "should_not_exist_xclade")]) == 1
    assert "inside the repo" in capsys.readouterr().err
    assert not (ROOT / "should_not_exist_xclade").exists()
    assert bi.main(args + ["--out", str(tmp_path)]) == 1
    assert "exists" in capsys.readouterr().err


# --- summarize -------------------------------------------------------------------------------------------------------------

QA, QB, QC = "Q59L12_PF11766_43_308", "J9VQE5_PF21671_20_199", "H2_PF21671_1_100"


def _hit(q, t, e, qs=1, qe=50, ql=100):
    return [q, t, "0.3", "50", str(qs), str(qe), "1", "50", str(e), "40", str(ql), "300"]


def _tables(inp):
    return [
        sm.read_tsv(inp / f)
        for f in ("query_domains.tsv", "target_proteins.tsv", "truth_pf21671.tsv")
    ]


def _decoy_ids(inp):
    return {
        ln[1:].split()[0]
        for ln in (inp / "decoys.faa").read_text().splitlines()
        if ln.startswith(">")
    }


@pytest.fixture
def fixture_dirs(tmp_path):
    inp, res = tmp_path / "in", tmp_path / "res"
    inp.mkdir()
    res.mkdir()
    _write_tsv(
        inp / "query_domains.tsv",
        ["query_id", "panel", "accession", "gene"],
        [
            [QA, "asco", "Q59L12", "ALS3"],
            [QB, "basidio", "J9VQE5", "CPL1"],
            [QC, "basidio", "H2", "PARALOG"],
        ],
    )
    _write_tsv(
        inp / "target_proteins.tsv",
        ["accession", "species", "length"],
        [[a, "Rtor_NP11", "300"] for a in ("M1", "M2", "M3")]
        + [[a, "Cneo_H99", "300"] for a in ("J9VQE5", "H2", "H3")],
    )
    _write_tsv(
        inp / "truth_pf21671.tsv",
        ["accession", "species", "pfam", "in_targets"],
        [
            ["M1", "Rtor_NP11", "PF21671", "yes"],
            ["M2", "Rtor_NP11", "PF21671", "yes"],
            ["J9VQE5", "Cneo_H99", "PF21671", "yes"],
            ["H2", "Cneo_H99", "PF21671", "yes"],
            ["Z9", "Rtor_NP11", "PF21671", "not_in_staged_proteome"],
        ],
    )
    (inp / "decoys.faa").write_text(
        "".join(f">dec0{k}_{q} decoy\nAAA\n" for q in (QA, QB, QC) for k in (0, 1))
    )
    rows = {
        "foldseek_queries": [
            _hit(QB, "M1", 1e-10),
            _hit(QB, "M3", 5e-3),
            _hit(QB, "J9VQE5", 1e-50),  # self: dropped
            _hit(QB, "H2", 1e-20),  # paralog that is itself a query: dropped
            _hit(QA, "H3", 0.05),
        ],
        "foldseek_decoys": [
            _hit("dec00_" + QA, "M2", 1e-4),
            _hit("dec01_" + QB, "J9VQE5", 1e-4),  # decoy hitting its own forward protein: kept
            _hit("dec00_" + QB, "M1", 0.05),
        ],
        "mmseqs_queries": [_hit(QB, "M1", 1e-4)],
        "mmseqs_decoys": [],
    }
    for name, rr in rows.items():
        with open(res / f"{name}.tsv", "w", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rr)
    flags = {
        a: {
            "accession": a,
            "sp_yes": "0" if a == "H3" else "1",
            "gpi_yes": "0",
            "disorder_longest_run": "0",
            "broad_flag": "0",
            "narrow_flag": "0",
        }
        for a in ("M1", "M2", "M3", "J9VQE5", "H2", "H3")
    }
    return inp, res, flags


def test_summarize_cutoffs_decoys_and_query_filter(fixture_dirs):
    inp, res, flags = fixture_dirs
    hits, decoys, summary, per_query = sm.summarize(*_tables(inp), flags, res, _decoy_ids(inp))
    s = {(r["method"], r["panel"], r["evalue_max"]): r for r in summary}
    b3 = s[("foldseek", "basidio", 1e-3)]
    assert (
        b3["targets_hit"] == 1
        and b3["pf21671_recovered_Rtor_NP11"] == "1/2"
        and b3["non_pf21671_hits_Rtor_NP11"] == 0
    )
    assert b3["pf21671_recovered_Cneo_H99"] == "0/0"  # every H99 member is a query: not scored
    assert b3["non_pf21671_hits_Cneo_H99"] == 0  # self and paralog-query hits dropped
    assert b3["decoy_hits"] == 1 and b3["decoys_with_hit"] == "1/4" and b3["decoy_targets_hit"] == 1
    b2 = s[("foldseek", "basidio", 1e-2)]
    assert (
        b2["targets_hit"] == 2
        and b2["non_pf21671_hits_Rtor_NP11"] == 1
        and b2["targets_hit_with_sp"] == 2
    )
    b1 = s[("foldseek", "basidio", 0.1)]
    assert b1["decoy_hits"] == 2 and b1["decoys_with_hit"] == "2/4" and b1["decoy_targets_hit"] == 2
    a3 = s[("foldseek", "asco", 1e-3)]
    assert a3["targets_hit"] == 0 and a3["decoy_hits"] == 1 and a3["decoys_with_hit"] == "1/2"
    assert (
        s[("foldseek", "asco", 0.1)]["targets_hit"] == 1
        and s[("foldseek", "asco", 0.1)]["targets_hit_with_sp"] == 0
    )
    assert s[("mmseqs", "basidio", 1e-3)]["pf21671_recovered_Rtor_NP11"] == "1/2"
    assert s[("mmseqs", "basidio", 1e-3)]["decoys_with_hit"] == "0/4"
    assert {h[4] for h in hits if h[0] == "foldseek"} == {"M1", "M3", "H3"}
    assert all(h[-1] == ("yes" if h[4] == "M1" else "no") for h in hits)
    assert len([d for d in decoys if d[0] == "foldseek"]) == 3
    pq = {(r[0], r[1], r[4]): r[5:] for r in per_query}
    assert pq[("foldseek", QB, 1e-2)] == [
        0,
        2,
        2,
        "0.0001",
    ]  # Cneo_H99, Rtor_NP11, n_decoys, best decoy E
    assert pq[("foldseek", QC, 1e-2)] == [0, 0, 2, "NA"]


def test_summarize_fails_on_unknown_target_or_decoy(fixture_dirs):
    inp, res, flags = fixture_dirs
    with open(res / "mmseqs_decoys.tsv", "w") as fh:
        fh.write("\t".join(_hit("dec00_" + QA, "NOPE", 1)) + "\n")
    with pytest.raises(sm.SummaryError, match="NOPE"):
        sm.summarize(*_tables(inp), flags, res, _decoy_ids(inp))
    with open(res / "mmseqs_decoys.tsv", "w") as fh:
        fh.write("\t".join(_hit("dec09_" + QA, "M1", 1)) + "\n")
    with pytest.raises(sm.SummaryError, match="not in decoys.faa"):
        sm.summarize(*_tables(inp), flags, res, _decoy_ids(inp))


def test_summarize_main_writes_tables(fixture_dirs, tmp_path):
    inp, res, flags = fixture_dirs
    fl = tmp_path / "flags.tsv"
    cols = ["accession", *sm.FLAG_COLS]
    _write_tsv(fl, cols, [[f[c] for c in cols] for f in flags.values()])
    out = tmp_path / "summary_out"
    args = ["--inputs", str(inp), "--results", str(res), "--flags", str(fl), "--out", str(out)]
    assert sm.main(args) == 0
    assert {p.name for p in out.iterdir()} == {
        "hits.tsv",
        "decoy_hits.tsv",
        "summary.tsv",
        "per_query.tsv",
    }
    header = (out / "hits.tsv").read_text().splitlines()[0].split("\t")
    assert header[header.index("alnlen") + 1 :][:4] == ["tstart", "tend", "qcov", "tcov"]
    assert len((out / "per_query.tsv").read_text().splitlines()[1].split("\t")) == 9
    assert sm.main(args) == 1
