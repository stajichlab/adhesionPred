"""Module hydrophobin_relaxed: the seven hydrophobin-class Pfam models with a relaxed full-sequence cutoff."""

import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules import hydrophobin_relaxed as hr
from cellsurface_sorting_hat.modules import pfam
from cellsurface_sorting_hat.modules.cli import main

HMM = """HMMER3/f [3.4 | Aug 2023]
NAME  Hydrophobin
ACC   PF01185.24
LENG  60
GA    2.60 -1000.00;
HMM          A
//
HMMER3/f [3.4 | Aug 2023]
NAME  Eas
ACC   PF22354.3
LENG  75
GA    2.60 -1000.00;
HMM          A
//
"""

EIGHT = (
    "MKLLVAAG"
    + "C"
    + "A" * 5
    + "CC"
    + "A" * 10
    + "C"
    + "A" * 8
    + "C"
    + "A" * 5
    + "CC"
    + "A" * 9
    + "CGGSTN"
)


def opts(nobias=True, cut_ga=True, query="hydrophobin_relaxed.hmm"):
    flags = (
        "hmmsearch " + ("--nobias " if nobias else "") + ("--cut_ga " if cut_ga else "") + "--noali"
    )
    return f"# Option settings:     {flags} {query} in.fasta\n# Query file:   {query}\n"


def row(target, query="Hydrophobin", qlen=60, seq=5.0, dom=-3.2, acc="PF01185.24"):
    return f"{target} - 100 {query} {acc} {qlen} 1e-3 {seq} 0.1 1 1 1e-3 1e-3 {dom} 0.0 1 60 1 60 1 60 0.9 -\n"


def domtbl(path, rows, **kw):
    path.write_text(opts(**kw) + "".join(rows) + "# [ok]\n")
    return path


def prot(pid, seq=EIGHT, state="ok"):
    return Protein(pid, seq, "sha", state, "", 0.0)


@pytest.fixture
def hmm(tmp_path):
    p = tmp_path / "hydrophobin_relaxed.hmm"
    p.write_text(HMM)
    return p


def test_parse_domtblout_reads_the_full_sequence_score_and_a_negative_domain_score(tmp_path):
    p = domtbl(tmp_path / "d.tbl", [row("p1", seq=5.0, dom=-3.2)])
    h = pfam.parse_domtblout(p)
    assert h[0]["seq_score"] == 5.0 and h[0]["score"] == -3.2


def test_model_info_reads_names_lengths_and_ga(hmm):
    info = hr.model_info(hmm)
    assert info["names"] == {"Hydrophobin": 60, "Eas": 75}
    assert info["ga"] == (2.6, -1000.0)


def test_model_info_refuses_models_with_different_ga(tmp_path):
    p = tmp_path / "m.hmm"
    p.write_text(HMM.replace("GA    2.60 -1000.00;", "GA    3.00 -1000.00;", 1))
    with pytest.raises(ValueError, match="GA"):
        hr.model_info(p)


def test_hit_needs_a_table_row_eight_cysteines_and_r0_called(hmm):
    info = hr.model_info(hmm)
    proteins = [
        prot("a"),
        prot("b"),
        prot("c"),
        prot("d", seq="MKC" * 2 + "A" * 30),
        prot("e"),
        prot("f", seq="", state="na_invalid"),
    ]
    hits = [
        {"target": t, "query": "Hydrophobin", "seq_score": 5.0, "score": -1.0}
        for t in ("a", "b", "c", "d")
    ]
    sp = {"a": "called", "b": "not_called", "d": "called", "e": "called"}  # c has no usable R0 row
    rows = {r["id"]: r for r in hr.relaxed_rows(proteins, hits, sp, info)}
    assert (rows["a"]["state"], rows["a"]["hit"]) == ("ok", "1")
    assert rows["b"]["hit"] == "0"
    assert (rows["c"]["state"], rows["c"]["hit"]) == (
        "ok",
        "",
    )  # would be a hit, R0 missing: not assessable
    assert rows["d"]["hit"] == "0"  # fewer than 8 cysteines
    assert rows["e"]["hit"] == "0"  # no table row
    assert rows["f"]["state"] == "na_invalid"
    assert rows["a"]["score"] == "5.0" and rows["a"]["n_cys"] == "8"


def test_check_table_refuses_a_mismatch(tmp_path, hmm):
    info = hr.model_info(hmm)
    good = domtbl(tmp_path / "g.tbl", [row("a")])
    hr.check_table(good, info, "nobias")  # no error
    with pytest.raises(ValueError, match="nobias"):
        hr.check_table(domtbl(tmp_path / "x.tbl", [row("a")], nobias=False), info, "nobias")
    with pytest.raises(ValueError, match="nobias"):
        hr.check_table(good, info, "default")
    with pytest.raises(ValueError, match="query file"):
        hr.check_table(domtbl(tmp_path / "y.tbl", [row("a")], query="other.hmm"), info, "nobias")
    with pytest.raises(ValueError, match="model"):
        hr.check_table(
            domtbl(tmp_path / "z.tbl", [row("a", query="Cerato-platanin")]), info, "nobias"
        )
    with pytest.raises(ValueError, match="length"):
        hr.check_table(domtbl(tmp_path / "w.tbl", [row("a", qlen=61)]), info, "nobias")
    with pytest.raises(ValueError, match="cut_ga"):
        hr.check_table(domtbl(tmp_path / "v.tbl", [row("a")], cut_ga=False), info, "nobias")


def test_params_hold_the_decisions_and_the_hmm_hash_enters_the_identity(tmp_path, hmm):
    info = hr.model_info(hmm)
    p = hr.module_params(
        info, "nobias", {"module": "step1_rule@R0", "params_hash": "1", "artefact_hash": "2"}
    )
    assert (
        p["cutoff_seq_bits"] == 2.6 and p["domain_ga"] == -1000.0 and p["search_option"] == "nobias"
    )
    assert p["min_cysteines"] == 8 and p["conditions"]["sp_module"]["params_hash"] == "1"
    assert p["models"] == ["Eas", "Hydrophobin"]


def test_cli_writes_the_module_and_stale_identity_follows_the_hmm_file(tmp_path, hmm):
    fasta = tmp_path / "p.faa"
    fasta.write_text(f">A one\n{EIGHT}\n>B two\n{EIGHT}\n>C three\nMKTAYIAKQRQ\n")
    res = tmp_path / "prediction_results.txt"
    res.write_text(
        "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 1\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
        "A x\tSP\t0.1\t0.9\t\nB y\tOTHER\t0.9\t0.1\t\nC z\tOTHER\t0.9\t0.1\t\n"
    )
    wd = tmp_path / "wd"
    common = ["--fasta", str(fasta), "--workdir", str(wd)]
    assert main(["signalp", *common, "--results", str(res), "--signalp-version", "6"]) == 0
    dom = domtbl(tmp_path / "d.tbl", [row("A"), row("B"), row("C", seq=4.0)])
    args = [
        "hydrophobin_relaxed",
        *common,
        "--domtbl",
        str(dom),
        "--hmm",
        str(hmm),
        "--search-option",
        "nobias",
        "--hmmer-version",
        "3.4",
    ]
    assert main(args) == 0
    with gzip.open(wd / "modules/hydrophobin_relaxed.tsv.gz", "rt") as fh:
        rows = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert (rows["A"]["hit"], rows["B"]["hit"], rows["C"]["hit"]) == ("1", "0", "0")
    rec1 = json.loads((wd / "modules/hydrophobin_relaxed.json").read_text())
    assert rec1["run_state"] == "ok"
    hmm.write_text(HMM.replace("2.60", "3.10"))
    dom2 = domtbl(tmp_path / "d2.tbl", [row("A")])
    assert (
        main(
            [
                "hydrophobin_relaxed",
                *common,
                "--domtbl",
                str(dom2),
                "--hmm",
                str(hmm),
                "--search-option",
                "nobias",
                "--hmmer-version",
                "3.4",
            ]
        )
        == 0
    )
    rec2 = json.loads((wd / "modules/hydrophobin_relaxed.json").read_text())
    assert (
        rec1["artefact_hash"] != rec2["artefact_hash"]
        and rec1["params_hash"] != rec2["params_hash"]
    )


def test_cli_needs_the_r0_table(tmp_path, hmm):
    fasta = tmp_path / "p.faa"
    fasta.write_text(f">A one\n{EIGHT}\n")
    dom = domtbl(tmp_path / "d.tbl", [row("A")])
    code = main(
        [
            "hydrophobin_relaxed",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--domtbl",
            str(dom),
            "--hmm",
            str(hmm),
            "--search-option",
            "nobias",
            "--hmmer-version",
            "3.4",
        ]
    )
    assert code == 2
