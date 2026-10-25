"""The ``disorder`` wrapper: per-residue AIUPred/IUPred3 text -> module table (SYNTHETIC rows only)."""

import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.modules import disorder
from cellsurface_sorting_hat.modules.cli import main

FASTA = ">P1 a\nMKSTA\n>P2 b\nMNLL\n>P3 c\nMKQ\n>BAD\nMK*T\n"
SEQ = {"P1": "MKSTA", "P2": "MNLL", "P3": "MKQ"}
BANNER = "#  AIUPred (synthetic banner)\n# Position\tResidue\tDisorder\n"


def rec(name, seq, scores, prefix="#>"):
    body = "".join(f"{i}\t{a}\t{s}\n" for i, (a, s) in enumerate(zip(seq, scores, strict=True), 1))
    return f"{prefix}{name}\n{body}\n"


AIUPRED = (
    BANNER + rec("P1", "MKSTA", [0.1, 0.6, 0.7, 0.8, 0.2]) + rec("P2", "MNLL", [0.5, 0.4, 0.5, 0.5])
)


def test_aiupred_layout_values_by_hand():
    out = disorder.parse_disorder_stream(AIUPRED.splitlines(), SEQ)
    assert out["P1"] == {"longest_run": 3, "frac": 0.6, "mean": pytest.approx(0.48), "n_res": 5}
    # P2: a score of exactly 0.5 counts as disordered; runs are [0.5] and [0.5, 0.5], so the longest run is 2 and 3 of 4 residues count
    assert out["P2"]["longest_run"] == 2 and out["P2"]["frac"] == 0.75


def test_iupred3_layout_and_one_shot_generator():
    text = rec("P1", "MKSTA", [0.9, 0.9, 0.1, 0.9, 0.9], prefix=">") + rec(
        "P3", "MKQ", [0.0, 0.0, 0.0], prefix=">"
    )
    out = disorder.parse_disorder_stream(
        (ln for ln in text.splitlines()), SEQ
    )  # a generator: one pass, no list needed
    assert (
        out["P1"]["longest_run"] == 2 and out["P3"]["longest_run"] == 0 and out["P3"]["frac"] == 0.0
    )


@pytest.mark.parametrize(
    "score, ok",
    [("1.04", True), ("-0.04", True), ("1.06", False), ("-0.06", False), ("nan", False)],
)
def test_score_slack(score, ok):
    text = f"#>P3\n1\tM\t{score}\n2\tK\t0.1\n3\tQ\t0.1\n"
    if ok:
        assert disorder.parse_disorder_stream(text.splitlines(), SEQ)["P3"]["n_res"] == 3
    else:
        with pytest.raises(disorder.DisorderFormatError, match="outside"):
            disorder.parse_disorder_stream(text.splitlines(), SEQ)


@pytest.mark.parametrize(
    "text, msg",
    [
        ("#>P3\n1\tM\t0.1\n3\tQ\t0.1\n", "out of order"),
        ("#>P3\n1\tM\t0.1\n2\tK\t0.1\n", "has 2 residues, the FASTA sequence has 3"),
        ("#>P3\n1\tM\t0.1\n2\tK\t0.1\n3\tW\t0.1\n", "differ from the FASTA"),
        (
            "#>P3\n1\tM\t0.1\n2\tK\t0.1\n3\tQ\t0.1\n#>P3\n1\tM\t0.1\n2\tK\t0.1\n3\tQ\t0.1\n",
            "duplicate record",
        ),
        ("#>P3\n#>P2\n1\tM\t0.1\n", "no residues"),
        ("1\tM\t0.1\n", "before any record"),
        ("#>P3\n1\tM\n", "pos aa score"),
        ("#>P3\nx\tM\t0.1\n", "not a number"),
        ("#>\n", "empty record name"),
        (BANNER, "no records"),
    ],
)
def test_stop_conditions(text, msg):
    with pytest.raises(disorder.DisorderFormatError, match=msg):
        disorder.parse_disorder_stream(text.splitlines(), SEQ)


def test_rows_include_error_for_a_missing_record_and_invalid_proteins(tmp_path):
    from cellsurface_sorting_hat.fasta import read_fasta

    f = tmp_path / "p.faa"
    f.write_text(FASTA)
    parsed = disorder.parse_disorder_stream(AIUPRED.splitlines(), SEQ)
    rows = {r["id"]: r for r in disorder.disorder_rows(read_fasta(f), parsed)}
    assert (
        rows["P1"]["longest_run"] == 3
        and rows["P1"]["n_res"] == 5
        and rows["P1"]["frac"] == "0.6000"
    )
    assert rows["P3"] == {"id": "P3", "state": "error"} and rows["BAD"]["state"] == "na_invalid"
    assert "call" not in rows["P1"]  # calls read longest_run through {test: ...}


def read_module(wd):
    with gzip.open(wd / "modules" / "disorder.tsv.gz", "rt") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_command_end_to_end_and_identity(tmp_path):
    fa, raw = tmp_path / "p.faa", tmp_path / "raw.txt"
    fa.write_text(FASTA)
    raw.write_text(AIUPRED)
    hashes = []
    for v, d in (("2.x", "a"), ("2.y", "b")):
        wd = tmp_path / d
        code = main(
            [
                "disorder",
                "--fasta",
                str(fa),
                "--workdir",
                str(wd),
                "--raw",
                str(raw),
                "--tool",
                "aiupred",
                "--tool-version",
                v,
            ]
        )
        assert code == 0
        record = json.loads((wd / "modules" / "disorder.json").read_text())
        assert (
            record["run_state"] == "partial" and record["params"]["residue_threshold"] == 0.5
        )  # P3 has no record
        hashes.append(record["artefact_hash"])
    assert hashes[0] != hashes[1]
    assert read_module(tmp_path / "a")["P1"]["longest_run"] == "3"


def test_command_stops_on_unknown_ids_and_on_a_residue_mismatch(tmp_path, capsys):
    fa, raw = tmp_path / "p.faa", tmp_path / "raw.txt"
    fa.write_text(FASTA)
    args = [
        "disorder",
        "--fasta",
        str(fa),
        "--workdir",
        str(tmp_path / "wd"),
        "--raw",
        str(raw),
        "--tool",
        "aiupred",
        "--tool-version",
        "2.x",
    ]
    raw.write_text(rec("OTHER", "MKSTA", [0.1] * 5))
    assert main(args) == 2 and "not in the FASTA" in capsys.readouterr().err
    raw.write_text(rec("P1", "MKSTW", [0.1] * 5))
    assert main(args) == 2 and "differ from the FASTA" in capsys.readouterr().err
