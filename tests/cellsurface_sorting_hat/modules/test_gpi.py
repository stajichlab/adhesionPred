"""The ``gpi`` wrapper: PredGPI gff3 -> module table (SYNTHETIC rows only)."""

import csv
import gzip

import pytest

from cellsurface_sorting_hat.fasta import read_fasta
from cellsurface_sorting_hat.modules import gpi
from cellsurface_sorting_hat.modules.cli import main

LONG = "M" + "A" * 59  # 60 residues
SHORT = "M" + "A" * 39  # 40 residues: never scored by PredGPI
FASTA = (
    f">P_gpi x\n{LONG}\n>P_chain y\n{LONG}\n>P_short z\n{SHORT}\n>P_none w\n{LONG}\n>BAD\nMK*T\n"
)


def row(pid, kind, start, end, score="1.0", source="PredGPI"):
    return (
        "\t".join([pid, source, kind, str(start), str(end), score, ".", ".", "evidence=x"]) + "\n"
    )


GOOD = (
    row("P_gpi", "GPI-anchor", 55, 55, "0.70")
    + row("P_chain", "Chain", 1, 60)
    + row("P_short", "Chain", 1, 40)
)


@pytest.fixture
def fasta(tmp_path):
    p = tmp_path / "p.faa"
    p.write_text(FASTA)
    return p


def lengths():
    return {p.id: len(p.sequence) for p in read_fasta_text()}


def read_fasta_text():
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "x.faa"
        f.write_text(FASTA)
        return read_fasta(f)


def test_parse_gives_kind_omega_and_confidence_class(tmp_path):
    f = tmp_path / "g.gff3"
    f.write_text(GOOD)
    out = gpi.parse_gpi(f, lengths())
    assert out["P_gpi"] == {"kind": "GPI-anchor", "omega": 55, "cls": "probable", "ub": "0.005"}
    assert out["P_chain"]["kind"] == "Chain" and out["P_chain"]["omega"] is None


@pytest.mark.parametrize(
    "score, cls, ub",
    [("1.0", "high", "0.0015"), ("0.70", "probable", "0.005"), ("0.55", "weak", "0.01")],
)
def test_confidence_classes(tmp_path, score, cls, ub):
    f = tmp_path / "g.gff3"
    f.write_text(row("P_gpi", "GPI-anchor", 50, 50, score))
    out = gpi.parse_gpi(f, lengths())["P_gpi"]
    assert (out["cls"], out["ub"]) == (cls, ub)


def test_rows_call_not_called_not_assessable_error_and_invalid(tmp_path):
    f = tmp_path / "g.gff3"
    f.write_text(GOOD)
    proteins = read_fasta_text()
    rows = {r["id"]: r for r in gpi.gpi_rows(proteins, gpi.parse_gpi(f, lengths()))}
    assert (
        rows["P_gpi"]["call"] == "called"
        and rows["P_gpi"]["omega_site"] == 55
        and rows["P_gpi"]["confidence_class"] == "probable"
    )
    assert rows["P_chain"]["call"] == "not_called" and rows["P_chain"]["omega_site"] == ""
    assert rows["P_short"] == {
        "id": "P_short",
        "state": "not_assessable",
    }  # a silent PredGPI "no" is not a not_called
    assert rows["P_none"] == {"id": "P_none", "state": "error"}  # the tool returned nothing for it
    assert rows["BAD"]["state"] == "na_invalid"


@pytest.mark.parametrize(
    "text, msg",
    [
        (row("P_gpi", "Region", 1, 60), "unknown row type"),
        (row("P_gpi", "GPI-anchor", 50, 51), "one site"),
        (row("P_gpi", "GPI-anchor", 70, 70), "one site"),
        (row("P_chain", "Chain", 1, 59), "does not span"),
        (row("P_gpi", "GPI-anchor", 50, 50, "0.9"), "unknown confidence score"),
        (row("P_gpi", "GPI-anchor", 50, 50, "x"), "not a number"),
        (row("P_gpi", "GPI-anchor", 50, 50, source="Other"), "expected PredGPI"),
        (row("P_gpi", "Chain", 1, 60) * 2, "duplicate ID"),
        (row("P_other", "Chain", 1, 60), "not a protein of the FASTA"),
        ("P_gpi\tPredGPI\tChain\t1\t60\n", "9 TAB separated"),
        ("", "no rows"),
    ],
)
def test_stop_conditions(tmp_path, text, msg):
    f = tmp_path / "g.gff3"
    f.write_text(text)
    with pytest.raises(gpi.GpiFormatError, match=msg):
        gpi.parse_gpi(f, lengths())


def read_module(wd):
    with gzip.open(wd / "modules" / "gpi.tsv.gz", "rt") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_command_writes_module_and_partial_run_state(tmp_path, fasta):
    g = tmp_path / "g.gff3"
    g.write_text(GOOD)
    wd = tmp_path / "wd"
    code = main(
        [
            "gpi",
            "--fasta",
            str(fasta),
            "--workdir",
            str(wd),
            "--gff3",
            str(g),
            "--predgpi-version",
            "202001",
        ]
    )
    assert code == 0
    rec = __import__("json").loads((wd / "modules" / "gpi.json").read_text())
    assert (
        rec["module"] == "gpi"
        and rec["run_state"] == "partial"
        and rec["tools"] == {"predgpi": "202001"}
    )  # P_none has no row
    assert read_module(wd)["P_gpi"]["call"] == "called"


def test_command_stops_on_ids_that_are_not_in_the_fasta(tmp_path, fasta, capsys):
    g = tmp_path / "g.gff3"
    g.write_text(row("OTHER", "Chain", 1, 60))
    code = main(
        [
            "gpi",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--gff3",
            str(g),
            "--predgpi-version",
            "202001",
        ]
    )
    assert code == 2 and "error" in capsys.readouterr().err


def test_identity_depends_on_the_tool_version(tmp_path, fasta):
    g = tmp_path / "g.gff3"
    g.write_text(GOOD)
    hashes = []
    for v, d in (("202001", "a"), ("202002", "b")):
        wd = tmp_path / d
        main(
            [
                "gpi",
                "--fasta",
                str(fasta),
                "--workdir",
                str(wd),
                "--gff3",
                str(g),
                "--predgpi-version",
                v,
            ]
        )
        hashes.append(
            __import__("json").loads((wd / "modules" / "gpi.json").read_text())["artefact_hash"]
        )
    assert hashes[0] != hashes[1]
