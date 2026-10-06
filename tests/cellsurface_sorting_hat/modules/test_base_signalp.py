import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.base import (
    ModuleSpec,
    artefact_hash,
    params_hash,
    sha256_file,
    write_module,
)
from cellsurface_sorting_hat.modules.signalp import (
    COLUMNS,
    SignalPFormatError,
    parse_signalp,
    signalp_rows,
)

SIGNALP = (
    "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 20260930100126\n"
    "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
    "CIMG_00013-t26_1-p1 | transcript=CIMG_00013-t26_1 | gene=CIMG_00013 | organism=Coccidioides_immitis_RS"
    "\tOTHER\t1.000000\t0.000000\t\n"
    "CIMG_04613-t26_1-p1 | gene=CIMG_04613 | protein_length=324\tSP\t0.000100\t0.999800\tCS pos: 19-20. Pr: 0.9\n"
    "XP_3\tLIPO\t0.1\t0.2\tCS pos: 20-21. Pr: 0.8\n"
)


def prot(pid, state="ok"):
    return Protein(pid, "MKT", "x" * 64, state, "", 0.0)


def test_parse_signalp_uses_the_first_token_of_the_header_as_id(tmp_path):
    path = tmp_path / "prediction_results.txt"
    path.write_text(SIGNALP)
    got = parse_signalp(path)
    assert set(got) == {"CIMG_00013-t26_1-p1", "CIMG_04613-t26_1-p1", "XP_3"}
    assert got["CIMG_04613-t26_1-p1"] == {"prediction": "SP", "sp_prob": 0.9998}


def test_only_sp_is_a_call_lipo_is_not(tmp_path):
    path = tmp_path / "p.txt"
    path.write_text(SIGNALP)
    rows = signalp_rows(
        [prot("CIMG_04613-t26_1-p1"), prot("XP_3"), prot("CIMG_00013-t26_1-p1")],
        parse_signalp(path),
    )
    assert [r["call"] for r in rows] == ["called", "not_called", "not_called"]
    assert rows[1]["prediction"] == "LIPO"


def test_missing_and_invalid_proteins_are_not_ok(tmp_path):
    path = tmp_path / "p.txt"
    path.write_text(SIGNALP)
    rows = signalp_rows([prot("ABSENT"), prot("BAD", state="na_invalid")], parse_signalp(path))
    assert [(r["id"], r["state"]) for r in rows] == [("ABSENT", "error"), ("BAD", "na_invalid")]


@pytest.mark.parametrize(
    "text,message",
    [
        ("", "no prediction rows"),
        ("A\tSP\t0.1\n", "at least 4"),
        ("A\tSP\t0.1\tx\n", "not a number"),
        ("A\tSP\t0.1\t0.9\nA\tSP\t0.1\t0.9\n", "duplicate ID"),
    ],
)
def test_bad_signalp_files_are_refused(tmp_path, text, message):
    path = tmp_path / "p.txt"
    path.write_text(text)
    with pytest.raises(SignalPFormatError, match=message):
        parse_signalp(path)


def test_hashes(tmp_path):
    a, b = tmp_path / "a.bin", tmp_path / "b.bin"
    a.write_bytes(b"1")
    b.write_bytes(b"2")
    assert sha256_file(a) != sha256_file(b)
    assert params_hash({"x": 1, "y": 2}) == params_hash({"y": 2, "x": 1})
    assert params_hash({"x": 1}) != params_hash({"x": 2})
    assert artefact_hash([a, b]) == artefact_hash([b, a])  # order does not matter
    assert artefact_hash([a]) != artefact_hash([b])


def test_write_module_writes_the_table_and_the_run_record(tmp_path):
    db = tmp_path / "db.hmm"
    db.write_bytes(b"hmm")
    spec = ModuleSpec("step1_rule@R0", "1", {"mode": "fast"}, (db,), {"signalp": "6.0h"})
    rows = [
        {"id": "A", "state": "ok", "call": "called", "sp_prob": "0.9", "prediction": "SP"},
        {"id": "B", "state": "na_invalid"},
    ]
    rec = write_module(tmp_path, spec, COLUMNS, rows)
    with gzip.open(tmp_path / "modules" / "step1_rule@R0.tsv.gz", "rt") as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    assert table[0]["call"] == "called" and table[1]["call"] == ""
    on_disk = json.loads((tmp_path / "modules" / "step1_rule@R0.json").read_text())
    assert on_disk == rec
    assert rec["params_hash"] == params_hash({"mode": "fast"})
    assert rec["artefact_hash"] == artefact_hash([db]) and rec["n_rows"] == 2


def test_write_module_refuses_a_duplicate_id(tmp_path):
    with pytest.raises(ValueError, match="duplicate row"):
        write_module(
            tmp_path,
            ModuleSpec("m", "1"),
            [],
            [{"id": "A", "state": "ok"}, {"id": "A", "state": "ok"}],
        )
