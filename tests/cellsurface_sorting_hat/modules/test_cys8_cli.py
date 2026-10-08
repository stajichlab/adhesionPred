"""The cys8 module command end to end (synthetic sequences)."""

import csv
import gzip
import json
from pathlib import Path

from cellsurface_sorting_hat.modules.cli import main

SPACING = Path(__file__).resolve().parents[3] / "data/sorting_hat/cys8_spacing.yaml"
SEQ = (
    "MKLLVAAG"
    + "C"
    + "A" * 10
    + "CC"
    + "A" * 11
    + "C"
    + "A" * 16
    + "C"
    + "A" * 8
    + "CC"
    + "A" * 10
    + "CGSTNPAD"
)


def _read(wd, name):
    with gzip.open(wd / "modules" / f"{name}.tsv.gz", "rt") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def _setup(tmp_path):
    fasta = tmp_path / "p.faa"
    fasta.write_text(f">XP_1 a\n{SEQ}\n>XP_2 b\n{SEQ}\n>XP_3 c\nMKTAYIAKQRQAAAAAA\n")
    res = tmp_path / "prediction_results.txt"
    res.write_text(
        "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 1\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
        "XP_1 x\tSP\t0.1\t0.9\t\nXP_2 y\tOTHER\t0.9\t0.1\t\nXP_3 z\tOTHER\t0.9\t0.1\t\n"
    )
    wd = tmp_path / "wd"
    args = ["--fasta", str(fasta), "--workdir", str(wd)]
    assert main(["signalp", *args, "--results", str(res), "--signalp-version", "6"]) == 0
    return args, wd


def test_cys8_command_writes_module_with_r0_identity(tmp_path, capsys):
    args, wd = _setup(tmp_path)
    assert main(["cys8", *args, "--spacing", str(SPACING)]) == 0
    rows = _read(wd, "cys8_pattern")
    assert (rows["XP_1"]["hit"], rows["XP_2"]["hit"], rows["XP_3"]["hit"]) == ("1", "0", "0")
    assert rows["XP_2"]["pattern_match"] == "1" and rows["XP_3"]["pattern_match"] == "0"
    run = json.loads((wd / "modules" / "cys8_pattern.json").read_text())
    cond = run["params"]["conditions"]["sp_module"]
    r0 = json.loads((wd / "modules" / "step1_rule@R0.json").read_text())
    assert cond["module"] == "step1_rule@R0" and cond["params_hash"] == r0["params_hash"]
    assert run["run_state"] == "ok"


def test_cys8_command_needs_the_r0_table(tmp_path, capsys):
    fasta = tmp_path / "p.faa"
    fasta.write_text(f">XP_1 a\n{SEQ}\n")
    code = main(
        [
            "cys8",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--spacing",
            str(SPACING),
        ]
    )
    assert code == 2
