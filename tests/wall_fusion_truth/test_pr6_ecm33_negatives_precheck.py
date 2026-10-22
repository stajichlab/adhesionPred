"""The pre-freeze check of the ecm33 hard negatives (fixtures; the real shared negative set for the expected members)."""

from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "analysis/wall_fusion_truth"
NEG = ROOT / "data/controls/hard-negatives-shared"


def load(name):
    spec = importlib.util.spec_from_file_location(name, WF / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


pre = load("pr6_ecm33_negatives_precheck")


def tbl_line(target, evalue, score, query="Q5AGC4"):
    cols = [
        target,
        "-",
        query,
        "-",
        f"{evalue:g}",
        f"{score}",
        "0.1",
        f"{evalue:g}",
        f"{score}",
        "0.1",
    ]
    cols += ["1.0", "1", "0", "0", "1", "1", "1", "0", "desc of " + target]
    return " ".join(cols)


def write_tbl(path, hits, ok=True):
    text = "#header\n" + "\n".join(tbl_line(*h) for h in hits) + "\n"
    text += "# [ok]\n" if ok else ""
    path.write_text(text)
    return path


CONTROLS = {
    a: {"accession": a, "gene": g, "species": sp, "stratum": st}
    for a, g, sp, st in [
        ("P38248", "ECM33", "Saccharomyces cerevisiae", "gpi_wall_enzyme"),
        ("Q4WNS8", "ecm33", "Aspergillus fumigatus", "gpi_wall_enzyme"),
        ("X1", "xyz1", "Candida albicans", "wall_hydrolase"),
        ("N1", "n1", "Candida albicans", "laccase_etc"),
        ("N2", "n2", "Candida albicans", "laccase_etc"),
    ]
}
LENGTHS = {"P38248": 429, "Q4WNS8": 400, "X1": 350, "N1": 500, "N2": 510}


def test_parse_tblout_needs_the_ok_trailer_and_full_lines(tmp_path):
    ok = write_tbl(tmp_path / "ok.tbl", [("P38248", 1e-80, 250.0)])
    assert pre.parse_tblout(ok) == [
        {"target": "P38248", "query": "Q5AGC4", "evalue": 1e-80, "score": 250.0}
    ]
    with pytest.raises(ValueError, match="trailer"):
        pre.parse_tblout(write_tbl(tmp_path / "bad.tbl", [("P38248", 1e-80, 250.0)], ok=False))
    short = tmp_path / "short.tbl"
    short.write_text("P38248 - Q5AGC4\n# [ok]\n")
    with pytest.raises(ValueError, match="short line"):
        pre.parse_tblout(short)


def test_the_table_the_proposed_list_and_the_edit_line():
    strict = [
        {"target": "P38248", "query": "Q", "evalue": 1e-80, "score": 250.0},
        {"target": "Q4WNS8", "query": "Q", "evalue": 1e-70, "score": 230.0},
        {
            "target": "X1",
            "query": "Q",
            "evalue": 1e-9,
            "score": 40.0,
        },  # a member under another name
    ]
    weak = strict + [
        {
            "target": "N1",
            "query": "Q",
            "evalue": 1e-4,
            "score": 20.0,
        },  # 1e-5 < E <= 1e-3: descriptive only
        {"target": "N2", "query": "Q", "evalue": 0.5, "score": 5.0},  # above 1e-3: dropped
    ]
    res = pre.build(strict, weak, CONTROLS, LENGTHS, ["P38248"])
    assert [r["accession"] for r in res["rows"]] == ["P38248", "Q4WNS8", "X1"]  # by E-value
    assert [r["listed_in_freeze"] for r in res["rows"]] == ["yes", "no", "no"]
    assert res["rows"][1]["gene"] == "ecm33" and res["rows"][1]["stratum"] == "gpi_wall_enzyme"
    assert res["rows"][2]["length"] == 350
    assert [r["accession"] for r in res["weak_rows"]] == ["N1"]  # not in the proposed list
    assert res["proposed"] == ["P38248", "Q4WNS8", "X1"] and res["new"] == ["Q4WNS8", "X1"]
    assert res["loud"] == []
    assert res["strata"] == {"gpi_wall_enzyme": ["P38248", "Q4WNS8"], "wall_hydrolase": ["X1"]}
    assert pre.json_edit(res["proposed"]) == '"negatives_ignored": ["P38248", "Q4WNS8", "X1"],'


def test_the_loud_lines_for_a_missing_expected_member_an_unknown_hit_and_an_unfound_listed_one():
    strict = [
        {"target": "P38248", "query": "Q", "evalue": 1e-80, "score": 250.0},
        {"target": "ZZ", "query": "Q", "evalue": 1e-20, "score": 90.0},  # not in controls.tsv
    ]
    res = pre.build(strict, strict, CONTROLS, LENGTHS, ["P38248", "N2"])
    text = "\n".join(res["loud"])
    assert "expected family member Q4WNS8 is NOT found" in text
    assert "hit accession ZZ is not in controls.tsv" in text
    assert "already listed in the freeze but not found by the search: N2" in text
    assert "P38248 is NOT found" not in text
    # the proposed list keeps what the freeze lists; it only adds what the search finds in controls.tsv
    assert res["proposed"] == ["N2", "P38248"]


def test_an_evalue_above_the_definition_is_not_a_member():
    strict = [{"target": "X1", "query": "Q", "evalue": 2e-5, "score": 30.0}]
    res = pre.build(strict, strict, CONTROLS, LENGTHS, [])
    assert res["rows"] == [] and [r["accession"] for r in res["weak_rows"]] == ["X1"]


def test_main_writes_the_files_without_sequences_and_refuses_a_used_directory(tmp_path, capsys):
    strict = write_tbl(
        tmp_path / "s.tbl", [("P38248", 1e-80, 250.0), ("Q4WNS8", 1e-70, 230.0), ("X1", 1e-9, 40.0)]
    )
    weak = write_tbl(tmp_path / "w.tbl", [("P38248", 1e-80, 250.0), ("N1", 1e-4, 20.0)])
    controls = tmp_path / "controls.tsv"
    with open(controls, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "gene", "species", "stratum"])
        for c in CONTROLS.values():
            w.writerow([c["accession"], c["gene"], c["species"], c["stratum"]])
    fasta = tmp_path / "neg.faa"
    fasta.write_text("".join(f">{a}\n{'ACDEFGHIKL' * (n // 10)}\n" for a, n in LENGTHS.items()))
    out = tmp_path / "out"
    args = [
        "summarize",
        "--strict",
        str(strict),
        "--weak",
        str(weak),
        "--controls",
        str(controls),
        "--negatives-fasta",
        str(fasta),
        "--already-ignored",
        "P38248",
        "--out",
        str(out),
    ]
    assert pre.main(args) == 0
    assert sorted(p.name for p in out.iterdir()) == [
        "negatives_family_table.tsv",
        "proposed_negatives_ignored.json",
        "summary.txt",
        "weak_E1e-3.tsv",
    ]
    assert (out / "proposed_negatives_ignored.json").read_text() == (
        '"negatives_ignored": ["P38248", "Q4WNS8", "X1"],\n'
    )
    summary = (out / "summary.txt").read_text()
    assert "NEW (not yet listed): Q4WNS8, X1" in summary and "LOUD" not in summary
    # no sequence in any output (the fixture sequence is ACDEFGHIKL repeated)
    assert all("ACDEFGHIKLACDEFGHIKL" not in p.read_text() for p in out.iterdir())
    tab = list(csv.DictReader(open(out / "negatives_family_table.tsv"), delimiter="\t"))
    assert [r["accession"] for r in tab] == ["P38248", "Q4WNS8", "X1"]
    assert pre.main(args) == 1  # a directory that is already filled is never overwritten
    assert "is not empty" in capsys.readouterr().err


def test_the_real_shared_negatives_hold_both_known_ecm33_members_and_the_fasta_matches():
    with open(NEG / "controls.tsv", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by = {r["accession"]: r for r in rows}
    assert len(rows) == 161
    for acc in pre.EXPECTED:
        assert by[acc]["stratum"] == "gpi_wall_enzyme" and by[acc]["label"] == "negative"
    lengths = pre.fasta_lengths(NEG / "sequences.faa")
    assert set(lengths) == set(
        by
    )  # the ids stay exactly as in the file (no #neg copies in the shared set)


SBATCH = (WF / "pr6_ecm33_negatives_precheck.sbatch").read_text()


def test_the_job_script_is_a_safe_draft():
    assert (
        subprocess.run(["bash", "-n", str(WF / "pr6_ecm33_negatives_precheck.sbatch")]).returncode
        == 0
    )
    head = "\n".join(SBATCH.splitlines()[:6])
    assert "DRAFT, NOT SUBMITTED" in head
    assert (
        "#SBATCH -p short" in SBATCH and "#SBATCH -c 2" in SBATCH and "#SBATCH --mem=2G" in SBATCH
    )
    assert "${SCRATCH:?" in SBATCH and "/tmp" not in SBATCH
    assert "module load hmmer/3.4" in SBATCH
    # the pin check, the refusal of an existing OUT, the sidecar check, the exact-token query extraction
    assert (
        'git -C "$PROJ_ROOT" diff --quiet "$CODE_COMMIT" HEAD' in SBATCH
        and "status --porcelain" in SBATCH
    )
    assert '[ ! -e "$OUT" ]' in SBATCH and "sidecar_ok" in SBATCH
    assert "name == q" in SBATCH and "exactly once" in SBATCH
    # the definition of the freeze (E 1e-5) and a separate descriptive pass (E 1e-3)
    assert "phmmer -E 1e-5 --domE 1e-5" in SBATCH and "phmmer -E 1e-3 --domE 1e-3" in SBATCH
    # the freeze file is neither read nor changed by this job
    assert "pr6_ecm33_freeze.json" not in SBATCH
    assert "pr6_ecm33_negatives_precheck.py" in SBATCH and "--already-ignored P38248" in SBATCH
