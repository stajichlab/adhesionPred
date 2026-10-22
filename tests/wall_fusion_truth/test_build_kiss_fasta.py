"""Tests of build_kiss_fasta.py (plan Task 5.2/5.3 input). Fixture files only; no JGI data."""

import csv
import gzip
import hashlib
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "build_kiss_fasta", ROOT / "analysis/wall_fusion_truth/build_kiss_fasta.py"
)
bkf = importlib.util.module_from_spec(spec)
sys.modules["build_kiss_fasta"] = bkf
spec.loader.exec_module(bkf)

PUB = "M" + "ACDEFGHIKL" * 6
LOC = "M" + "LKIHGFEDCA" * 6
REP = "M" + "KLMNPQRSTV" * 6
NEG = "M" + "VTSRQPNMLK" * 6


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def write_tsv(path, header, rows):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def make(tmp_path, neg_accession="N00001", tamper_local=False):
    cdir = tmp_path / "set"
    cdir.mkdir()
    cols = ["accession", "stratum", "cluster_id", "jgi_accept", "sequence_sha256", "protein_id_jgi"]
    write_tsv(
        cdir / "controls.tsv",
        cols,
        [["P00001", "wall", "KISS-wall-1", "row_verified", sha(PUB), "1_Neucr2"]],
    )
    write_tsv(
        cdir / "controls_local_sequence.tsv",
        cols,
        [
            [
                "jgi:Neucr2:44",
                "adhesion",
                "KISS-adhesion-2",
                "portal_validated",
                sha(LOC + "A" if tamper_local else LOC),
                "44_Neucr2",
            ]
        ],
    )
    (cdir / "sequences.faa").write_text(f">P00001\n{PUB}\n")
    write_tsv(
        cdir / "family_reps.tsv",
        ["accession", "slice", "kiss_cluster_number", "sequence_sha256"],
        [["P00009", "wall", "7", sha(REP)]],
    )
    (cdir / "family_rep_sequences.faa").write_text(f">P00009\n{REP}\n")
    ndir = tmp_path / "negset"
    ndir.mkdir()
    write_tsv(ndir / "controls.tsv", ["accession", "label"], [[neg_accession, "negative"]])
    (ndir / "sequences.faa").write_text(f">{neg_accession}\n{NEG}\n")
    jdir = tmp_path / "jgi"
    f = jdir / "Neucr2" / "Mycocosm" / "Proteins" / "Neucr2.fasta.gz"
    f.parent.mkdir(parents=True)
    with gzip.open(f, "wt") as fh:
        fh.write(f">jgi|Neucr2|44|NCU00044T0\n{LOC}*\n>jgi|Neucr2|45|NCU00045T0\n{REP}*\n")
    table = tmp_path / "jgi_files.tsv"
    write_tsv(
        table,
        ["portal", "organism", "file", "role", "n_proteins", "sha256"],
        [["Neucr2", "x", f.name, "use", 2, hashlib.sha256(f.read_bytes()).hexdigest()]],
    )
    return cdir, ndir, jdir, table


def test_public_and_local_sequences_are_joined(tmp_path):
    cdir, _, jdir, table = make(tmp_path)
    seqs, members = bkf.gather(cdir, jdir, table, [], False)
    assert seqs == {"P00001": PUB, "jgi:Neucr2:44": LOC}
    assert {m["accession"]: m["class"] for m in members} == {
        "P00001": "positive",
        "jgi:Neucr2:44": "positive",
    }


def test_family_reps_and_negatives_only_when_asked_for(tmp_path):
    cdir, ndir, jdir, table = make(tmp_path)
    seqs, members = bkf.gather(cdir, jdir, table, [ndir], True)
    cls = {m["accession"]: m["class"] for m in members}
    assert cls["P00009"] == "family_rep" and cls["N00001"] == "negative"
    assert members[-1]["slice"] == "" and members[-1]["source"] == "negset"


def test_a_sequence_that_does_not_match_its_sha256_stops_everything(tmp_path):
    cdir, _, jdir, table = make(tmp_path, tamper_local=True)
    with pytest.raises(ValueError, match="differs from sequence_sha256"):
        bkf.gather(cdir, jdir, table, [], False)
    (tmp_path / "again").mkdir()
    cdir2, _, jdir2, table2 = make(tmp_path / "again")
    (cdir2 / "sequences.faa").write_text(f">P00001\n{PUB}K\n")
    with pytest.raises(ValueError, match="sequences.faa differs"):
        bkf.gather(cdir2, jdir2, table2, [], False)


def test_a_changed_jgi_file_is_refused_by_the_file_table(tmp_path):
    cdir, _, jdir, table = make(tmp_path)
    f = next(jdir.rglob("Neucr2.fasta.gz"))
    with gzip.open(f, "wt") as fh:
        fh.write(f">jgi|Neucr2|44|NCU00044T0\n{LOC}*\n")
    with pytest.raises(ValueError, match="sha256"):
        bkf.gather(cdir, jdir, table, [], False)


def test_refuses_an_output_directory_inside_the_repository(tmp_path):
    cdir, _, jdir, table = make(tmp_path)
    inside = ROOT / "analysis" / "wall_fusion_truth" / "should_not_exist"
    with pytest.raises(SystemExit, match="inside the repository"):
        bkf.main(
            [
                "--jgi-dir",
                str(jdir),
                "--out-dir",
                str(inside),
                "--control-dir",
                str(cdir),
                "--jgi-table",
                str(table),
            ]
        )
    assert not inside.exists()
    with pytest.raises(SystemExit, match="inside the repository"):
        bkf.refuse_repo_path(ROOT)


def test_a_protein_in_two_classes_is_an_error_unless_keep_both(tmp_path):
    cdir, ndir, jdir, table = make(tmp_path, neg_accession="P00001")
    with pytest.raises(bkf.OverlapError, match="P00001"):
        bkf.gather(cdir, jdir, table, [ndir], False)
    seqs, members = bkf.gather(cdir, jdir, table, [ndir], False, on_overlap="keep-both")
    assert seqs["P00001#neg"] == NEG
    neg = next(m for m in members if m["accession"] == "P00001#neg")
    assert neg["class"] == "negative" and neg["overlap_with"] == "P00001"


def test_outputs_are_written_outside_the_repo_with_checksums(tmp_path):
    cdir, ndir, jdir, table = make(tmp_path)
    out = tmp_path / "out"
    assert (
        bkf.main(
            [
                "--jgi-dir",
                str(jdir),
                "--out-dir",
                str(out),
                "--control-dir",
                str(cdir),
                "--jgi-table",
                str(table),
                "--negative-set",
                str(ndir),
            ]
        )
        == 0
    )
    fa = out / "kiss_sequences.faa"
    assert (out / "kiss_sequences.faa.sha256").read_text().strip() == hashlib.sha256(
        fa.read_bytes()
    ).hexdigest()
    assert LOC in fa.read_text()  # the JGI sequence goes to the local file, never into the repo
    assert (out / "members.tsv").read_text().splitlines()[0].startswith("accession\tclass")
