import gzip
import hashlib
import shutil
import subprocess

import pytest

from cellsurface_sorting_hat.fasta import NA_INVALID, OK, FastaError, read_fasta


def _write(tmp_path, text, name="p.faa"):
    path = tmp_path / name
    path.write_bytes(text.encode())
    return path


def test_reads_proteins_and_hashes_the_sequence(tmp_path):
    p = read_fasta(_write(tmp_path, ">A first protein\nMKT\nAYI\n>B\nMKTAYI\n"))
    assert [x.id for x in p] == ["A", "B"]
    assert p[0].sequence == "MKTAYI"
    assert p[0].sha256 == hashlib.sha256(b"MKTAYI").hexdigest()
    assert p[0].sha256 == p[1].sha256  # identical sequences, different IDs, are allowed


def test_trailing_stop_is_stripped_silently(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\nMKTAYI*\n>B\nMKT\n"))
    assert (p[0].sequence, p[0].state) == ("MKTAYI", OK)


def test_internal_stop_is_invalid(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\nMKT*AYI\n>B\nMKT\n"))
    assert p[0].state == NA_INVALID and "stop" in p[0].note


def test_non_residue_characters_are_invalid(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\nMKT-AY1\n>B\nMKT\n"))
    assert p[0].state == NA_INVALID and "-" in p[0].note


def test_ambiguous_residues_are_allowed_and_counted(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\nMKXXAYIB\n>B\nMKT\n"))
    assert p[0].state == OK
    assert p[0].ambiguous_fraction == pytest.approx(3 / 8)


def test_empty_sequence_is_invalid(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\n>B\nMKT\n"))
    assert p[0].state == NA_INVALID


def test_crlf_and_lowercase(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\r\nmktayi\r\n>B\r\nMKT\r\n"))
    assert p[0].sequence == "MKTAYI"


def test_duplicate_ids_stop_the_run(tmp_path):
    with pytest.raises(FastaError, match="duplicate ID: A"):
        read_fasta(_write(tmp_path, ">A\nMKT\n>A\nMKS\n"))


def test_empty_file_and_no_valid_protein_stop_the_run(tmp_path):
    with pytest.raises(FastaError, match="no records"):
        read_fasta(_write(tmp_path, ""))
    with pytest.raises(FastaError, match="no valid proteins"):
        read_fasta(_write(tmp_path, ">A\nMK*T\n"))


def test_gzip(tmp_path):
    path = tmp_path / "p.faa.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(">A\nMKT\n")
    assert read_fasta(path)[0].sequence == "MKT"


@pytest.mark.skipif(shutil.which("zstd") is None, reason="zstd is not installed")
def test_zstd(tmp_path):
    plain = _write(tmp_path, ">A\nMKT\n")
    subprocess.run(["zstd", "-q", str(plain), "-o", str(tmp_path / "p.faa.zst")], check=True)
    assert read_fasta(tmp_path / "p.faa.zst")[0].sequence == "MKT"


def test_trailing_stop_is_flagged(tmp_path):
    p = read_fasta(_write(tmp_path, ">A\nMKT*\n>B\nMKT\n"))
    assert (p[0].trailing_stop, p[1].trailing_stop) == (True, False)


def test_sequence_lines_before_the_first_header_stop_the_run(tmp_path):
    with pytest.raises(FastaError, match="before the first header"):
        read_fasta(_write(tmp_path, "MKT\n>A\nMKT\n"))


def test_a_byte_order_mark_is_ignored(tmp_path):
    path = tmp_path / "p.faa"
    path.write_bytes(b"\xef\xbb\xbf>A\nMKT\n")
    assert read_fasta(path)[0].id == "A"


def test_bytes_that_are_not_utf8_stop_the_run_with_the_file_name(tmp_path):
    path = tmp_path / "p.faa"
    path.write_bytes(b">A \xe9\nMKT\n")
    with pytest.raises(FastaError, match="p.faa: cannot be read as text"):
        read_fasta(path)


def test_a_corrupt_gz_names_the_file(tmp_path):
    path = _write(tmp_path, "not gzip data", name="p.faa.gz")
    with pytest.raises(FastaError, match="p.faa.gz"):
        read_fasta(path)


def test_a_missing_zstd_binary_names_the_file(tmp_path, monkeypatch):
    path = _write(tmp_path, "x", name="p.faa.zst")

    def boom(*a, **k):
        raise FileNotFoundError("zstd")

    monkeypatch.setattr(subprocess, "run", boom)
    with pytest.raises(FastaError, match="p.faa.zst"):
        read_fasta(path)


def test_a_failing_zstd_names_the_file(tmp_path, monkeypatch):
    path = _write(tmp_path, "x", name="p.faa.zst")

    def boom(*a, **k):
        raise subprocess.CalledProcessError(1, "zstd")

    monkeypatch.setattr(subprocess, "run", boom)
    with pytest.raises(FastaError, match="p.faa.zst"):
        read_fasta(path)
