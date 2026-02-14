import tempfile
from pathlib import Path

from adhesion_predict.io import find_fasta_files, process_fasta_file


class TestFindFastaFiles:
    def test_find_fasta_files_faa(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "test.faa"
            test_file.write_text(">seq1\nMSEQ\n")

            result = find_fasta_files(tmpdir)
            assert len(result) == 1
            assert result[0].name == "test.faa"

    def test_find_fasta_files_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "sequences.fasta"
            test_file.write_text(">seq1\nMSEQ\n")

            result = find_fasta_files(tmpdir)
            assert len(result) == 1
            assert result[0].name == "sequences.fasta"

    def test_find_fasta_files_pep(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "proteins.pep"
            test_file.write_text(">seq1\nMSEQ\n")

            result = find_fasta_files(tmpdir)
            assert len(result) == 1

    def test_find_fasta_files_uppercase_ext(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "test.FAA"
            test_file.write_text(">seq1\nMSEQ\n")

            result = find_fasta_files(tmpdir)
            assert len(result) == 1

    def test_find_fasta_files_nested(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            subdir = Path(tmpdir) / "subdir"
            subdir.mkdir()
            test_file = subdir / "nested.faa"
            test_file.write_text(">seq1\nMSEQ\n")

            result = find_fasta_files(tmpdir)
            assert len(result) == 1

    def test_find_fasta_files_multiple_extensions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "test1.faa").write_text(">seq1\nMSEQ\n")
            (Path(tmpdir) / "test2.fasta").write_text(">seq2\nMSEQ\n")
            (Path(tmpdir) / "test3.pep").write_text(">seq3\nMSEQ\n")

            result = find_fasta_files(tmpdir)
            assert len(result) == 3

    def test_find_fasta_files_no_matches(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "test.txt").write_text("not a fasta")

            result = find_fasta_files(tmpdir)
            assert len(result) == 0

    def test_find_fasta_files_empty_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = find_fasta_files(tmpdir)
            assert len(result) == 0


class TestProcessFastaFile:
    def test_process_fasta_file_basic(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "test.faa"
            test_file.write_text(">seq1\nMSEQ\n>seq2\nMPQRS\n")

            result = process_fasta_file(test_file)
            assert len(result) == 2
            assert result[0]["id"] == "seq1"
            assert result[0]["sequence"] == "MSEQ"
            assert result[1]["id"] == "seq2"
            assert result[1]["sequence"] == "MPQRS"

    def test_process_fasta_file_single_sequence(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "test.faa"
            test_file.write_text(">protein1\nMKLLIA\n")

            result = process_fasta_file(test_file)
            assert len(result) == 1
            assert result[0]["id"] == "protein1"
            assert result[0]["sequence"] == "MKLLIA"

    def test_process_fasta_file_empty_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "empty.faa"
            test_file.write_text("")

            result = process_fasta_file(test_file)
            assert len(result) == 0

    def test_process_fasta_file_invalid_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "invalid.faa"
            test_file.write_text("not valid fasta content")

            result = process_fasta_file(test_file)
            assert len(result) == 0

    def test_process_fasta_file_with_existing_test_file(self):
        test_file = Path("tests/input_tests/Saccharomyces.pep")
        if test_file.exists():
            result = process_fasta_file(test_file)
            assert isinstance(result, list)
