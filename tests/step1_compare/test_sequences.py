import gzip

import manifest
import pytest
import seqhash
import sequences
import truth_table
from conftest import load_script

SGD_FASTA = (
    ">YAL001C TFC3 SGDID:S000000001, Chr I from 151006-147594, Genome Release 64-5-1\n"
    "MVLTIYPDELVQ\nIVSDKI*\n"
    ">YAL002W VPS8 SGDID:S000000002, Chr I\nMEQNGLDHDSRS\n"
)
CGD_FASTA = ">C1_00010W_A\nMSTQKA\n>C1_00020C_A\nMKKLLV\n>C1_00020C_A\nMKKLLV\n"
POMBASE_FASTA = ">SPAC1002.01.1:pep mrx11|component\nMSKF\n>SPBC21H7.03c.1:pep x\nMAAA\n"
UNIPROT_FASTA = ">sp|P22146|GAS1_YEAST 1,3-beta-glucanosyltransferase\nMLFKS\n>tr|Q4WXC4|Q4WXC4_ASPFU CspA\nMKVA\n"


def _write(tmp_path, name, text, compress=False):
    path = tmp_path / name
    data = text.encode()
    path.write_bytes(gzip.compress(data) if compress else data)
    return path


def _row(gene_id, synonym1="", label="P-ext", source="X"):
    return {
        "source_id": source,
        "gene_id": gene_id,
        "symbol": gene_id,
        "synonym1": synonym1,
        "label": label,
    }


def test_sgd_mapping_reads_sgdid(tmp_path):
    index, dups = sequences.index_fasta(_write(tmp_path, "s.fa.gz", SGD_FASTA, True), "sgd")
    assert set(index) == {"S000000001", "S000000002"} and dups == 0
    assert index["S000000001"] == ("YAL001C", "MVLTIYPDELVQIVSDKI*")


def test_cgd_mapping_uses_first_synonym_and_merges_identical_duplicates(tmp_path):
    index, dups = sequences.index_fasta(_write(tmp_path, "c.fa", CGD_FASTA), "cgd")
    assert dups == 1
    matched, unmatched = sequences.attach(
        [_row("CAL0000173921", "C1_00010W_A"), _row("CAL0000000001", "snR5a")], index, "cgd"
    )
    assert [m["fasta_id"] for m in matched] == ["C1_00010W_A"]
    assert [u["gene_id"] for u in unmatched] == ["CAL0000000001"]


def test_duplicate_key_with_different_sequence_raises(tmp_path):
    bad = _write(tmp_path, "c.fa", ">C1_00010W_A\nMSTQKA\n>C1_00010W_A\nMSTQKV\n")
    with pytest.raises(sequences.MappingError, match="two different sequences"):
        sequences.index_fasta(bad, "cgd")


def test_pombase_suffix_and_case_insensitive_match(tmp_path):
    index, _ = sequences.index_fasta(_write(tmp_path, "p.fa", POMBASE_FASTA), "pombase")
    matched, unmatched = sequences.attach(
        [_row("SPAC1002.01"), _row("spbc21h7.03C"), _row("SPAC1002.02")], index, "pombase"
    )
    assert [m["fasta_id"] for m in matched] == ["SPAC1002.01.1:pep", "SPBC21H7.03c.1:pep"]
    assert [u["gene_id"] for u in unmatched] == ["SPAC1002.02"]


def test_uniprot_isoform_and_version_differences(tmp_path):
    index, _ = sequences.index_fasta(_write(tmp_path, "u.fa", UNIPROT_FASTA), "uniprot")
    matched, _ = sequences.attach([_row("P22146-2"), _row(" q4wxc4 ")], index, "uniprot")
    assert [m["fasta_id"] for m in matched] == ["sp|P22146|GAS1_YEAST", "tr|Q4WXC4|Q4WXC4_ASPFU"]


def test_clean_and_hash_ignore_case_whitespace_and_stop():
    assert seqhash.clean("mk j*\nA") == "MKLA"
    assert seqhash.seq_sha256("MKLA") == seqhash.seq_sha256("mkja*")


def test_clean_matches_package_clean_sequence():
    io_mod = pytest.importorskip("surface_glyco.io")
    for seq in ("MKJA*", "MSTQKA", "JJ*"):
        assert seqhash.clean(seq) == io_mod.clean_sequence(seq)


def test_attach_script_reports_unmatched_by_label(tmp_path):
    attach = load_script("02_attach_sequences")
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    _write(input_dir, "c.fa", CGD_FASTA)
    species = [{"source_id": "X", "fasta_file": "c.fa", "id_mapping": "cgd"}]
    manifest_rows = [{"file": "c.fa", "mode": "record", "sha256": ""}]
    truth = [
        _row("CAL1", "C1_00010W_A", "P-ext"),
        _row("CAL2", "tRNA1", "N-int"),
        _row("CAL3", "C1_00020C_A", "unlabelled"),
    ]
    seq_rows, unmatched, counts = attach.run(truth, species, manifest_rows, input_dir, tmp_path)
    assert counts[0]["matched"] == "2" and counts[0]["unmatched_n_int"] == "1"
    assert counts[0]["fasta_duplicate_records"] == "1"
    assert counts[0]["fasta_sha256"] == manifest.sha256_file(input_dir / "c.fa")
    assert seq_rows[0]["seq_sha256"] == seqhash.seq_sha256("MSTQKA")
    assert truth_table.read_tsv(tmp_path / "unmatched_ids.tsv")[0]["gene_id"] == "CAL2"


def test_attach_stops_when_a_source_matches_nothing(tmp_path):
    attach = load_script("02_attach_sequences")
    _write(tmp_path, "c.fa", CGD_FASTA)
    species = [{"source_id": "X", "fasta_file": "c.fa", "id_mapping": "cgd"}]
    manifest_rows = [{"file": "c.fa", "mode": "record", "sha256": ""}]
    truth = [_row("CAL9", "C9_99999W_A", "N-sec")]
    with pytest.raises(sequences.MappingError, match="no gene matched"):
        attach.run(truth, species, manifest_rows, tmp_path, tmp_path)
    assert not (tmp_path / "truth_sequences.tsv.gz").exists()


def test_attach_stops_when_a_positive_or_ambiguous_gene_is_unmatched(tmp_path):
    attach = load_script("02_attach_sequences")
    _write(tmp_path, "c.fa", CGD_FASTA)
    species = [{"source_id": "X", "fasta_file": "c.fa", "id_mapping": "cgd"}]
    manifest_rows = [{"file": "c.fa", "mode": "record", "sha256": ""}]
    truth = [_row("CAL1", "C1_00010W_A", "N-sec"), _row("CAL5", "C5_00000W_A", "ambiguous")]
    with pytest.raises(sequences.MappingError, match="ambiguous gene CAL5"):
        attach.run(truth, species, manifest_rows, tmp_path, tmp_path)
    assert truth_table.read_tsv(tmp_path / "unmatched_ids.tsv")[0]["gene_id"] == "CAL5"


def _manifest_row():
    return {c: "" for c in manifest.MANIFEST_COLUMNS} | {"file": "c.fa", "mode": "record"}


def _setup_main(tmp_path, truth):
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    _write(input_dir, "c.fa", CGD_FASTA)
    species = tmp_path / "species.tsv"
    species.write_text("source_id\tfasta_file\tid_mapping\nX\tc.fa\tcgd\n")
    manifest_tsv = tmp_path / "manifest.tsv"
    manifest.write_manifest(manifest_tsv, [_manifest_row()])
    truth_table.write_tsv(tmp_path / "truth_set.tsv.gz", truth_table.TRUTH_COLUMNS, truth)
    return [
        "--species", str(species), "--manifest", str(manifest_tsv),
        "--input-dir", str(input_dir), "--work-dir", str(tmp_path),
    ]  # fmt: skip


def _full_row(gene_id, synonym1, label):
    return {c: "" for c in truth_table.TRUTH_COLUMNS} | _row(gene_id, synonym1, label)


def test_main_writes_outputs_and_no_temp_files(tmp_path, capsys):
    attach = load_script("02_attach_sequences")
    argv = _setup_main(tmp_path, [_full_row("CAL1", "C1_00010W_A", "P-ext")])
    assert attach.main(argv) == 0
    names = sorted(p.name for p in tmp_path.iterdir() if p.name.startswith((".tmp", "truth_seq")))
    assert names == ["truth_sequences.tsv.gz"]
    assert "matched=1" in capsys.readouterr().out


def test_main_unknown_source_stops_with_valid_ids(tmp_path, capsys):
    attach = load_script("02_attach_sequences")
    argv = _setup_main(tmp_path, [_full_row("CAL1", "C1_00010W_A", "P-ext")])
    assert attach.main([*argv, "--sources", "Nope"]) == 2
    err = capsys.readouterr().err
    assert err.startswith("STOP:") and "valid ids: X" in err
    assert not (tmp_path / "truth_sequences.tsv.gz").exists()


def test_main_stop_keeps_previous_truth_sequences(tmp_path, capsys):
    attach = load_script("02_attach_sequences")
    argv = _setup_main(tmp_path, [_full_row("CAL1", "C1_00010W_A", "P-ext")])
    assert attach.main(argv) == 0
    before = (tmp_path / "truth_sequences.tsv.gz").read_bytes()
    truth_table.write_tsv(
        tmp_path / "truth_set.tsv.gz",
        truth_table.TRUTH_COLUMNS,
        [_full_row("CAL5", "C5_00000W_A", "P-ext")],
    )
    assert attach.main(argv) == 2
    assert "STOP:" in capsys.readouterr().err
    assert (tmp_path / "truth_sequences.tsv.gz").read_bytes() == before
    assert not list(tmp_path.glob(".tmp.*"))


def test_main_missing_manifest_row_stops(tmp_path, capsys):
    attach = load_script("02_attach_sequences")
    argv = _setup_main(tmp_path, [_full_row("CAL1", "C1_00010W_A", "P-ext")])
    manifest.write_manifest(tmp_path / "manifest.tsv", [])
    assert attach.main(argv) == 2
    assert "no manifest row" in capsys.readouterr().err
