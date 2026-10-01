import gzip
import json

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


def _setup_main(tmp_path, truth, all_sources=True):
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    _write(input_dir, "c.fa", CGD_FASTA)
    species = tmp_path / "species.tsv"
    species.write_text("source_id\tfasta_file\tid_mapping\nX\tc.fa\tcgd\n")
    manifest_tsv = tmp_path / "manifest.tsv"
    manifest.write_manifest(manifest_tsv, [_manifest_row()])
    truth_table.write_tsv(tmp_path / "truth_set.tsv.gz", truth_table.TRUTH_COLUMNS, truth)
    (tmp_path / "extract_log.json").write_text(json.dumps({"all_sources": all_sources}))
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


def test_main_stop_deletes_previous_truth_sequences(tmp_path, capsys):
    attach = load_script("02_attach_sequences")
    argv = _setup_main(tmp_path, [_full_row("CAL1", "C1_00010W_A", "P-ext")])
    assert attach.main(argv) == 0
    assert (tmp_path / "truth_sequences.tsv.gz").exists()
    truth_table.write_tsv(
        tmp_path / "truth_set.tsv.gz",
        truth_table.TRUTH_COLUMNS,
        [_full_row("CAL5", "C5_00000W_A", "P-ext")],
    )
    assert attach.main(argv) == 2
    assert "STOP:" in capsys.readouterr().err
    assert not (tmp_path / "truth_sequences.tsv.gz").exists()
    assert truth_table.read_tsv(tmp_path / "unmatched_ids.tsv")[0]["gene_id"] == "CAL5"
    assert json.loads((tmp_path / "sequence_run.json").read_text())["truth_set_sha256"]
    assert not list(tmp_path.glob(".tmp.*"))


def test_main_missing_manifest_row_stops(tmp_path, capsys):
    attach = load_script("02_attach_sequences")
    argv = _setup_main(tmp_path, [_full_row("CAL1", "C1_00010W_A", "P-ext")])
    manifest.write_manifest(tmp_path / "manifest.tsv", [])
    assert attach.main(argv) == 2
    assert "no manifest row" in capsys.readouterr().err


# ---- review round 1 ----


def _run(tmp_path, fasta, truth, mapping="cgd", sources=("X",)):
    attach = load_script("02_attach_sequences")
    input_dir = tmp_path / "in"
    input_dir.mkdir(exist_ok=True)
    species, manifest_rows = [], []
    for name, text in fasta.items():
        _write(input_dir, name, text)
        manifest_rows.append({"file": name, "mode": "record", "sha256": ""})
    for source_id, name in zip(sources, fasta, strict=False):
        species.append({"source_id": source_id, "fasta_file": name, "id_mapping": mapping})
    return attach, lambda **kw: attach.run(truth, species, manifest_rows, input_dir, tmp_path, **kw)


def test_all_outputs_columns_and_content(tmp_path):
    attach, go = _run(
        tmp_path,
        {"c.fa": ">C1_00010W_A\nmstq\nka*\n>C1_00020C_A\nMKKLLV\n"},
        [_row("CAL1", "C1_00010W_A", "P-ext"), _row("CAL2", "tRNA1", "N-int")],
    )
    go(truth_set_sha256="abc")
    rows = truth_table.read_tsv(tmp_path / "truth_sequences.tsv.gz")
    with gzip.open(tmp_path / "truth_sequences.tsv.gz", "rt") as handle:
        assert handle.readline().rstrip("\n").split("\t") == [
            "source_id", "gene_id", "label", "fasta_id", "length", "seq_sha256", "sequence",
        ]  # fmt: skip
    assert len(rows) == 1 and rows[0]["sequence"] == "MSTQKA" and "*" not in rows[0]["sequence"]
    assert rows[0]["length"] == str(len("MSTQKA"))
    assert rows[0]["seq_sha256"] == seqhash.seq_sha256("MSTQKA")
    header = (tmp_path / "unmatched_ids.tsv").read_text().splitlines()[0]
    assert header.split("\t") == [
        "source_id", "gene_id", "symbol", "synonym1", "label", "reason",
    ]  # fmt: skip
    assert truth_table.read_tsv(tmp_path / "unmatched_ids.tsv")[0]["reason"] == "no_fasta_record"
    header = (tmp_path / "sequence_counts.tsv").read_text().splitlines()[0]
    assert header.split("\t") == list(attach.SEQ_COUNT_COLUMNS)
    assert attach.SEQ_COUNT_COLUMNS[:3] == ("source_id", "truth_set_sha256", "fasta_file")
    assert truth_table.read_tsv(tmp_path / "sequence_counts.tsv")[0]["truth_set_sha256"] == "abc"
    assert not list(tmp_path.glob(".tmp.*"))


def test_empty_cleaned_sequence_is_unmatched_and_stops_for_p_ext(tmp_path):
    attach, go = _run(
        tmp_path,
        {"c.fa": ">C1_00010W_A\n*\n>C1_00020C_A\nMK\n"},
        [_row("CAL1", "C1_00010W_A", "P-ext"), _row("CAL2", "C1_00020C_A", "N-int")],
    )
    with pytest.raises(sequences.MappingError, match="P-ext gene CAL1"):
        go()
    reason = truth_table.read_tsv(tmp_path / "unmatched_ids.tsv")[0]["reason"]
    assert reason == "empty_sequence"


def test_empty_cleaned_sequence_for_n_int_is_only_unmatched(tmp_path):
    attach, go = _run(
        tmp_path,
        {"c.fa": ">C1_00010W_A\n*\n>C1_00020C_A\nMK\n"},
        [_row("CAL1", "C1_00010W_A", "N-int"), _row("CAL2", "C1_00020C_A", "N-int")],
    )
    seq_rows, unmatched, counts = go()
    assert [r["gene_id"] for r in seq_rows] == ["CAL2"]
    assert unmatched[0]["reason"] == "empty_sequence" and counts[0]["unmatched_n_int"] == "1"


def test_non_ascii_sequence_stops_with_diagnostics(tmp_path):
    attach, go = _run(
        tmp_path,
        {"c.fa": ">C1_00010W_A\nMK\u00c9A\n>C1_00020C_A\nMK\n"},
        [_row("CAL1", "C1_00010W_A", "N-int"), _row("CAL2", "C1_00020C_A", "N-int")],
    )
    (tmp_path / "truth_sequences.tsv.gz").write_bytes(b"old")
    with pytest.raises(sequences.MappingError, match="non-ASCII"):
        go()
    assert (tmp_path / "unmatched_ids.tsv").exists()
    assert (tmp_path / "sequence_counts.tsv").exists()
    assert not (tmp_path / "truth_sequences.tsv.gz").exists()


def test_gene_side_key_collision_raises_naming_both_genes():
    index = {"C1_00010W_A": ("C1_00010W_A", "MK")}
    rows = [_row("CAL1", "C1_00010W_A"), _row("CAL2", " c1_00010w_a ")]
    with pytest.raises(sequences.MappingError, match="CAL1 and CAL2"):
        sequences.attach(rows, index, "cgd")
    with pytest.raises(sequences.MappingError, match="CAL1 occurs twice"):
        sequences.attach([_row("CAL1", "A"), _row("CAL1", "A")], index, "cgd")
    with pytest.raises(sequences.MappingError, match="CAL1 occurs twice"):
        sequences.attach([_row("CAL1", "A"), _row("CAL1", "B")], index, "cgd")


def test_empty_gene_keys_do_not_collide():
    rows = [_row("CAL1", ""), _row("CAL2", "")]
    matched, unmatched = sequences.attach(rows, {}, "cgd")
    assert matched == [] and len(unmatched) == 2


def test_fasta_side_uniprot_isoform_is_stripped_but_not_for_cgd_or_pombase(tmp_path):
    index, _ = sequences.index_fasta(
        _write(tmp_path, "u.fa", ">sp|P22146-2|GAS1_YEAST x\nMK\n"), "uniprot"
    )
    assert set(index) == {"P22146"}
    index, _ = sequences.index_fasta(_write(tmp_path, "c.fa", ">C1_00010W-2\nMK\n"), "cgd")
    assert set(index) == {"C1_00010W-2"}
    index, _ = sequences.index_fasta(
        _write(tmp_path, "p.fa", ">SPAC1002.01-2.1:pep x\nMK\n"), "pombase"
    )
    assert set(index) == {"SPAC1002.01-2"}
    assert sequences.normalize_id("ab-2", "cgd") == "AB-2"
    assert sequences.normalize_id("ab-2", "pombase") == "AB-2"
    assert sequences.normalize_id("ab-2", "uniprot") == "AB"


def test_duplicates_differing_only_in_case_or_stop_are_merged(tmp_path):
    index, dups = sequences.index_fasta(
        _write(tmp_path, "c.fa", ">K\nMSTQKA\n>K\nmstqka*\n>K\nMST QKA\n"), "cgd"
    )
    assert dups == 2 and len(index) == 1
    with pytest.raises(sequences.MappingError, match="two different sequences"):
        sequences.index_fasta(_write(tmp_path, "d.fa", ">K\nMSTQKA\n>K\nMSTQKV*\n"), "cgd")


def test_zero_match_guard_is_per_source(tmp_path):
    attach, go = _run(
        tmp_path,
        {"a.fa": CGD_FASTA, "b.fa": CGD_FASTA},
        [
            _row("CAL1", "C1_00010W_A", "N-sec", "A"),
            _row("CAL2", "C9_99999W_A", "N-sec", "B"),
        ],
        sources=("A", "B"),
    )
    with pytest.raises(sequences.MappingError, match="B: no gene matched"):
        go()


@pytest.mark.parametrize("label", ["P-ext", "ambiguous"])
def test_each_unmatched_positive_label_stops_alone(tmp_path, label):
    attach, go = _run(
        tmp_path,
        {"c.fa": CGD_FASTA},
        [_row("CAL1", "C1_00010W_A", "N-sec"), _row("CAL5", "C5_00000W_A", label)],
    )
    with pytest.raises(sequences.MappingError, match=f"{label} gene CAL5"):
        go()


def test_partial_run_is_recorded_and_partial_truth_set_is_refused(tmp_path, capsys):
    attach = load_script("02_attach_sequences")
    argv = _setup_main(tmp_path, [_full_row("CAL1", "C1_00010W_A", "P-ext")], all_sources=False)
    assert attach.main(argv) == 2
    err = capsys.readouterr().err
    assert err.startswith("STOP:") and "all_sources" in err and "--allow-partial-truth-set" in err
    assert not (tmp_path / "truth_sequences.tsv.gz").exists()
    assert attach.main([*argv, "--allow-partial-truth-set"]) == 0
    log = json.loads((tmp_path / "sequence_run.json").read_text())
    assert log["all_sources"] is True and log["sources"] == ["X"]
    assert log["truth_set_sha256"] == manifest.sha256_file(tmp_path / "truth_set.tsv.gz")
    assert log["fasta_sha256"] == {"X": manifest.sha256_file(tmp_path / "in" / "c.fa")}
    assert (
        log["git_commit"]
        and log["python"]
        and log["arguments"] == [*argv, "--allow-partial-truth-set"]
    )


def test_missing_extract_log_is_refused(tmp_path, capsys):
    attach = load_script("02_attach_sequences")
    argv = _setup_main(tmp_path, [_full_row("CAL1", "C1_00010W_A", "P-ext")])
    (tmp_path / "extract_log.json").unlink()
    assert attach.main(argv) == 2
    assert "extract_log.json" in capsys.readouterr().err


def test_sources_subset_marks_run_as_not_all_sources(tmp_path):
    attach = load_script("02_attach_sequences")
    argv = _setup_main(tmp_path, [_full_row("CAL1", "C1_00010W_A", "P-ext")])
    species = tmp_path / "species.tsv"
    species.write_text("source_id\tfasta_file\tid_mapping\nX\tc.fa\tcgd\nY\tc.fa\tcgd\n")
    assert attach.main([*argv, "--sources", "X"]) == 0
    log = json.loads((tmp_path / "sequence_run.json").read_text())
    assert log["sources"] == ["X"] and log["all_sources"] is False
    assert attach.main(argv) == 0
    assert json.loads((tmp_path / "sequence_run.json").read_text())["all_sources"] is True


def test_run_json_is_byte_stable(tmp_path):
    attach = load_script("02_attach_sequences")
    argv = _setup_main(tmp_path, [_full_row("CAL1", "C1_00010W_A", "P-ext")])
    assert attach.main(argv) == 0
    first = (tmp_path / "sequence_run.json").read_bytes()
    assert attach.main(argv) == 0
    assert (tmp_path / "sequence_run.json").read_bytes() == first
