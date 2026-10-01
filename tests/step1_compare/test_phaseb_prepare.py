import gzip
import json

import pytest
import seqhash
import seqsets
import truth_table
from conftest import load_script

LONG = "MKV" + "ST" * 600  # 1,203 aa: has a C-terminal window
SHARED = "MKTAYIAKQRQISFVKSHFSRQ"


def _truth_rows():
    rows = []
    for source_id, gene_id, seq, label in (
        ("Scer_SGD", "S1", SHARED, "P-ext"),
        ("Scer_SGD", "S2", LONG, "N-sec"),
        ("Calb_CGD", "C1", SHARED, "N-int"),  # same sequence as S1 in another source
    ):
        rows.append(
            {
                "source_id": source_id,
                "gene_id": gene_id,
                "label": label,
                "fasta_id": gene_id,
                "length": str(len(seq)),
                "seq_sha256": seqhash.seq_sha256(seq),
                "sequence": seq,
            }
        )
    return rows


def _phaseb_inputs(tmp_path, truth_rows=None):
    work = tmp_path / "work"
    downloads = work / "downloads"
    downloads.mkdir(parents=True)
    truth_table.write_tsv(
        work / "truth_sequences.tsv.gz", list(_truth_rows()[0]), truth_rows or _truth_rows()
    )
    with gzip.open(work / "keyword_sequences.fasta.gz", "wt") as handle:
        handle.write(f">sp|P11111|KW1_YEAST desc\n{SHARED}\n>tr|Q22222|KW2_CANAL\nMSSPLLA*\n")
    (downloads / "prot.fasta").write_text(
        f">YAL001C desc\n{SHARED.lower()}\n>YAL002W\nMKLLV\n>YAL003W\n\n"
    )
    site = tmp_path / "site"
    site.mkdir()
    (site / "rs.fasta").write_text(">CIMG_1-t1-p1 | gene=CIMG_1\nMPPPP\n")
    for name in ("sequence_run.json", "keyword_tier_run.json"):
        (work / name).write_text(json.dumps({"all_sources": True}))
    sets = [
        {"set_id": "truth", "kind": "truth", "location": "truth_sequences.tsv.gz", "note": ""},
        {"set_id": "kw", "kind": "keyword", "location": "keyword_sequences.fasta.gz", "note": ""},
        {"set_id": "Scer_proteome", "kind": "download", "location": "prot.fasta", "note": ""},
        {"set_id": "Cimm", "kind": "site", "location": "cocci:rs.fasta", "note": ""},
    ]
    sets_path = tmp_path / "sets.tsv"
    truth_table.write_tsv(sets_path, seqsets.SET_COLUMNS, sets)
    return work, downloads, site, sets_path


def test_prepare_dedupes_by_hash_across_sets(tmp_path, monkeypatch):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 0
    out = work / "phaseb"
    members = truth_table.read_tsv(out / "sequence_members.tsv.gz")
    unique = truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    # 3 truth + 2 keyword + 2 proteome (one empty record skipped) + 1 site = 8 members
    assert len(members) == 8
    # SHARED (5 members), LONG, MSSPLLA, MKLLV, MPPPP
    assert len(unique) == 5
    assert [r["seq_sha256"] for r in unique] == sorted(r["seq_sha256"] for r in unique)
    assert [r["row"] for r in unique] == ["0", "1", "2", "3", "4"]
    shared = seqhash.seq_sha256(SHARED)
    assert sum(m["seq_sha256"] == shared for m in members) == 4
    kw = {m["gene_id"] for m in members if m["set_id"] == "kw"}
    assert kw == {"P11111", "Q22222"}
    cimm = [m for m in members if m["set_id"] == "Cimm"]
    assert cimm[0]["gene_id"] == "CIMG_1-t1-p1" and cimm[0]["source_id"] == "Cimm"
    log = json.loads((out / "prepare_run.json").read_text())
    assert log["all_sources"] is True
    assert log["inputs"]["Scer_proteome"]["empty_records"] == 1
    assert log["over_max_residues"] == 1


def test_cterm_row_numbers_only_long_sequences(tmp_path):
    seqs = {seqhash.seq_sha256(s): s for s in (SHARED, LONG, "A" * 1023, "A" * 1022)}
    rows = seqsets.unique_rows(seqs)
    assert len(rows) == 4
    long_rows = [r for r in rows if int(r["length"]) > 1022]
    assert sorted(r["cterm_row"] for r in long_rows) == ["0", "1"]
    assert all(r["cterm_row"] == "" for r in rows if int(r["length"]) <= 1022)


def test_prepare_stops_on_a_non_esm_character(tmp_path, monkeypatch, capsys):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    (downloads / "prot.fasta").write_text(">YAL009W\nMK-LV\n")
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 2
    assert "not ESM-2 residues" in capsys.readouterr().err
    assert not (work / "phaseb").exists()


def test_prepare_stops_on_a_stale_truth_hash(tmp_path, monkeypatch, capsys):
    rows = _truth_rows()
    rows[0]["seq_sha256"] = "0" * 64
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path, rows)
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 2
    assert "re-run 02_attach_sequences.py" in capsys.readouterr().err


def test_prepare_refuses_a_partial_truth_set(tmp_path, monkeypatch, capsys):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    (work / "sequence_run.json").write_text(json.dumps({"all_sources": False}))
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 2
    assert "02_attach_sequences.py" in capsys.readouterr().err
    assert prepare.main(argv + ["--allow-partial-truth-set"]) == 0


def test_prepare_stops_on_a_missing_input(tmp_path, monkeypatch, capsys):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    (downloads / "prot.fasta").unlink()
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 2
    assert "not found" in capsys.readouterr().err


def test_same_gene_with_two_sequences_stops(tmp_path):
    path = tmp_path / "dup.fasta"
    path.write_text(">G1\nMKV\n>G1\nMKL\n")
    with pytest.raises(seqsets.SequenceSetError, match="two different sequences"):
        seqsets.fasta_members("x", path, "download")


def test_unique_fasta_headers_are_hashes(tmp_path, monkeypatch):
    work, downloads, site, sets_path = _phaseb_inputs(tmp_path)
    prepare = load_script("05_prepare_sequences")
    monkeypatch.setattr(prepare.paths, "site_value", lambda key: str(site))
    argv = ["--sets", str(sets_path), "--work-dir", str(work), "--input-dir", str(downloads)]
    assert prepare.main(argv) == 0
    text = gzip.open(work / "phaseb" / "unique_sequences.fasta.gz", "rt").read()
    headers = [line[1:] for line in text.splitlines() if line.startswith(">")]
    unique = truth_table.read_tsv(work / "phaseb" / "unique_sequences.tsv.gz")
    assert headers == [r["seq_sha256"] for r in unique]
    first = (work / "phaseb" / "unique_sequences.fasta.gz").read_bytes()
    assert prepare.main(argv) == 0
    assert (work / "phaseb" / "unique_sequences.fasta.gz").read_bytes() == first  # byte-stable


def test_real_sets_file_is_valid():
    import paths

    rows = seqsets.read_sets(paths.STEP1_DIR / "sequence_sets.tsv")
    assert [r["set_id"] for r in rows][:2] == ["truth", "uniprot_kw"]
    assert sum(r["kind"] == "download" for r in rows) == 8
    assert any(r["set_id"] == "Cimm_RS_proteome" for r in rows)


def test_truth_row_with_a_non_esm_character_stops():
    seq = "MK-LV"  # seqhash.clean keeps '-', sanitize_sequence would delete it
    row = {
        "source_id": "Scer_SGD",
        "gene_id": "S9",
        "seq_sha256": seqhash.seq_sha256(seq),
        "sequence": seq,
    }
    with pytest.raises(seqsets.SequenceSetError, match="not ESM-2 residues"):
        seqsets.truth_members([row])


def test_keyword_isoform_collision_names_the_cause(tmp_path):
    path = tmp_path / "kw.fasta"
    path.write_text(">sp|P11111|A_YEAST\nMKV\n>sp|P11111-2|A_YEAST\nMKVLL\n")
    with pytest.raises(seqsets.SequenceSetError, match="isoform suffix"):
        seqsets.fasta_members("uniprot_kw", path, "keyword")
