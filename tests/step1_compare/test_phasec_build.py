"""08_build_eval_tables.py: class mapping on real columns, dedupe, precedence (C-6), checks."""

import json

import pytest

np = pytest.importorskip("numpy")

import dedupe  # noqa: E402
import phasec_fixture as pf  # noqa: E402
import truth_table  # noqa: E402
from conftest import load_phasec, load_script  # noqa: E402


@pytest.fixture
def fx(tmp_path):
    return pf.make_work(tmp_path)


def _run(fx):
    return load_phasec("08_build_eval_tables").main(pf.build_argv(fx))


def _table(fx):
    return {
        r["seq_sha256"]: r
        for r in truth_table.read_tsv(fx["work"] / "phasec" / "eval_table.tsv.gz")
    }


def test_fixture_columns_equal_the_phase_b_constants():
    import keyword_tier
    import seqsets

    f07 = load_script("07_build_features")
    assert pf.MEMBER_FEATURE_COLUMNS == f07.MEMBER_FEATURE_COLUMNS
    assert pf.UNIQUE_FEATURE_COLUMNS == f07.UNIQUE_FEATURE_COLUMNS
    assert pf.UNIQUE_COLUMNS == seqsets.UNIQUE_COLUMNS
    assert pf.MEMBER_COLUMNS == seqsets.MEMBER_COLUMNS
    assert pf.KEYWORD_COLUMNS == keyword_tier.KEYWORD_COLUMNS


def test_build_writes_all_outputs_with_recorded_hashes(fx):
    import manifest

    assert _run(fx) == 0
    out = fx["work"] / "phasec"
    log = json.loads((out / "build_run.json").read_text())
    names = load_phasec("08_build_eval_tables").OUTPUT_NAMES
    assert sorted(p.name for p in out.iterdir()) == sorted(names)
    for name in names[:-1]:
        assert log["outputs_sha256"][name] == manifest.sha256_file(out / name)
    assert log["all_sources"] is True and log["truth_set_sha256"] == pf.TRUTH_SHA


def test_alternate_files_excluded(fx):
    assert _run(fx) == 0
    table = _table(fx)
    assert fx["named"]["ALT_ONLY"] not in table
    assert not any("Spom_SCHPO-mod" in r["source_ids"] for r in table.values())
    # the Spom genes are in the table once, from Spom_PomBase only
    spom = [r for r in table.values() if r["source_ids"] == "Spom_PomBase"]
    assert len(spom) == 6 + 8 + 6 + 1


def test_unlabelled_genes_are_not_in_the_table(fx):
    assert _run(fx) == 0
    assert not any("unlabelled" in r["label"] for r in _table(fx).values())


def test_tc_row_yields_to_go_label(fx):
    assert _run(fx) == 0
    table = _table(fx)
    n = fx["named"]
    assert table[n["SHARED_NSEC"]]["class"] == "neg" and table[n["SHARED_NSEC"]]["origin"] == "go"
    assert table[n["SHARED_UNRES"]]["class"] == "excluded"
    assert table[n["SHARED_UNRES"]]["stratum"] == "pm-unresolved"
    assert table[n["SHARED_AMBIG"]]["class"] == "excluded"
    assert table[n["SHARED_AMBIG"]]["stratum"] == "ambiguous"
    assert table[n["SHARED_POS"]]["origin"] == "go"
    log = truth_table.read_tsv(fx["work"] / "phasec" / "eval_dedupe_log.tsv")
    wins = [r for r in log if r["reason"] == "go_label_wins"]
    assert {r["seq_sha256"] for r in wins} == {
        n[k] for k in ("SHARED_NSEC", "SHARED_UNRES", "SHARED_AMBIG", "SHARED_POS")
    }
    run = json.loads((fx["work"] / "phasec" / "build_run.json").read_text())
    assert run["tc_rows_dropped_by_go_class"] == {"excluded": 2, "neg": 1, "pos": 1}
    tc = [r for r in table.values() if r["origin"] == "tc"]
    assert len(tc) == 16 + 1  # 16 T-c rows with own sequences, one row for the duplicate pair
    assert table[n["TC_DUP"]]["gene_ids"] == "QDUP01,QDUP02"
    assert all(r["class"] == "pos" for r in tc)


def test_dedupe_drops_a_hash_present_in_both_classes():
    base = {"label": "P-ext", "subset": "wall", "stratum": "wall", "d8_class": "",
            "homology_only": "no", "internal_evidence_htp_only": "no", "species": "s",
            "role": "train", "clade": "c", "taxon_id": "1", "length": "100", "emb_row": "0",
            "emb_cterm_row": ""}  # fmt: skip
    pos = {**base, "seq_sha256": "a", "class": "pos", "source_id": "S", "gene_id": "g1"}
    neg = {**base, "seq_sha256": "a", "class": "neg", "source_id": "S", "gene_id": "g2",
           "label": "N-sec", "subset": "", "stratum": "N-sec"}  # fmt: skip
    keep = {**base, "seq_sha256": "b", "class": "pos", "source_id": "S", "gene_id": "g3"}
    keep2 = {**keep, "source_id": "T", "gene_id": "g4", "homology_only": "yes"}
    rows, log = dedupe.merge_go([pos, neg, keep, keep2])
    assert set(rows) == {"b"}
    assert rows["b"]["source_ids"] == "S,T" and rows["b"]["homology_only"] == "no"
    assert [(r["gene_id"], r["reason"]) for r in log] == [
        ("g1", "both_classes"),
        ("g2", "both_classes"),
    ]
    excl = {
        **keep,
        "gene_id": "g5",
        "class": "excluded",
        "label": "ambiguous",
        "stratum": "ambiguous",
    }
    rows, log = dedupe.merge_go([keep, excl])
    assert rows["b"]["class"] == "excluded" and {r["reason"] for r in log} == {"class_and_excluded"}


def test_truth_set_has_no_iea(fx):
    path = fx["work"] / "truth_set_triaged.tsv.gz"
    rows = truth_table.read_tsv(path)
    row = next(r for r in rows if r["label"] == "N-sec")
    row["evidence_codes"] = "IEA"
    truth_table.write_tsv(path, list(rows[0]), rows)
    pf.rewrite_json(fx["work"] / "phaseb" / "features_run.json", input_sha256={
        **json.loads((fx["work"] / "phaseb" / "features_run.json").read_text())["input_sha256"],
        "truth_set_triaged.tsv.gz": pf._sha(path)})  # fmt: skip
    m = load_phasec("08_build_eval_tables")
    assert m.iea_problems(rows) == [
        f"{row['source_id']}:{row['gene_id']} has label N-sec with IEA evidence only"
    ]
    assert _run(fx) == 2
    assert not (fx["work"] / "phasec").exists()


def test_iea_in_an_evidence_column_is_found():
    m = load_phasec("08_build_eval_tables")
    row = {"source_id": "S", "gene_id": "g", "label": "P-ext", "evidence_codes": "IDA,IEA",
           "surface_evidence": "IDA,IEA", "internal_evidence": "", "secretory_evidence": ""}  # fmt: skip
    assert m.iea_problems([row]) == ["S:g has IEA in surface_evidence"]


def test_tc_taxon_table_is_complete(fx, capsys):
    rows = truth_table.read_tsv(fx["clades"])
    truth_table.write_tsv(
        fx["clades"], list(rows[0]), [r for r in rows if r["taxon_id"] != "246410"]
    )
    assert _run(fx) == 2
    assert "tc_taxon_clades.tsv has no row for taxon_id 246410" in capsys.readouterr().err


def test_committed_tc_taxon_table_covers_the_phase_a_taxa():
    # The seven taxon_ids of the Phase A keyword_tier.tsv.gz (3,092 rows, counted 2026-10-01).
    m = load_phasec("08_build_eval_tables")
    import evalio

    clades = m.read_tc_clades(evalio.PHASEC_DIR / "tc_taxon_clades.tsv")
    assert clades == {
        "237561": "Saccharomycotina", "246410": "Eurotiomycetes", "284593": "Saccharomycotina",
        "330879": "Eurotiomycetes", "443226": "Eurotiomycetes", "498019": "Saccharomycotina",
        "559292": "Saccharomycotina",
    }  # fmt: skip
    species = truth_table.read_tsv(m.paths.STEP1_DIR / "species.tsv")
    assert set(clades.values()) <= {s["in_clade"] for s in species} | {"other"}


def test_stale_input_stops(fx, capsys):
    # M3: features and embeddings from different unique sequence sets
    pf.rewrite_json(
        fx["work"] / "phaseb" / "emb" / "embedding_run.json", unique_sequences_sha256="0" * 64
    )
    assert _run(fx) == 2
    assert "differs from embedding_run.json unique_sequences_sha256" in capsys.readouterr().err
    assert not (fx["work"] / "phasec").exists()


def test_changed_unique_sequences_file_stops(fx, capsys):
    path = fx["work"] / "phaseb" / "unique_sequences.tsv.gz"
    path.write_bytes(path.read_bytes() + b"\n")
    assert _run(fx) == 2
    assert "differs from the file 07 read" in capsys.readouterr().err


def test_truth_hash_mismatch_stops(fx, capsys):
    pf.rewrite_json(fx["work"] / "keyword_tier_run.json", truth_set_sha256="9" * 64)
    assert _run(fx) == 2
    assert "truth_set_sha256 differs" in capsys.readouterr().err


def test_literature_rows(fx):
    assert _run(fx) == 0
    lit = {
        r["gene"]: r for r in truth_table.read_tsv(fx["work"] / "phasec" / "eval_literature.tsv")
    }
    assert set(lit) == {"LIT1", "LIT2", "LIT3"}
    assert [lit[g]["literature_positive"] for g in ("LIT1", "LIT2", "LIT3")] == ["yes", "yes", "no"]
    run = json.loads((fx["work"] / "phasec" / "build_run.json").read_text())
    assert run["literature_no_accession"] == ["LIT4"] and run["literature_positives"] == 2


def test_fasta_holds_table_and_literature_sequences(fx):
    assert _run(fx) == 0
    lines = pf.gz_lines(fx["work"] / "phasec" / "eval_sequences.fasta.gz")
    heads = [x[1:] for x in lines if x.startswith(">")]
    table = _table(fx)
    lit = truth_table.read_tsv(fx["work"] / "phasec" / "eval_literature.tsv")
    assert heads == sorted(set(table) | {r["seq_sha256"] for r in lit})
    assert lines[lines.index(">" + heads[0]) + 1] == fx["seqs"][heads[0]]


def test_missing_feature_value_stops(fx, capsys):
    # Review Focus 3: a table hash with an empty sp_prob must stop, not become NaN
    path = fx["work"] / "phaseb" / "features_unique.tsv.gz"
    rows = truth_table.read_tsv(path)
    target = fx["named"]["SHARED_POS"]
    for r in rows:
        if r["seq_sha256"] == target:
            r["sp_prob"] = ""
    truth_table.write_tsv(path, pf.UNIQUE_FEATURE_COLUMNS, rows)
    assert _run(fx) == 2
    assert f"{target} sp_prob: '' is not a number" in capsys.readouterr().err


def test_identical_sequence_in_two_sources_is_one_row(fx):
    # Review Focus 4: one row; both sources listed
    path = fx["work"] / "phaseb" / "features.tsv.gz"
    rows = truth_table.read_tsv(path)
    afum = next(r for r in rows if r["source_id"] == "Afum_ASPFU" and r["label"] == "N-sec")
    h99 = next(r for r in rows if r["source_id"] == "Cneo_H99_GOA" and r["label"] == "N-sec")
    for c in ("seq_sha256", "length", "emb_row", "emb_cterm_row"):
        h99[c] = afum[c]
    truth_table.write_tsv(path, pf.MEMBER_FEATURE_COLUMNS, rows)
    assert _run(fx) == 0
    row = _table(fx)[afum["seq_sha256"]]
    assert row["source_ids"] == "Afum_ASPFU,Cneo_H99_GOA"
    assert row["clades"] == "Basidiomycota,Eurotiomycetes"
