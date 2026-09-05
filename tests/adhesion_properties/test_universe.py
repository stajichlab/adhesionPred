"""Tests for analysis/adhesion_properties/universe.py."""

import random

from universe import (
    build_protein_universe,
    load_taxonomy_by_locustag,
    read_fai_ids,
    read_result_ids,
    sample_background_ids,
)

_SUMMARY_HEADER = (
    "locustag,asmid,species_name,phylum,subphylum,class,order,family,genus,"
    "total_proteins,adhesion_count,adhesion_fraction,mean_adhesion_prob,median_adhesion_prob\n"
)


def test_read_result_ids(tmp_path):
    result = tmp_path / "sp.adhesion_predict.csv"
    result.write_text(
        "id,prediction,probability_adhesion\n"
        "LOC1_000001-T1,Adhesion,0.9\n"
        "LOC1_000002-T1,Adhesion,0.7\n"
    )
    assert read_result_ids(result) == [("LOC1_000001-T1", 0.9), ("LOC1_000002-T1", 0.7)]


def test_read_fai_ids(tmp_path):
    fai = tmp_path / "sp.proteins.fa.fai"
    fai.write_text("LOC1_000001-T1\t100\t10\t60\t61\nLOC1_000002-T1\t100\t10\t60\t61\n")
    assert read_fai_ids(fai) == ["LOC1_000001-T1", "LOC1_000002-T1"]


def test_load_taxonomy_by_locustag(tmp_path):
    summary = tmp_path / "species_adhesion_summary.csv"
    summary.write_text(
        _SUMMARY_HEADER + "LOC1,ASM1,Sp one,Ascomycota,Pezizomycotina,Sordariomycetes,Hypocreales,"
        "Nectriaceae,Fusarium,100,5,0.05,0.8,0.8\n"
    )
    lookup = load_taxonomy_by_locustag(summary)
    assert lookup["LOC1"]["phylum"] == "Ascomycota"
    assert lookup["LOC1"]["genus"] == "Fusarium"


def test_sample_background_ids_returns_all_when_scarce():
    rng = random.Random(1)
    result = sample_background_ids(["A", "B", "C"], {"A"}, n_wanted=10, rng=rng)
    assert sorted(result) == ["B", "C"]


def test_sample_background_ids_caps_at_n_wanted():
    rng = random.Random(1)
    all_ids = [f"P{i}" for i in range(100)]
    result = sample_background_ids(all_ids, set(), n_wanted=10, rng=rng)
    assert len(result) == 10
    assert set(result).issubset(set(all_ids))


def test_build_protein_universe_matched_species(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_one.proteins.fa.fai").write_text(
        "".join(f"LOC1_{i:06d}-T1\t100\t10\t60\t61\n" for i in range(1, 11))
    )
    (results_dir / "Sp_one.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\n"
        "LOC1_000001-T1,Adhesion,0.9\n"
        "LOC1_000002-T1,Adhesion,0.8\n"
    )
    summary = tmp_path / "species_adhesion_summary.csv"
    summary.write_text(
        _SUMMARY_HEADER + "LOC1,ASM1,Sp one,Ascomycota,Pezizomycotina,Sordariomycetes,Hypocreales,"
        "Nectriaceae,Fusarium,10,2,0.2,0.85,0.85\n"
    )

    adhesion_rows, background_rows = build_protein_universe(results_dir, input_dir, summary, seed=1)

    assert len(adhesion_rows) == 2
    assert {r["protein_id"] for r in adhesion_rows} == {"LOC1_000001-T1", "LOC1_000002-T1"}
    assert adhesion_rows[0]["phylum"] == "Ascomycota"
    assert len(background_rows) == 2
    assert all(r["protein_id"] not in {"LOC1_000001-T1", "LOC1_000002-T1"} for r in background_rows)


def test_build_protein_universe_skips_unmatched_species(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_two.proteins.fa.fai").write_text("LOC2_000001-T1\t100\t10\t60\t61\n")
    (results_dir / "Sp_two.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\nLOC2_000001-T1,Adhesion,0.9\n"
    )
    summary = tmp_path / "species_adhesion_summary.csv"
    summary.write_text(_SUMMARY_HEADER)  # no rows — LOC2 not present

    adhesion_rows, background_rows = build_protein_universe(results_dir, input_dir, summary, seed=1)
    assert adhesion_rows == []
    assert background_rows == []
