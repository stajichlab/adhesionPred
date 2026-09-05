"""Tests for analysis/kingdom_survey/join.py."""

from join import (
    build_species_table,
    count_fai_lines,
    count_result_rows,
    find_missing_results,
    first_id_from_fai,
    first_id_from_result_csv,
    load_samples_taxonomy,
    locustag_from_id,
    probability_stats,
)


def test_locustag_from_id():
    assert locustag_from_id("F07B100A_000481-T1") == "F07B100A"


def test_first_id_from_fai(tmp_path):
    fai = tmp_path / "sp.proteins.fa.fai"
    fai.write_text("LOC1_000001-T1\t100\t10\t60\t61\nLOC1_000002-T1\t100\t150\t60\t61\n")
    assert first_id_from_fai(fai) == "LOC1_000001-T1"


def test_count_fai_lines(tmp_path):
    fai = tmp_path / "sp.proteins.fa.fai"
    fai.write_text("a\tb\tc\td\te\n" * 5)
    assert count_fai_lines(fai) == 5


def test_first_id_from_result_csv(tmp_path):
    result = tmp_path / "sp.adhesion_predict.csv"
    result.write_text("id,prediction,probability_adhesion\nLOC1_000002-T1,Adhesion,0.9\n")
    assert first_id_from_result_csv(result) == "LOC1_000002-T1"


def test_count_result_rows(tmp_path):
    result = tmp_path / "sp.adhesion_predict.csv"
    result.write_text(
        "id,prediction,probability_adhesion\n"
        "LOC1_000002-T1,Adhesion,0.9\n"
        "LOC1_000003-T1,Adhesion,0.6\n"
    )
    assert count_result_rows(result) == 2


def test_probability_stats(tmp_path):
    result = tmp_path / "sp.adhesion_predict.csv"
    result.write_text(
        "id,prediction,probability_adhesion\n"
        "LOC1_000002-T1,Adhesion,0.9\n"
        "LOC1_000003-T1,Adhesion,0.7\n"
    )
    mean, median = probability_stats(result)
    assert mean == 0.8
    assert median == 0.8


def test_load_samples_taxonomy(tmp_path):
    samples = tmp_path / "samples.csv"
    samples.write_text(
        "ASMID,SPECIESIN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,"
        "CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,LOCUSTAG\n"
        "ASM1,Sp one,strain1,PRJ1,1,dikarya,Ascomycota,Pezizomycotina,"
        "Sordariomycetes,,Hypocreales,Nectriaceae,Fusarium,Fusarium one,LOC1\n"
    )
    taxonomy = load_samples_taxonomy(samples)
    assert taxonomy["LOC1"]["phylum"] == "Ascomycota"
    assert taxonomy["LOC1"]["genus"] == "Fusarium"


def test_build_species_table_matched(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_one.proteins.fa.fai").write_text(
        "LOC1_000001-T1\t100\t10\t60\t61\n"
        "LOC1_000002-T1\t100\t10\t60\t61\n"
        "LOC1_000003-T1\t100\t10\t60\t61\n"
        "LOC1_000004-T1\t100\t10\t60\t61\n"
    )
    (results_dir / "Sp_one.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\nLOC1_000002-T1,Adhesion,0.9\n"
    )
    taxonomy = {
        "LOC1": {
            "asmid": "ASM1",
            "phylum": "Ascomycota",
            "subphylum": "Pezizomycotina",
            "class": "Sordariomycetes",
            "order": "Hypocreales",
            "family": "Nectriaceae",
            "genus": "Fusarium",
            "species_name": "Fusarium one",
        }
    }
    matched, unmatched, mismatched = build_species_table(results_dir, input_dir, taxonomy)
    assert len(matched) == 1
    assert unmatched == []
    assert mismatched == []
    row = matched[0]
    assert row["locustag"] == "LOC1"
    assert row["total_proteins"] == 4
    assert row["adhesion_count"] == 1
    assert row["adhesion_fraction"] == 0.25


def test_build_species_table_mismatched_locustag(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_two.proteins.fa.fai").write_text("LOC2_000001-T1\t100\t10\t60\t61\n")
    (results_dir / "Sp_two.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\nXP_999999.1,Adhesion,0.9\n"
    )
    taxonomy = {
        "LOC2": {
            "asmid": "ASM2",
            "phylum": "Ascomycota",
            "subphylum": None,
            "class": None,
            "order": None,
            "family": None,
            "genus": None,
            "species_name": "Sp two",
        }
    }
    matched, unmatched, mismatched = build_species_table(results_dir, input_dir, taxonomy)
    assert matched == []
    assert len(mismatched) == 1
    assert mismatched[0]["fai_locustag"] == "LOC2"
    assert mismatched[0]["result_locustag"] == "XP"


def test_build_species_table_unmatched_locustag(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_three.proteins.fa.fai").write_text("LOC3_000001-T1\t100\t10\t60\t61\n")
    (results_dir / "Sp_three.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\nLOC3_000001-T1,Adhesion,0.9\n"
    )
    matched, unmatched, mismatched = build_species_table(results_dir, input_dir, {})
    assert matched == []
    assert len(unmatched) == 1
    assert "LOC3" in unmatched[0]["reason"]


def test_find_missing_results(tmp_path):
    input_dir = tmp_path / "input"
    results_dir = tmp_path / "results"
    input_dir.mkdir()
    results_dir.mkdir()
    (input_dir / "Sp_one.proteins.fa.fai").write_text("LOC1_000001-T1\t100\t10\t60\t61\n")
    (input_dir / "Sp_two.proteins.fa.fai").write_text("LOC2_000001-T1\t100\t10\t60\t61\n")
    (results_dir / "Sp_one.adhesion_predict.csv").write_text(
        "id,prediction,probability_adhesion\nLOC1_000001-T1,Adhesion,0.9\n"
    )
    missing = find_missing_results(results_dir, input_dir)
    assert len(missing) == 1
    assert missing[0]["stem"] == "Sp_two"
