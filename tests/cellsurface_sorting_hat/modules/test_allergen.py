import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.allergen import (
    allergen_name,
    allergen_rows,
    best_hit_other_species,
    build_allergen_fasta,
    check_blast_ids,
    lso_report,
    parse_blast,
    parse_blast_hits,
    read_meta,
    recall_by_rule,
    species_code,
)


def prot(pid, state="ok"):
    return Protein(pid, "MKT", "x" * 64, state, "", 0.0)


SEQ1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"  # 33 residues
SEQ2 = "MNLLPQWERTYIPASDFGHKLCVNMQRSTWYAAA"  # 34 residues
IUIS = (
    "AllergenID\tIsoName\tName\tSequence\n"
    f"11\tAsp f 1.0101\tAsp f 1\t{SEQ1[:20]} {SEQ1[20:]}\n"
    f"12\tAsp f 2.0101\tAsp f 2\t{SEQ2}\n"
    "13\tNoSeq\tNo seq\t\n"
    "14\tEpi p 1.0101\tEpi p 1\tN-TERMINAL PEPTIDE: ADGIVAVELDTY >INTERNAL PEPTIDE: RGSFXK\n"
    "15\tPep 1.0101\tPep 1\tADGIVAVELDTY\n"
)
ALLERGENS = (
    "AllergenID\tName\tSpecies\tTaxOrder\tExposure\tAllergenicity\n"
    "11\tAsp f 1\tAspergillus fumigatus\tEurotiales\tAirway\tIgE binding in 75% of 40 sera\n"
    "12\tAsp f 2\tAspergillus fumigatus\tEurotiales\tAirway\t\n"
)


def test_allergen_fasta_is_clean_and_has_a_meta_table(tmp_path):
    src, al, out = tmp_path / "i.tsv", tmp_path / "a.tsv", tmp_path / "a.faa"
    src.write_text(IUIS)
    al.write_text(ALLERGENS)
    n, skipped = build_allergen_fasta(src, out, al)
    assert n == 2
    assert (
        out.read_text() == f">Asp_f_1.0101|11\n{SEQ1}\n>Asp_f_2.0101|12\n{SEQ2}\n"
    )  # spaces removed
    assert skipped == [
        ("Epi_p_1.0101|14", "not a protein sequence"),
        ("Pep_1.0101|15", "peptide fragment of 12 residues"),
    ]
    meta = read_meta(str(out) + ".meta.tsv")
    assert (
        meta["Asp_f_1.0101|11"]["exposure"] == "Airway"
        and "75% of 40 sera" in meta["Asp_f_1.0101|11"]["evidence"]
    )
    assert (
        meta["Asp_f_2.0101|12"]["evidence"] == ""
        and meta["Asp_f_2.0101|12"]["species"] == "Aspergillus fumigatus"
    )


def test_duplicate_allergen_ids_are_refused(tmp_path):
    src = tmp_path / "i.tsv"
    src.write_text(
        f"AllergenID\tIsoName\tName\tSequence\n1\tA 1\tA 1\t{SEQ1}\n1\tA 1\tA 1\t{SEQ2}\n"
    )
    with pytest.raises(ValueError, match="duplicate allergen ID"):
        build_allergen_fasta(src, tmp_path / "a.faa")


BLAST = (
    "P1\tAsp_f_1.0101|11\t99.0\t100\t120\t125\t200\t1e-50\n"
    "P1\tAsp_f_2.0101|12\t45.0\t90\t120\t300\t60\t1e-5\n"
    "P2\tAsp_f_2.0101|12\t38.0\t85\t200\t300\t50\t1e-3\n"
)


def test_blast_best_hit_is_by_bit_score_and_fields_are_derived(tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text(BLAST)
    best = parse_blast(path)
    meta = {
        "Asp_f_1.0101|11": {
            "species": "Aspergillus fumigatus",
            "exposure": "Airway",
            "evidence": "75% of 40",
        }
    }
    rows = {
        r["id"]: r
        for r in allergen_rows(
            [prot("P1"), prot("P2"), prot("P3"), prot("BAD", "na_invalid")], best, meta
        )
    }
    assert rows["P1"]["identity"] == "99.00" and rows["P1"]["allergen_name"] == "Asp_f_1.0101"
    assert rows["P1"]["coverage"] == "80.0" and rows["P1"]["exposure"] == "Airway"
    assert (rows["P2"]["identity"], rows["P2"]["aligned_length"], rows["P2"]["coverage"]) == (
        "38.00",
        "85",
        "28.3",
    )
    assert (rows["P3"]["identity"], rows["P3"]["aligned_length"], rows["P3"]["allergen_name"]) == (
        "0",
        "0",
        "",
    )
    assert rows["BAD"]["state"] == "na_invalid"


def test_blast_line_with_the_wrong_field_count_is_refused(tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text("P1\tS\t99.0\n")
    with pytest.raises(ValueError, match="expected 8 fields"):
        parse_blast(path)


def test_a_blast_table_from_another_proteome_is_refused():
    best = {"X9": {"subject": "S", "identity": 90.0, "length": 100, "qlen": 100, "slen": 100}}
    with pytest.raises(ValueError, match="not in the FASTA"):
        check_blast_ids(best, ["P1", "P2"])
    check_blast_ids({"P1": best["X9"]}, ["P1", "P2"])  # no error


def test_names_and_species_codes():
    assert allergen_name("Asp_f_1.0101|11") == "Asp_f_1.0101"
    assert species_code("Asp_f_1.0101|11") == "Asp_f"
    assert species_code("Cand_a_3.0101|5") == "Cand_a"


def h(q, s, ident, length, qlen=100, slen=100, bits=100.0):
    return {
        "query": q,
        "subject": s,
        "identity": ident,
        "length": length,
        "qlen": qlen,
        "slen": slen,
        "bitscore": bits,
    }


IDS = ["Asp_f_1.0101|1", "Asp_f_2.0101|2", "Asp_n_1.0101|3", "Alt_a_1.0101|4"]


def test_other_species_hits_exclude_the_same_species():
    hits = [
        h(IDS[0], IDS[1], 99, 100, bits=300),
        h(IDS[0], IDS[2], 60, 90, bits=100),
        h(IDS[0], IDS[3], 40, 85, bits=50),
    ]
    best = best_hit_other_species(hits, IDS)
    assert best[IDS[0]]["subject"] == IDS[2]  # the 99% hit is another allergen of the same species


def test_recall_by_rule_uses_the_rules_of_the_engine_and_the_fasta_denominator():
    best = {IDS[0]: h(IDS[0], IDS[2], 60, 90), IDS[1]: h(IDS[1], IDS[3], 72, 85)}
    rules = [("similarity", 35.0, 80, 0.0), ("homolog", 70.0, 0, 80.0)]
    got = {r["rule"]: (r["recovered"], r["n"]) for r in recall_by_rule(IDS, best, rules)}
    assert got == {
        "similarity": (2, 4),
        "homolog": (1, 4),
    }  # IDS[2], IDS[3] have no hit: counted as missed


def test_lso_report_refuses_a_sequence_with_no_blast_line(tmp_path):
    path = tmp_path / "b.tsv"
    lines = [f"{i}\t{i}\t100.0\t100\t100\t100\t500\t0" for i in IDS[:3]]  # IDS[3] absent
    path.write_text("\n".join(lines) + "\n")
    assert len(parse_blast_hits(path)) == 3
    with pytest.raises(ValueError, match="no BLAST line"):
        lso_report(path, IDS)
    lines.append(f"{IDS[3]}\t{IDS[3]}\t100.0\t100\t100\t100\t500\t0")
    path.write_text("\n".join(lines) + "\n")
    rep = lso_report(path, IDS)
    assert (rep["n_sequences"], rep["n_species"]) == (4, 3)
    assert [r["recovered"] for r in rep["recall"]] == [0, 0]


# ---- additions beyond the brief (review lessons from Tasks 1 to 3) ----


def test_allergen_rows_refuses_a_blast_table_of_another_proteome():
    best = {"X9": {"subject": "S", "identity": 90.0, "length": 100, "qlen": 100, "slen": 100}}
    with pytest.raises(ValueError, match=r"1 BLAST query ID.*'X9'"):
        allergen_rows([prot("P1")], best)


def test_lso_report_refuses_ids_that_are_not_in_the_fasta(tmp_path):
    path = tmp_path / "b.tsv"
    lines = [f"{i}\t{i}\t100.0\t100\t100\t100\t500\t0" for i in IDS]
    path.write_text("\n".join(lines) + "\nZZ\tZZ\t100.0\t100\t100\t100\t500\t0\n")
    with pytest.raises(ValueError, match="not in the FASTA"):
        lso_report(path, IDS)
    lines.append(f"{IDS[0]}\tZZ\t50.0\t100\t100\t100\t50\t0")  # unknown subject
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="not in the FASTA"):
        lso_report(path, IDS)


@pytest.mark.parametrize(
    "line",
    [
        "P1\tS\tabc\t100\t100\t100\t50\t0",  # not a number
        "P1\tS\tnan\t100\t100\t100\t50\t0",
        "P1\tS\tinf\t100\t100\t100\t50\t0",
        "P1\tS\t-1.0\t100\t100\t100\t50\t0",  # negative identity
        "P1\tS\t101.0\t100\t100\t100\t50\t0",  # identity above 100
        "P1\tS\t90.0\t-5\t100\t100\t50\t0",  # negative length
        "P1\tS\t90.0\t1.5\t100\t100\t50\t0",  # length is an integer
        "P1\tS\t90.0\t100\t0\t100\t50\t0",  # zero query length
        "P1\tS\t90.0\t100\t100\t0\t50\t0",  # zero subject length
        "P1\tS\t90.0\t100\t100\t100\tnan\t0",
        "P1\tS\t90.0\t100\t100\t100\t50\tinf",
    ],
)
def test_bad_blast_numbers_name_path_and_line(tmp_path, line):
    path = tmp_path / "b.tsv"
    good = "P0\tS\t90.0\t100\t100\t100\t50\t0\n"
    path.write_text(good + line + "\n")
    with pytest.raises(ValueError, match=r"b\.tsv:2:"):
        parse_blast(path)
    with pytest.raises(ValueError, match=r"b\.tsv:2:"):
        parse_blast_hits(path)


def test_field_count_error_names_the_line_in_both_parsers(tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text("P1\tS\t90.0\t100\t100\t100\t50\t0\nP1\tS\n")
    for fn in (parse_blast, parse_blast_hits):
        with pytest.raises(ValueError, match=r"b\.tsv:2: expected 8 fields"):
            fn(path)


def test_blast_evalue_zero_and_blank_lines_are_accepted(tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text("\nP1\tS\t90.0\t100\t100\t100\t50\t0.0\n\n")
    assert parse_blast(path)["P1"]["evalue"] == 0.0
    assert len(parse_blast_hits(path)) == 1


def test_best_hit_ties_on_bit_score_go_to_higher_identity(tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text(
        "P1\tA\t50.0\t100\t100\t100\t80\t0\nP1\tB\t60.0\t100\t100\t100\t80\t0\n"
        "P1\tC\t99.0\t100\t100\t100\t70\t0\n"
    )
    assert parse_blast(path)["P1"]["subject"] == "B"  # not the first line, not the last
    hits = parse_blast_hits(path)
    ids = ["P1", "A_x|1", "B_y|2"]
    assert best_hit_other_species(hits, ids) == {}  # P1 is not an allergen ID


def test_coverage_is_capped_at_100_and_is_of_the_allergen_length():
    best = {"P1": {"subject": "S_s_1|1", "identity": 50.0, "length": 150, "qlen": 500, "slen": 100}}
    row = allergen_rows([prot("P1")], best)[0]
    assert row["coverage"] == "100.0"  # not 30.0 (query length) and not 150.0 (uncapped)
    assert row["aligned_length"] == "150"


def test_the_best_hit_is_not_the_first_or_the_highest_identity(tmp_path):
    path = tmp_path / "b.tsv"
    path.write_text(
        "P1\tA\t99.0\t100\t100\t100\t50\t0\nP1\tB\t40.0\t100\t100\t100\t90\t0\n"
        "P1\tC\t70.0\t100\t100\t100\t60\t0\n"
    )
    assert parse_blast(path)["P1"]["subject"] == "B"  # highest bit score


def test_a_meta_value_is_not_corrupted_and_tabs_do_not_split_columns(tmp_path):
    src, al, out = tmp_path / "i.tsv", tmp_path / "a.tsv", tmp_path / "a.faa"
    src.write_text(f"AllergenID\tIsoName\tName\tSequence\n11\tAsp f 1.0101\tAsp f 1\t{SEQ1}\n")
    al.write_text(
        "AllergenID\tName\tSpecies\tTaxOrder\tExposure\tAllergenicity\n"
        "11\tAsp f 1\tCinnamomum nanum\tEurotiales\tAirway\tbinds nan  IgE\n"
    )
    build_allergen_fasta(src, out, al)
    m = read_meta(str(out) + ".meta.tsv")["Asp_f_1.0101|11"]
    assert m["species"] == "Cinnamomum nanum"  # "nan" inside a word stays
    assert m["evidence"] == "binds nan IgE"


def test_a_pandas_nan_cell_is_an_empty_value(tmp_path):
    src, al, out = tmp_path / "i.tsv", tmp_path / "a.tsv", tmp_path / "a.faa"
    src.write_text(f"AllergenID\tIsoName\tName\tSequence\n11\tAsp f 1.0101\tAsp f 1\t{SEQ1}\n")
    al.write_text(
        "AllergenID\tName\tSpecies\tTaxOrder\tExposure\tAllergenicity\n"
        "11\tAsp f 1\tAspergillus fumigatus\tnan\tAirway\tnan\n"
    )
    build_allergen_fasta(src, out, al)
    m = read_meta(str(out) + ".meta.tsv")["Asp_f_1.0101|11"]
    assert (m["tax_order"], m["evidence"]) == ("", "")


def test_a_table_without_the_needed_columns_is_refused(tmp_path):
    src, al = tmp_path / "i.tsv", tmp_path / "a.tsv"
    src.write_text("Id\tSeq\n1\t" + SEQ1 + "\n")
    with pytest.raises(ValueError, match=r"i\.tsv.*AllergenID"):
        build_allergen_fasta(src, tmp_path / "a.faa")
    src.write_text(f"AllergenID\tIsoName\tName\tSequence\n1\tA 1\tA 1\t{SEQ1}\n")
    al.write_text("Id\tName\n1\tx\n")
    with pytest.raises(ValueError, match=r"a\.tsv.*AllergenID"):
        build_allergen_fasta(src, tmp_path / "a.faa", al)


def test_fragment_boundary_is_30_residues(tmp_path):
    src, out = tmp_path / "i.tsv", tmp_path / "a.faa"
    s30, s29 = "A" * 30, "A" * 29
    src.write_text(f"AllergenID\tIsoName\tName\tSequence\n1\tX 1\tX 1\t{s30}\n2\tX 2\tX 2\t{s29}\n")
    n, skipped = build_allergen_fasta(src, out)
    assert n == 1 and skipped == [("X_2|2", "peptide fragment of 29 residues")]


def test_similarity_rule_boundaries():
    ids = ["A_a_1|1", "B_b_1|2"]
    sim = [("similarity", 35.0, 80, 0.0)]

    def run(ident, length):
        best = {ids[0]: h(ids[0], ids[1], ident, length)}
        return recall_by_rule(ids, best, sim)[0]["recovered"]

    assert run(35.0, 80) == 1  # both limits are inclusive
    assert run(34.99, 80) == 0
    assert run(35.0, 79) == 0
    assert run(100.0, 79) == 0  # no sliding window, no coverage shortcut


def test_homolog_rule_boundaries():
    ids = ["A_a_1|1", "B_b_1|2"]
    rule = [("homolog", 70.0, 0, 80.0)]

    def run(ident, length, slen=100):
        best = {ids[0]: h(ids[0], ids[1], ident, length, slen=slen)}
        return recall_by_rule(ids, best, rule)[0]["recovered"]

    assert run(70.0, 80) == 1  # identity 70 and coverage 80 are inclusive
    assert run(69.99, 80) == 0
    assert run(70.0, 79) == 0
    assert run(70.0, 40, slen=50) == 1  # coverage uses the allergen length, not 100
    assert run(70.0, 130, slen=100) == 1  # coverage above 100 still passes


def test_recall_ignores_a_hit_of_a_sequence_that_is_not_in_ids():
    best = {"ZZ": h("ZZ", IDS[2], 99, 100)}
    got = recall_by_rule(IDS, best, [("similarity", 35.0, 80, 0.0)])
    assert got[0]["recovered"] == 0 and got[0]["n"] == 4


def test_lso_report_counts_species_by_name_not_by_sequence(tmp_path):
    ids = ["Asp_f_1.0101|1", "Asp_f_2.0101|2", "Alt_a_1.0101|3"]
    lines = [f"{i}\t{i}\t100.0\t100\t100\t100\t500\t0" for i in ids]
    lines.append("Asp_f_1.0101|1\tAlt_a_1.0101|3\t80.0\t100\t100\t100\t200\t0")
    lines.append("Asp_f_1.0101|1\tAsp_f_2.0101|2\t99.0\t100\t100\t100\t300\t0")
    path = tmp_path / "b.tsv"
    path.write_text("\n".join(lines) + "\n")
    rep = lso_report(path, ids)
    assert (rep["n_sequences"], rep["n_species"]) == (3, 2)
    got = {r["rule"]: r["recovered"] for r in rep["recall"]}
    assert got == {"iuis_allergen_similarity": 1, "iuis_allergen_homolog": 1}
