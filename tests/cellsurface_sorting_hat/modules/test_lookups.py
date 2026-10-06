import gzip

import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.lookups import (
    antigen_rows,
    cys_rows,
    expression_rows,
    gene_of_ranking_id,
    load_protein_map,
    ranking_by_gene,
    read_table,
    tm_rows,
)

RANKING = (
    "protein\trank\tpercentile\tantigenicity\tspecificity\tprevalence\tmax_fungal_crossreact_pid\n"
    "CIMG_04613-t26_1-p1\t651\t7.12\t2.5\t0.0\t0.9201\t0.0\n"
    "CIMG_04613-t26_2-p1\t900\t9.9\t2.0\t0.0\t0.9201\t0.0\n"
    "CIMG_09560-t26_1-p1\t986\t10.79\t0.1\t1.0\t0.9877\t69.4\n"
)
PMAP = "protein_id\tgene_id\tproduct\tlength\nXP_1\tCIMG_04613\tp\t324\nXP_2\tCIMG_09560\tp\t100\nXP_9\tCIMG_99999\tp\t10\n"


def prot(pid, state="ok"):
    return Protein(pid, "MKT", "x" * 64, state, "", 0.0)


def test_gene_of_ranking_id():
    assert gene_of_ranking_id("CIMG_04613-t26_1-p1") == "CIMG_04613"


def test_ranking_keeps_the_best_row_per_gene_and_counts_genes_with_several(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(RANKING)
    by_gene, several = ranking_by_gene(path)
    assert by_gene["CIMG_04613"]["rank"] == "651" and several == 1


def test_antigen_states_and_fields(tmp_path):
    (tmp_path / "r.tsv").write_text(RANKING)
    (tmp_path / "m.tsv").write_text(PMAP)
    by_gene, _ = ranking_by_gene(tmp_path / "r.tsv")
    pmap = load_protein_map(tmp_path / "m.tsv")
    taxa = {
        "XP_1": 246410,
        "XP_2": 246410,
        "XP_9": 246410,
        "XP_X": 246410,
        "AF1": 746128,
        "BAD": 246410,
    }
    rows = antigen_rows(
        [prot(p, "na_invalid" if p == "BAD" else "ok") for p in taxa], taxa, pmap, by_gene, {246410}
    )
    got = {r["id"]: r for r in rows}
    assert (
        got["XP_1"]["state"] == "ok"
        and got["XP_1"]["percentile"] == "7.12"
        and got["XP_1"]["max_crossreact"] == "0.0"
    )
    assert got["XP_2"]["percentile"] == "10.79"
    assert got["XP_9"]["state"] == "not_in_reference"  # gene not in the ranking
    assert got["XP_X"]["state"] == "not_in_reference"  # protein not in the map
    assert got["AF1"]["state"] == "not_applicable"  # another taxon
    assert got["BAD"]["state"] == "na_invalid"


def test_read_table_counts_repeated_keys(tmp_path):
    path = tmp_path / "t.tsv.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("k\tv\na\t1\na\t2\nb\t3\n")
    rows, repeated = read_table(path, "k")
    assert rows["a"]["v"] == "1" and repeated == 1
    with pytest.raises(ValueError, match="missing column"):
        read_table(path, "nope")


def test_cys_expression_and_tm_rows():
    taxa = {"A": 246410, "B": 246410, "C": 746128}
    cys = cys_rows(
        [prot(p) for p in "ABC"],
        taxa,
        {"A": {"tier": "cys_rich_sp_unassigned", "cys_frac": "0.1"}},
        {246410},
    )
    assert [(r["id"], r["state"]) for r in cys] == [
        ("A", "ok"),
        ("B", "not_in_reference"),
        ("C", "not_applicable"),
    ]
    assert cys[0]["tier"] == "cys_rich_sp_unassigned"
    expr = expression_rows(
        [prot("A"), prot("C")],
        taxa,
        {"A": "CIMG_1"},
        {"CIMG_1": {"log2fc_48h": "9.09", "padj_48h": "1e-5", "log2fc_8d": "8.0"}},
        {246410},
    )
    assert expr[0]["log2fc"] == "9.09" and expr[1]["state"] == "not_applicable"
    tm = tm_rows([prot("A"), prot("B")], {"A": {"pred_hel": "7", "topology": "o10-32i"}})
    assert tm[0]["n_tm"] == "7" and tm[1]["state"] == "error"


def test_tm_rows_count_helices_after_the_signal_peptide_window():
    from cellsurface_sorting_hat.modules.lookups import tm_rows

    tmhmm = {
        "SP_ONLY": {"pred_hel": "1", "topology": "o10-32i"},
        "RECEPTOR": {
            "pred_hel": "7",
            "topology": "i40-62o70-92i100-122o130-152i160-182o190-212i220-242o",
        },
        "MIXED": {"pred_hel": "2", "topology": "o7-25i40-62o"},
        "NONE": {"pred_hel": "0", "topology": "o"},
    }
    got = {r["id"]: (r["n_tm"], r["n_tm_mature"]) for r in tm_rows([prot(k) for k in tmhmm], tmhmm)}
    assert got == {"SP_ONLY": ("1", 0), "RECEPTOR": ("7", 7), "MIXED": ("2", 1), "NONE": ("0", 0)}


# ---- hardening tests (beyond the brief) ----
def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_read_table_short_row_names_path_and_line(tmp_path):
    path = write(tmp_path, "t.tsv", "k\tv\na\t1\nb\n")
    with pytest.raises(ValueError, match=r"t\.tsv:3"):
        read_table(path, "k")


def test_read_table_long_row_is_refused(tmp_path):
    path = write(tmp_path, "t.tsv", "k\tv\na\t1\t2\n")
    with pytest.raises(ValueError, match=r"t\.tsv:2"):
        read_table(path, "k")


@pytest.mark.parametrize("bad", ["nan", "inf", "-inf", "x", "-1"])
def test_read_table_refuses_bad_nonnegative_numbers(tmp_path, bad):
    path = write(tmp_path, "t.tsv", f"k\tv\na\t1\nb\t{bad}\n")
    with pytest.raises(ValueError, match=r"t\.tsv:3.*'v'"):
        read_table(path, "k", numeric={"v": "nonneg"})


def test_read_table_nonneg_accepts_zero_and_finite_accepts_negative(tmp_path):
    path = write(tmp_path, "t.tsv", "k\tv\tw\na\t0\t-2.5\n")
    rows, _ = read_table(path, "k", numeric={"v": "nonneg", "w": "finite"})
    assert rows["a"]["w"] == "-2.5"


@pytest.mark.parametrize("bad", ["nan", "inf", "x"])
def test_read_table_finite_refuses_non_numbers(tmp_path, bad):
    path = write(tmp_path, "t.tsv", f"k\tw\na\t{bad}\n")
    with pytest.raises(ValueError, match="'w'"):
        read_table(path, "k", numeric={"w": "finite"})


def test_read_table_blank_allowed_only_for_or_blank_modes(tmp_path):
    path = write(tmp_path, "t.tsv", "k\tv\tw\na\t\t\n")
    rows, _ = read_table(path, "k", numeric={"v": "nonneg_or_blank", "w": "finite_or_blank"})
    assert rows["a"]["v"] == ""
    with pytest.raises(ValueError, match="'v'"):
        read_table(path, "k", numeric={"v": "nonneg"})


def test_read_table_required_columns(tmp_path):
    path = write(tmp_path, "t.tsv", "k\tv\na\t1\n")
    with pytest.raises(ValueError, match="missing column.*'zz'"):
        read_table(path, "k", required=("v", "zz"))


def test_protein_map_requires_gene_id_column(tmp_path):
    path = write(tmp_path, "m.tsv", "protein_id\tproduct\nXP_1\tp\n")
    with pytest.raises(ValueError, match="gene_id"):
        load_protein_map(path)


@pytest.mark.parametrize(
    "row",
    [
        "X-t1_1-p1\tnan\t1\t1\t0\t0.5\t0\n",
        "X-t1_1-p1\t5\tinf\t1\t0\t0.5\t0\n",
        "X-t1_1-p1\t5\t1\t-1\t0\t0.5\t0\n",
        "X-t1_1-p1\t5\t1\t1\t0\t0.5\t-3\n",
        "X-t1_1-p1\tfive\t1\t1\t0\t0.5\t0\n",
        "X-t1_1-p1\t5\t1\t1\n",
    ],
)
def test_ranking_refuses_bad_values_with_path_and_line(tmp_path, row):
    path = write(
        tmp_path,
        "r.tsv",
        "protein\trank\tpercentile\tantigenicity\tspecificity\tprevalence\tmax_fungal_crossreact_pid\n"
        + row,
    )
    with pytest.raises(ValueError, match=r"r\.tsv:2"):
        ranking_by_gene(path)


def test_ranking_tie_keeps_first_row_and_lowest_rank_wins(tmp_path):
    head = "protein\trank\tpercentile\tantigenicity\tspecificity\tprevalence\tmax_fungal_crossreact_pid\n"
    path = write(
        tmp_path,
        "r.tsv",
        head
        + "G-t1_1-p1\t5\t1\t1\t0\t0\t0\n"
        + "G-t1_2-p1\t5\t2\t1\t0\t0\t0\n"
        + "G-t1_3-p1\t4\t3\t1\t0\t0\t0\n",
    )
    by_gene, several = ranking_by_gene(path)
    assert by_gene["G"]["percentile"] == "3" and several == 1
    path2 = write(
        tmp_path,
        "r2.tsv",
        head + "G-t1_1-p1\t5\t1\t1\t0\t0\t0\n" + "G-t1_2-p1\t5\t2\t1\t0\t0\t0\n",
    )
    assert ranking_by_gene(path2)[0]["G"]["percentile"] == "1"


def test_ranking_rank_must_be_an_integer(tmp_path):
    path = write(
        tmp_path,
        "r.tsv",
        "protein\trank\tpercentile\tantigenicity\tspecificity\tprevalence\tmax_fungal_crossreact_pid\n"
        "G-t1_1-p1\t5.5\t1\t1\t0\t0\t0\n",
    )
    with pytest.raises(ValueError, match=r"r\.tsv:2.*rank"):
        ranking_by_gene(path)


def test_antigen_rows_not_applicable_is_exact_match_on_taxon_id():
    pmap = {"XP_1": "G"}
    by_gene = {
        "G": {
            "rank": "1",
            "percentile": "1",
            "antigenicity": "1",
            "specificity": "1",
            "prevalence": "1",
            "max_fungal_crossreact_pid": "0",
        }
    }
    taxa = {"XP_1": 246411}  # a neighbour of 246410 is not the same taxon
    got = antigen_rows([prot("XP_1")], taxa, pmap, by_gene, {246410})
    assert got == [{"id": "XP_1", "state": "not_applicable"}]
    # a taxon of the same RS strain, passed in the set, is applicable
    assert antigen_rows([prot("XP_1")], taxa, pmap, by_gene, {246410, 246411})[0]["state"] == "ok"


def test_applicability_is_checked_before_the_table_lookup_for_every_lookup():
    taxa = {"A": 5}
    assert (
        cys_rows([prot("A")], taxa, {"A": {"tier": "t", "cys_frac": "0.1"}}, {1})[0]["state"]
        == "not_applicable"
    )
    assert (
        expression_rows([prot("A")], taxa, {"A": "G"}, {"G": {}}, {1})[0]["state"]
        == "not_applicable"
    )
    assert (
        antigen_rows([prot("A")], taxa, {"A": "G"}, {"G": {}}, {1})[0]["state"] == "not_applicable"
    )


def test_empty_applicable_taxa_makes_every_row_not_applicable():
    taxa = {"A": 246410}
    assert cys_rows([prot("A")], taxa, {"A": {"tier": "t", "cys_frac": "0"}}, set())[0][
        "state"
    ] == ("not_applicable")


def test_cimg_id_falls_back_to_gene_of_ranking_id():
    by_gene = {
        "CIMG_1": {
            "rank": "3",
            "percentile": "1",
            "antigenicity": "1",
            "specificity": "1",
            "prevalence": "1",
            "max_fungal_crossreact_pid": "0",
        }
    }
    taxa = {"CIMG_1-t26_1-p1": 246410}
    got = antigen_rows([prot("CIMG_1-t26_1-p1")], taxa, {}, by_gene, {246410})
    assert got[0]["state"] == "ok" and got[0]["rank"] == "3"
    assert got[0]["idmap_method"] == "gene_best_transcript"


def test_expression_reads_8d_and_padj_and_missing_gene_is_not_in_reference():
    taxa = {"A": 1, "B": 1}
    got = expression_rows(
        [prot("A"), prot("B")],
        taxa,
        {"A": "G"},
        {"G": {"log2fc_48h": "-0.5", "padj_48h": "0.32", "log2fc_8d": "1.7"}},
        {1},
    )
    assert (got[0]["log2fc"], got[0]["padj"], got[0]["log2fc_8d"]) == ("-0.5", "0.32", "1.7")
    assert got[1]["state"] == "not_in_reference"


def test_invalid_protein_is_na_invalid_even_when_taxon_not_applicable():
    taxa = {"A": 9}
    for rows in (
        cys_rows([prot("A", "na_invalid")], taxa, {}, {1}),
        expression_rows([prot("A", "na_invalid")], taxa, {}, {}, {1}),
        antigen_rows([prot("A", "na_invalid")], taxa, {}, {}, {1}),
        tm_rows([prot("A", "na_invalid")], {}),
    ):
        assert rows == [{"id": "A", "state": "na_invalid"}]


@pytest.mark.parametrize("fn", ["cys", "expr", "antigen"])
def test_taxa_must_cover_the_fasta(fn):
    def run(taxa):
        ps = [prot("A")]
        if fn == "cys":
            return cys_rows(ps, taxa, {}, {1})
        if fn == "expr":
            return expression_rows(ps, taxa, {}, {}, {1})
        return antigen_rows(ps, taxa, {}, {}, {1})

    with pytest.raises(ValueError, match=r"no taxon.*'A'"):
        run({})


def test_tm_refuses_ids_not_in_the_fasta():
    with pytest.raises(ValueError, match=r"2 .*not in the FASTA.*'X1'"):
        tm_rows([prot("A")], {"A": {"pred_hel": "0", "topology": "o"}, "X1": {}, "X2": {}})


@pytest.mark.parametrize("bad", ["nan", "-1", "2.5", "", "x"])
def test_tm_refuses_bad_helix_count(bad):
    with pytest.raises(ValueError, match=r"pred_hel.*'A'"):
        tm_rows([prot("A")], {"A": {"pred_hel": bad, "topology": "o"}})


def test_tm_signal_window_boundary_is_strictly_after_35():
    tmhmm = {
        "AT35": {"pred_hel": "1", "topology": "o35-57i"},
        "AT36": {"pred_hel": "1", "topology": "o36-58i"},
        "AT1": {"pred_hel": "1", "topology": "o1-20i"},
    }
    got = {r["id"]: r["n_tm_mature"] for r in tm_rows([prot(k) for k in tmhmm], tmhmm)}
    assert got == {"AT35": 0, "AT36": 1, "AT1": 0}


def test_tm_row_for_protein_missing_from_tmhmm_is_error_and_topology_kept():
    got = tm_rows([prot("A")], {"A": {"pred_hel": "1", "topology": "o40-62i"}})
    assert got[0]["topology"] == "o40-62i"
