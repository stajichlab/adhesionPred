"""Extra label tests built from small per-test GAF files (the golden fixture is untouched)."""

import gaf
import labels
import pytest
import truth_table

TAXON = "taxon:559292"
WALL, EXT, PM = "GO:0009277", "GO:0005576", "GO:0005886"
CYTO, NUC, MITO, PERIPH = "GO:0005829", "GO:0005634", "GO:0005739", "GO:0071944"
ER = "GO:0005783"


def _build(tmp_path, ontology, annotations):
    """annotations: list of (gene_id, term, evidence). The symbol is gene_id.lower()."""
    lines = ["!gaf-version: 2.2", ""]
    for gene, term, ev in annotations:
        sym = gene.lower()
        cols = ["SGD", gene, sym, "located_in", term, "PMID:1", ev, "", "C", sym]
        cols += [f"{sym}|{gene}-syn", "protein", TAXON, "20260101", "SGD", "", ""]
        lines.append("\t".join(cols))
    path = tmp_path / "m.gaf"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    info = truth_table.SourceInfo("Fix", "Fix", "559292", "Sacc", "train", "m.gaf", "s", "d", "o")
    filtered = gaf.filter_gaf(path, ontology)
    return truth_table.build_truth_rows(filtered, ontology, info)


def _by_gene(rows):
    return {r["gene_id"]: r for r in rows}


def test_homology_only_is_label_change_not_presence_of_homology_code(tmp_path, mini_ontology):
    rows = _by_gene(
        _build(
            tmp_path,
            mini_ontology,
            [
                ("N1", NUC, "IDA"),  # nucleus IDA + cytosol IBA: N-int with or without IBA
                ("N1", CYTO, "IBA"),
                ("W1", WALL, "IDA"),  # wall IDA + cytosol IBA: label changes
                ("W1", CYTO, "IBA"),
            ],
        )
    )
    n1 = rows["N1"]
    assert "IBA" in n1["evidence_codes"].split(",")
    assert (n1["label"], n1["label_no_homology"], n1["homology_only"]) == ("N-int", "N-int", "no")
    assert truth_table.in_direct_stratum(n1, "N-int")
    w1 = rows["W1"]
    assert (w1["label"], w1["label_no_homology"], w1["homology_only"]) == (
        "ambiguous",
        "P-ext",
        "yes",
    )


def test_label_no_homology_uses_no_homology_policy_not_experimental(tmp_path, mini_ontology):
    rows = _by_gene(_build(tmp_path, mini_ontology, [("T1", WALL, "TAS")]))
    t1 = rows["T1"]
    assert (t1["label"], t1["label_no_homology"], t1["label_experimental"]) == (
        "P-ext",
        "P-ext",
        "unlabelled",
    )
    assert t1["homology_only"] == "no"


def test_subset_and_pm_candidate_ignore_iea_terms(tmp_path, mini_ontology):
    rows = _by_gene(
        _build(
            tmp_path,
            mini_ontology,
            [("S1", EXT, "IDA"), ("S1", WALL, "IEA"), ("S1", PM, "IEA")],
        )
    )
    s1 = rows["S1"]
    assert s1["label"] == "P-ext"
    assert s1["subset"] == "extracellular-only"
    assert s1["stratum"] == "extracellular-only"
    assert s1["pm_candidate"] == "no"


def test_htp_only_ignores_iea_internal_rows(tmp_path, mini_ontology):
    rows = _by_gene(
        _build(
            tmp_path,
            mini_ontology,
            [
                ("H1", WALL, "IDA"),  # internal evidence: HDA plus an IEA row
                ("H1", MITO, "HDA"),
                ("H1", CYTO, "IEA"),
                ("H2", WALL, "IDA"),  # internal evidence: IDA and HDA
                ("H2", MITO, "HDA"),
                ("H2", CYTO, "IDA"),
            ],
        )
    )
    assert (rows["H1"]["label"], rows["H1"]["internal_evidence"]) == ("ambiguous", "HDA")
    assert rows["H1"]["internal_evidence_htp_only"] == "yes"
    assert (rows["H2"]["label"], rows["H2"]["internal_evidence"]) == ("ambiguous", "HDA,IDA")
    assert rows["H2"]["internal_evidence_htp_only"] == "no"


@pytest.mark.parametrize("code", ["HDA", "HMP", "HEP", "HGI", "HTP"])
def test_each_high_throughput_code_counts(code):
    assert labels.htp_only({code})


def test_non_high_throughput_code_does_not_count():
    assert not labels.htp_only({"IDA"})


@pytest.mark.parametrize("periph_evidence", ["IEA", "IDA"])
def test_cell_periphery_blocks_n_int(tmp_path, mini_ontology, periph_evidence):
    rows = _by_gene(
        _build(tmp_path, mini_ontology, [("P1", CYTO, "IDA"), ("P1", PERIPH, periph_evidence)])
    )
    assert rows["P1"]["label"] == "unlabelled"


def test_cell_periphery_blocks_n_int_in_classify():
    assert labels.classify(frozenset({CYTO}), frozenset({CYTO, "GO:0071944"})) == "unlabelled"


def test_rows_are_sorted_by_gene_id(tmp_path, mini_ontology):
    rows = _build(
        tmp_path,
        mini_ontology,
        [("Z9", WALL, "IDA"), ("M5", NUC, "IDA"), ("A1", ER, "IDA")],
    )
    assert [r["gene_id"] for r in rows] == ["A1", "M5", "Z9"]


def test_row_content_columns(tmp_path, mini_ontology):
    rows = _by_gene(
        _build(
            tmp_path,
            mini_ontology,
            [("C1", PM, "IDA"), ("C1", NUC, "IBA"), ("C1", ER, "IEA")],
        )
    )
    c1 = rows["C1"]
    assert c1["label"] == "N-sec"
    assert c1["tier"] == "T-a"
    assert c1["symbol"] == "c1"
    assert c1["synonym1"] == "c1"
    assert c1["evidence_codes"] == "IBA,IDA,IEA"
    assert c1["secretory_evidence"] == "IDA"
    assert c1["internal_evidence"] == "IBA"
    assert c1["surface_evidence"] == ""
    assert (c1["source_id"], c1["role"], c1["obo_sha256"]) == ("Fix", "train", "o")
