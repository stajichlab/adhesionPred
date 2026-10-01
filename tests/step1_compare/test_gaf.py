import gzip

import gaf
import pytest
from gaf_fixture import EXT, OBSOLETE, WALL, write_golden_gaf


def test_gaf_extract_golden(tmp_path, mini_ontology):
    path = write_golden_gaf(tmp_path / "golden.gaf")
    out = gaf.filter_gaf(path, mini_ontology)
    assert out.primary_db == "SGD"
    assert dict(out.dropped) == {
        "other_db": 5,
        "object_type": 2,
        "other_aspect": 6,
        "not_qualifier": 1,
        "obsolete_term": 1,
    }
    assert len(out.cc_rows) == 35
    assert len(out.all_aspect_triples) == 43
    genes = {r.gene_id for r in out.cc_rows}
    assert {"CPX-1", "R01", "P99999"}.isdisjoint(genes)
    assert len(genes) == 20
    assert out.unknown_term_rows == 0


def test_taxon_filter_accepts_ncbitaxon_prefix(tmp_path, mini_ontology):
    path = write_golden_gaf(tmp_path / "golden.gaf")
    out = gaf.filter_gaf(path, mini_ontology, taxon_id="559292")
    genes = {r.gene_id for r in out.cc_rows}
    assert "G21" not in genes and "G22" in genes
    assert out.dropped["taxon"] == 2
    # The taxon filter runs before triples are recorded: the 2 G21 triples are not counted.
    assert len(out.all_aspect_triples) == 41


def test_gzip_detected_by_magic_not_suffix(tmp_path, mini_ontology):
    plain = write_golden_gaf(tmp_path / "golden.gaf")
    packed = tmp_path / "golden_no_suffix"
    packed.write_bytes(gzip.compress(plain.read_bytes()))
    assert len(gaf.filter_gaf(packed, mini_ontology).cc_rows) == 35
    assert gaf.header_value(packed, "date-generated") == "2026-05-21T09:00"


def test_short_row_and_empty_file_raise(tmp_path, mini_ontology):
    bad = tmp_path / "bad.gaf"
    bad.write_text("!gaf-version: 2.2\nSGD\tS1\tX\n")
    with pytest.raises(gaf.GafFormatError, match="3 columns"):
        gaf.filter_gaf(bad, mini_ontology)
    empty = tmp_path / "empty.gaf"
    empty.write_text("!gaf-version: 2.2\n")
    with pytest.raises(gaf.GafFormatError, match="no data rows"):
        gaf.filter_gaf(empty, mini_ontology)


def test_evidence_code_sets_match_spec():
    assert gaf.HOMOLOGY_CODES == {"IBA", "IBD", "IKR", "IRD", "ISS", "ISO", "ISA", "ISM", "RCA"}
    assert "IEA" not in gaf.EXPERIMENTAL_CODES | gaf.HOMOLOGY_CODES


# ---- final review item 8: filter order and edge cases (rows built here, not in ROWS) ----


def _gaf(path, rows):
    """rows: (db, gene, qualifier, term, aspect, object_type, taxon); evidence is IDA."""
    lines = ["!gaf-version: 2.2"]
    for db, gene, qual, term, aspect, otype, taxon in rows:
        cols = [db, gene, gene, qual, term, "PMID:1", "IDA", "", aspect, "", "", otype, taxon]
        lines.append("\t".join(cols))
    path.write_text("\n".join(lines) + "\n")
    return path


T = "taxon:559292"
OK = [("SGD", f"S{i}", "located_in", WALL, "C", "protein", T) for i in range(3)]


def test_non_primary_db_is_checked_before_object_type(tmp_path, mini_ontology):
    rows = OK + [("ComplexPortal", "CPX-9", "part_of", WALL, "C", "protein_complex", T)]
    out = gaf.filter_gaf(_gaf(tmp_path / "a.gaf", rows), mini_ontology)
    assert dict(out.dropped) == {"other_db": 1}


def test_aspect_is_checked_before_not_and_obsolete_and_not_before_obsolete(tmp_path, mini_ontology):
    rows = OK + [
        ("SGD", "N1", "NOT|enables", WALL, "F", "protein", T),  # NOT and not C
        ("SGD", "O1", "involved_in", OBSOLETE, "P", "protein", T),  # obsolete and not C
        ("SGD", "B1", "NOT|located_in", OBSOLETE, "C", "protein", T),  # NOT and obsolete
    ]
    out = gaf.filter_gaf(_gaf(tmp_path / "b.gaf", rows), mini_ontology)
    assert dict(out.dropped) == {"other_aspect": 2, "not_qualifier": 1}
    assert len(out.cc_rows) == 3


def test_taxon_filtered_rows_give_no_triples(tmp_path, mini_ontology):
    rows = OK + [("SGD", "X1", "enables", EXT, "F", "protein", "taxon:4932")]
    path = _gaf(tmp_path / "c.gaf", rows)
    assert len(gaf.filter_gaf(path, mini_ontology).all_aspect_triples) == 4
    out = gaf.filter_gaf(path, mini_ontology, taxon_id="559292")
    assert out.dropped["taxon"] == 1
    assert len(out.all_aspect_triples) == 3


def test_taxon_match_is_exact_not_a_substring(tmp_path, mini_ontology):
    assert not gaf.taxon_matches("taxon:2375610", "237561")
    assert not gaf.taxon_matches("taxon:1237561", "237561")
    assert gaf.taxon_matches("taxon:237561|taxon:5476", "237561")
    rows = [("CGD", f"C{i}", "located_in", WALL, "C", "protein", "taxon:237561") for i in range(2)]
    rows.append(("CGD", "C9", "located_in", WALL, "C", "protein", "taxon:2375610"))
    out = gaf.filter_gaf(_gaf(tmp_path / "d.gaf", rows), mini_ontology, taxon_id="237561")
    assert out.dropped["taxon"] == 1
    assert {r.gene_id for r in out.cc_rows} == {"C0", "C1"}


def test_header_value_is_empty_for_a_missing_key(tmp_path):
    path = write_golden_gaf(tmp_path / "golden.gaf")
    assert gaf.header_value(path, "no-such-key") == ""


def test_short_row_and_empty_file_errors_name_the_file(tmp_path, mini_ontology):
    bad = tmp_path / "short_one.gaf"
    bad.write_text("!gaf-version: 2.2\nSGD\tS1\tX\n")
    with pytest.raises(gaf.GafFormatError) as info:
        gaf.filter_gaf(bad, mini_ontology)
    assert str(bad) in str(info.value)
    empty = tmp_path / "empty_one.gaf"
    empty.write_text("!gaf-version: 2.2\n")
    with pytest.raises(gaf.GafFormatError) as info:
        gaf.filter_gaf(empty, mini_ontology)
    assert str(empty) in str(info.value)


def test_experimental_codes_and_protein_types_are_exact():
    assert gaf.EXPERIMENTAL_CODES == {
        "EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP",
    }  # fmt: skip
    assert gaf.PROTEIN_TYPES == {"protein", "gene_product"}


def test_term_unknown_to_the_ontology_is_kept_and_counted(tmp_path, mini_ontology):
    rows = OK + [("SGD", "U1", "located_in", "GO:9999999", "C", "protein", T)]
    out = gaf.filter_gaf(_gaf(tmp_path / "e.gaf", rows), mini_ontology)
    assert out.unknown_term_rows == 1
    assert "U1" in {r.gene_id for r in out.cc_rows}
