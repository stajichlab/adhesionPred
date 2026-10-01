import gzip

import gaf
import pytest
from gaf_fixture import write_golden_gaf


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
