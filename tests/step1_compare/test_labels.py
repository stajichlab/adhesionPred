import gaf
import labels
import pytest
import truth_table
from gaf_fixture import EXPECTED_LABELS, write_golden_gaf

WALL, EXT, PM, CYTO, NUC = "GO:0005618", "GO:0005576", "GO:0005886", "GO:0005829", "GO:0005634"
ER, MEM = "GO:0012505", "GO:0016020"


@pytest.mark.parametrize(
    ("label_terms", "any_terms", "expected"),
    [
        ({WALL}, {WALL}, "P-ext"),
        ({EXT}, {EXT, CYTO}, "P-ext"),  # internal term at IEA level only
        ({WALL, CYTO}, {WALL, CYTO}, "ambiguous"),
        ({NUC}, {NUC}, "N-int"),
        ({NUC}, {NUC, MEM}, "unlabelled"),  # membrane at any level blocks N-int
        ({NUC}, {NUC, WALL}, "unlabelled"),
        ({PM}, {PM}, "N-sec"),
        ({PM, NUC}, {PM, NUC}, "N-sec"),
        ({ER}, {ER, EXT}, "unlabelled"),  # surface at any level blocks N-sec
        (set(), {WALL, CYTO}, "unlabelled"),  # wall and cytosol at IEA level only
        (set(), set(), "unlabelled"),
    ],
)
def test_label_rules_truth_table(label_terms, any_terms, expected):
    assert labels.classify(frozenset(label_terms), frozenset(any_terms)) == expected


def test_subset_and_pm_candidate():
    assert labels.subset_of("P-ext", {WALL, EXT}) == "wall"
    assert labels.subset_of("P-ext", {EXT}) == "extracellular-only"
    assert labels.subset_of("N-sec", {PM}) == ""
    assert labels.is_pm_candidate("P-ext", {WALL, PM})
    assert not labels.is_pm_candidate("N-sec", {PM})


def test_policies():
    assert labels.policy_accepts("non_iea", "IBA") and not labels.policy_accepts("non_iea", "IEA")
    assert not labels.policy_accepts("no_homology", "ISS") and labels.policy_accepts(
        "no_homology", "TAS"
    )
    assert labels.policy_accepts("experimental", "HDA") and not labels.policy_accepts(
        "experimental", "TAS"
    )
    with pytest.raises(ValueError):
        labels.policy_accepts("all", "IDA")


def _rows(tmp_path, ontology):
    filtered = gaf.filter_gaf(write_golden_gaf(tmp_path / "g.gaf"), ontology)
    info = truth_table.SourceInfo(
        "Fix",
        "Fixture yeast",
        "559292",
        "Saccharomycotina",
        "train",
        "g.gaf",
        "sha",
        "2026-05-21",
        "obo",
    )
    return {r["gene_id"]: r for r in truth_table.build_truth_rows(filtered, ontology, info)}


def test_golden_labels(tmp_path, mini_ontology):
    rows = _rows(tmp_path, mini_ontology)
    assert {g: (r["label"], r["subset"]) for g, r in rows.items()} == EXPECTED_LABELS
    assert [g for g, r in rows.items() if r["pm_candidate"] == "yes"] == ["G03"]


def test_conflicting_evidence_is_stored_not_collapsed(tmp_path, mini_ontology):
    g10 = _rows(tmp_path, mini_ontology)["G10"]
    assert (g10["label"], g10["label_no_homology"], g10["homology_only"]) == (
        "ambiguous",
        "P-ext",
        "yes",
    )
    assert g10["surface_evidence"] == "IBA,IDA"
    assert g10["internal_evidence"] == "IBA"
    g11 = _rows(tmp_path, mini_ontology)["G11"]
    assert (g11["label"], g11["label_no_homology"], g11["label_experimental"]) == (
        "P-ext",
        "unlabelled",
        "unlabelled",
    )


def test_truth_set_has_no_iea(tmp_path, mini_ontology):
    rows = _rows(tmp_path, mini_ontology)
    support = {
        "P-ext": "surface_evidence",
        "N-int": "internal_evidence",
        "N-sec": "secretory_evidence",
    }
    for r in rows.values():
        for column in ("surface_evidence", "internal_evidence", "secretory_evidence"):
            assert "IEA" not in r[column].split(",")
        if r["label"] in support:
            assert r[support[r["label"]]], f"{r['gene_id']} label has no non-IEA support"
    assert rows["G09"]["label"] == "unlabelled"  # wall + cytosol only at IEA level
    assert rows["G15"]["label"] == "unlabelled"  # IEA-only gene


def test_htp_only_rule():
    assert labels.htp_only({"HDA"}) and labels.htp_only({"HDA", "HTP"})
    assert not labels.htp_only({"HDA", "IDA"}) and not labels.htp_only(set())


def test_direct_stratum_is_an_intersection(tmp_path, mini_ontology):
    rows = _rows(tmp_path, mini_ontology)
    direct = sorted(g for g, r in rows.items() if truth_table.in_direct_stratum(r, "P-ext"))
    assert direct == ["G01", "G02", "G03", "G21", "G22"]
    # G10 (wall IDA + cytosol IBA) turns P-ext when labels are recomputed without homology
    # codes; the intersection rule keeps it out of the direct P-ext stratum.
    assert rows["G10"]["label_no_homology"] == "P-ext"
    assert not truth_table.in_direct_stratum(rows["G10"], "P-ext")
    assert not truth_table.in_direct_stratum(rows["G10"], "ambiguous")


def test_internal_evidence_htp_only_column(tmp_path, mini_ontology):
    rows = _rows(tmp_path, mini_ontology)
    g14, g04, g10 = rows["G14"], rows["G04"], rows["G10"]
    assert (g14["label"], g14["internal_evidence_htp_only"]) == ("ambiguous", "yes")  # mito HDA
    assert (g04["label"], g04["internal_evidence_htp_only"]) == ("ambiguous", "no")  # cytosol IDA
    assert g10["internal_evidence_htp_only"] == "no"  # cytosol IBA only
