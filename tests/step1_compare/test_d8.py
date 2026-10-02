import http.client
import json
import urllib.error
import urllib.parse

import d8_triage
import pytest
import truth_table
from conftest import load_script


def entry(acc, reviewed, gpi_eco=(), tm=(), sgd=None):
    features = []
    if gpi_eco:
        features.append(
            {
                "type": "Lipidation",
                "description": "GPI-anchor amidated asparagine",
                "evidences": [{"evidenceCode": c} for c in gpi_eco],
            }
        )
    for code in tm:
        features.append(
            {
                "type": "Transmembrane",
                "description": "Helical",
                "evidences": [{"evidenceCode": code}],
            }
        )
    xrefs = [{"database": "SGD", "id": sgd}] if sgd else []
    kind = "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)"
    return {
        "primaryAccession": acc,
        "entryType": kind,
        "features": features,
        "uniProtKBCrossReferences": xrefs,
    }


GAS1 = entry("P22146", True, gpi_eco=["ECO:0000269"], sgd="S000004924")
YPS1 = entry("P32329", True, gpi_eco=["ECO:0000255"], sgd="S000004110")
MSB2 = entry("P32334", True, tm=["ECO:0000255"], sgd="S000003246")
TREMBL_GPI = entry("A0A000", False, gpi_eco=["ECO:0000269"], sgd="S000000009")
BARE = entry("Q00001", True, sgd="S000000010")


def test_parse_entry_records_eco_and_review_state():
    ev = d8_triage.parse_entry(MSB2)
    assert (ev.accession, ev.reviewed, ev.tm_count, ev.tm_eco) == (
        "P32334",
        True,
        1,
        ["ECO:0000255"],
    )
    assert ev.xrefs == {"SGD": {"S000003246"}}
    assert d8_triage.parse_entry(GAS1).curated_gpi
    assert not d8_triage.parse_entry(YPS1).curated_gpi  # predicted GPI site
    assert not d8_triage.parse_entry(TREMBL_GPI).curated_gpi  # not reviewed


def test_classify_pm():
    p = d8_triage.parse_entry
    assert d8_triage.classify_pm([p(GAS1)], False)[0] == "P-gpi"
    assert d8_triage.classify_pm([p(MSB2)], False) == ("PM-TM", "P32334 1 TM ECO:0000255")
    assert d8_triage.classify_pm([p(YPS1)], False)[0] == "pm-unresolved"
    assert d8_triage.classify_pm([p(YPS1)], True)[0] == "P-gpi"  # literature row
    assert d8_triage.classify_pm([], False) == ("pm-unresolved", "no UniProt entry found")


def test_query_urls():
    terms = d8_triage.query_terms([{"gene_id": "S000003246"}], "sgd")
    assert terms == {"xref:sgd-S000003246": "S000003246"}
    assert d8_triage.query_terms([{"gene_id": "Q4WXC4"}], "uniprot") == {
        "accession:Q4WXC4": "Q4WXC4"
    }
    urls = d8_triage.search_urls([f"accession:P{i:05d}" for i in range(120)])
    assert len(urls) == 3
    query = urllib.parse.parse_qs(urllib.parse.urlparse(urls[0]).query)
    assert query["fields"] == [d8_triage.UNIPROT_FIELDS] and query["format"] == ["json"]
    assert "ft_lipid:GPI-anchor" in urllib.parse.unquote_plus(d8_triage.organism_gpi_url("559292"))


def test_triage_source_and_apply(tmp_path):
    triage = load_script("03_triage_pm")
    rows = [
        {
            "source_id": "Scer",
            "gene_id": "S000003246",
            "symbol": "MSB2",
            "label": "P-ext",
            "pm_candidate": "yes",
            "stratum": "wall",
        },
        {
            "source_id": "Scer",
            "gene_id": "S000004110",
            "symbol": "YPS1",
            "label": "P-ext",
            "pm_candidate": "yes",
            "stratum": "wall",
        },
        {
            "source_id": "Scer",
            "gene_id": "S000004924",
            "symbol": "GAS1",
            "label": "ambiguous",
            "pm_candidate": "no",
            "stratum": "ambiguous",
        },
    ]
    sp = {"source_id": "Scer", "id_mapping": "sgd", "taxon_id": "559292"}

    def fake_fetch(url, tag):
        if "organism_id" in urllib.parse.unquote_plus(url):
            return [GAS1, YPS1], "2026_03"
        return [MSB2, YPS1], "2026_03"

    t, outside, counts = triage.triage_source(sp, rows, set(), fake_fetch)
    assert {r["symbol"]: r["d8_class"] for r in t} == {"MSB2": "PM-TM", "YPS1": "pm-unresolved"}
    assert [o["symbol"] for o in outside] == ["GAS1"]
    assert (counts["pm_tm"], counts["pm_unresolved"], counts["p_gpi"]) == ("1", "1", "0")
    assert counts["uniprot_release"] == "2026_03"
    triaged = triage.apply_triage(rows, t)
    assert [r["stratum"] for r in triaged] == ["PM-TM", "pm-unresolved", "ambiguous"]
    assert triaged[0]["label"] == "P-ext"  # label kept; stratum carries the D8 class


def _work(tmp_path):
    truth = [
        {
            "source_id": "Scer",
            "gene_id": "S000003246",
            "symbol": "MSB2",
            "label": "P-ext",
            "pm_candidate": "yes",
            "stratum": "extracellular-only",
        },
    ]
    truth_table.write_tsv(tmp_path / "truth_set.tsv.gz", list(truth[0]), truth)
    species = [{"source_id": "Scer", "id_mapping": "sgd", "taxon_id": "559292"}]
    truth_table.write_tsv(tmp_path / "species.tsv", list(species[0]), species)
    (tmp_path / "curated_gpi.tsv").write_text("source_id\tgene_id\tsymbol\tpmid\tnote\n")
    # Hardening: main() refuses a truth table without a full-run marker.
    (tmp_path / "extract_log.json").write_text(json.dumps({"all_sources": True}))
    return [
        "--species",
        str(tmp_path / "species.tsv"),
        "--work-dir",
        str(tmp_path),
        "--curated-gpi",
        str(tmp_path / "curated_gpi.tsv"),
    ]


def test_no_candidate_matched_stops_and_writes_nothing(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    assert triage.main(argv, fetch=lambda url, tag: ([], "2026_03")) == 2
    assert "none of 1 candidates matched" in capsys.readouterr().err
    assert not (tmp_path / "truth_set_triaged.tsv.gz").exists()
    assert not (tmp_path / "d8_triage.tsv").exists()


def test_http_error_stops(tmp_path):
    triage = load_script("03_triage_pm")

    def broken(url, tag):
        raise urllib.error.HTTPError(url, 500, "Server Error", {}, None)

    assert triage.main(_work(tmp_path), fetch=broken) == 2
    assert not (tmp_path / "truth_set_triaged.tsv.gz").exists()


def test_entry_without_expected_fields_raises():
    with pytest.raises(d8_triage.UniprotError, match="primaryAccession"):
        d8_triage.parse_entry({"features": []})


# ---- hardening tests (not in the brief) ----

OUTPUTS = (
    "d8_triage.tsv",
    "d8_gpi_outside_pext.tsv",
    "d8_counts.tsv",
    "truth_set_triaged.tsv.gz",
    "d8_run.json",
)


def good_fetch(url, tag):
    if "organism_id" in urllib.parse.unquote_plus(url):
        return [GAS1], "2026_03"
    return [MSB2], "2026_03"


def test_successful_run_writes_all_outputs_and_run_json(tmp_path):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    assert triage.main(argv, fetch=good_fetch) == 0
    for name in OUTPUTS:
        assert (tmp_path / name).exists(), name
    assert not list(tmp_path.glob(".tmp.*"))
    run = json.loads((tmp_path / "d8_run.json").read_text())
    assert run["sources"] == ["Scer"] and run["all_sources"] is True
    assert run["uniprot_release"] == {"Scer": "2026_03"}
    assert len(run["truth_set_sha256"]) == 64
    assert run["arguments"] == argv and run["python"]
    assert "unknown" == run["git_commit"] or len(run["git_commit"]) >= 7
    assert "time" not in "".join(run)  # no timestamp keys
    rows = truth_table.read_tsv(tmp_path / "truth_set_triaged.tsv.gz")
    assert rows[0]["stratum"] == "PM-TM" and rows[0]["d8_class"] == "PM-TM"


@pytest.mark.parametrize(
    "exc",
    [
        urllib.error.HTTPError("u", 500, "Server Error", {}, None),
        urllib.error.URLError("no route"),
        OSError("connection reset"),
        TimeoutError("timed out"),
        http.client.IncompleteRead(b"x"),
        json.JSONDecodeError("bad", "x", 0),
    ],
)
def test_request_errors_stop_with_exit_2(tmp_path, capsys, exc):
    triage = load_script("03_triage_pm")

    def broken(url, tag):
        raise exc

    assert triage.main(_work(tmp_path), fetch=broken) == 2
    assert capsys.readouterr().err.startswith("STOP: ")
    assert not [n for n in OUTPUTS if (tmp_path / n).exists()]


def test_response_without_results_stops(tmp_path, capsys):
    triage = load_script("03_triage_pm")

    class Opener:
        def __call__(self, request, timeout=0):
            import io

            class R(io.BytesIO):
                headers = {"x-uniprot-release": "2026_03"}

            return R(b'{"messages": ["bad query"]}')

    raw = tmp_path / "raw"
    import manifest

    def fetch(url, tag):
        orig = manifest.uniprot_pages

        def pages(u):
            return orig(u, Opener())

        manifest.uniprot_pages = pages
        try:
            return triage.uniprot_fetch(url, raw, tag)
        finally:
            manifest.uniprot_pages = orig

    assert triage.main(_work(tmp_path), fetch=fetch) == 2
    assert "no 'results'" in capsys.readouterr().err
    assert not [n for n in OUTPUTS if (tmp_path / n).exists()]


def test_entry_without_fields_stops_in_main(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    assert triage.main(_work(tmp_path), fetch=lambda u, t: ([{"features": []}], "r")) == 2
    assert "primaryAccession" in capsys.readouterr().err


def test_xref_without_database_is_uniprot_error():
    bad = entry("P1", True)
    bad["uniProtKBCrossReferences"] = [{"id": "S1"}]
    with pytest.raises(d8_triage.UniprotError, match="cross-reference"):
        d8_triage.parse_entry(bad)


def test_failure_leaves_previous_outputs_untouched(tmp_path):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    assert triage.main(argv, fetch=good_fetch) == 0
    before = {n: (tmp_path / n).read_bytes() for n in OUTPUTS}
    assert triage.main(argv, fetch=lambda u, t: ([], "2026_03")) == 2
    assert {n: (tmp_path / n).read_bytes() for n in OUTPUTS} == before
    assert not list(tmp_path.glob(".tmp.*"))


def test_partial_truth_set_is_refused(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    (tmp_path / "extract_log.json").write_text(json.dumps({"all_sources": False}))
    assert triage.main(argv, fetch=good_fetch) == 2
    assert "all_sources" in capsys.readouterr().err
    assert not (tmp_path / "d8_triage.tsv").exists()
    assert triage.main([*argv, "--allow-partial-truth-set"], fetch=good_fetch) == 0


def test_missing_extract_log_is_refused(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    (tmp_path / "extract_log.json").unlink()
    assert triage.main(argv, fetch=good_fetch) == 2
    assert "STOP:" in capsys.readouterr().err
    assert not (tmp_path / "d8_triage.tsv").exists()


def test_unknown_or_empty_sources_stop(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    assert triage.main([*argv, "--sources", "Nope"], fetch=good_fetch) == 2
    err = capsys.readouterr().err
    assert "unknown source Nope" in err and "Scer" in err
    assert triage.main([*argv, "--sources"], fetch=good_fetch) == 2
    assert "no sources selected" in capsys.readouterr().err
    assert not (tmp_path / "d8_triage.tsv").exists()


def test_sources_subset_marks_run_partial(tmp_path):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    species = [
        {"source_id": "Scer", "id_mapping": "sgd", "taxon_id": "559292"},
        {"source_id": "Other", "id_mapping": "sgd", "taxon_id": "1"},
    ]
    truth_table.write_tsv(tmp_path / "species.tsv", list(species[0]), species)
    assert triage.main([*argv, "--sources", "Scer"], fetch=good_fetch) == 0
    run = json.loads((tmp_path / "d8_run.json").read_text())
    assert run["sources"] == ["Scer"] and run["all_sources"] is False


def test_literature_row_makes_p_gpi_and_missing_curated_file_stops(tmp_path, capsys):
    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    (tmp_path / "curated_gpi.tsv").write_text(
        "source_id\tgene_id\tsymbol\tpmid\tnote\nScer\tS000003246\tMSB2\t1\tx\n"
    )
    assert triage.main(argv, fetch=good_fetch) == 0
    assert truth_table.read_tsv(tmp_path / "d8_triage.tsv")[0]["d8_class"] == "P-gpi"
    (tmp_path / "curated_gpi.tsv").unlink()
    assert triage.main(argv, fetch=good_fetch) == 2
    assert "STOP:" in capsys.readouterr().err


def test_curated_gpi_needs_eco_269_and_review():
    p = d8_triage.parse_entry
    assert p(entry("P1", True, gpi_eco=["ECO:0000269", "ECO:0000255"])).curated_gpi
    assert not p(entry("P2", True, gpi_eco=["ECO:0000255", "ECO:0000250"])).curated_gpi
    assert not p(entry("P3", False, gpi_eco=["ECO:0000269"])).curated_gpi


def test_pm_tm_needs_a_tm_feature():
    p = d8_triage.parse_entry
    assert d8_triage.classify_pm([p(BARE)], False) == (
        "pm-unresolved",
        "no curated GPI evidence and no TM feature",
    )
    assert (
        d8_triage.classify_pm([p(entry("Q2", True, tm=["ECO:0000255", "ECO:0000269"]))], False)[0]
        == "PM-TM"
    )


def test_d8_run_json_records_hash_of_current_truth_table(tmp_path):
    import manifest

    triage = load_script("03_triage_pm")
    argv = _work(tmp_path)
    assert triage.main(argv, fetch=good_fetch) == 0
    run = json.loads((tmp_path / "d8_run.json").read_text())
    assert run["truth_set_sha256"] == manifest.sha256_file(tmp_path / "truth_set.tsv.gz")


# ---- final review item 5: lipidation filter and the extra d8_counts columns ----


def _with_features(acc, features, sgd=None):
    e = entry(acc, True, sgd=sgd)
    e["features"] = features
    return e


ECO269 = [{"evidenceCode": "ECO:0000269"}]


def test_non_lipidation_feature_with_gpi_in_description_does_not_count():
    feature = {"type": "Region", "description": "GPI-anchor signal", "evidences": ECO269}
    ev = d8_triage.parse_entry(_with_features("P1", [feature]))
    assert ev.gpi_feature_count == 0 and ev.gpi_eco == [] and not ev.curated_gpi


def test_lipidation_not_starting_with_gpi_anchor_does_not_count():
    feature = {"type": "Lipidation", "description": "N-myristoyl glycine; GPI", "evidences": ECO269}
    ev = d8_triage.parse_entry(_with_features("P1", [feature]))
    assert ev.gpi_feature_count == 0 and ev.gpi_eco == [] and not ev.curated_gpi


def test_parser_counts_gpi_features_and_those_without_evidence():
    features = [
        {"type": "Lipidation", "description": "GPI-anchor amidated serine", "evidences": ECO269},
        {"type": "Lipidation", "description": "GPI-anchor amidated glycine"},
    ]
    ev = d8_triage.parse_entry(_with_features("P1", features))
    assert ev.gpi_feature_count == 2
    assert ev.gpi_features_without_eco == 1
    assert ev.curated_gpi  # classification unchanged: one feature has ECO:0000269
    assert d8_triage.parse_entry(GAS1).gpi_feature_count == 1
    assert d8_triage.parse_entry(GAS1).gpi_features_without_eco == 0


def test_d8_counts_report_no_uniprot_entry_and_gpi_feature_no_evidence():
    triage = load_script("03_triage_pm")
    rows = [
        {"source_id": "Scer", "gene_id": g, "symbol": s, "label": "P-ext", "pm_candidate": "yes"}
        for g, s in (("S1", "NOEV"), ("S2", "MISSING"), ("S000003246", "MSB2"))
    ]
    no_ev = _with_features(
        "P9", [{"type": "Lipidation", "description": "GPI-anchor amidated serine"}], sgd="S1"
    )
    sp = {"source_id": "Scer", "id_mapping": "sgd", "taxon_id": "559292"}

    def fetch(url, tag):
        if "organism_id" in urllib.parse.unquote_plus(url):
            return [], "2026_03"
        return [no_ev, MSB2], "2026_03"

    t, _, counts = triage.triage_source(sp, rows, set(), fetch)
    assert counts["no_uniprot_entry"] == "1"  # S2 only
    assert counts["gpi_feature_no_evidence"] == "1"  # P9 only
    assert {r["symbol"]: r["d8_class"] for r in t}["NOEV"] == "pm-unresolved"
    assert set(triage.D8_COUNT_COLUMNS) >= {"no_uniprot_entry", "gpi_feature_no_evidence"}
