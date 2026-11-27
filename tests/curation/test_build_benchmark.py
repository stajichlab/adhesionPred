"""Tests for analysis/curation/build_benchmark.py (the frozen adhesion benchmark builder).

Fixtures are synthetic: no real accession, PMID or curated row is used, so a change to the curated
tables cannot make these tests pass or fail for the wrong reason.
"""

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis/curation/build_benchmark.py"
_spec = importlib.util.spec_from_file_location("build_benchmark", SCRIPT)
bb = importlib.util.module_from_spec(_spec)
sys.modules["build_benchmark"] = bb
_spec.loader.exec_module(bb)


def _source_row(**over):
    row = {c: "" for c in bb.REQUIRED_SOURCE_COLUMNS}
    row.update(
        accession="X00001",
        accession_db="UniProtKB",
        gene="TST1",
        protein_name="Test adhesin",
        organism="Testomyces fictus",
        taxon_id="999999",
        clade_class="Testomycetes",
        label="adhesin",
        mechanism="domain_binding",
        evidence_code="ko_assay",
        assay_summary="Deletion reduced adhesion to plastic by 60%.",
        pmids="11111111",
        uniprot_length="400",
        verified="pmid+acc",
        notes="",
    )
    row.update(over)
    return row


def _curated_row(**over):
    row = {
        "accession": "C00001",
        "gene": "CUR1",
        "protein_name": "Curated protein",
        "organism": "Testomyces fictus",
        "taxon_id": "999999",
        "cls": "adhesin",
        "evidence_level": "E1",
        "family": "PA14",
        "evidence_summary": "binds fibronectin",
        "evidence_codes": "literature",
        "pmids": "22222222",
        "length": "500",
        "signal_peptide": "yes",
        "gpi_anchor": "yes",
        "surface": "yes",
        "moonlighting": "",
    }
    row.update(over)
    return row


def _write_tsv(path, rows, columns=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(columns or rows[0])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


# --- tier rules on literature-curation rows ---------------------------------------------------


def test_adhesin_with_direct_evidence_is_a_positive():
    tier, indirect, reason = bb.classify_source_row(_source_row(), "f")
    assert (tier, indirect, reason) == ("positive", "", "")


@pytest.mark.parametrize("code", ["ko_assay", "het_expr", "binding"])
def test_each_direct_evidence_code_is_accepted(code):
    tier, _, _ = bb.classify_source_row(_source_row(evidence_code=code), "f")
    assert tier == "positive"


def test_adhesin_without_direct_evidence_is_excluded_with_a_reason():
    tier, _, reason = bb.classify_source_row(_source_row(evidence_code="expression_induced"), "f")
    assert tier == "excluded"
    assert "direct evidence" in reason


def test_moonlighting_adhesin_goes_to_its_own_tier_not_the_positives():
    # It has real binding evidence but no secretory architecture, so no layout method can predict it.
    tier, _, reason = bb.classify_source_row(
        _source_row(mechanism="moonlighting", evidence_code="binding"), "f"
    )
    assert (tier, reason) == ("positive_moonlighting", "")


def test_hard_negative_needs_a_demonstrated_function():
    ok, _, _ = bb.classify_source_row(
        _source_row(label="hard_negative", evidence_code="function_shown"), "f"
    )
    bad, _, reason = bb.classify_source_row(
        _source_row(label="hard_negative", evidence_code=""), "f"
    )
    assert ok == "hard_negative"
    assert bad == "excluded" and "function_shown" in reason


def test_indirect_negative_is_flagged_so_it_is_never_read_as_untested():
    tier, indirect, _ = bb.classify_source_row(
        _source_row(
            label="hard_negative",
            evidence_code="function_shown",
            notes="adhesion phenotype on deletion but indirect mechanism; owner decision",
        ),
        "f",
    )
    assert (tier, indirect) == ("hard_negative", "yes")


@pytest.mark.parametrize("field", ["accession", "pmids"])
def test_a_row_without_an_accession_or_a_pmid_is_excluded(field):
    tier, _, reason = bb.classify_source_row(_source_row(**{field: ""}), "f")
    assert tier == "excluded"
    assert reason


def test_unknown_label_is_excluded_rather_than_guessed():
    tier, _, reason = bb.classify_source_row(_source_row(label="probable_adhesin"), "f")
    assert tier == "excluded" and "label" in reason


# --- tier rules on the existing curated table -------------------------------------------------


def test_only_e1_curated_adhesins_become_positives():
    assert bb.classify_curated_row(_curated_row(evidence_level="E1"))[0] == "positive"
    for level in ("E2", "E3", ""):
        tier, _, reason = bb.classify_curated_row(_curated_row(evidence_level=level))
        assert tier == "excluded", level
        assert "not direct" in reason


def test_curated_adhesin_without_secretory_architecture_is_not_a_plain_positive():
    tier, _, _ = bb.classify_curated_row(_curated_row(signal_peptide="no", gpi_anchor="no"))
    assert tier == "positive_moonlighting"


def test_indirect_surface_class_becomes_a_flagged_hard_negative():
    tier, indirect, _ = bb.classify_curated_row(
        _curated_row(cls="surface_other_adhesion_phenotype", evidence_level="")
    )
    assert (tier, indirect) == ("hard_negative", "yes")


def test_the_surface_filter_keeps_transcription_factors_out_of_the_negatives():
    # A negative set of "surface proteins" holding cytoplasmic regulators would make C1 too easy.
    tier, _, reason = bb.classify_curated_row(
        _curated_row(
            cls="surface_other_adhesion_phenotype",
            evidence_level="",
            signal_peptide="no",
            gpi_anchor="no",
        )
    )
    assert tier == "excluded"
    assert "not secreted" in reason


def test_indirect_regulator_class_is_not_used_at_all():
    tier, _, reason = bb.classify_curated_row(_curated_row(cls="indirect_regulator"))
    assert tier == "excluded" and "not used" in reason


# --- merging ----------------------------------------------------------------------------------


def test_merge_prefers_a_classified_row_over_an_excluded_duplicate():
    rows = [
        {"accession": "A1", "tier": "excluded", "source": "curated:adhesins.tsv"},
        {"accession": "A1", "tier": "positive", "source": "curation:basidio"},
    ]
    assert bb.merge(rows)[0]["tier"] == "positive"
    assert len(bb.merge(rows)) == 1


def test_merge_raises_on_a_real_tier_conflict_instead_of_choosing():
    rows = [
        {"accession": "A1", "tier": "positive", "source": "curation:one"},
        {"accession": "A1", "tier": "hard_negative", "source": "curation:two"},
    ]
    with pytest.raises(bb.BuildError, match="tier conflict"):
        bb.merge(rows)


# --- schema and end to end --------------------------------------------------------------------


def test_a_source_file_missing_a_column_is_refused(tmp_path):
    rows = [_source_row()]
    del rows[0]["pmids"]
    path = tmp_path / "bad.tsv"
    _write_tsv(path, rows)
    with pytest.raises(bb.BuildError, match="missing columns"):
        bb.from_sources(tmp_path)


def test_empty_sources_directory_is_refused(tmp_path):
    with pytest.raises(bb.BuildError, match="no source"):
        bb.from_sources(tmp_path)


def test_build_writes_the_frozen_outputs_and_counts(tmp_path):
    curated = tmp_path / "adhesins.tsv"
    _write_tsv(
        curated,
        [
            _curated_row(accession="C1"),
            _curated_row(accession="C2", cls="hard_negative", evidence_level="N1"),
            _curated_row(accession="C3", cls="indirect_regulator"),
        ],
    )
    sources = tmp_path / "sources"
    _write_tsv(
        sources / "testclade.tsv",
        [
            _source_row(accession="S1"),
            _source_row(accession="S2", mechanism="moonlighting", evidence_code="binding"),
            _source_row(accession="S3", label="hard_negative", evidence_code="function_shown"),
        ],
        columns=bb.REQUIRED_SOURCE_COLUMNS,
    )
    out = tmp_path / "adhesion_vtest"
    clade_map = tmp_path / "clade_map.tsv"
    _write_clade_map(clade_map, "999999\tTestomycetes\tTestomyces fictus\n")

    rc = bb.main(
        [
            "--version",
            "adhesion_vtest",
            "--out",
            str(out),
            "--curated",
            str(curated),
            "--sources",
            str(sources),
            "--clade-map",
            str(clade_map),
        ]
    )
    assert rc == 0

    rows = {r["accession"]: r for r in csv.DictReader(open(out / "proteins.tsv"), delimiter="\t")}
    assert rows["C1"]["tier"] == "positive"
    assert rows["C2"]["tier"] == "hard_negative"
    assert rows["C3"]["tier"] == "excluded"
    assert rows["S1"]["tier"] == "positive"
    assert rows["S2"]["tier"] == "positive_moonlighting"
    assert rows["S3"]["tier"] == "hard_negative"
    # Nothing is dropped silently: an excluded row keeps its reason.
    assert rows["C3"]["exclude_reason"]

    version = (out / "VERSION").read_text()
    assert "adhesion_vtest" in version
    assert "tier positive\t2" in version
    assert "sha256" in version

    counts = (out / "counts.tsv").read_text()
    # C1 and S1 are the positives, C2 and S3 the hard negatives; all four share the fixture taxon.
    assert "Testomycetes\t2\t2" in counts

    sums = (out / "SHA256SUMS").read_text()
    for name in ("proteins.tsv", "VERSION", "counts.tsv"):
        assert name in sums


def test_build_refuses_to_overwrite_a_frozen_directory(tmp_path, capsys):
    curated = tmp_path / "adhesins.tsv"
    _write_tsv(curated, [_curated_row()])
    sources = tmp_path / "sources"
    _write_tsv(sources / "c.tsv", [_source_row()], columns=bb.REQUIRED_SOURCE_COLUMNS)
    out = tmp_path / "frozen"
    out.mkdir()
    clade_map = tmp_path / "clade_map.tsv"
    _write_clade_map(clade_map, "999999\tTestomycetes\tTestomyces fictus\n")

    assert (
        bb.main(
            [
                "--version",
                "v",
                "--out",
                str(out),
                "--curated",
                str(curated),
                "--sources",
                str(sources),
                "--clade-map",
                str(clade_map),
            ]
        )
        == 1
    )
    assert "exists" in capsys.readouterr().err


def test_build_refuses_a_set_with_no_positives(tmp_path, capsys):
    curated = tmp_path / "adhesins.tsv"
    _write_tsv(curated, [_curated_row(cls="indirect_regulator")])
    sources = tmp_path / "sources"
    _write_tsv(
        sources / "c.tsv",
        [_source_row(label="hard_negative", evidence_code="function_shown")],
        columns=bb.REQUIRED_SOURCE_COLUMNS,
    )
    clade_map = tmp_path / "clade_map.tsv"
    _write_clade_map(clade_map, "999999\tTestomycetes\tTestomyces fictus\n")
    assert (
        bb.main(
            [
                "--version",
                "v",
                "--out",
                str(tmp_path / "o"),
                "--curated",
                str(curated),
                "--sources",
                str(sources),
                "--clade-map",
                str(clade_map),
            ]
        )
        == 1
    )
    assert "no positives" in capsys.readouterr().err


def test_dry_run_writes_nothing(tmp_path):
    curated = tmp_path / "adhesins.tsv"
    _write_tsv(curated, [_curated_row()])
    sources = tmp_path / "sources"
    _write_tsv(sources / "c.tsv", [_source_row()], columns=bb.REQUIRED_SOURCE_COLUMNS)
    out = tmp_path / "o"
    clade_map = tmp_path / "clade_map.tsv"
    _write_clade_map(clade_map, "999999\tTestomycetes\tTestomyces fictus\n")
    assert (
        bb.main(
            [
                "--version",
                "v",
                "--out",
                str(out),
                "--curated",
                str(curated),
                "--sources",
                str(sources),
                "--clade-map",
                str(clade_map),
                "--dry-run",
            ]
        )
        == 0
    )
    assert not out.exists()


def test_looks_secreted_accepts_either_route():
    assert bb.looks_secreted("yes", "no")
    assert bb.looks_secreted("no", "yes")
    assert not bb.looks_secreted("no", "no")
    assert not bb.looks_secreted("", "")


# --- clade map --------------------------------------------------------------------------------


def _write_clade_map(path, body):
    path.write_text("# comment line is ignored\ntaxon_id\tclade_class\torganism\n" + body)


def test_clade_map_reads_rows_and_ignores_comments(tmp_path):
    path = tmp_path / "clade_map.tsv"
    _write_clade_map(path, "999999\tTestomycetes\tTestomyces fictus\n")
    assert bb.read_clade_map(path) == {"999999": "Testomycetes"}


def test_clade_map_refuses_a_taxon_mapped_to_two_classes(tmp_path):
    path = tmp_path / "clade_map.tsv"
    _write_clade_map(path, "1\tA\tx\n1\tB\ty\n")
    with pytest.raises(bb.BuildError, match="two classes"):
        bb.read_clade_map(path)


def test_clade_map_refuses_an_incomplete_row(tmp_path):
    path = tmp_path / "clade_map.tsv"
    _write_clade_map(path, "1\t\tx\n")
    with pytest.raises(bb.BuildError, match="no taxon_id or no clade_class"):
        bb.read_clade_map(path)


def test_apply_clades_fills_a_missing_class_from_the_map():
    rows = [{"clade_class": "", "taxon_id": "999999", "tier": "positive", "exclude_reason": ""}]
    assert bb.apply_clades(rows, {"999999": "Testomycetes"})[0]["clade_class"] == "Testomycetes"


def test_apply_clades_keeps_a_class_the_curator_already_set():
    rows = [
        {"clade_class": "Stated", "taxon_id": "999999", "tier": "positive", "exclude_reason": ""}
    ]
    assert bb.apply_clades(rows, {"999999": "Other"})[0]["clade_class"] == "Stated"


def test_apply_clades_excludes_an_unmapped_taxon_rather_than_guessing():
    # Guessing would put the protein in the wrong leave-one-clade-out fold.
    rows = [{"clade_class": "", "taxon_id": "123", "tier": "positive", "exclude_reason": ""}]
    out = bb.apply_clades(rows, {"999999": "Testomycetes"})[0]
    assert out["tier"] == "excluded"
    assert "clade map" in out["exclude_reason"]


def test_apply_clades_leaves_already_excluded_rows_alone():
    rows = [{"clade_class": "", "taxon_id": "123", "tier": "excluded", "exclude_reason": "earlier"}]
    assert bb.apply_clades(rows, {})[0]["exclude_reason"] == "earlier"


def test_the_committed_clade_map_is_valid_and_covers_the_curated_table():
    clade_map = bb.read_clade_map(bb.ROOT / "data/benchmark/clade_map.tsv")
    assert len(clade_map) >= 20
    rows = bb.from_curated(bb.ROOT / "data/curated/adhesins/adhesins.tsv")
    scored = [r for r in rows if r["tier"] != "excluded"]
    real = {r["taxon_id"] for r in scored if r["taxon_id"] not in ("", "0")}
    unmapped = sorted(real - set(clade_map))
    # A taxon that carries a scored protein must be classified on purpose, not dropped by default.
    assert unmapped == [], unmapped
