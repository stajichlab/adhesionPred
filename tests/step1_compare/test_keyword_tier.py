import gzip
import json

import keyword_tier
import manifest
import pytest
import seqhash
from conftest import load_script

SEEDS_TSV = (
    "# comment line\n"
    "gene\tuniprot_query\tspecies\torder\tfamily\tclass\tevidence_level\tmoonlighting\tpmids"
    "\tevidence_summary\n"
    "SOWgp58\taccession:Q8NK60\tCoccidioides immitis\tOnygenales\tf\tadhesin\tE1\tno\t1\ts\n"
    "YPS3\t\tHistoplasma capsulatum\tOnygenales\tf\tadhesin\tE2\tno\t1\tunresolved\n"
)


def srow(acc, taxon="330879", source="uniprot_surface_kw", gene="g"):
    return {"accession": acc, "gene": gene, "genome": "G", "taxon_id": taxon, "source": source}


def test_read_seed_accessions_skips_comments_and_blank(tmp_path):
    path = tmp_path / "seeds.tsv"
    path.write_text(SEEDS_TSV)
    assert keyword_tier.read_seed_accessions(path) == {"Q8NK60": "SOWgp58"}


def test_tc_excludes_test_proteins():
    surface = [
        srow("Q8NK60", taxon="246410"),  # literature accession
        srow("O14000", taxon="284812"),  # S. pombe
        srow("Q4WAAA"),  # held-out accession (A. fumigatus truth gene)
        srow("Q4WBBB"),  # no sequence found
        srow("A0ACG8", taxon="443226"),  # same sequence as Q8NK60 under another accession
        srow("Q4WCCC"),  # same sequence as a held-out truth protein
        srow("Q4WDDD"),  # kept
        srow("P50142", source="curated_literature"),  # not T-c: ignored
    ]
    seqs = {
        "Q8NK60": "MKSOW",
        "O14000": "MPOM",
        "Q4WAAA": "MAAA",
        "A0ACG8": "mksow*",
        "Q4WCCC": "MCCC",
        "Q4WDDD": "MDDD",
        "P50142": "MHSP",
    }
    kept, removed = keyword_tier.build_keyword_tier(
        surface,
        seqs,
        {"Q8NK60": "SOWgp58"},
        {"Q4WAAA": "Afum_ASPFU:Q4WAAA"},
        {seqhash.seq_sha256("MCCC"): "Afum_ASPFU:Q4WZZZ"},
    )
    assert [r["accession"] for r in kept] == ["Q4WDDD"]
    assert {r["accession"]: r["reason"] for r in removed} == {
        "Q8NK60": "literature_accession",
        "O14000": "spombe_taxon",
        "Q4WAAA": "heldout_accession",
        "Q4WBBB": "no_sequence",
        "A0ACG8": "literature_hash",
        "Q4WCCC": "heldout_hash",
    }
    assert kept[0]["seq_sha256"] == seqhash.seq_sha256("MDDD") and kept[0]["tier"] == "T-c"


def test_heldout_sets_skip_training_and_unlabelled_rows():
    build = load_script("04_build_keyword_tier")
    species = [
        {"source_id": "Scer_SGD", "role": "train", "id_mapping": "sgd"},
        {"source_id": "Afum_ASPFU", "role": "test_clade", "id_mapping": "uniprot"},
        {"source_id": "Spom_PomBase", "role": "test_species", "id_mapping": "pombase"},
    ]
    seq_rows = [
        {"source_id": "Scer_SGD", "gene_id": "S1", "label": "P-ext", "seq_sha256": "h1"},
        {"source_id": "Afum_ASPFU", "gene_id": "Q4W1", "label": "N-sec", "seq_sha256": "h2"},
        {"source_id": "Afum_ASPFU", "gene_id": "Q4W2", "label": "unlabelled", "seq_sha256": "h3"},
        {"source_id": "Spom_PomBase", "gene_id": "SPAC1.01", "label": "P-ext", "seq_sha256": "h4"},
    ]
    accessions, hashes = build.heldout_sets(species, seq_rows)
    assert accessions == {"Q4W1": "Afum_ASPFU:Q4W1"}
    assert hashes == {"h2": "Afum_ASPFU:Q4W1", "h4": "Spom_PomBase:SPAC1.01"}


def test_cached_sequences_without_sidecar_are_refused(tmp_path, capsys):
    build = load_script("04_build_keyword_tier")
    (tmp_path / "keyword_sequences.fasta.gz").write_bytes(gzip.compress(b">sp|P1|X\nMK\n"))
    surface = tmp_path / "surface.tsv"
    surface.write_text("accession\tgene\tgenome\ttaxon_id\tsource\n")
    seeds = tmp_path / "seeds.tsv"
    seeds.write_text(SEEDS_TSV)
    argv = ["--surface", str(surface), "--seeds", str(seeds), "--work-dir", str(tmp_path)]
    assert build.main(argv) == 2
    assert "keyword_sequences.json missing" in capsys.readouterr().err


def test_check_cached_accepts_matching_sidecar_and_rejects_changed_file(tmp_path):
    build = load_script("04_build_keyword_tier")
    seq_path = tmp_path / "keyword_sequences.fasta.gz"
    seq_path.write_bytes(gzip.compress(b">sp|P1|X\nMK\n"))
    meta = {"sha256": manifest.sha256_file(seq_path), "uniprot_release": "2026_03"}
    (tmp_path / "keyword_sequences.json").write_text(json.dumps(meta))
    assert build.check_cached(seq_path)["uniprot_release"] == "2026_03"
    seq_path.write_bytes(gzip.compress(b">sp|P1|X\nMA\n"))
    with pytest.raises(manifest.DownloadError, match="does not match"):
        build.check_cached(seq_path)


# ---- hardening tests (added; not in the brief) ----

SPECIES_TSV = (
    "source_id\trole\tid_mapping\n" "Scer_SGD\ttrain\tsgd\n" "Afum_ASPFU\ttest_clade\tuniprot\n"
)
SURFACE_HEADER = "accession\tgene\tgenome\ttaxon_id\tsource\n"
SURFACE_ROWS = (
    "Q8NK60\tSOW\tG\t246410\tuniprot_surface_kw\n"
    "O14000\tPOM\tG\t4896\tuniprot_surface_kw\n"
    "Q4WAAA\tA\tG\t330879\tuniprot_surface_kw\n"
    "Q4WCCC\tC\tG\t330879\tuniprot_surface_kw\n"
    "Q4WDDD\tD\tG\t330879\tuniprot_surface_kw\n"
)
FASTA = (
    ">sp|Q8NK60|X\nMKSOW\n>sp|O14000|X\nMPOM\n>sp|Q4WAAA|X\nMAAA\n"
    ">sp|Q4WCCC|X\nMCCC\n>sp|Q4WDDD|X\nMDDD\n"
)
TRUTH_SEQ = (
    "source_id\tgene_id\tlabel\tseq_sha256\n"
    f"Scer_SGD\tS1\tP-ext\t{seqhash.seq_sha256('MDDD')}\n"
    f"Afum_ASPFU\tQ4WZZZ\tN-sec\t{seqhash.seq_sha256('MCCC')}\n"
    "Afum_ASPFU\tQ4WAAA\tP-ext\th-aaa\n"
    "Afum_ASPFU\tQ4WUNL\tunlabelled\th-unl\n"
)


def make_work(tmp_path, extract_all=True, sequence_all=True, with_logs=True):
    work = tmp_path / "work"
    work.mkdir()
    seq_path = work / "keyword_sequences.fasta.gz"
    seq_path.write_bytes(gzip.compress(FASTA.encode()))
    meta = {"sha256": manifest.sha256_file(seq_path), "uniprot_release": "2026_03"}
    (work / "keyword_sequences.json").write_text(json.dumps(meta))
    with gzip.open(work / "truth_sequences.tsv.gz", "wt") as handle:
        handle.write(TRUTH_SEQ)
    if with_logs:
        (work / "extract_log.json").write_text(json.dumps({"all_sources": extract_all}))
        (work / "sequence_run.json").write_text(json.dumps({"all_sources": sequence_all}))
    (tmp_path / "surface.tsv").write_text(SURFACE_HEADER + SURFACE_ROWS)
    (tmp_path / "seeds.tsv").write_text(SEEDS_TSV)
    (tmp_path / "species.tsv").write_text(SPECIES_TSV)
    return work


def argv_for(tmp_path, work, *extra):
    return [
        "--surface", str(tmp_path / "surface.tsv"),
        "--seeds", str(tmp_path / "seeds.tsv"),
        "--species", str(tmp_path / "species.tsv"),
        "--work-dir", str(work),
        *extra,
    ]  # fmt: skip


def read_removed(work):
    import truth_table

    return {r["accession"]: r for r in truth_table.read_tsv(work / "keyword_tier_removed.tsv")}


def test_main_end_to_end_writes_tier_removal_log_and_run_json(tmp_path, capsys):
    build = load_script("04_build_keyword_tier")
    import truth_table

    work = make_work(tmp_path)
    assert build.main(argv_for(tmp_path, work)) == 0
    kept = truth_table.read_tsv(work / "keyword_tier.tsv.gz")
    assert [r["accession"] for r in kept] == ["Q4WDDD"]
    removed = read_removed(work)
    assert {a: r["reason"] for a, r in removed.items()} == {
        "Q8NK60": "literature_accession",
        "O14000": "spombe_taxon",
        "Q4WAAA": "heldout_accession",
        "Q4WCCC": "heldout_hash",
    }
    assert removed["Q4WCCC"]["matched"] == "Afum_ASPFU:Q4WZZZ"
    assert len(removed) == 4  # one row per removed protein
    run = json.loads((work / "keyword_tier_run.json").read_text())
    assert run["uniprot_release"] == "2026_03"
    assert run["kept"] == 1 and run["removed"] == 4
    assert set(run["inputs_sha256"]) == {
        "surface",
        "seeds",
        "species",
        "truth_sequences",
        "keyword_sequences",
    }
    assert not [p for p in work.iterdir() if p.name.startswith(".tmp")]
    assert "kept=1 removed=4" in capsys.readouterr().out


def test_training_source_hash_does_not_remove_a_protein(tmp_path):
    # Scer_SGD (train) has a gene with the sequence of Q4WDDD: Q4WDDD must stay.
    build = load_script("04_build_keyword_tier")
    work = make_work(tmp_path)
    assert build.main(argv_for(tmp_path, work)) == 0
    assert "Q4WDDD" not in read_removed(work)


def test_run_json_is_reproducible_without_timestamps(tmp_path):
    build = load_script("04_build_keyword_tier")
    work = make_work(tmp_path)
    assert build.main(argv_for(tmp_path, work)) == 0
    first = (work / "keyword_tier_run.json").read_bytes()
    assert build.main(argv_for(tmp_path, work)) == 0
    assert (work / "keyword_tier_run.json").read_bytes() == first


@pytest.mark.parametrize("taxon", ["284812", "4896"])
def test_both_spombe_taxa_are_removed(taxon):
    kept, removed = keyword_tier.build_keyword_tier(
        [srow("O1", taxon=taxon)], {"O1": "MK"}, {}, {}, {}
    )
    assert not kept and removed[0]["reason"] == "spombe_taxon"


@pytest.mark.parametrize(
    ("extract_all", "sequence_all", "with_logs"),
    [(False, True, True), (True, False, True), (True, True, False)],
)
def test_partial_or_unlogged_truth_tables_are_refused(
    tmp_path, capsys, extract_all, sequence_all, with_logs
):
    build = load_script("04_build_keyword_tier")
    work = make_work(tmp_path, extract_all, sequence_all, with_logs)
    assert build.main(argv_for(tmp_path, work)) == 2
    err = capsys.readouterr().err
    assert err.startswith("STOP:")
    assert "all_sources" in err or "cannot read" in err
    assert not (work / "keyword_tier.tsv.gz").exists()
    assert not (work / "keyword_tier_run.json").exists()
    assert build.main(argv_for(tmp_path, work, "--allow-partial-truth-set")) == 0


def test_stop_leaves_previous_outputs_untouched(tmp_path, capsys):
    build = load_script("04_build_keyword_tier")
    work = make_work(tmp_path)
    assert build.main(argv_for(tmp_path, work)) == 0
    before = {n: (work / n).read_bytes() for n in build.OUTPUT_NAMES}
    (work / "extract_log.json").write_text(json.dumps({"all_sources": False}))
    assert build.main(argv_for(tmp_path, work)) == 2
    assert {n: (work / n).read_bytes() for n in build.OUTPUT_NAMES} == before


def test_missing_truth_sequences_stop_without_outputs(tmp_path, capsys):
    build = load_script("04_build_keyword_tier")
    work = make_work(tmp_path)
    (work / "truth_sequences.tsv.gz").unlink()
    assert build.main(argv_for(tmp_path, work)) == 2
    assert capsys.readouterr().err.startswith("STOP:")
    assert not [p for p in work.iterdir() if p.name.startswith("keyword_tier")]


def test_cache_with_changed_sha_is_refused_by_main(tmp_path, capsys):
    build = load_script("04_build_keyword_tier")
    work = make_work(tmp_path)
    (work / "keyword_sequences.fasta.gz").write_bytes(gzip.compress(b">sp|P1|X\nMA\n"))
    assert build.main(argv_for(tmp_path, work)) == 2
    assert "does not match" in capsys.readouterr().err


def test_cache_with_empty_release_is_refused(tmp_path):
    build = load_script("04_build_keyword_tier")
    seq_path = tmp_path / "keyword_sequences.fasta.gz"
    seq_path.write_bytes(gzip.compress(b">sp|P1|X\nMK\n"))
    meta = {"sha256": manifest.sha256_file(seq_path), "uniprot_release": ""}
    (tmp_path / "keyword_sequences.json").write_text(json.dumps(meta))
    with pytest.raises(manifest.DownloadError, match="does not match"):
        build.check_cached(seq_path)


class FakeResponse:
    def __init__(self, body, headers):
        self._body = body
        self.headers = headers

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_fetch_sequences_writes_file_and_sidecar_with_release(tmp_path):
    build = load_script("04_build_keyword_tier")

    def opener(request, timeout):
        return FakeResponse(b">sp|P1|X\nMK\n", {"X-UniProt-Release": "2026_03"})

    dest = tmp_path / "keyword_sequences.fasta.gz"
    meta = build.fetch_sequences(["P1", "P2"], dest, opener=opener)
    assert meta["uniprot_release"] == "2026_03" and meta["accessions"] == "2"
    assert build.check_cached(dest)["sha256"] == manifest.sha256_file(dest)
    assert not list(tmp_path.glob("*.part"))


def test_network_error_during_fetch_is_a_stop_with_no_files(tmp_path, capsys, monkeypatch):
    import urllib.error

    build = load_script("04_build_keyword_tier")
    work = make_work(tmp_path)
    (work / "keyword_sequences.fasta.gz").unlink()
    (work / "keyword_sequences.json").unlink()

    def boom(url, opener=None):
        raise urllib.error.URLError("no route")

    monkeypatch.setattr(manifest, "http_get", boom)
    assert build.main(argv_for(tmp_path, work, "--fetch")) == 2
    assert capsys.readouterr().err.startswith("STOP:")
    assert not (work / "keyword_sequences.fasta.gz").exists()
    assert not [p for p in work.iterdir() if p.name.endswith(".part")]
    assert not (work / "keyword_tier.tsv.gz").exists()
