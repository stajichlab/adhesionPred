"""Tests for analysis/curation/cluster_benchmark.py.

No network and no MMseqs2: fetching is driven through an injected `fetch`, and the clustering step
is tested on a cluster map rather than by running the binary.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis/curation/cluster_benchmark.py"
_spec = importlib.util.spec_from_file_location("cluster_benchmark", SCRIPT)
cb = importlib.util.module_from_spec(_spec)
sys.modules["cluster_benchmark"] = cb
_spec.loader.exec_module(cb)

GROUPS = {
    "Saccharomycetes": "Saccharomycotina",
    "Sordariomycetes": "Pezizomycotina",
    "Leotiomycetes": "Pezizomycotina",
}


def _row(acc, clade="Saccharomycetes", tier="positive"):
    return {"accession": acc, "clade_class": clade, "tier": tier, "accession_db": "UniProtKB"}


# --- FASTA ------------------------------------------------------------------------------------


def test_parse_fasta_reads_multiline_records():
    text = ">A one\nMKT\nLLV\n>B\nMSTQ\n"
    assert cb.parse_fasta(text) == {"A": "MKTLLV", "B": "MSTQ"}


def test_parse_fasta_of_empty_text_is_empty():
    assert cb.parse_fasta("") == {}


def test_parse_fasta_takes_the_first_token_as_the_id():
    assert list(cb.parse_fasta(">sp|P1|NAME_X desc\nMK\n")) == ["sp|P1|NAME_X"]


# --- fetching ---------------------------------------------------------------------------------


def test_fetch_sequence_uses_uniprot_for_a_uniprot_accession():
    seen = []

    def fetch(url):
        seen.append(url)
        return ">P12345\nMKTL\n"

    assert cb.fetch_sequence("P12345", "UniProtKB", fetch) == "MKTL"
    assert "rest.uniprot.org" in seen[0]


def test_fetch_sequence_falls_back_to_ncbi_when_uniprot_has_no_record():
    seen = []

    def fetch(url):
        seen.append(url)
        if "uniprot" in url:
            raise cb.FetchError("404")
        return ">XP_1.1\nMSTQ\n"

    assert cb.fetch_sequence("XP_000001.1", "RefSeq", fetch) == "MSTQ"
    assert "eutils" in seen[-1]


def test_fetch_sequence_returns_empty_when_both_sources_fail():
    def fetch(url):
        raise cb.FetchError("nope")

    assert cb.fetch_sequence("ZZZ", "UniProtKB", fetch) == ""


def test_http_get_retries_then_succeeds(monkeypatch):
    calls = {"n": 0}

    class Resp:
        def read(self):
            return b"ok"

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def urlopen(req, timeout=0):
        calls["n"] += 1
        if calls["n"] < 3:
            raise TimeoutError("slow")
        return Resp()

    monkeypatch.setattr(cb.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(cb.time, "sleep", lambda s: None)
    assert cb.http_get("http://x") == "ok"
    assert calls["n"] == 3


def test_http_get_does_not_retry_a_404(monkeypatch):
    calls = {"n": 0}

    def urlopen(req, timeout=0):
        calls["n"] += 1
        raise cb.urllib.error.HTTPError("http://x", 404, "gone", {}, None)

    monkeypatch.setattr(cb.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(cb.time, "sleep", lambda s: None)
    with pytest.raises(cb.FetchError):
        cb.http_get("http://x")
    assert calls["n"] == 1


def test_read_existing_makes_the_fetch_step_resumable(tmp_path):
    faa = tmp_path / "sequences.faa"
    faa.write_text(">A\nMK\n")
    assert cb.read_existing(faa) == {"A": "MK"}
    assert cb.read_existing(tmp_path / "absent.faa") == {}


# --- folds on clusters ------------------------------------------------------------------------


def test_paralogues_in_one_cluster_count_once():
    # The point of the whole step: 6 FLO-like paralogues are not 6 independent observations.
    rows = [_row(f"P{i}") for i in range(6)] + [
        _row(f"N{i}", tier="hard_negative") for i in range(6)
    ]
    clusters = {f"P{i}": "P0" for i in range(6)}
    clusters.update({f"N{i}": f"N{i}" for i in range(6)})
    folds, _ = cb.fold_rows_from_clusters(rows, clusters, GROUPS, 5, 5)
    (fold,) = folds
    assert fold["n_positive"] == 1
    assert fold["n_hard_negative"] == 6
    assert fold["scored"] == "no"  # one positive cluster cannot carry a fold


def test_a_clade_with_enough_clusters_on_both_sides_is_scored():
    rows = [_row(f"P{i}") for i in range(5)] + [
        _row(f"N{i}", tier="hard_negative") for i in range(5)
    ]
    clusters = {r["accession"]: r["accession"] for r in rows}
    (fold,) = cb.fold_rows_from_clusters(rows, clusters, GROUPS, 5, 5)[0]
    assert fold["scored"] == "yes" and fold["level"] == "class"


def test_thin_clades_pool_into_their_group_on_cluster_counts():
    rows = [_row(f"P{i}", "Sordariomycetes") for i in range(3)] + [
        _row(f"Q{i}", "Leotiomycetes") for i in range(3)
    ]
    rows += [_row(f"N{i}", "Sordariomycetes", "hard_negative") for i in range(6)]
    clusters = {r["accession"]: r["accession"] for r in rows}
    (fold,) = cb.fold_rows_from_clusters(rows, clusters, GROUPS, 5, 5)[0]
    assert fold["fold"] == "Pezizomycotina"
    assert fold["n_positive"] == 6 and fold["scored"] == "yes"


def test_a_protein_without_a_sequence_is_reported_not_counted():
    rows = [_row("P1"), _row("P2"), _row("N1", tier="hard_negative")]
    folds, unclustered = cb.fold_rows_from_clusters(rows, {"P1": "P1", "N1": "N1"}, GROUPS, 1, 1)
    assert unclustered == ["P2"]
    assert folds[0]["n_positive"] == 1


def test_moonlighting_and_excluded_rows_are_not_clustered_into_a_fold():
    rows = [
        _row("P1"),
        _row("M1", tier="positive_moonlighting"),
        _row("X1", tier="excluded"),
        _row("N1", tier="hard_negative"),
    ]
    clusters = {r["accession"]: r["accession"] for r in rows}
    folds, _ = cb.fold_rows_from_clusters(rows, clusters, GROUPS, 1, 1)
    assert folds[0]["n_positive"] == 1 and folds[0]["n_hard_negative"] == 1


def test_an_unmapped_class_is_an_error():
    rows = [_row("P1", "Unknownomycetes")]
    with pytest.raises(cb.BuildError, match="no group"):
        cb.fold_rows_from_clusters(rows, {"P1": "P1"}, GROUPS, 1, 1)


# --- mmseqs wrapper ---------------------------------------------------------------------------


def test_run_mmseqs_raises_with_the_stderr_when_the_binary_fails(tmp_path, monkeypatch):
    class Proc:
        returncode = 1
        stderr = "boom"
        stdout = ""

    monkeypatch.setattr(cb.subprocess, "run", lambda *a, **k: Proc())
    with pytest.raises(cb.BuildError, match="boom"):
        cb.run_mmseqs("mmseqs", tmp_path / "x.faa", tmp_path / "tmp", 0.3, 0.5)


def test_run_mmseqs_raises_when_the_cluster_table_is_absent(tmp_path, monkeypatch):
    class Proc:
        returncode = 0
        stderr = ""
        stdout = ""

    monkeypatch.setattr(cb.subprocess, "run", lambda *a, **k: Proc())
    with pytest.raises(cb.BuildError, match="was not written"):
        cb.run_mmseqs("mmseqs", tmp_path / "x.faa", tmp_path / "tmp", 0.3, 0.5)


def test_run_mmseqs_reads_the_cluster_table(tmp_path, monkeypatch):
    tmpdir = tmp_path / "tmp"
    tmpdir.mkdir()

    class Proc:
        returncode = 0
        stderr = ""
        stdout = ""

    def fake_run(cmd, **kwargs):
        (tmpdir / "clu_cluster.tsv").write_text("R1\tR1\nR1\tM2\nR3\tR3\n")
        return Proc()

    monkeypatch.setattr(cb.subprocess, "run", fake_run)
    assert cb.run_mmseqs("mmseqs", tmp_path / "x.faa", tmpdir, 0.3, 0.5) == {
        "R1": "R1",
        "M2": "R1",
        "R3": "R3",
    }


def test_cluster_step_refuses_to_run_before_the_fetch_step(tmp_path, capsys):
    (tmp_path / "proteins.tsv").write_text("accession\ttier\n")
    rc = cb.main(["cluster", "--benchmark", str(tmp_path), "--tmpdir", str(tmp_path / "t")])
    assert rc == 1
    assert "fetch step first" in capsys.readouterr().err
