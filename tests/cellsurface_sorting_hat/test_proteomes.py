import json

import pytest

from cellsurface_sorting_hat.proteomes import ProteomeError, write_provenance


def test_provenance_checks_the_count_and_records_the_digest(tmp_path):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKT\n>B\nMKT\n>C\nMKS\n")
    rec = write_provenance(
        tmp_path / "prov.json", "toy", "https://example.org/x", fasta, "2026-10-05", 2, 5
    )
    assert (rec["n_proteins"], rec["n_unique_sequences"], rec["n_invalid"]) == (3, 2, 0)
    assert json.loads((tmp_path / "prov.json").read_text())["sha256"] == rec["sha256"]
    with pytest.raises(ProteomeError, match="expected 10 to 20"):
        write_provenance(tmp_path / "p2.json", "toy", "x", fasta, "2026-10-05", 10, 20)


@pytest.fixture
def fasta(tmp_path):
    f = tmp_path / "p.faa"
    f.write_text(">A\nMKT\n>B\nMK*T\n>C\nMKS\n")
    return f


def test_count_bounds_are_inclusive(fasta, tmp_path):
    write_provenance(tmp_path / "lo.json", "t", "s", fasta, "2026-10-05", 3, 9)
    write_provenance(tmp_path / "hi.json", "t", "s", fasta, "2026-10-05", 1, 3)
    for lo, hi in ((4, 9), (1, 2)):
        with pytest.raises(ProteomeError, match="3 proteins"):
            write_provenance(tmp_path / "x.json", "t", "s", fasta, "2026-10-05", lo, hi)


def test_invalid_proteins_are_counted_and_ids_recorded(fasta, tmp_path):
    rec = write_provenance(tmp_path / "p.json", "t", "s", fasta, "2026-10-05", 1, 9)
    assert rec["n_invalid"] == 1
    assert rec["first_ids"] == ["A", "B", "C"]


@pytest.mark.parametrize(
    "lo,hi",
    [
        (-1, 5),
        (1, -5),
        (float("nan"), 5),
        (1, float("inf")),
        (5, 1),
        (True, 5),
        ("1", 5),
        (None, 5),
    ],
)
def test_bad_bounds_are_refused_before_any_output(fasta, tmp_path, lo, hi):
    out = tmp_path / "out" / "p.json"
    with pytest.raises(ProteomeError, match="expected_m"):
        write_provenance(out, "t", "s", fasta, "2026-10-05", lo, hi)
    assert not out.parent.exists()


@pytest.mark.parametrize("field", ["name", "source", "retrieved"])
@pytest.mark.parametrize("bad", ["", "  ", None])
def test_empty_text_fields_are_refused(fasta, tmp_path, field, bad):
    args = {"name": "t", "source": "s", "retrieved": "2026-10-05"}
    args[field] = bad
    with pytest.raises(ProteomeError, match=field):
        write_provenance(
            tmp_path / "p.json", args["name"], args["source"], fasta, args["retrieved"], 1, 9
        )
    assert not (tmp_path / "p.json").exists()


@pytest.mark.parametrize("bad", ["05/10/2026", "2026-13-01", "2026-10", "today"])
def test_retrieved_must_be_an_iso_date(fasta, tmp_path, bad):
    with pytest.raises(ProteomeError, match="retrieved"):
        write_provenance(tmp_path / "p.json", "t", "s", fasta, bad, 1, 9)
    assert not (tmp_path / "p.json").exists()


def test_unreadable_fasta_names_the_path_and_writes_nothing(tmp_path):
    bad = tmp_path / "empty.faa"
    bad.write_text("")
    with pytest.raises(ProteomeError, match="empty.faa"):
        write_provenance(tmp_path / "p.json", "t", "s", bad, "2026-10-05", 0, 9)
    with pytest.raises(ProteomeError, match="missing.faa"):
        write_provenance(
            tmp_path / "p.json", "t", "s", tmp_path / "missing.faa", "2026-10-05", 0, 9
        )
    assert not (tmp_path / "p.json").exists()


def test_digest_is_of_the_file_bytes_and_path_is_absolute(fasta, tmp_path):
    import hashlib

    rec = write_provenance(tmp_path / "p.json", "t", "s", fasta, "2026-10-05", 1, 9)
    assert rec["sha256"] == hashlib.sha256(fasta.read_bytes()).hexdigest()
    assert rec["fasta"] == str(fasta.resolve())
    assert rec["n_unique_sequences"] == 3
