"""Reading result CSVs."""

import pytest

from surface_glyco.results import read_all, read_called

NEW = (
    "id,prediction,surface_glycoprotein_score\n"
    "A,surface_glycoprotein,0.9\nB,other,0.1\nC,surface_glycoprotein,0.6\n"
)


def _write(tmp_path, text):
    p = tmp_path / "r.csv"
    p.write_text(text)
    return p


def test_only_called_rows_are_returned(tmp_path):
    assert read_called(_write(tmp_path, NEW)) == [("A", 0.9), ("C", 0.6)]


def test_header_only_file_is_empty_not_an_error(tmp_path):
    assert read_called(_write(tmp_path, "id,prediction,surface_glycoprotein_score\n")) == []


def test_ids_with_commas_survive(tmp_path):
    text = 'id,prediction,surface_glycoprotein_score\n"sp|P1,x",surface_glycoprotein,0.5\n'
    assert read_called(_write(tmp_path, text)) == [("sp|P1,x", 0.5)]


def test_old_schema_is_rejected_with_the_file_name(tmp_path):
    p = _write(tmp_path, "id,prediction,probability_adhesion\nA,Adhesion,0.9\n")
    with pytest.raises(ValueError, match="r.csv"):
        read_called(p)


def test_read_all_returns_every_row(tmp_path):
    assert read_all(_write(tmp_path, NEW))[1] == ("B", "other", 0.1)
