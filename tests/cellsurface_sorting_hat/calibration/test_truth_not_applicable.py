"""Task 3.2 (D5, D18, FD-A2): ``calibrate truth`` accepts the value ``not_applicable``."""

import csv
import gzip

from test_truth_call_status import args, call_file, setup

from cellsurface_sorting_hat.calibration.cli import main
from cellsurface_sorting_hat.call_status import load_call_source
from cellsurface_sorting_hat.engine import load_config


def _make_two_not_applicable(tmp_path, positives=1):
    """Turn the call value of the first ``positives`` positives and of one negative to NA."""
    path = tmp_path / "c.tsv.gz"
    with gzip.open(path, "rt", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    for row in rows[1:]:
        if row[0] in {f"P{k}" for k in range(positives)} | {"N0"}:
            row[3] = "not_applicable"
    with gzip.open(path, "wt", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(rows)


def test_a_table_with_not_applicable_is_accepted_and_the_count_is_reported(tmp_path):
    setup(tmp_path, n_pos=40, n_neg=40, called_pos=30, called_neg=2)
    _make_two_not_applicable(tmp_path, positives=1)
    assert main(args(tmp_path)) == 0
    entry = load_call_source(call_file(tmp_path), load_config()).entries[0]
    notes = entry.measure["notes"]
    assert "not assessable: 2 (not applicable: 2)" in notes  # one positive and one negative
    assert entry.measure["n_pos"] == 39 and entry.measure["n_neg"] == 39


def test_a_not_applicable_positive_counts_as_missed_in_the_pessimistic_bound(tmp_path):
    setup(tmp_path, n_pos=40, n_neg=40, called_pos=30, called_neg=2)
    _make_two_not_applicable(tmp_path, positives=1)
    assert main(args(tmp_path)) == 0
    notes = load_call_source(call_file(tmp_path), load_config()).entries[0].measure["notes"]
    assert "sensitivity if not assessable positives count as missed: 0.725" in notes  # 29 of 40


def test_the_notes_are_unchanged_when_no_value_is_not_applicable(tmp_path):
    setup(tmp_path, n_pos=40, n_neg=40, called_pos=30, called_neg=2)
    assert main(args(tmp_path)) == 0
    notes = load_call_source(call_file(tmp_path), load_config()).entries[0].measure["notes"]
    assert "not assessable: 0;" in notes and "not applicable" not in notes
