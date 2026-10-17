"""Task 3.1 (D16): ``calibrate truth`` records the basis of its negatives; the default is ``assumed``."""

import pytest
from test_truth_call_status import args, call_file, setup

from cellsurface_sorting_hat.calibration.cli import main
from cellsurface_sorting_hat.call_status import load_call_source
from cellsurface_sorting_hat.engine import load_config

XFAIL_31 = pytest.mark.xfail(strict=True, reason="Task 3.1: calibrate truth --negative-basis")


def _entry(tmp_path):
    return load_call_source(call_file(tmp_path), load_config()).entries[0]


def _tight(tmp_path):
    # 200 positives and 200 negatives, all in their own cluster: narrow intervals
    setup(tmp_path, n_pos=200, n_neg=200, called_pos=190, called_neg=2)


@XFAIL_31
def test_curated_negatives_can_reach_estimated(tmp_path):
    _tight(tmp_path)
    assert main(args(tmp_path, "--negative-basis", "curated")) == 0
    entry = _entry(tmp_path)
    assert entry.measure["negative_basis"] == {"all": {"basis": "curated", "n": 200}}
    assert entry.status == "estimated"


@XFAIL_31
def test_the_default_basis_is_assumed_and_caps_the_status_at_smoke(tmp_path):
    _tight(tmp_path)
    assert main(args(tmp_path)) == 0
    entry = _entry(tmp_path)
    assert entry.measure["negative_basis"] == {"all": {"basis": "assumed", "n": 200}}
    assert entry.status == "smoke"


@XFAIL_31
def test_a_basis_other_than_curated_or_assumed_is_refused(tmp_path, capsys):
    _tight(tmp_path)
    with pytest.raises(SystemExit):
        main(args(tmp_path, "--negative-basis", "plausible"))
    assert "invalid choice: 'plausible'" in capsys.readouterr().err
