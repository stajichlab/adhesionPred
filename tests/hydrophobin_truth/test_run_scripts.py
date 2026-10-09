"""Site rules for the one-off run drivers in analysis/hydrophobin_truth/run/ (the lint tests of
scripts/sorting_hat do not look here)."""

import re
import subprocess
from pathlib import Path

import pytest

RUN = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/run"
SCRIPTS = sorted(p for p in RUN.iterdir() if p.suffix in (".sh", ".sbatch"))


def problems(text):
    out = []
    if "BASH_SOURCE" in text:
        out.append("BASH_SOURCE")
    if "set -euo pipefail" not in text:
        out.append("no set -euo pipefail")
    if re.search(r"(^|[^A-Za-z_])/tmp\b", text):
        out.append("/tmp")
    if re.search(r"\$\{?SCRATCH\b(?!:\?)", text) and "SCRATCH:?" not in text:
        out.append("SCRATCH used without :?")
    if re.search(r"(^|[ (=])python3?( |$)", text, re.M):
        out.append("python not called as /usr/bin/python3.12")
    return out


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_drivers_follow_the_site_rules(script):
    assert subprocess.run(["bash", "-n", str(script)], capture_output=True).returncode == 0
    assert problems(script.read_text()) == []


def test_the_check_fails_on_each_rule():
    assert "BASH_SOURCE" in problems("set -euo pipefail\nx=$(dirname ${BASH_SOURCE[0]})")
    assert "no set -euo pipefail" in problems("echo hi")
    assert "/tmp" in problems("set -euo pipefail\ncd /tmp/x")
    assert "SCRATCH used without :?" in problems("set -euo pipefail\ncd $SCRATCH/x")
    assert problems("set -euo pipefail\ncd ${SCRATCH:?}/x") == []
    assert "python not called as /usr/bin/python3.12" in problems("set -euo pipefail\npython3 a.py")
