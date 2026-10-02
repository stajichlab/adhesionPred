"""12_report.py: every number comes from metrics.json or findings.json; findings are quoted."""

import json
import re

import pytest

pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import phasec_fixture as pf  # noqa: E402
from conftest import load_phasec  # noqa: E402


@pytest.fixture
def reported(phasec_chain, tmp_path):
    fx = pf.copy_work(phasec_chain, tmp_path)
    assert load_phasec("12_report").main(["--work-dir", str(fx["work"])]) == 0
    out = fx["work"] / "phasec"
    return (
        fx,
        (out / "report.md").read_text(),
        json.loads((out / "metrics.json").read_text()),
        json.loads((out / "findings.json").read_text()),
    )


def test_report_numbers_come_from_metrics(reported):
    _, text, m, f = reported
    rep = load_phasec("12_report")
    allowed = rep.allowed_numbers(m, f)
    assert rep.NUMBER.findall(text)  # the check sees numbers at all
    assert rep.unsupported_numbers(text, allowed) == []
    first = re.search(r"\| (\d\.\d{3}) \[", text).group(1)
    altered = text.replace(f"| {first} [", "| 0.4242 [", 1)
    assert rep.unsupported_numbers(altered, allowed) == ["0.4242"]
    assert rep.unsupported_numbers("negative -0.4242 here", allowed) == ["-0.4242"]
    assert rep.unsupported_numbers("ruling C-11, model M35, set S2-x, R2", allowed) == []


def test_report_stops_on_a_number_without_source(phasec_chain, tmp_path, monkeypatch, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    rep = load_phasec("12_report")
    real = rep.render
    monkeypatch.setattr(rep, "render", lambda m, f: real(m, f) + "\nRecall is 0.4242.\n")
    assert rep.main(["--work-dir", str(fx["work"])]) == 2
    assert "0.4242" in capsys.readouterr().err
    assert not (fx["work"] / "phasec" / "report.md").exists()


def test_report_quotes_findings_and_labels(reported):
    _, text, m, f = reported
    holds = "yes" if f["a_b1_not_saturated"]["holds"] else "no"
    assert f"holds = {holds}." in text
    # ruling C-13: no headline ML candidate; (b) and (c) say "any" and list holds_for
    b = ", ".join(f["b_ml_beats_b1_and_r2_on_nsec_s1"]["holds_for"]) or "no ML candidate"
    want = "- (b) Any ML candidate beats B1 and R2 on N-sec FPR at the recall of R2 (S1:all): "
    assert f"{want}holds for {b}." in text
    c = ", ".join(f["c_same_under_s2"]["holds_for"]) or "no ML candidate"
    assert f"- (c) Any ML candidate, the same on every S2 test set: holds for {c}." in text
    assert "What the data cannot show" in text and "assumed, not measured" in text
    for name, ts in m["test_sets"].items():
        assert f"## {name} ({ts['label']})" in text
    # every metric row of a test set carries the label of that test set
    section = text.split("## S1:all (")[1].split("\n## ")[0]
    label = m["test_sets"]["S1:all"]["label"]
    rows = [r for r in section.splitlines() if r.startswith("| ") and "[" in r]
    assert rows and all(f"| {label} |" in r for r in rows)


def test_report_states_the_caveats(reported):
    # review M-3, M-4, M-7: fixed caveat lines, not free text about accuracy
    _, text, m, _ = reported
    rep = load_phasec("12_report")
    caveats = text.split("## Caveats\n")[1].split("\n## ")[0]
    assert len(rep.CAVEATS) == 3
    for line in rep.CAVEATS:
        assert f"- {line}" in caveats
    for review in ("M-3", "M-4", "M-7"):
        assert review in caveats
    floor = m["settings"]["estimate_min_direct_positives"]
    assert f"at least {floor} direct-evidence positives" in text


def test_literature_section_reports_recall_only(reported):
    # review I-1: the literature set has no negatives; spec 3.2 asks for recall only
    _, text, m, _ = reported
    lit = m["test_sets"]["S3-Eurotiomycetes:literature"]
    section = text.split(f"## S3-Eurotiomycetes:literature ({lit['label']})")[1].split("\n## ")[0]
    assert "| Candidate | Variant | Label | recall |" in section
    for word in ("precision", "pr_auc", "roc_auc", "fpr", "prevalence"):
        assert word not in section, word
    # a set with negatives still shows the precision columns
    s1 = text.split("## S1:all (")[1].split("\n## ")[0]
    assert "| precision |" in s1 and "Precision at R2 recall" in s1


def test_stale_input_stops(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    path = fx["work"] / "phasec" / "metrics.json"
    path.write_text(path.read_text().replace('"ci_level": 95', '"ci_level": 90'))
    assert load_phasec("12_report").main(["--work-dir", str(fx["work"])]) == 2
    assert "metrics.json differs from the SHA-256 in evaluate_run.json" in capsys.readouterr().err


def test_null_metric_prints_n_a(phasec_chain, tmp_path):
    # Task 9 fix round: by-construction null metrics print n/a, never a number
    import copy

    rep = load_phasec("12_report")
    assert rep.ci({"value": None, "lo": None, "hi": None}) == "n/a"
    assert rep.fmt(None) == "n/a"
    fx = pf.copy_work(phasec_chain, tmp_path)
    out = fx["work"] / "phasec"
    m = json.loads((out / "metrics.json").read_text())
    f = json.loads((out / "findings.json").read_text())
    m = copy.deepcopy(m)
    ts = m["test_sets"]["S1:all"]
    cand = m["settings"]["candidates"][0]
    cell = ts["truth"]["direct"]["metrics"]["all"]["V-go"][cand]
    cell["precision"] = {"value": None, "lo": None, "hi": None}
    text = rep.render(m, f)
    section = text.split("## S1:all (")[1].split("\n## ")[0]
    row = next(r for r in section.splitlines() if r.startswith(f"| {cand} | V-go |") and "n/a" in r)
    assert "n/a" in row
    assert rep.unsupported_numbers(text, rep.allowed_numbers(m, f)) == []
