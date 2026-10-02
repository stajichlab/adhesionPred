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
    lab_b = m["test_sets"]["S1:all"]["label"]
    want = (
        "- (b) Any ML candidate beats B1 and R2 on N-sec FPR at the recall of R2 "
        f"(S1:all, {lab_b}): "
    )
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


# ---- fix round 1 -------------------------------------------------------------------------


def _sections(text):
    """{heading: lines} for every '## ' section."""
    out, name = {}, None
    for line in text.splitlines():
        if line.startswith("## "):
            name = line
            out[name] = []
        elif name:
            out[name].append(line)
    return out


def _is_header_or_rule(lines, i):
    return lines[i].startswith("|") and i + 1 < len(lines) and lines[i + 1].startswith("|---")


def _tables(lines):
    """[(preceding non-blank line, header cells, rows)] for the tables of a section."""
    out, i = [], 0
    while i < len(lines):
        if _is_header_or_rule(lines, i):
            before = next((x for x in reversed(lines[:i]) if x.strip()), "")
            cells = lambda ln: [c.strip() for c in ln.strip().strip("|").split(" | ")]  # noqa: E731
            j = i + 2
            rows = []
            while j < len(lines) and lines[j].startswith("|"):
                rows.append(cells(lines[j]))
                j += 1
            out.append((before, cells(lines[i]), rows))
            i = j
        else:
            i += 1
    return out


def test_label_beside_every_number_of_a_test_set(reported):
    _, text, m, _ = reported
    rep = load_phasec("12_report")
    secs = _sections(text)
    for name, ts in m["test_sets"].items():
        key = f"## {name} ({ts['label']})"
        lines = secs[key]
        for i, line in enumerate(lines):
            if _is_header_or_rule(lines, i) or not rep.NUMBER.findall(line):
                continue
            assert ts["label"] in line, (name, line)
    fl = secs["## Findings (quoted from findings.json)"]
    labels = {ts["label"] for ts in m["test_sets"].values()}
    for i, line in enumerate(fl):
        if _is_header_or_rule(fl, i) or not rep.NUMBER.findall(line):
            continue
        if line.startswith("The difference is"):  # sign sentence: the interval level only
            continue
        assert any(lab in line for lab in labels), line


def test_findings_b_and_c_print_sign_sentence_and_tables(reported):
    _, text, m, f = reported
    rep = load_phasec("12_report")
    sentence = rep.SIGN.format(level=rep.fmt(m["settings"]["ci_level"]))
    assert text.count(sentence) == 2
    assert "entirely above zero" in sentence
    secs = _sections(text)["## Findings (quoted from findings.json)"]
    tabs = _tables(secs)
    assert len(tabs) == 2
    assert tabs[1][1] == rep.FINDING_HEADER
    c = f["c_same_under_s2"]
    assert len(tabs[1][2]) == sum(len(e["candidates"]) for e in c["test_sets"].values())
    assert f"Found {rep.fmt(c['n_s2'])} of {rep.fmt(c['expected_s2'])} expected S2" in text


def test_interval_entirely_below_zero_prints_sentence_and_row(reported):
    import copy

    _, _, m, f = reported
    rep = load_phasec("12_report")
    f2 = copy.deepcopy(f)
    name, entry = next(iter(f2["c_same_under_s2"]["test_sets"].items()))
    cand, e = next(iter(entry["candidates"].items()))
    e["B1"].update(lo=-0.5, value=-0.3, hi=-0.1)
    text = rep.render(m, f2)
    assert "A value below zero means the ML candidate has the higher N-sec false-positive" in text
    assert "(the baseline is better)" in text
    findings = text.split("## Findings (quoted from findings.json)")[1].split("\n## ")[0]
    row = next(
        r for r in findings.splitlines() if r.startswith(f"| {name} |") and f"| {cand} |" in r
    )
    assert "-0.300 [-0.500, -0.100]" in row and row.rstrip().endswith("| no |")


def test_decision_table_lists_every_candidate_and_glossary(reported):
    _, text, m, _ = reported
    rep = load_phasec("12_report")
    secs = _sections(text)
    # the decision table follows the caveats and precedes the estimate table
    order = [
        h for h in secs if h in ("## Caveats", "## Decision table", "## Estimate or smoke test")
    ]
    assert order == ["## Caveats", "## Decision table", "## Estimate or smoke test"]
    ((_, header, rows),) = _tables(secs["## Decision table"])
    assert header == ["Test set", "Label", "Candidate", "recall", "FPR", "N-sec FPR", "ROC-AUC"]
    cands = m["settings"]["candidates"]
    assert len(rows) == len(m["test_sets"]) * len(cands)
    for name, ts in m["test_sets"].items():
        got = [r[2] for r in rows if r[0] == name]
        assert got == cands and all(r[1] == ts["label"] for r in rows if r[0] == name)
    assert "best candidate" in text and "names no best candidate" in text
    gl = "\n".join(secs["## Glossary"])
    wanted = (
        "B0 B1 R0 R1 R2 M8 M35 M8-C M35-C H V-go V-kw S1 S2 S3 FULL N-int N-sec PM-TM pos neg "
        "oof in_sample final estimate"
    ).split() + ["smoke test"]
    for code in wanted:
        assert f"- {code}: " in gl, code
    assert [c for c, _ in rep.GLOSSARY if c in wanted] == wanted or set(wanted) <= {
        c for c, _ in rep.GLOSSARY
    }
    # every code in a table header is in the glossary
    code_re = re.compile(r"^(?:[A-Z]+[0-9]*(?:-[A-Za-z]+)?|oof|final|in_sample)$")
    defined = {c for c, _ in rep.GLOSSARY}
    for lines in secs.values():
        for _, header, _ in _tables(lines):
            for cell_ in header:
                for tok in re.split(r"[\s/,]+", cell_):
                    if code_re.match(tok) and tok != "ID":
                        assert tok in defined, (tok, header)


def _expected(m, name, truth, stratum, v, c, key):
    block = m["test_sets"][name]["truth"][truth]
    mm = block["metrics"][stratum][v][c].get(key)
    rep = load_phasec("12_report")
    return rep.ci(mm)


def check_cells(text, m):
    """Mismatches between report cells and metrics.json (cell-level check, review M-7)."""
    rep = load_phasec("12_report")
    bad = []
    secs = _sections(text)
    for name, ts in m["test_sets"].items():
        lit = ts["kind"] == "literature"
        heads = rep.LITERATURE_HEADLINE if lit else rep.HEADLINE
        for tb, header, rows in _tables(secs[f"## {name} ({ts['label']})"]):
            truth = next((t for t in ("direct", "all") if tb.startswith(f"Truth `{t}`")), None)
            if truth:
                assert header[3:] == list(heads)
                for r in rows:
                    for k, got in zip(heads, r[3:], strict=True):
                        if got != _expected(m, name, truth, "all", r[1], r[0], k):
                            bad.append((name, truth, r[0], r[1], k))
            elif tb.startswith("Strata"):
                cols = header[5:]
                for r in rows:
                    for k, got in zip(cols, r[5:], strict=True):
                        if got != _expected(m, name, "direct", r[0], "V-go", r[1], k):
                            bad.append((name, "strata", r[0], r[1], k))
    ((_, _, drows),) = _tables(secs["## Decision table"])
    for r in drows:
        for k, got in zip(("recall", "fpr", None, "roc_auc"), r[3:], strict=True):
            stratum = "N-sec" if k is None else "all"
            key = "fpr" if k is None else k
            if got != _expected(m, r[0], "direct", stratum, "V-go", r[2], key):
                bad.append((r[0], "decision", r[2], key, stratum))
    return bad


def _mutated(m, fn):
    import copy

    m2 = copy.deepcopy(m)
    for ts in m2["test_sets"].values():
        for truth in ts["truth"].values():
            for stratum in truth["metrics"].values():
                fn(stratum)
    return m2


def test_cells_equal_metrics_json_and_mutations_are_caught(reported):
    _, text, m, f = reported
    rep = load_phasec("12_report")
    assert check_cells(text, m) == []
    n_cells = sum(1 for _ in rep.NUMBER.finditer(text))
    assert n_cells > 100

    def swap_variants(stratum):
        stratum["V-go"], stratum["V-kw"] = stratum["V-kw"], stratum["V-go"]

    def swap_columns(stratum):
        for per in stratum.values():
            for mm in per.values():
                if "recall" in mm and "fpr" in mm:
                    mm["recall"], mm["fpr"] = mm["fpr"], mm["recall"]

    def shift(stratum):
        for per in stratum.values():
            for mm in per.values():
                for d in mm.values():
                    if d.get("value") is not None:
                        d["value"] += 0.0123

    def swap_bounds(stratum):
        for per in stratum.values():
            for mm in per.values():
                for d in mm.values():
                    d["lo"], d["hi"] = d["hi"], d["lo"]

    for fn in (swap_variants, swap_columns, shift, swap_bounds):
        mutated_text = rep.render(_mutated(m, fn), f)
        assert check_cells(mutated_text, m), fn.__name__


def test_named_panel_agreement_context_and_prevalence_text(reported):
    _, text, m, _ = reported
    # context keys in words, species counts absent from metrics.json
    assert "- S1 fold 0 (key S1|0): " in text
    assert "FULL model (whole training pool) (key FULL|0)" in text
    assert "The counts per species (C. immitis, C. posadasii) are not in metrics.json." in text
    # prevalence columns
    assert "| 1% prevalence | 3% prevalence | 5% prevalence | 10% prevalence |" in text
    # agreement: score-source counts
    assert "| oof | final | in_sample |" in text
    counts = m["agreement"]["score_source_counts"]["proteomes"]
    set_id = next(iter(counts))
    row = next(r for r in text.splitlines() if r.startswith(f"| {set_id} |"))
    rep = load_phasec("12_report")
    for src in ("oof", "final", "in_sample"):
        want = rep.fmt(counts[set_id][src]) if src in counts[set_id] else "none"
        assert row.split(" | ")[7 + ("oof", "final", "in_sample").index(src)].strip(" |") == want
    assert "mixes out-of-fold and final-model scores" in text
    # named panel: V-go and V-kw both printed for the found proteins
    assert "Calls (V-go / V-kw):" in text and "Scores (V-go / V-kw):" in text
    found = next(p for p in m["named_panel"] if p["found"])
    c0 = m["settings"]["candidates"][0]
    want = f"{rep.fmt(found['calls'][f'{c0}|V-go'])} / {rep.fmt(found['calls'][f'{c0}|V-kw'])}"
    assert f"| {found['name']} | {want} " in text
    # plain words
    for gone in ("Phase C spec", "ruling C-", "precision of the estimate", "1 of 1"):
        assert gone not in text, gone
    for gone in (r"\bR-B\b", r"\bQ8\b"):
        assert not re.search(gone, text), gone


def test_stop_removes_an_older_report(phasec_chain, tmp_path, capsys):
    fx = pf.copy_work(phasec_chain, tmp_path)
    rep = load_phasec("12_report")
    out = fx["work"] / "phasec"
    assert rep.main(["--work-dir", str(fx["work"])]) == 0
    assert (out / "report.md").exists() and (out / "report_run.json").exists()
    path = out / "metrics.json"
    path.write_text(path.read_text().replace('"ci_level": 95', '"ci_level": 90'))
    assert rep.main(["--work-dir", str(fx["work"])]) == 2
    assert not (out / "report.md").exists() and not (out / "report_run.json").exists()
    assert "STOP" in capsys.readouterr().err


def test_null_in_a_hand_edited_input_is_a_stop(phasec_chain, tmp_path, capsys):
    import hashlib

    fx = pf.copy_work(phasec_chain, tmp_path)
    rep = load_phasec("12_report")
    out = fx["work"] / "phasec"
    assert rep.main(["--work-dir", str(fx["work"])]) == 0
    m = json.loads((out / "metrics.json").read_text())
    m["test_sets"]["S1:all"]["truth"]["direct"]["n"] = None
    (out / "metrics.json").write_text(json.dumps(m))
    run = json.loads((out / "evaluate_run.json").read_text())
    run["outputs_sha256"]["metrics.json"] = hashlib.sha256(
        (out / "metrics.json").read_bytes()
    ).hexdigest()
    (out / "evaluate_run.json").write_text(json.dumps(run))
    assert rep.main(["--work-dir", str(fx["work"])]) == 2
    err = capsys.readouterr().err
    assert err.startswith("STOP: an input has a null or a wrong type (TypeError")
    assert not (out / "report.md").exists()
