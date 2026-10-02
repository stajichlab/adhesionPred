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
    # the list of cell-checked tables comes from the constant that the cell tests use
    assert "; ".join(rep.CELL_CHECKED) in caveats
    assert "checked only by the typo guard" in caveats
    assert set(CHECKERS) == set(rep.CELL_CHECKED)
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
        if line.startswith("A comparator counts as beaten"):  # rule: findings.json fraction
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
    extra = "g t T-c P-gpi homology_only direct all Platt".split() + [
        "Platt scaling",
        "Brier score",
        "half-width",
        "paired interval",
        "S1:all",
        "S2-<source>:<source>",
        "S3-<clade>:clade",
        "literature",
        "ML",
        "FPR",
        "ROC-AUC",
        "PR-AUC",
        "J",
        "C",
        "grid edge",
    ]
    extra.remove("Platt")
    codes = [c for c, _ in rep.GLOSSARY]
    assert len(codes) == len(set(codes))
    assert set(codes) == set(wanted) | set(extra)  # exact content
    for code in extra:
        assert f"- {code}: " in gl, code
    # a jargon term used outside the glossary is defined in it
    outside = "\n".join(ln for h, lines in secs.items() if h != "## Glossary" for ln in lines)
    for pattern, term in (
        (r"\bPlatt\b", "Platt scaling"),
        (r"\bBrier\b", "Brier score"),
        (r"half-width", "half-width"),
        (r"\bpaired\b", "paired interval"),
        (r"\bT-c\b", "T-c"),
        (r"\bP-gpi\b", "P-gpi"),
        (r"homology_only", "homology_only"),
    ):
        if re.search(pattern, outside):
            assert f"- {term}: " in gl, term
    assert rep.OPERATING_POINT in text
    assert text.count(rep.OPERATING_POINT) == 1
    assert rep.ZERO_WIDTH in text
    # every code in a table header is in the glossary
    code_re = re.compile(r"^(?:[A-Z]+[0-9]*(?:-[A-Za-z]+)?|oof|final|in_sample)$")
    defined = {c for c, _ in rep.GLOSSARY}
    for lines in secs.values():
        for _, header, _ in _tables(lines):
            for cell_ in header:
                for tok in re.split(r"[\s/,]+", cell_):
                    if code_re.match(tok) and tok != "ID":
                        assert tok in defined, (tok, header)


# ---- cell-level checks (review M-7) ----
# The expected text is built here from the raw values of metrics.json and findings.json with
# its own formatter. It does not call the report's fmt() or ci().


def f3(x):
    return "n/a" if x is None else f"{x:.3f}"


def n3(x):
    return f"{x:,}"


def iv(d):
    if d is None or d.get("value") is None:
        return "n/a"
    return f"{d['value']:.3f} [{d['lo']:.3f}, {d['hi']:.3f}]"


def ivn(d, n_res):
    """iv() plus the defined-resample count when it is below the number of resamples."""
    text = iv(d)
    if text != "n/a" and d.get("n_defined") is not None and d["n_defined"] < n_res:
        text += f" (n_defined {d['n_defined']:,} of {n_res:,})"
    return text


def yn(x):
    return "yes" if x else "no"


def _truth_block(m, name, truth="direct"):
    return m["test_sets"][name]["truth"][truth]


def _met(block, stratum, v, c, k):
    return block["metrics"].get(stratum, {}).get(v, {}).get(c, {}).get(k)


def _sec(text, name, ts):
    return _sections(text)[f"## {name} ({ts['label']})"]


def check_decision(text, m, f):
    bad = []
    ((_, header, rows),) = _tables(_sections(text)["## Decision table"])
    assert header == ["Test set", "Label", "Candidate", "recall", "FPR", "N-sec FPR", "ROC-AUC"]
    for r in rows:
        b = _truth_block(m, r[0])
        want = [
            iv(_met(b, "all", "V-go", r[2], "recall")),
            iv(_met(b, "all", "V-go", r[2], "fpr")),
            iv(_met(b, "N-sec", "V-go", r[2], "fpr")),
            iv(_met(b, "all", "V-go", r[2], "roc_auc")),
        ]
        if r[3:] != want:
            bad.append(("decision", r[0], r[2]))
    return bad


def check_metric_tables(text, m, f):
    bad = []
    head = [
        "recall",
        "precision",
        "fpr",
        "roc_auc",
        "pr_auc",
        "precision_at_recall_0.8",
        "precision_at_recall_0.9",
        "recall_at_fpr_0.01",
    ]
    for name, ts in m["test_sets"].items():
        heads = ["recall"] if ts["kind"] == "literature" else head
        for tb, header, rows in _tables(_sec(text, name, ts)):
            truth = next((t for t in ("direct", "all") if tb.startswith(f"Truth `{t}`")), None)
            if not truth:
                continue
            assert header == ["Candidate", "Variant", "Label", *heads]
            b = _truth_block(m, name, truth)
            for r in rows:
                want = [iv(_met(b, "all", r[1], r[0], k)) for k in heads]
                if r[3:] != want:
                    bad.append(("metric", name, truth, r[0], r[1]))
    return bad


def check_strata(text, m, f):
    bad = []
    for name, ts in m["test_sets"].items():
        cols = ["recall"] if ts["kind"] == "literature" else ["recall", "fpr", "roc_auc"]
        for tb, header, rows in _tables(_sec(text, name, ts)):
            if not tb.startswith("Strata"):
                continue
            assert header[5:] == cols
            b = _truth_block(m, name)
            for r in rows:
                n = b["n"][r[0]]
                want = [n3(n["pos"]), n3(n["neg"])] + [
                    iv(_met(b, r[0], "V-go", r[1], k)) for k in cols
                ]
                if r[3:] != want:
                    bad.append(("strata", name, r[0], r[1]))
    return bad


def check_variant_effect(text, m, f):
    bad = []
    for name, ts in m["test_sets"].items():
        for tb, header, rows in _tables(_sec(text, name, ts)):
            if not tb.startswith("Variant effect"):
                continue
            cols = header[2:]
            assert cols == (["recall"] if ts["kind"] == "literature" else
                            ["recall", "fpr", "roc_auc", "pr_auc"])  # fmt: skip
            b = _truth_block(m, name)
            for r in rows:
                eff = b["variant_effect"]["all"][r[0]]
                if r[2:] != [iv(eff.get(k)) for k in cols]:
                    bad.append(("variant_effect", name, r[0]))
    return bad


def check_vs_rule(text, m, f):
    bad = []
    for name, ts in m["test_sets"].items():
        for tb, _, rows in _tables(_sec(text, name, ts)):
            if not tb.startswith("ML against the rule"):
                continue
            b = _truth_block(m, name)
            for r in rows:
                vr = b["vs_rule"][r[1]][r[0]]
                if r[3:] != [iv(vr["precision_at_rule_recall"]), iv(vr["recall_at_rule_fpr"])]:
                    bad.append(("vs_rule", name, r[0], r[1]))
    return bad


def check_prevalence(text, m, f):
    bad = []
    for name, ts in m["test_sets"].items():
        for tb, header, rows in _tables(_sec(text, name, ts)):
            if not tb.startswith("Precision at assumed prevalence"):
                continue
            b = _truth_block(m, name)
            for r in rows:
                pv = b["prevalence"]["V-go"][r[0]]
                keys = list(pv)
                assert header[2:] == [f"{float(k) * 100:g}% prevalence" for k in keys]
                if r[2:] != [iv(pv[k]) for k in keys]:
                    bad.append(("prevalence", name, r[0]))
    return bad


def check_brier(text, m, f):
    bad = []
    for name, ts in m["test_sets"].items():
        for _, header, rows in _tables(_sec(text, name, ts)):
            if header == ["Candidate", "Variant", "Label", "Brier"]:
                for r in rows:
                    if r[3] != f3(ts["calibration"][r[1]][r[0]]["brier"]):
                        bad.append(("brier", name, r[0], r[1]))
    return bad


def check_reliability(text, m, f):
    bad = []
    for name, ts in m["test_sets"].items():
        for tb, header, rows in _tables(_sec(text, name, ts)):
            if not tb.startswith("Reliability"):
                continue
            for r in rows:
                bins = ts["calibration"][r[1]][r[0]]["bins"]
                if header[3:] != [f"{f3(x['lo'])} to {f3(x['hi'])}" for x in bins]:
                    bad.append(("reliability header", name, r[0], r[1]))
                if r[3:] != [f"{f3(x['frac_pos'])} ({n3(x['n'])})" for x in bins]:
                    bad.append(("reliability", name, r[0], r[1]))
    return bad


def check_agreement(text, m, f):
    bad = []
    ag = m["agreement"]
    for sec_name, lines in _sections(text).items():
        if not sec_name.startswith("## Agreement"):
            continue
        for _, header, rows in _tables(lines):
            assert header[3:] == ["both", "rule only", "ML only", "neither", "oof", "final",
                                  "in_sample"]  # fmt: skip
            proteome = header[0] == "Set"
            for r in rows:
                grp = ag["proteomes"] if proteome else ag["truth"]
                k = grp[r[0]][r[2]][r[1]]
                have = ag["score_source_counts"]["proteomes" if proteome else "truth"][r[0]]
                want = [n3(k["rule1_ml1"]), n3(k["rule1_ml0"]), n3(k["rule0_ml1"]),
                        n3(k["rule0_ml0"])] + [n3(have[s]) if s in have else "none"
                                               for s in ("oof", "final", "in_sample")]  # fmt: skip
                if r[3:] != want:
                    bad.append(("agreement", r[0], r[1], r[2]))
    return bad


def check_findings(text, m, f):
    bad = []
    tabs = _tables(_sections(text)["## Findings (quoted from findings.json)"])
    b = f["b_ml_beats_b1_and_r2_on_nsec_s1"]
    wanted = [(b["test_set"], c, e) for c, e in b["candidates"].items()]
    for name, entry in f["c_same_under_s2"]["test_sets"].items():
        wanted += [(name, c, e) for c, e in entry["candidates"].items()]
    rows = [r for _, _, tab in tabs for r in tab]
    assert len(rows) == len(wanted)
    for r, (name, c, e) in zip(rows, wanted, strict=True):
        fp = _truth_block(m, name)["nsec_fpr_at_rule_recall"]["V-go"]["fpr"]
        n = m["settings"]["n_resamples"]
        want = [name, m["test_sets"][name]["label"], c, ivn(fp["B1"], n), ivn(fp["R2"], n),
                ivn(fp[c], n), ivn(e["B1"], n), ivn(e["R2"], n), yn(e["beats_B1"]),
                yn(e["beats_R2"]), yn(e["holds"])]  # fmt: skip
        if r != want:
            bad.append(("findings", name, c))
    return bad


def check_named_panel(text, m, f):
    bad = []
    cands = m["settings"]["candidates"]
    found = {p["name"]: p for p in m["named_panel"] if p["found"]}
    for tb, header, rows in _tables(_sections(text)["## Named panel"]):
        field = {"Calls (V-go / V-kw):": "calls", "Scores (V-go / V-kw):": "scores"}.get(tb)
        if not field:
            continue
        assert header == ["Protein", *cands]
        assert [r[0] for r in rows] == list(found)
        for r in rows:
            p = found[r[0]]
            want = []
            for c in cands:
                pair = []
                for v in ("V-go", "V-kw"):
                    x = p[field].get(f"{c}|{v}")
                    pair.append(f3(x) if field == "scores" else ("n/a" if x is None else str(x)))
                want.append(" / ".join(pair))
            if r[1:] != want:
                bad.append(("panel", field, r[0]))
    return bad


CHECKERS = {
    "decision table": check_decision,
    "metric tables of each test set (direct and all truth)": check_metric_tables,
    "strata tables": check_strata,
    "variant-effect tables": check_variant_effect,
    "ML-against-rule tables": check_vs_rule,
    "prevalence tables": check_prevalence,
    "Brier tables": check_brier,
    "reliability tables": check_reliability,
    "agreement tables": check_agreement,
    "findings (b) and (c) tables": check_findings,
    "named-panel call and score tables": check_named_panel,
    "fitted-settings table": None,  # set below
}


def check_fitted(text, m, f):
    """Final review I-1, I-3: one row per unit and candidate, values from metrics.json."""
    bad = []
    ((tb, header, rows),) = _tables(_sections(text)["## Fitted settings"])
    assert tb.startswith("Fitted settings per unit")
    assert header == ["Unit", "Variant", "Candidate", "g", "t", "J", "C", "H variant",
                      "Threshold", "Platt a", "Platt b", "Convergence warnings", "At grid edge"]  # fmt: skip
    fs = m["fitted_settings"]
    assert len(rows) == sum(len(per) for per in fs.values())
    for r in rows:
        e = fs[f"{r[0]}|{r[1]}"][r[2]]
        cw = e.get("convergence_warnings")
        edges = [k for k in ("C", "g", "t") if e["at_grid_edge"].get(k)]
        want = [e["g"] if e.get("g") else "n/a", f3(e.get("t")), f3(e.get("j")), f3(e.get("C")),
                e.get("h_variant") or "n/a", f3(e.get("threshold")), f3(e.get("platt_a")),
                f3(e.get("platt_b")), "n/a" if cw is None else n3(cw),
                ", ".join(edges) or "no"]  # fmt: skip
        if r[3:] != want:
            bad.append(("fitted", r[0], r[1], r[2]))
    return bad


CHECKERS["fitted-settings table"] = check_fitted


@pytest.mark.parametrize("name", list(CHECKERS))
def test_cells_equal_json_values(reported, name):
    _, text, m, f = reported
    assert CHECKERS[name](text, m, f) == []


def test_checker_names_equal_the_report_constant():
    assert list(CHECKERS) == list(load_phasec("12_report").CELL_CHECKED)


def test_named_panel_agreement_context_and_prevalence_text(reported):
    _, text, m, _ = reported
    # context keys in words, species counts absent from metrics.json
    assert "- S1 fold 0 (key S1|0): " in text
    assert "FULL model (whole training pool) (key FULL|0)" in text
    assert "Each count covers the C. immitis and C. posadasii rows together." in text
    assert "not in metrics.json" not in text
    assert "1 rows" not in text and "(key S1|4): 1 row." in text
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


# ---- final fix wave (final review I-1, I-2, I-3) ----


def test_fitted_settings_section_and_edge_line(reported):
    _, text, m, _ = reported
    rep = load_phasec("12_report")
    sec = "\n".join(_sections(text)["## Fitted settings"])
    ec = m["grid_edge_counts"]
    line = (
        "Settings at a grid edge: "
        + "; ".join(f"{k} {ec[k]['at_edge']:,} of {ec[k]['n']:,}" for k in ("C", "g", "t"))
        + "."
    )
    assert line in sec
    assert rep.EDGE_SENTENCE in sec
    assert "the data may prefer a value outside the grid; the report cannot say." in sec
    # the counts are the counts of the flags
    fs = m["fitted_settings"]
    for k in ("C", "g", "t"):
        flags = [e["at_grid_edge"][k] for per in fs.values() for e in per.values()
                 if k in e["at_grid_edge"]]  # fmt: skip
        assert ec[k] == {"at_edge": sum(flags), "n": len(flags)}, k
    # H's variant per unit is in the table
    rows = _tables(_sections(text)["## Fitted settings"])[0][2]
    h = [r for r in rows if r[2] == "H"]
    assert h and len(h) == len(fs)
    assert all(r[7] == fs[f"{r[0]}|{r[1]}"]["H"]["h_variant"] for r in h)


def test_rule_subset_sentence_beside_the_decision_table(reported):
    _, text, _, _ = reported
    rep = load_phasec("12_report")
    dec = "\n".join(_sections(text)["## Decision table"])
    assert rep.RULE_SUBSET in dec and text.count(rep.RULE_SUBSET) == 1
    assert "see Fitted settings" in rep.RULE_SUBSET and "## Fitted settings" in text


def test_edge_flag_mutation_is_seen_by_the_cell_check(reported):
    import copy

    _, text, m, f = reported
    rep = load_phasec("12_report")
    m2 = copy.deepcopy(m)
    unit = next(iter(m2["fitted_settings"]))
    e = m2["fitted_settings"][unit]["M8"]
    e["at_grid_edge"]["C"] = not e["at_grid_edge"]["C"]
    assert check_fitted(rep.render(m2, f), m, f) != []


def _finding_row(text, cand):
    sec = _sections(text)["## Findings (quoted from findings.json)"]
    return next(r for r in sec if r.startswith("| S1:all |") and f"| {cand} |" in r)


def test_finding_rows_print_n_defined_below_the_resample_count(reported):
    # final review M-3: a hand edit puts one difference at 90% defined resamples
    import copy

    _, _, m, f = reported
    rep = load_phasec("12_report")
    n = m["settings"]["n_resamples"]
    f2 = copy.deepcopy(f)
    b = f2["b_ml_beats_b1_and_r2_on_nsec_s1"]
    cand, e = next(iter(b["candidates"].items()))
    e["B1"].update(value=0.2, lo=0.05, hi=0.3, n_defined=int(0.9 * n))
    text = rep.render(m, f2)
    row = _finding_row(text, cand)
    assert f"0.200 [0.050, 0.300] (n_defined {int(0.9 * n):,} of {n:,})" in row
    e["B1"]["n_defined"] = n
    text = rep.render(m, f2)
    row = _finding_row(text, cand)
    assert "0.200 [0.050, 0.300] |" in row and "0.200 [0.050, 0.300] (n_defined" not in row
    rule = rep.DEFINED_RULE.format(fraction=f"{b['min_defined_fraction']:.3f}")
    assert text.count(rule) == 2
