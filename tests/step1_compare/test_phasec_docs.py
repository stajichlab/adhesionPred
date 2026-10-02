"""COLUMNS.md and README.md document the Phase C outputs exactly as the code writes them.

Same rules as test_phaseb_docs.py: the set of names in a table section must equal the code
tuple; the bullet list after a `**Keys...:**` marker must equal the keys of a real output.
Grids and counts in the text (C grid, t grid, count floor, number of S2 test sets) must equal
the code constants.
"""

import gzip
import json
import re

import pytest

pytest.importorskip("numpy")
pytest.importorskip("sklearn")

import dedupe  # noqa: E402
import paths  # noqa: E402
import splits  # noqa: E402
from conftest import load_phasec  # noqa: E402
from test_phaseb_docs import COLUMNS, README, key_names, section, table_names  # noqa: E402

SCRIPTS = ("08_build_eval_tables", "09_make_splits", "10_fit_and_score", "11_evaluate", "12_report")
EXTRA_HEADINGS = ("phasec/logs/mmseqs.log", "phasec/logs/wall.<job>.txt")
RUN_JSON_KEYS = "**Keys of `input_sha256`:**"


def _json(work, name):
    return json.loads((work / "phasec" / name).read_text())


def _bullet(heading, key):
    lines = section(heading).splitlines()
    return next(x for x in lines if x.startswith(f"- `{key}` "))


def test_phasec_outputs_have_columns_md_headings():
    names = [f"phasec/{n}" for s in SCRIPTS for n in load_phasec(s).OUTPUT_NAMES]
    names += list(EXTRA_HEADINGS)
    headings = [line[3:].split()[0] for line in COLUMNS.splitlines() if line.startswith("## ")]
    missing = [n for n in names if n not in headings]
    assert not missing, f"COLUMNS.md has no heading for {missing}"


def test_table_columns_equal_the_code_constants():
    s08 = load_phasec("08_build_eval_tables")
    s10 = load_phasec("10_fit_and_score")
    expected = {
        "phasec/eval_table.tsv.gz": dedupe.TABLE_COLUMNS,
        "phasec/eval_literature.tsv": s08.LITERATURE_COLUMNS,
        "phasec/eval_dedupe_log.tsv": dedupe.LOG_COLUMNS,
        "phasec/clusters.tsv.gz": splits.CLUSTER_COLUMNS,
        "phasec/split_members.tsv.gz": splits.MEMBER_COLUMNS,
        "phasec/tc_removed.tsv": splits.REMOVED_COLUMNS,
        "phasec/max_identity.tsv.gz": splits.IDENTITY_COLUMNS,
        "phasec/scores.tsv.gz": s10.SCORE_COLUMNS,
    }
    for heading, columns in expected.items():
        assert table_names(heading) == set(columns), heading


def test_documented_columns_equal_the_written_headers(phasec_chain):
    """The code constants are the writers' constants only if the files carry them."""
    w = phasec_chain["work"] / "phasec"

    def header(name):
        raw = (w / name).read_bytes()
        text = gzip.decompress(raw).decode() if name.endswith(".gz") else raw.decode()
        return text.splitlines()[0].split("\t")

    for name in ("eval_table.tsv.gz", "eval_literature.tsv", "eval_dedupe_log.tsv",
                 "clusters.tsv.gz", "split_members.tsv.gz", "tc_removed.tsv",
                 "max_identity.tsv.gz", "scores.tsv.gz"):  # fmt: skip
        assert set(header(name)) == table_names(f"phasec/{name}"), name


def test_proteome_calls_columns_equal_the_code(phasec_chain):
    s11 = load_phasec("11_evaluate")
    patterns = {"call_<candidate>_<variant>", "score_<candidate>_<variant>"}
    got = table_names("phasec/proteome_calls.tsv.gz")
    assert got == set(s11.PROTEOME_BASE) | patterns
    m = _json(phasec_chain["work"], "metrics.json")
    cols = s11.proteome_columns(tuple(m["settings"]["candidates"]))
    raw = (phasec_chain["work"] / "phasec" / "proteome_calls.tsv.gz").read_bytes()
    header = gzip.decompress(raw).decode().splitlines()[0].split("\t")
    assert header == list(cols)
    dynamic = set(cols) - set(s11.PROTEOME_BASE)
    assert {c.split("_")[0] + "_<candidate>_<variant>" for c in dynamic} == patterns


def test_run_json_keys_equal_the_code(phasec_chain, tmp_path):
    import phasec_fixture as pf

    fx = pf.copy_work(phasec_chain, tmp_path)
    assert load_phasec("12_report").main(["--work-dir", str(fx["work"])]) == 0
    w = fx["work"]
    build = _json(w, "build_run.json")
    assert key_names("phasec/build_run.json") == set(build)
    assert key_names("phasec/build_run.json", RUN_JSON_KEYS) == set(build["input_sha256"])
    splits_run = _json(w, "splits_run.json")
    assert key_names("phasec/splits_run.json") == set(splits_run)
    assert key_names("phasec/splits_run.json", RUN_JSON_KEYS) == set(splits_run["input_sha256"])
    scores = _json(w, "scores_run.json")
    assert key_names("phasec/scores_run.json") == set(scores)
    assert key_names("phasec/scores_run.json", RUN_JSON_KEYS) == set(scores["input_sha256"])
    unit = scores["units"]["S1|0|V-go"]
    rule = "**Keys of a rule object in `units`:**"
    assert key_names("phasec/scores_run.json", rule) == set(unit["R2"]) == set(unit["R0"])
    lr = "**Keys of a logistic-regression object in `units`:**"
    assert key_names("phasec/scores_run.json", lr) == set(unit["M8"]) == set(unit["H"])
    assert set(unit["B0"]) == set(unit["B1"]) == set(unit["M8"])
    evaluate = _json(w, "evaluate_run.json")
    assert key_names("phasec/evaluate_run.json") == set(evaluate)
    assert key_names("phasec/evaluate_run.json", RUN_JSON_KEYS) == set(evaluate["input_sha256"])
    report = _json(w, "report_run.json")
    assert key_names("phasec/report_run.json") == set(report)
    assert key_names("phasec/report_run.json", RUN_JSON_KEYS) == set(report["input_sha256"])


def test_run_json_documents_every_outputs_sha256_name(phasec_chain, tmp_path):
    import phasec_fixture as pf

    fx = pf.copy_work(phasec_chain, tmp_path)
    assert load_phasec("12_report").main(["--work-dir", str(fx["work"])]) == 0
    marker = "**Keys of `outputs_sha256`:**"
    for name in ("build_run", "splits_run", "scores_run", "evaluate_run", "report_run"):
        got = key_names(f"phasec/{name}.json", marker)
        assert got == set(_json(fx["work"], f"{name}.json")["outputs_sha256"]), name


def test_metrics_and_findings_keys_equal_the_code(phasec_chain):
    m = _json(phasec_chain["work"], "metrics.json")
    f = _json(phasec_chain["work"], "findings.json")
    assert key_names("phasec/metrics.json") == set(m)
    assert key_names("phasec/findings.json") == set(f)
    s2 = m["test_sets"]["S2-Spom_PomBase:Spom_PomBase"]
    s1 = m["test_sets"]["S1:all"]
    assert key_names("phasec/metrics.json", "**Keys of a test set object:**") == set(s2)
    assert set(s1) == set(s2) - {"calibration"}  # calibration only for S2 test sets
    assert key_names("phasec/metrics.json", "**Keys of a truth object:**") == set(
        s1["truth"]["direct"]
    )
    assert key_names("phasec/metrics.json", "**Keys of `settings`:**") == set(m["settings"])
    assert key_names("phasec/metrics.json", "**Keys of `agreement`:**") == set(m["agreement"])
    panel = {bool(e["found"]): e for e in m["named_panel"]}
    assert set(panel) == {True, False}, "the fixture needs a found and a missing panel entry"
    assert key_names("phasec/metrics.json", "**Keys of a named panel object:**") == set(
        panel[False]
    )
    added = key_names("phasec/metrics.json", "**Keys added when `found` is `true`:**")
    assert added == set(panel[True]) - set(panel[False])
    assert f"`{m['schema']}`" in section("phasec/metrics.json")
    assert f"`{f['schema']}`" in section("phasec/findings.json")
    for key, marker in (
        ("a_b1_not_saturated", "**Keys of `a_b1_not_saturated`:**"),
        ("b_ml_beats_b1_and_r2_on_nsec_s1", "**Keys of `b_ml_beats_b1_and_r2_on_nsec_s1`:**"),
        ("c_same_under_s2", "**Keys of `c_same_under_s2`:**"),
    ):
        assert key_names("phasec/findings.json", marker) == set(f[key]), key


def test_fitted_settings_keys_equal_the_code(phasec_chain):
    # final review I-1, I-3: metrics.json copies the fitted settings with grid-edge flags
    m = _json(phasec_chain["work"], "metrics.json")
    units = _json(phasec_chain["work"], "scores_run.json")["units"]
    fs = m["fitted_settings"]
    assert set(fs) == set(units)
    heading = "phasec/metrics.json"
    rule = key_names(heading, "**Keys of a rule object in `fitted_settings`:**")
    lr = key_names(heading, "**Keys of a logistic-regression object in `fitted_settings`:**")
    for unit, per in fs.items():
        assert set(per) == set(units[unit])
        for c, e in per.items():
            assert set(e) == (rule if c in ("R0", "R1", "R2") else lr), (unit, c)
            src = {k: v for k, v in units[unit][c].items() if k != "rule_grid"}
            assert {k: v for k, v in e.items() if k != "at_grid_edge"} == src, (unit, c)
    edge = {k for per in fs.values() for e in per.values() for k in e["at_grid_edge"]}
    marker = "**Keys of `at_grid_edge` (only the parameters that are not null for the candidate):**"
    assert key_names(heading, marker) == edge
    assert key_names(heading, "**Keys of `grid_edge_counts`:**") == set(m["grid_edge_counts"])


def test_at_grid_edge_values_in_columns_md_equal_the_code():
    import models
    import rules

    s11 = load_phasec("11_evaluate")
    assert s11.GRID_EDGES == {
        "C": (models.C_GRID[0], models.C_GRID[-1]),
        "g": (rules.G_VALUES[0], rules.G_VALUES[-1]),
        "t": (rules.T_VALUES[0], rules.T_VALUES[-1]),
    }
    c_line = _bullet("phasec/metrics.json", "C")
    assert f"{models.C_GRID[0]:g} or {models.C_GRID[-1]:g}" in c_line
    g_line = _bullet("phasec/metrics.json", "g")
    assert f"`{rules.G_VALUES[0]}` or `{rules.G_VALUES[-1]}`" in g_line
    t_line = _bullet("phasec/metrics.json", "t")
    assert f"{rules.T_VALUES[0]:g} or {rules.T_VALUES[-1]:g}" in t_line


def test_strata_in_columns_md_equal_the_code(phasec_chain):
    # Task 12 Important 1: the strata of metrics.json, set equality with STRATA
    s11 = load_phasec("11_evaluate")
    marker = (
        "**Strata (keys of `n` and of `metrics` of a truth object; `STRATA` in `11_evaluate.py`):**"
    )
    assert key_names("phasec/metrics.json", marker) == set(s11.STRATA)
    m = _json(phasec_chain["work"], "metrics.json")
    for name, ts in m["test_sets"].items():
        got = set(ts["truth"]["direct"]["n"])
        want = set(s11.STRATA) - ({"identity_below_0.3"} if name.startswith("S1:") else set())
        assert got == want == set(ts["truth"]["direct"]["metrics"]), name
    long_line = _bullet("phasec/metrics.json", "long")
    assert f"{s11.LONG_CUTOFF:,} aa" in long_line


def test_context_key_in_columns_md_equals_the_code(phasec_chain):
    s11 = load_phasec("11_evaluate")
    m = _json(phasec_chain["work"], "metrics.json")
    assert key_names("phasec/metrics.json", "**Keys of `context`:**") == set(m["context"])
    line = _bullet("phasec/metrics.json", "onygenales_tc_rows_in_vkw_training")
    for taxon in s11.ONYGENALES_TC_TAXA:
        assert taxon in line


def test_findings_candidate_keys_equal_the_code(phasec_chain):
    f = _json(phasec_chain["work"], "findings.json")
    one = next(iter(f["b_ml_beats_b1_and_r2_on_nsec_s1"]["candidates"].values()))
    marker = "**Keys of a candidate object in `b_ml_beats_b1_and_r2_on_nsec_s1`:**"
    assert key_names("phasec/findings.json", marker) == set(one)


def test_readme_names_the_phase_c_commands_and_rules():
    for text in ("phasec/08_build_eval_tables.py", '"$S1/phasec/c1_evaluate.sh"', "phasec/12_report.py",
                 "C1_STEPS", "STEP1_MMSEQS", "Illegal", "tc_taxon_clades.tsv", "wall_seconds"):  # fmt: skip
        assert text in README, text
    job = (paths.STEP1_DIR / "phasec" / "c1_evaluate.sh").read_text()
    assert "C1_STEPS" in job and "wall_seconds" in job
    assert '`input_sha256["unique_sequences.tsv.gz"]` must equal `embedding_run.json`' in README
    rule = "A job that runs an AVX2 tool must request a node feature that has AVX2"
    assert rule in " ".join(README.split()) and rule in " ".join(job.split())


def test_readme_job_resources_equal_the_job_script():
    job = (paths.STEP1_DIR / "phasec" / "c1_evaluate.sh").read_text()
    flat = " ".join(README.split())
    for sbatch, shown in (("-c 16", "-c 16"), ("--mem=32G", "--mem=32G"),
                          ("--time=4:00:00", "--time=4:00:00"),
                          ("-p epyc", "partition epyc"), ("--constraint=ryzen", "--constraint=ryzen")):  # fmt: skip
        assert f"#SBATCH {sbatch}" in job, sbatch
        assert shown in flat, shown
    assert "review estimates, not measurements" in flat
    assert "No Phase C job has been run" in flat


def test_readme_names_the_09_log_and_the_part_list():
    assert "phasec/logs/mmseqs.log" in README
    assert "`test,test_lit`" in section("phasec/scores.tsv.gz")


def test_phase_c_documents_do_not_use_the_forbidden_script_path_variable():
    """The Phase A README names the variable once, to forbid it; Phase C text must not."""
    forbidden = "BASH" + "_SOURCE"  # split so that this file does not contain the word itself
    columns = COLUMNS.split("# Phase C outputs (evaluation)", 1)
    readme = README.split("## Phase C: rule versus ML evaluation", 1)
    assert len(columns) == 2 and len(readme) == 2, "the Phase C sections are missing"
    assert forbidden not in columns[1] and forbidden not in readme[1]


def test_grids_in_columns_md_equal_the_code():
    import findings
    import models
    import rules

    c_line = _bullet("phasec/scores_run.json", "c_grid")
    values = c_line.split(":", 1)[1].split("(")[0].split(",")
    assert tuple(float(x) for x in values) == models.C_GRID  # ruling C-12
    t_line = _bullet("phasec/scores_run.json", "t")
    lo, hi = (float(x) for x in t_line.split("cut, ")[1].split(" in steps")[0].split(" to "))
    assert (lo, hi) == (rules.T_VALUES[0], rules.T_VALUES[-1])
    step = float(re.search(r"in steps of ([0-9.]+)", t_line).group(1))
    gaps = {round(b - a, 10) for a, b in zip(rules.T_VALUES, rules.T_VALUES[1:], strict=False)}
    assert gaps == {step}
    floor = _bullet("phasec/metrics.json", "floor_met")
    assert f"at least {findings.MIN_DIRECT_POSITIVES} " in floor
    label = _bullet("phasec/metrics.json", "label")
    assert "`floor_met`" in label
    settings = _bullet("phasec/metrics.json", "estimate_min_direct_positives")
    assert f"{findings.MIN_DIRECT_POSITIVES}" in settings


def test_s2_count_and_resample_count_in_columns_md_equal_the_code(phasec_chain):
    import bootstrap
    import findings

    m = _json(phasec_chain["work"], "metrics.json")
    n_s2 = _bullet("phasec/findings.json", "n_s2")
    expected = _bullet("phasec/findings.json", "expected_s2")
    assert f"{findings.N_S2} S2" in expected
    assert "`n_s2`" in n_s2
    assert f"{m['settings']['seed']}" in _bullet("phasec/scores_run.json", "seed")
    assert m["settings"]["n_resamples"] == 50  # the fixture run; the default is the constant
    assert f"{bootstrap.N_RESAMPLES:,} " in section("phasec/metrics.json")


# --- parser self-tests: a stale or missing name must make the comparison fail ---


def test_table_parser_sees_a_missing_row_and_a_stale_name():
    row = next(ln for ln in COLUMNS.splitlines() if ln.startswith("| seq_sha256 | Hash of the seq"))
    heading = "phasec/eval_table.tsv.gz"
    assert "seq_sha256" not in table_names(heading, COLUMNS.replace(row + "\n", "", 1))
    stale = COLUMNS.replace(row, row + "\n| old_column | Gone. |", 1)
    assert table_names(heading, stale) != set(dedupe.TABLE_COLUMNS)
    assert table_names(heading) == set(dedupe.TABLE_COLUMNS)


def test_key_parser_sees_a_missing_key_and_an_extra_key():
    line = next(ln for ln in COLUMNS.splitlines() if ln.startswith("- `labelled_genes_without"))
    heading = "phasec/build_run.json"
    assert "labelled_genes_without_sequence" in key_names(heading)
    missing = COLUMNS.replace(line + "\n", "", 1)
    assert "labelled_genes_without_sequence" not in key_names(heading, text=missing)
    extra = COLUMNS.replace(line, line + "\n- `not_written_by_the_code`", 1)
    assert "not_written_by_the_code" in key_names(heading, text=extra)
