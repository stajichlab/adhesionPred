"""Task 3.2: the Calls table of the report has a fourth value column."""

from test_cli import run_cli


def test_the_calls_table_counts_not_applicable_apart(tmp_path, write_module, nodes_dmp):
    code, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    assert code == 0
    text = (out / "report.md").read_text()
    assert "| call | variant | called | not_called | not_assessable | not_applicable |" in text
    # the toy run: antigen_lookup is not applicable for ENZ1 and STAR1 (valid) and BAD1 is invalid
    row = next(ln for ln in text.splitlines() if ln.startswith("| cocci_specificity_rank_top15 |"))
    assert row == "| cocci_specificity_rank_top15 | - | 2 | 0 | 1 | 2 |"
