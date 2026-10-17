"""scripts/sorting_hat/refresh_3_10/compare_old_new_rule.py on a synthetic out_final (no HPCC data)."""

import csv
import gzip
import importlib.util
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[2] / "scripts/sorting_hat/refresh_3_10/compare_old_new_rule.py"
)
spec = importlib.util.spec_from_file_location("compare_old_new_rule", SCRIPT)
cmp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cmp)

# protein: (signal, repeat, broad, narrow, new, surface)
C, N, A, X = "called", "not_called", "not_assessable", "not_applicable"
PROTEINS = {
    "both": (C, C, C, C, C, C),  # called by both rules
    "broad_only": (C, N, C, N, N, C),  # the old rule calls it, the new rule does not
    "no_sp": (N, N, C, N, N, N),
    "unknown_sp": (A, C, N, N, A, A),  # not_assessable: neither called nor not_called
    "na_leaf": (C, X, N, N, N, N),  # not_applicable repeat acts as not_called
    "na_call": (X, N, N, N, X, X),  # not_applicable result
}
FAMILIES = {"both": "PF1", "broad_only": "PF2,PF3", "no_sp": "PF2"}
CALLS = [
    "signal_peptide_protein",
    "tandem_repeat_protein",
    "wall_family_domain",
    "wall_family_narrow_domain",
    "cell_wall_adhesion_candidate",
    "surface_attachment_candidate",
]


def _write(path, header, rows):
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def make_out(tmp_path, long=True, wide=False, evidence=True, drop=None):
    out = tmp_path / "Prot_A" / "out_final"
    out.mkdir(parents=True)
    rows, wide_rows = [], []
    for pid, vals in PROTEINS.items():
        wrow = [pid]
        for call, value in zip(CALLS, vals, strict=True):
            variant = "R0" if call in CALLS[0:1] + CALLS[4:] else ""
            if call != drop:
                rows.append([pid, call, variant, value, "unvalidated", "", "", ""])
            wrow += [value, "unvalidated"]
        wide_rows.append(wrow)
    if long:
        _write(
            out / "calls.long.tsv.gz",
            cmp.LONG_NEEDED + ["status", "status_basis", "other_basis", "derived_from"],
            rows,
        )
    if wide:
        head = ["protein"]
        for call in CALLS:
            name = f"{call}[R0]" if call in CALLS[0:1] + CALLS[4:] else call
            head += [name, name + "_status"]
        _write(out / "calls.wide.tsv.gz", head, wide_rows)
    if evidence:
        ev = [[p, "pfam_adhesion", "families", f] for p, f in FAMILIES.items()]
        ev.append(["both", "tm", "n_tm", "0"])
        _write(out / "evidence.tsv.gz", ["protein", "module", "field", "value"], ev)
    return out


def run(out, tmp_path, capsys, *extra):
    res = tmp_path / "counts.tsv"
    code = cmp.main([str(out), "--out", str(res), *extra])
    return code, res, capsys.readouterr()


def read_rows(path):
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def test_counts_and_not_values_apart(tmp_path, capsys):
    out = make_out(tmp_path)
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    code, res, _ = run(out, tmp_path, capsys)
    assert code == 0
    (row,) = read_rows(res)
    assert row["proteome"] == "Prot_A" and row["n_proteins"] == "6"
    assert (row["new_called"], row["new_not_called"]) == ("1", "3")
    assert (row["new_not_assessable"], row["new_not_applicable"]) == ("1", "1")
    # old rule: both, broad_only called; unknown_sp unknown; na_leaf not_called (X acts as N)
    assert (row["old_called"], row["old_not_assessable"], row["old_not_applicable"]) == (
        "2",
        "1",
        "0",
    )
    assert (row["old_only"], row["new_only"], row["rule_mismatch"]) == ("1", "0", "0")
    assert row["surface_called"] == "2" and row["surface_not_applicable"] == "1"
    assert res.read_text().startswith("#")
    assert before == {p.name: p.read_bytes() for p in out.iterdir()}  # input untouched


def test_family_split(tmp_path, capsys):
    out = make_out(tmp_path)
    _, res, _ = run(out, tmp_path, capsys)
    fam = read_rows(res.with_name("counts.tsv.by_family.tsv"))
    assert {(r["family"], r["n_proteins"]) for r in fam} == {("PF2", "1"), ("PF3", "1")}


def test_wide_only_gives_same_counts_and_no_evidence_skips_split(tmp_path, capsys):
    out = make_out(tmp_path, long=False, wide=True, evidence=False)
    code, res, captured = run(out, tmp_path, capsys)
    assert code == 0
    (row,) = read_rows(res)
    assert (row["new_called"], row["old_called"], row["old_only"]) == ("1", "2", "1")
    assert "family split skipped" in captured.err
    assert not res.with_name("counts.tsv.by_family.tsv").exists()


def test_refuses_missing_column_and_missing_call(tmp_path, capsys):
    out = make_out(tmp_path)
    _write(out / "calls.long.tsv.gz", ["protein", "call", "value"], [["a", "x", "called"]])
    code, _, captured = run(out, tmp_path, capsys)
    assert code == 2 and "missing column(s) variant" in captured.err
    out2 = make_out(tmp_path / "second", drop="wall_family_domain")
    code, _, captured = run(out2, tmp_path, capsys)
    assert code == 2 and "missing call(s) wall_family_domain" in captured.err


def test_refuses_bad_value(tmp_path, capsys):
    out = make_out(tmp_path)
    _write(out / "calls.long.tsv.gz", cmp.LONG_NEEDED, [["a", "x", "", "maybe"]])
    code, _, captured = run(out, tmp_path, capsys)
    assert code == 2 and "not a call value" in captured.err


def test_missing_directory(tmp_path, capsys):
    code, _, captured = run(tmp_path / "nope", tmp_path, capsys)
    assert code == 2 and "not a directory" in captured.err
    with pytest.raises(SystemExit):
        cmp.main([str(tmp_path), "--out", str(tmp_path / "o.tsv"), "--label", "x", str(tmp_path)])
