"""Task 3.6 (D6): the engine compares the R0 table with the one that a Pfam module was built against."""

import json

from test_cli import VALID, build
from test_narrow_call import DOM, OPTIONS, fam_row, write_table

from cellsurface_sorting_hat.cli import main

R0 = "step1_rule@R0"
R0_CALLS = dict(zip(VALID, ["called", "not_called", "called", "called"], strict=True))


def digest(values):
    from cellsurface_sorting_hat.modules.base import values_digest

    return values_digest(values)


def test_the_digest_depends_on_ids_and_values_not_on_order():
    a = digest({"A": "called", "B": "not_called"})
    assert a == digest({"B": "not_called", "A": "called"})
    assert a != digest({"A": "not_called", "B": "called"})  # same ids, other values
    assert a != digest({"A": "called", "C": "not_called"})  # same values, other ids


def test_the_pfam_command_records_a_digest_of_the_condition_table(tmp_path):
    from cellsurface_sorting_hat.modules.cli import main as module_main

    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKTAYIAKQRQISFVKSHFSRQ\n>B\nMNLLPQWERTYIPASDFG\n")
    domtbl = tmp_path / "d.domtbl"
    domtbl.write_text(OPTIONS + DOM.format(t="A", n="X", a="PF00001") + "# [ok]\n")
    table = write_table(tmp_path / "t.tsv", [fam_row("PF00001")])
    wd = tmp_path / "wd"
    from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module

    write_module(
        wd,
        ModuleSpec(R0, "1"),
        ["call"],
        [
            {"id": "A", "state": "ok", "call": "called"},
            {"id": "B", "state": "ok", "call": "not_called"},
        ],
    )
    argv = ["pfam", "--fasta", str(fasta), "--workdir", str(wd), "--domtbl", str(domtbl)]
    argv += ["--family-table", str(table), "--pfam-release", "38.2", "--pfam-sha256", "ab"]
    argv += ["--hmmer-version", "3.4", "--sp-module", R0]
    assert module_main(argv) == 0
    record = json.loads((wd / "modules" / "pfam_adhesion.json").read_text())
    cond = record["params"]["conditions"]["sp_module"]
    assert cond["module"] == R0 and "values_sha256" not in cond and "column" not in cond
    snap = record["condition_snapshots"]["sp_module"]
    assert snap["column"] == "call"
    assert snap["values_sha256"] == digest({"A": "called", "B": "not_called"})


def _run(tmp_path, write_module, nodes_dmp, cond):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    rows = [{"id": i, "state": "ok", "hit": "0", "hit_narrow": "0"} for i in VALID]
    meta = {}
    if cond is not None:
        ident = {k: v for k, v in cond.items() if k not in ("column", "values_sha256")}
        snap = {k: v for k, v in cond.items() if k in ("column", "values_sha256")}
        meta = {
            "params": {"conditions": {"sp_module": ident}},
            "condition_snapshots": {"sp_module": snap},
        }
    write_module(wd, "pfam_adhesion", rows, meta=meta)
    out = tmp_path / "out"
    argv = ["--fasta", str(fasta), "--taxon-map", str(taxon_map), "--taxdump", str(nodes_dmp)]
    argv += ["--workdir", str(wd), "--out", str(out)]
    assert main(argv) == 0
    return out, json.loads((out / "run.json").read_text()), (out / "report.md").read_text()


def _cond(values, **over):
    cond = {"module": R0, "column": "call", "params_hash": "p", "artefact_hash": "a"}
    cond["values_sha256"] = digest(values)
    cond.update(over)
    return cond


def test_a_pfam_module_built_against_another_r0_table_gives_a_warning(
    tmp_path, write_module, nodes_dmp
):
    other = {**R0_CALLS, "ENZ1": "called"}  # table B differs in one value
    _, run_json, report = _run(tmp_path, write_module, nodes_dmp, _cond(other))
    assert "**WARNING: " in report
    warning = next(ln for ln in report.splitlines() if "WARNING" in ln and "pfam_adhesion" in ln)
    assert R0 in warning
    checks = run_json["condition_checks"]
    assert checks == [
        {"module": "pfam_adhesion", "condition": "sp_module", "source": R0, "result": "mismatch"}
    ]


def test_the_same_r0_table_gives_no_warning(tmp_path, write_module, nodes_dmp):
    _, run_json, report = _run(tmp_path, write_module, nodes_dmp, _cond(R0_CALLS))
    assert not [ln for ln in report.splitlines() if "WARNING" in ln and "pfam_adhesion" in ln]
    assert run_json["condition_checks"][0]["result"] == "match"


def test_an_old_record_without_a_digest_gives_a_note_not_an_error(
    tmp_path, write_module, nodes_dmp
):
    old = {"module": R0, "params_hash": "p", "artefact_hash": "a"}
    _, run_json, report = _run(tmp_path, write_module, nodes_dmp, old)
    assert run_json["condition_checks"][0]["result"] == "no_digest"
    assert "## Notes" in report and "pfam_adhesion" in report.split("## Notes")[1]
    assert not [ln for ln in report.splitlines() if "WARNING" in ln and "pfam_adhesion" in ln]


def test_a_module_without_conditions_is_not_checked(tmp_path, write_module, nodes_dmp):
    _, run_json, report = _run(tmp_path, write_module, nodes_dmp, None)
    assert run_json["condition_checks"] == []
    assert "## Notes" not in report


def test_run_json_always_has_the_field(tmp_path, write_module, nodes_dmp):
    _, run_json, _ = _run(tmp_path, write_module, nodes_dmp, _cond(R0_CALLS))
    assert "condition_checks" in run_json
