"""Task 3.5 (D6, D13): the derived marker in calls.long and the leakage exception in the basis."""

import gzip

from test_cli import VALID, build, read_long
from test_status_loader_rules import measure

from cellsurface_sorting_hat.cli import main

EXCEPTION = {
    "reason": "leakage unmeasured (SignalP training overlap), uncapped by owner decision",
    "signer": "owner",
    "date": "2026-10-06",
}


def _run(tmp_path, write_module, nodes_dmp, exception=True):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    notes = "call=signal_peptide_protein; leakage: unknown (not measured)"
    m = measure(notes=notes, **({"leakage_exception": EXCEPTION} if exception else {}))
    rows = [
        {"id": i, "state": "ok", "call": v}
        for i, v in zip(VALID, ["called", "not_called", "called", "called"], strict=True)
    ]
    write_module(
        wd,
        "step1_rule@R0",
        rows,
        status=[{"taxa": [40], "status": "estimated", "source": "phasec", "measure": m}],
    )
    out = tmp_path / "out"
    argv = ["--fasta", str(fasta), "--taxon-map", str(taxon_map), "--taxdump", str(nodes_dmp)]
    argv += ["--workdir", str(wd), "--out", str(out)]
    return main(argv), out


def test_a_record_whose_status_comes_from_the_r0_entry_has_the_exception_text(
    tmp_path, write_module, nodes_dmp
):
    code, out = _run(tmp_path, write_module, nodes_dmp)
    assert code == 0
    got = read_long(out)
    r0 = got[("ENZ1", "signal_peptide_protein", "R0")]
    assert r0["status"] == "estimated"
    assert "uncapped by owner decision" in r0["status_basis"]
    leaf = got[("ENZ1", "tandem_repeat_protein", "")]
    assert "uncapped" not in leaf["status_basis"]


def test_a_composite_that_reads_r0_shows_the_exception_too(tmp_path, write_module, nodes_dmp):
    _, out = _run(tmp_path, write_module, nodes_dmp)
    composite = read_long(out)[("ENZ1", "serodiagnostic_marker_candidate", "R0")]
    assert "uncapped by owner decision" in composite["status_basis"]


def test_calls_long_has_a_derived_from_column_and_calls_wide_does_not(
    tmp_path, write_module, nodes_dmp
):
    code, out = _run(tmp_path, write_module, nodes_dmp)
    assert code == 0
    got = read_long(out)
    assert got[("SOW1", "cell_wall_adhesion_candidate", "R0")]["derived_from"] != ""
    assert got[("SOW1", "signal_peptide_protein", "R0")]["derived_from"] == ""
    with gzip.open(out / "calls.wide.tsv.gz", "rt") as fh:
        assert "derived_from" not in fh.readline()
