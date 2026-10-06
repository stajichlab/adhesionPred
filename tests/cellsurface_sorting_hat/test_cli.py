"""End to end: FASTA + module tables + status sources -> calls, wide table, report."""

import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.cli import main

FASTA = (
    ">SOW1 SOWgp-like\nMKTAYIAKQRQISFVKSHFSRQ\n"
    ">ENZ1 enzyme\nMNLLPQWERTYIPASDFG\n"
    ">DUP1 same sequence as SOW1\nMKTAYIAKQRQISFVKSHFSRQ\n"
    ">BAD1 internal stop\nMKT*AYI\n"
    ">STAR1 CRLF, lowercase, trailing stop\nmktayhhhh*\n"
)
TAXA = {"SOW1": 41, "DUP1": 41, "ENZ1": 40, "BAD1": 40, "STAR1": 40}
VALID = ["SOW1", "ENZ1", "DUP1", "STAR1"]


def build(tmp_path, write_module):
    fasta = tmp_path / "p.faa"
    fasta.write_bytes(FASTA.encode())
    taxon_map = tmp_path / "taxa.tsv"
    taxon_map.write_text("".join(f"{k}\t{v}\n" for k, v in TAXA.items()))
    wd = tmp_path / "wd"

    def call(values):
        return [{"id": i, "state": "ok", "call": v} for i, v in zip(VALID, values, strict=True)]

    write_module(
        wd,
        "step1_rule@R0",
        call(["called", "not_called", "called", "called"]),
        status=[{"taxa": [40], "status": "estimated", "source": "phasec"}],
    )
    write_module(wd, "repeat02", call(["called", "not_called", "called", "not_called"]))
    write_module(wd, "repeat14", call(["not_called"] * 4))
    write_module(wd, "pfam_adhesion", [{"id": i, "state": "ok", "hit": "0"} for i in VALID])
    write_module(wd, "pfam_allergen", [{"id": i, "state": "ok", "hit": "0"} for i in VALID])
    write_module(
        wd,
        "antigen_lookup",
        [
            {"id": "SOW1", "state": "ok", "percentile": "7.1"},
            {"id": "DUP1", "state": "ok", "percentile": "7.1"},
            {"id": "ENZ1", "state": "not_applicable", "percentile": ""},
            {"id": "STAR1", "state": "not_applicable", "percentile": ""},
        ],
    )
    write_module(
        wd,
        "allergen_homology",
        [
            {"id": "SOW1", "state": "ok", "identity": "0", "coverage": "0"},
            {"id": "ENZ1", "state": "ok", "identity": "82", "coverage": "90"},
            {"id": "DUP1", "state": "ok", "identity": "0", "coverage": "0"},
            {"id": "STAR1", "state": "ok", "identity": "40", "coverage": "85"},
        ],
    )
    return fasta, taxon_map, wd


def run_cli(tmp_path, write_module, nodes_dmp, extra=None, **over):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    out = tmp_path / "out"
    argv = [
        "--fasta",
        str(fasta),
        "--taxon-map",
        str(taxon_map),
        "--taxdump",
        str(nodes_dmp),
        "--workdir",
        str(wd),
        "--out",
        str(out),
    ] + (extra or [])
    return main(argv), out, wd


def read_long(out):
    with gzip.open(out / "calls.long.tsv.gz", "rt") as fh:
        return {
            (r["protein"], r["call"], r["variant"]): r for r in csv.DictReader(fh, delimiter="\t")
        }


def test_golden_calls(tmp_path, write_module, nodes_dmp):
    code, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    assert code == 0
    got = read_long(out)

    def v(p, c, var=""):
        r = got[(p, c, var)]
        return r["value"], r["status"]

    # SOW1: Coccidioides (taxon 41), not covered by the step 1 estimate (tested on taxon 40)
    assert v("SOW1", "surface_glycoprotein", "R0") == ("called", "unvalidated")
    assert v("SOW1", "adhesion_repeat") == ("called", "unvalidated")
    assert v("SOW1", "cell_wall_adhesion_candidate", "R0") == ("called", "unvalidated")
    assert v("SOW1", "antigen_candidate") == ("called", "unvalidated")
    assert v("SOW1", "antigen_candidate_surface", "R0") == ("called", "unvalidated")
    assert v("SOW1", "allergen_candidate")[0] == "not_called"
    assert v("SOW1", "other_not_surface", "R0")[0] == "not_called"
    assert v("SOW1", "other_surface_no_mechanism", "R0")[0] == "not_called"
    assert got[("SOW1", "surface_glycoprotein", "R0")]["status_basis"] == (
        "step1_rule@R0:taxon not tested"
    )
    # ENZ1: A. fumigatus (40), step 1 estimated; antigen not applicable; allergen homology called
    assert v("ENZ1", "surface_glycoprotein", "R0") == ("not_called", "estimated")
    assert v("ENZ1", "antigen_candidate")[0] == "not_assessable"
    assert v("ENZ1", "antigen_candidate_surface", "R0") == ("not_called", "estimated")
    assert v("ENZ1", "allergen_candidate")[0] == "called"
    assert v("ENZ1", "other_not_surface", "R0")[0] == "not_called"
    # STAR1: surface, no mechanism called, antigen left out of the basis
    r = got[("STAR1", "other_surface_no_mechanism", "R0")]
    assert (r["value"], r["other_basis"]) == ("called", "antigen_candidate")
    # BAD1: invalid protein, every call unknown
    assert {r["value"] for k, r in got.items() if k[0] == "BAD1"} == {"not_assessable"}


def test_identical_sequences_with_different_ids_get_identical_calls(
    tmp_path, write_module, nodes_dmp
):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    got = read_long(out)
    keys = {(c, var) for (p, c, var) in got if p == "SOW1"}
    assert keys == {(c, var) for (p, c, var) in got if p == "DUP1"}
    for c, var in keys:
        a, b = got[("SOW1", c, var)], got[("DUP1", c, var)]
        assert (a["value"], a["status"], a["other_basis"]) == (
            b["value"],
            b["status"],
            b["other_basis"],
        )


def test_crlf_lowercase_trailing_stop_fasta_is_read(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    fasta.write_bytes(FASTA.replace("\n", "\r\n").encode())
    out = tmp_path / "out"
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon-map",
            str(taxon_map),
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(wd),
            "--out",
            str(out),
        ]
    )
    assert code == 0
    assert read_long(out)[("STAR1", "surface_glycoprotein", "R0")]["value"] == "called"


def test_unknown_taxon_stops_the_run(tmp_path, write_module, nodes_dmp, capsys):
    fasta, _, wd = build(tmp_path, write_module)
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon",
            "999",
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(wd),
            "--out",
            str(tmp_path / "out"),
        ]
    )
    assert code == 2
    assert "taxon 999" in capsys.readouterr().err


def test_a_protein_with_no_taxon_stops_the_run(tmp_path, write_module, nodes_dmp, capsys):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_text("SOW1\t41\n")  # the others have no taxon and no --taxon
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon-map",
            str(taxon_map),
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(wd),
            "--out",
            str(tmp_path / "out"),
        ]
    )
    assert code == 2
    assert "no taxon for protein" in capsys.readouterr().err


def test_taxon_map_overrides_taxon(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_text("SOW1\t41\n")
    out = tmp_path / "out"
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon",
            "40",
            "--taxon-map",
            str(taxon_map),
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(wd),
            "--out",
            str(out),
        ]
    )
    assert code == 0
    got = read_long(out)
    assert got[("SOW1", "surface_glycoprotein", "R0")]["status"] == "unvalidated"  # taxon 41
    assert got[("STAR1", "surface_glycoprotein", "R0")]["status"] == "estimated"  # --taxon 40


def test_partial_module_output_is_reported_and_missing_proteins_are_unknown(
    tmp_path, write_module, nodes_dmp
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(wd, "repeat14", [{"id": "SOW1", "state": "ok", "call": "not_called"}])
    out = tmp_path / "out"
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon-map",
            str(taxon_map),
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(wd),
            "--out",
            str(out),
        ]
    )
    assert code == 0
    assert "| repeat14 | partial |" in (out / "report.md").read_text()
    got = read_long(out)
    # ENZ1: repeat02 false and repeat14 missing -> unknown, not false
    assert got[("ENZ1", "adhesion_repeat", "")]["value"] == "not_assessable"
    # SOW1: repeat02 true wins over anything
    assert got[("SOW1", "adhesion_repeat", "")]["value"] == "called"


def test_module_with_run_state_unavailable_is_treated_as_absent(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(wd, "step1_rule@R2", [], meta={"run_state": "unavailable"})
    out = tmp_path / "out"
    assert (
        main(
            [
                "--fasta",
                str(fasta),
                "--taxon-map",
                str(taxon_map),
                "--taxdump",
                str(nodes_dmp),
                "--workdir",
                str(wd),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    report = (out / "report.md").read_text()
    assert "| step1_rule@R2 | unavailable |" in report
    assert "Unavailable step 1 variants:" in report and "step1_rule@R2" in report
    assert not any(k[2] == "R2" for k in read_long(out))


def test_no_module_at_all_still_writes_a_report_with_a_warning(tmp_path, nodes_dmp):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKTAYI\n")
    (tmp_path / "wd").mkdir()
    out = tmp_path / "out"
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon",
            "40",
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(tmp_path / "wd"),
            "--out",
            str(out),
        ]
    )
    assert code == 0
    report = (out / "report.md").read_text()
    assert "WARNING: The default gate step1_rule@R0 is not available" in report
    got = read_long(out)
    assert {r["value"] for r in got.values()} == {"not_assessable"}
    assert not any(k[2] for k in got)  # no per-variant records without a step 1 variant


def test_wide_table_and_run_json(tmp_path, write_module, nodes_dmp):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    with gzip.open(out / "calls.wide.tsv.gz", "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["protein"] for r in rows] == ["SOW1", "ENZ1", "DUP1", "BAD1", "STAR1"]
    star = rows[-1]
    assert star["other_surface_no_mechanism[R0]"] == "called"
    assert star["other_surface_no_mechanism[R0]_basis"] == "antigen_candidate"
    assert star["surface_glycoprotein[R0]_status"] == "estimated"
    run = json.loads((out / "run.json").read_text())
    assert run["n_proteins"] == 5 and run["n_invalid"] == 1
    assert run["taxa"] == {"40": 3, "41": 2}


def test_report_has_the_known_limits_and_counts(tmp_path, write_module, nodes_dmp):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    report = (out / "report.md").read_text()
    assert "## Known limits" in report
    assert 'means "SignalP calls a signal peptide"' in report
    assert "| surface_glycoprotein | R0 | 3 | 1 | 1 |" in report
    assert "- antigen_candidate: " in report


def _run(tmp_path, nodes_dmp, fasta, taxon_map, wd, out="out"):
    argv = [
        "--fasta",
        str(fasta),
        "--taxdump",
        str(nodes_dmp),
        "--workdir",
        str(wd),
        "--out",
        str(tmp_path / out),
    ]
    if taxon_map is not None:
        argv += ["--taxon-map", str(taxon_map)]
    return main(argv), tmp_path / out


def test_evidence_and_protein_tables(tmp_path, write_module, nodes_dmp):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    with gzip.open(out / "evidence.tsv.gz", "rt") as fh:
        ev = {
            (r["protein"], r["module"], r["field"]): r["value"]
            for r in csv.DictReader(fh, delimiter="\t")
        }
    assert ev[("SOW1", "antigen_lookup", "percentile")] == "7.1"
    assert ev[("ENZ1", "allergen_homology", "identity")] == "82"
    assert ("ENZ1", "antigen_lookup", "percentile") not in ev  # state not_applicable
    with gzip.open(out / "proteins.tsv.gz", "rt") as fh:
        prot = {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert prot["SOW1"]["sha256"] == prot["DUP1"]["sha256"]
    assert prot["STAR1"]["trailing_stop"] == "1" and prot["SOW1"]["trailing_stop"] == "0"
    assert prot["BAD1"]["state"] == "na_invalid" and prot["SOW1"]["taxon"] == "41"


def test_report_shows_thresholds_invalid_proteins_and_trailing_stops(
    tmp_path, write_module, nodes_dmp
):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    report = (out / "report.md").read_text()
    assert "antigen_percentile_max = 15" in report and "allergen_identity_min = 70" in report
    assert "- BAD1: internal stop codon" in report
    assert "1 had a trailing `*`" in report
    assert "| antigen_lookup | not_applicable | 2 |" in report


def test_identical_sequences_with_different_module_rows_are_flagged(
    tmp_path, write_module, nodes_dmp
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(
        wd,
        "repeat02",
        [
            {"id": i, "state": "ok", "call": "called" if i == "SOW1" else "not_called"}
            for i in VALID
        ],
    )
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert (
        "Module repeat02 gives different rows to 1 group(s) of identical sequences"
        in (out / "report.md").read_text()
    )


def test_a_bad_call_value_is_unknown_and_counted(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(
        wd,
        "repeat14",
        [
            {"id": i, "state": "ok", "call": "Called" if i == "SOW1" else "not_called"}
            for i in VALID
        ],
    )
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert "| repeat14 | bad_value | 1 |" in (out / "report.md").read_text()
    got = read_long(out)
    assert got[("SOW1", "adhesion_repeat", "")]["value"] == "called"  # repeat02 is true
    assert got[("STAR1", "adhesion_repeat", "")]["value"] == "not_called"


def test_a_module_whose_ids_match_no_protein_is_an_error(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(
        wd, "repeat14", [{"id": f"sp|{i}|X", "state": "ok", "call": "not_called"} for i in VALID]
    )
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    report = (out / "report.md").read_text()
    assert "| repeat14 | error |" in report
    assert "no FASTA ID matches the 4 row(s)" in report
    assert "WARNING: Module repeat14 is in state error" in report


def test_module_table_without_a_state_column_stops_the_run(
    tmp_path, write_module, nodes_dmp, capsys
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    with gzip.open(wd / "modules" / "repeat14.tsv.gz", "wt") as fh:
        fh.write("id\tcall\nSOW1\tcalled\n")
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert "repeat14.tsv.gz: missing column 'state'" in capsys.readouterr().err


def test_duplicate_rows_in_a_module_table_stop_the_run(tmp_path, write_module, nodes_dmp, capsys):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    with gzip.open(wd / "modules" / "repeat14.tsv.gz", "wt") as fh:
        fh.write("id\tstate\tcall\nSOW1\tok\tcalled\nSOW1\tok\tnot_called\n")
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert "duplicate row for ID 'SOW1'" in capsys.readouterr().err


def test_a_run_record_without_a_table_is_listed_with_its_state(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    (wd / "modules" / "step1_rule@R1.json").write_text(json.dumps({"run_state": "unavailable"}))
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert "| step1_rule@R1 | unavailable |" in (out / "report.md").read_text()


def test_absent_required_modules_are_listed(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    (wd / "modules" / "pfam_adhesion.tsv.gz").unlink()
    (wd / "modules" / "pfam_adhesion.json").unlink()
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert (
        "Modules the rules need and that were not found: pfam_adhesion"
        in (out / "report.md").read_text()
    )
    assert read_long(out)[("SOW1", "adhesion_domain", "")]["value"] == "not_assessable"


@pytest.mark.parametrize(
    "text,message",
    [
        ("SOW1 41\n", "taxa.tsv:1: expected"),
        ("protein\ttaxon\n", "taxa.tsv:1: taxon ID is not an integer"),
    ],
)
def test_bad_taxon_map_lines_name_the_line(
    tmp_path, write_module, nodes_dmp, capsys, text, message
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_text(text)
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert message in capsys.readouterr().err


def test_taxon_map_ids_not_in_the_fasta_are_reported(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    with open(taxon_map, "a") as fh:
        fh.write("GHOST\t40\n")
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert "1 ID(s) in --taxon-map are not in the FASTA" in (out / "report.md").read_text()


def test_a_broken_status_source_stops_the_run_and_names_the_file(
    tmp_path, write_module, nodes_dmp, capsys
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    (wd / "status" / "step1_rule@R0.json").write_text(json.dumps({"module": "x"}))
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert "step1_rule@R0.json: not a valid status source" in capsys.readouterr().err


def test_a_program_error_is_not_hidden_as_a_user_error(
    tmp_path, write_module, nodes_dmp, monkeypatch
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    import cellsurface_sorting_hat.cli as cli

    def boom(*a, **k):
        raise KeyError("bug")

    monkeypatch.setattr(cli, "evaluate", boom)
    with pytest.raises(KeyError):
        _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)


def _module_bytes(wd, name, data):
    with gzip.open(wd / "modules" / f"{name}.tsv.gz", "wb") as fh:
        fh.write(data)


@pytest.mark.parametrize(
    "data,message",
    [
        (
            b"id\tstate\tcall\nSOW1\tok\tcalled\textra\n",
            "repeat14.tsv.gz:2: 4 fields, the header has 3",
        ),
        (b"id\tstate\tcall\nSOW1\tok\n", "repeat14.tsv.gz:2: 2 fields"),
    ],
)
def test_a_module_row_with_the_wrong_number_of_fields_names_file_and_line(
    tmp_path, write_module, nodes_dmp, capsys, data, message
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    _module_bytes(wd, "repeat14", data)
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2 and message in capsys.readouterr().err


def test_a_truncated_module_table_names_the_file(tmp_path, write_module, nodes_dmp, capsys):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    path = wd / "modules" / "repeat14.tsv.gz"
    path.write_bytes(path.read_bytes()[:20])
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2 and "repeat14.tsv.gz: cannot be read" in capsys.readouterr().err


@pytest.mark.parametrize("text", ["{not json", "[]"])
def test_a_bad_run_record_names_the_file(tmp_path, write_module, nodes_dmp, capsys, text):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    (wd / "modules" / "repeat14.json").write_text(text)
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2 and "repeat14.json:" in capsys.readouterr().err


def test_a_failed_run_leaves_no_partial_output(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    _module_bytes(wd, "repeat14", b"id\tcall\nSOW1\tcalled\n")  # no state column
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2 and not out.exists()


def test_byte_order_marks_in_the_taxon_map_and_module_table_are_ignored(
    tmp_path, write_module, nodes_dmp
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_bytes(b"\xef\xbb\xbf" + taxon_map.read_bytes())
    path = wd / "modules" / "repeat14.tsv.gz"
    with gzip.open(path, "rb") as fh:
        raw = fh.read()
    _module_bytes(wd, "repeat14", b"\xef\xbb\xbf" + raw)
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert read_long(out)[("SOW1", "surface_glycoprotein", "R0")]["value"] == "called"


def test_a_repeated_id_in_the_taxon_map_names_both_lines(tmp_path, write_module, nodes_dmp, capsys):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    taxon_map.write_text("SOW1\t41\nENZ1\t40\nSOW1\t40\n")
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert "taxa.tsv:3: ID 'SOW1' already on line 1" in capsys.readouterr().err


def test_an_unknown_taxon_stops_the_run_even_with_no_status_file(tmp_path, nodes_dmp, capsys):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKTAYI\n")
    (tmp_path / "wd").mkdir()
    code = main(
        [
            "--fasta",
            str(fasta),
            "--taxon",
            "999",
            "--taxdump",
            str(nodes_dmp),
            "--workdir",
            str(tmp_path / "wd"),
            "--out",
            str(tmp_path / "out"),
        ]
    )
    assert code == 2 and "taxon 999" in capsys.readouterr().err


def test_calibration_is_shown_when_the_status_source_has_a_measure(
    tmp_path, write_module, nodes_dmp
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    status = json.loads((wd / "status" / "step1_rule@R0.json").read_text())
    status["entries"][0]["measure"] = {
        "calibration_set": "S1:all",
        "n_pos": 232,
        "n_neg": 4244,
        "sensitivity": {"value": 0.603, "lo": 0.55, "hi": 0.65},
        "specificity": {"value": 0.963, "lo": 0.95, "hi": 0.97},
    }
    (wd / "status" / "step1_rule@R0.json").write_text(json.dumps(status))
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    report = (out / "report.md").read_text()
    assert "## Module calibration" in report
    assert (
        "| step1_rule@R0 | 40 | estimated | S1:all | 232 | 4244 | 0.603 [0.550, 0.650] | 0.963 [0.950, 0.970] |"
        in report
    )
    # taxon 41 is not covered by the entry; a module with no status file says so
    assert (
        "| step1_rule@R0 | 41 | unvalidated | - | - | - | not measured | not measured |" in report
    )
    assert "| repeat02 | 40 | unvalidated | - | - | - | not measured | not measured |" in report
    run = json.loads((out / "run.json").read_text())
    row = next(c for c in run["calibration"] if c["module"] == "step1_rule@R0" and c["taxon"] == 40)
    assert row["n_pos"] == 232 and row["matched_taxon"] == 40


@pytest.mark.parametrize("content", [">A\nMKT\n>A\nMKT\n", "", "MKT\n"])
def test_a_fasta_error_names_the_fasta_file(tmp_path, write_module, nodes_dmp, capsys, content):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    fasta.write_text(content)
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2
    assert str(fasta) in capsys.readouterr().err


@pytest.mark.parametrize(
    "mutate",
    [
        lambda d: d["entries"][0].update(taxa=["abc"]),
        lambda d: d.pop("version"),
        lambda d: d["entries"][0].update(taxa=5),
    ],
)
def test_status_source_errors_name_the_file(tmp_path, write_module, nodes_dmp, capsys, mutate):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    path = wd / "status" / "step1_rule@R0.json"
    data = json.loads(path.read_text())
    mutate(data)
    path.write_text(json.dumps(data))
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2 and "step1_rule@R0.json: not a valid status source" in capsys.readouterr().err


def test_status_source_with_bad_json_names_the_file(tmp_path, write_module, nodes_dmp, capsys):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    (wd / "status" / "step1_rule@R0.json").write_text("{not json")
    code, _ = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 2 and "step1_rule@R0.json: not a valid status source" in capsys.readouterr().err


def test_a_bad_hit_value_is_bad_value(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    write_module(wd, "pfam_adhesion", [{"id": i, "state": "ok", "hit": "yes"} for i in VALID])
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert "| pfam_adhesion | bad_value | 4 |" in (out / "report.md").read_text()


def test_calibration_table_has_rows_even_with_no_module_and_no_status(tmp_path, nodes_dmp):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKTAYI\n")
    (tmp_path / "wd").mkdir()
    code, out = _run(tmp_path, nodes_dmp, fasta, None, tmp_path / "wd")
    assert code == 2  # no taxon given
    argv = [
        "--fasta",
        str(fasta),
        "--taxon",
        "40",
        "--taxdump",
        str(nodes_dmp),
        "--workdir",
        str(tmp_path / "wd"),
        "--out",
        str(out),
    ]
    assert main(argv) == 0
    report = (out / "report.md").read_text()
    assert "## Module calibration" in report
    assert (
        "| step1_rule@R0 | 40 | unvalidated | - | - | - | not measured | not measured |" in report
    )


def test_module_rows_with_ids_not_in_the_fasta_are_counted_and_change_no_call(
    tmp_path, write_module, nodes_dmp
):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    code, base = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd, out="base")
    rows = [{"id": i, "state": "ok", "call": "not_called"} for i in VALID]
    rows += [{"id": f"GHOST{n}", "state": "ok", "call": "called"} for n in range(3)]
    write_module(wd, "repeat14", rows)
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert (
        "Module repeat14 has 3 row(s) whose ID is not in the FASTA"
        in (out / "report.md").read_text()
    )
    assert json.loads((out / "run.json").read_text())["unmatched_module_ids"] == {"repeat14": 3}
    assert read_long(out) == read_long(base)


def test_a_bad_value_gives_unknown_not_false(tmp_path, write_module, nodes_dmp):
    fasta, taxon_map, wd = build(tmp_path, write_module)
    # ENZ1: repeat02 not_called, repeat14 bad value -> the OR cannot be decided
    write_module(
        wd,
        "repeat14",
        [{"id": i, "state": "ok", "call": "bogus" if i == "ENZ1" else "not_called"} for i in VALID],
    )
    code, out = _run(tmp_path, nodes_dmp, fasta, taxon_map, wd)
    assert code == 0
    assert read_long(out)[("ENZ1", "adhesion_repeat", "")]["value"] == "not_assessable"


def test_every_table_output_has_a_matching_sha256_sidecar(tmp_path, write_module, nodes_dmp):
    import hashlib

    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    for name in ("calls.long", "calls.wide", "evidence", "proteins"):
        path = out / f"{name}.tsv.gz"
        side = out / f"{name}.tsv.gz.sha256"
        assert side.read_text().strip() == hashlib.sha256(path.read_bytes()).hexdigest()
