"""The module command end to end on tiny inputs."""

import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.modules.cli import main
from cellsurface_sorting_hat.modules.pfam import FAMILY_COLUMNS

FASTA = ">XP_1 a\nMKTAYIAKQRQ\n>XP_2 b\nMNLLPQWERT\n>BAD\nMK*T\n"

OPTIONS = "# Option settings:     hmmsearch --cut_ga --cpu 2 --noali fam.hmm in.fasta\n"


def write_domtbl(path, body):
    path.write_text(OPTIONS + body + "# [ok]\n")


def read(workdir, name):
    with gzip.open(workdir / "modules" / f"{name}.tsv.gz", "rt") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "p.faa"
    path.write_text(FASTA)
    return path


def test_signalp_command_writes_the_r0_module(tmp_path, fasta, capsys):
    res = tmp_path / "prediction_results.txt"
    res.write_text(
        "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 1\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\nXP_1 x\tSP\t0.1\t0.9\t\nXP_2 y\tOTHER\t0.9\t0.1\t\n"
    )
    wd = tmp_path / "wd"
    code = main(
        [
            "signalp",
            "--fasta",
            str(fasta),
            "--workdir",
            str(wd),
            "--results",
            str(res),
            "--signalp-version",
            "6.0h-gpu",
        ]
    )
    assert code == 0 and json.loads(capsys.readouterr().out)["module"] == "step1_rule@R0"
    rows = read(wd, "step1_rule@R0")
    assert (rows["XP_1"]["call"], rows["XP_2"]["call"], rows["BAD"]["state"]) == (
        "called",
        "not_called",
        "na_invalid",
    )
    run = json.loads((wd / "modules" / "step1_rule@R0.json").read_text())
    assert run["tools"] == {"signalp": "6.0h-gpu"} and run["params"]["rule"] == "R0"


def test_pfam_command_writes_both_modules_and_uses_the_recorded_pfam_digest(tmp_path, fasta):
    table = tmp_path / "f.tsv"
    row = dict.fromkeys(FAMILY_COLUMNS, "")
    row.update(
        pfam_acc="PF05730",
        name="CFEM",
        module="pfam_adhesion",
        **{"class": "2b-i"},
        active="yes",
        active_by="owner",
        active_date="2026-10-05",
    )
    row2 = dict(row, pfam_acc="PF16541", name="AltA1", module="pfam_allergen")
    table.write_text(
        "\n".join(
            "\t".join(r[c] for c in FAMILY_COLUMNS)
            for r in [dict(zip(FAMILY_COLUMNS, FAMILY_COLUMNS, strict=True)), row, row2]
        )
        + "\n"
    )
    dom = tmp_path / "d.domtbl"
    write_domtbl(
        dom, "XP_1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    wd = tmp_path / "wd"
    assert (
        main(
            [
                "pfam",
                "--fasta",
                str(fasta),
                "--workdir",
                str(wd),
                "--domtbl",
                str(dom),
                "--family-table",
                str(table),
                "--hmmer-version",
                "3.4",
                "--pfam-release",
                "38.2",
                "--pfam-sha256",
                "abc",
            ]
        )
        == 0
    )
    assert read(wd, "pfam_adhesion")["XP_1"]["hit"] == "1"
    assert read(wd, "pfam_allergen")["XP_1"]["hit"] == "0"
    rec = json.loads((wd / "modules" / "pfam_adhesion.json").read_text())
    assert rec["artefact_hash"].startswith("abc:") and rec["params"]["families"] == ["PF05730"]


def test_repeat_command(tmp_path, fasta):
    tab = tmp_path / "r.tsv"
    tab.write_text(
        "protein\trep_period\trep_n_copies\trep_coverage\nXP_1\t5\t3.0\t0.5\nXP_2\t0\t0\t0\n"
    )
    wd = tmp_path / "wd"
    assert main(["repeat02", "--fasta", str(fasta), "--workdir", str(wd), "--table", str(tab)]) == 0
    assert read(wd, "repeat02")["XP_1"]["call"] == "called"


def test_allergen_command(tmp_path, fasta):
    blast = tmp_path / "b.tsv"
    blast.write_text("XP_2\tAsp_f_1.0101|11\t82.0\t90\t120\t100\t150\t1e-20\n")
    ref = tmp_path / "a.faa"
    ref.write_text(">Asp_f_1.0101|11\nMKT\n")
    wd = tmp_path / "wd"
    assert (
        main(
            [
                "allergen",
                "--fasta",
                str(fasta),
                "--workdir",
                str(wd),
                "--blast",
                str(blast),
                "--allergen-fasta",
                str(ref),
                "--blast-version",
                "2.16.0+",
            ]
        )
        == 0
    )
    r = read(wd, "allergen_homology")["XP_2"]
    assert (r["identity"], r["coverage"], r["allergen_name"]) == ("82.00", "90.0", "Asp_f_1.0101")


def test_antigen_command_needs_a_taxon_and_marks_other_taxa_not_applicable(tmp_path, fasta, capsys):
    rank = tmp_path / "r.tsv"
    rank.write_text(
        "protein\trank\tpercentile\tantigenicity\tspecificity\tprevalence\tmax_fungal_crossreact_pid\n"
        "CIMG_1-t26_1-p1\t5\t0.05\t2\t1\t1\t0\n"
    )
    pmap = tmp_path / "m.tsv"
    pmap.write_text("protein_id\tgene_id\tproduct\tlength\nXP_1\tCIMG_1\tp\t11\n")
    wd = tmp_path / "wd"
    base = [
        "antigen",
        "--fasta",
        str(fasta),
        "--workdir",
        str(wd),
        "--ranking",
        str(rank),
        "--protein-map",
        str(pmap),
    ]
    assert main(base) == 2 and "give --taxon or --taxon-map" in capsys.readouterr().err
    assert main(base + ["--taxon", "746128"]) == 0
    assert read(wd, "antigen_lookup")["XP_1"]["state"] == "not_applicable"
    assert main(base + ["--taxon", "246410"]) == 0
    assert read(wd, "antigen_lookup")["XP_1"]["percentile"] == "0.05"


def test_tm_command(tmp_path, fasta):
    tab = tmp_path / "t.tsv"
    tab.write_text(
        "protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nXP_1\t11\t0\t0\t0\to\nXP_2\t10\t20\t5\t1\ti5-27o\n"
    )
    wd = tmp_path / "wd"
    assert main(["tm", "--fasta", str(fasta), "--workdir", str(wd), "--table", str(tab)]) == 0
    assert read(wd, "tm")["XP_2"]["n_tm"] == "1"


def test_input_errors_exit_2(tmp_path, fasta, capsys):
    assert (
        main(
            [
                "signalp",
                "--fasta",
                str(fasta),
                "--workdir",
                str(tmp_path),
                "--results",
                str(tmp_path / "none"),
                "--signalp-version",
                "6",
            ]
        )
        == 2
    )
    assert "cellsurface_sorting_hat_module: error" in capsys.readouterr().err


def test_a_result_file_whose_ids_do_not_match_the_fasta_is_refused(tmp_path, fasta, capsys):
    res = tmp_path / "prediction_results.txt"
    res.write_text(
        "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 1\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\nsp|Q1|X\tSP\t0.1\t0.9\t\n"
    )
    code = main(
        [
            "signalp",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--results",
            str(res),
            "--signalp-version",
            "6.0h-gpu",
        ]
    )
    assert code == 2 and "no protein has a result" in capsys.readouterr().err
    assert not (tmp_path / "wd" / "modules").exists()  # nothing was written


def test_the_signalp_module_identity_does_not_depend_on_the_proteome(tmp_path, fasta):
    """A measurement of rule R0 belongs to the tool version, so it must stay valid for another proteome."""
    other = tmp_path / "o.faa"
    other.write_text(">Z1\nMKTAYI\n")
    ids = {}
    for name, f, pid in (("a", fasta, "XP_1"), ("b", other, "Z1")):
        res = tmp_path / f"{name}.txt"
        res.write_text(
            f"# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 1\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n{pid}\tSP\t0.1\t0.9\t\n"
        )
        wd = tmp_path / name
        assert (
            main(
                [
                    "signalp",
                    "--fasta",
                    str(f),
                    "--workdir",
                    str(wd),
                    "--results",
                    str(res),
                    "--signalp-version",
                    "6.0h-gpu",
                ]
            )
            == 0
        )
        ids[name] = json.loads((wd / "modules" / "step1_rule@R0.json").read_text())
    assert (
        ids["a"]["artefact_hash"] == ids["b"]["artefact_hash"]
        and ids["a"]["params_hash"] == ids["b"]["params_hash"]
    )
    changed = tmp_path / "c"
    assert (
        main(
            [
                "signalp",
                "--fasta",
                str(fasta),
                "--workdir",
                str(changed),
                "--results",
                str(tmp_path / "a.txt"),
                "--signalp-version",
                "6.0i-gpu",
            ]
        )
        == 0
    )
    assert (
        json.loads((changed / "modules" / "step1_rule@R0.json").read_text())["artefact_hash"]
        != ids["a"]["artefact_hash"]
    )


def test_pfam_command_applies_the_no_tm_condition_from_the_tm_table(tmp_path, fasta):
    table = tmp_path / "f.tsv"
    row = dict.fromkeys(FAMILY_COLUMNS, "")
    row.update(
        pfam_acc="PF05730",
        name="CFEM",
        module="pfam_adhesion",
        **{"class": "2b-i"},
        second_condition="no_tm",
        active="yes",
        active_by="owner",
        active_date="2026-10-05",
    )
    header = dict(zip(FAMILY_COLUMNS, FAMILY_COLUMNS, strict=True))
    table.write_text(
        "\n".join("\t".join(r[c] for c in FAMILY_COLUMNS) for r in [header, row]) + "\n"
    )
    dom = tmp_path / "d.domtbl"
    write_domtbl(
        dom,
        "XP_1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
        "XP_2 - 10 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n",
    )
    tm = tmp_path / "t.tsv"
    tm.write_text(
        "protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nXP_1\t11\t25\t20\t1\ti5-27o\nXP_2\t10\t150\t20\t7\ti5-27o40-62i70-92o\n"
    )
    wd = tmp_path / "wd"
    assert main(["tm", "--fasta", str(fasta), "--workdir", str(wd), "--table", str(tm)]) == 0
    assert (
        main(
            [
                "pfam",
                "--fasta",
                str(fasta),
                "--workdir",
                str(wd),
                "--domtbl",
                str(dom),
                "--family-table",
                str(table),
                "--hmmer-version",
                "3.4",
                "--pfam-release",
                "38.2",
                "--pfam-sha256",
                "abc",
                "--tm-module",
                "tm",
            ]
        )
        == 0
    )
    rows = read(wd, "pfam_adhesion")
    assert (rows["XP_1"]["hit"], rows["XP_2"]["hit"]) == (
        "1",
        "0",
    )  # XP_1: only the signal peptide read as a helix; XP_2: helices after residue 35


def _family_table(tmp_path, active):
    table = tmp_path / "f.tsv"
    base = dict.fromkeys(FAMILY_COLUMNS, "")
    base.update(
        pfam_acc="PF05730", name="CFEM", module="pfam_adhesion", **{"class": "2b-i"}, active="no"
    )
    if active:
        base.update(active="yes", active_by="owner", active_date="2026-10-05")
    header = dict(zip(FAMILY_COLUMNS, FAMILY_COLUMNS, strict=True))
    table.write_text(
        "\n".join("\t".join(r[c] for c in FAMILY_COLUMNS) for r in [header, base]) + "\n"
    )
    return table


def _pfam_args(tmp_path, fasta, table, dom, wd):
    return [
        "pfam",
        "--fasta",
        str(fasta),
        "--workdir",
        str(wd),
        "--domtbl",
        str(dom),
        "--family-table",
        str(table),
        "--hmmer-version",
        "3.4",
        "--pfam-release",
        "38.2",
        "--pfam-sha256",
        "abc",
    ]


def test_pfam_command_writes_unavailable_modules_when_no_family_is_active(tmp_path, fasta):
    dom = tmp_path / "d.domtbl"
    write_domtbl(
        dom, "XP_1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    wd = tmp_path / "wd"
    assert main(_pfam_args(tmp_path, fasta, _family_table(tmp_path, active=False), dom, wd)) == 0
    rec = json.loads((wd / "modules" / "pfam_adhesion.json").read_text())
    assert rec["run_state"] == "unavailable" and "no active family" in rec["note"]
    assert {r["state"] for r in read(wd, "pfam_adhesion").values() if r["id"] != "BAD"} == {
        "unavailable"
    }


def test_pfam_command_refuses_a_domain_table_from_another_proteome(tmp_path, fasta, capsys):
    dom = tmp_path / "d.domtbl"
    write_domtbl(
        dom, "OTHER1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    code = main(
        _pfam_args(tmp_path, fasta, _family_table(tmp_path, active=True), dom, tmp_path / "wd")
    )
    assert code == 2 and "not in the FASTA" in capsys.readouterr().err


def test_pfam_command_refuses_an_unfinished_domain_table(tmp_path, fasta, capsys):
    dom = tmp_path / "d.domtbl"
    dom.write_text(
        OPTIONS
        + "XP_1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    code = main(
        _pfam_args(tmp_path, fasta, _family_table(tmp_path, active=True), dom, tmp_path / "wd")
    )
    assert code == 2 and "[ok]" in capsys.readouterr().err


def test_allergen_command_refuses_a_blast_table_from_another_proteome(tmp_path, fasta, capsys):
    blast = tmp_path / "b.tsv"
    blast.write_text("OTHER\tAsp_f_1.0101|11\t82.0\t90\t120\t100\t150\t1e-20\n")
    ref = tmp_path / "a.faa"
    ref.write_text(">Asp_f_1.0101|11\nMKT\n")
    code = main(
        [
            "allergen",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--blast",
            str(blast),
            "--allergen-fasta",
            str(ref),
            "--blast-version",
            "2.14.0+",
        ]
    )
    assert code == 2 and "not in the FASTA" in capsys.readouterr().err


def test_partial_results_set_the_run_state_partial(tmp_path, fasta):
    tab = tmp_path / "r.tsv"
    tab.write_text("protein\trep_period\trep_n_copies\trep_coverage\nXP_1\t5\t3.0\t0.5\n")
    wd = tmp_path / "wd"
    long_fasta = tmp_path / "long.faa"
    long_fasta.write_text(">XP_1\n" + "M" * 100 + "\n>XP_2\n" + "M" * 100 + "\n")
    assert (
        main(["repeat02", "--fasta", str(long_fasta), "--workdir", str(wd), "--table", str(tab)])
        == 0
    )
    rec = json.loads((wd / "modules" / "repeat02.json").read_text())
    assert rec["run_state"] == "partial" and "1 protein(s) have no result" in rec["note"]


def test_a_lookup_that_matches_no_applicable_protein_is_refused(tmp_path, fasta, capsys):
    cand = tmp_path / "c.tsv"
    cand.write_text("protein_id\ttier\tcys_frac\nCIMG_9\tcys_rich_sp_unassigned\t0.1\n")
    code = main(
        [
            "cys",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--taxon",
            "246410",
            "--candidates",
            str(cand),
        ]
    )
    assert code == 2 and "no protein has a result" in capsys.readouterr().err


def test_signalp_file_from_a_non_eukarya_run_is_refused(tmp_path, fasta, capsys):
    res = tmp_path / "p.txt"
    res.write_text(
        "# SignalP-6.0\tOrganism: Other\tTimestamp: 1\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\nXP_1\tSP\t0.1\t0.9\t\n"
    )
    code = main(
        [
            "signalp",
            "--fasta",
            str(fasta),
            "--workdir",
            str(tmp_path / "wd"),
            "--results",
            str(res),
            "--signalp-version",
            "6",
        ]
    )
    assert code == 2 and "not a Eukarya run" in capsys.readouterr().err


def _signalp_file(tmp_path, ids):
    res = tmp_path / "sp.txt"
    body = "".join(f"{i}\tSP\t0.1\t0.9\t\n" for i in ids)
    res.write_text(
        "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 1\n# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
        + body
    )
    return res


def test_a_table_with_some_ids_outside_the_fasta_is_refused_and_names_the_file(
    tmp_path, fasta, capsys
):
    res = _signalp_file(tmp_path, ["XP_1", "STRAY"])
    args = ["--fasta", str(fasta), "--workdir", str(tmp_path / "wd")]
    assert main(["signalp", *args, "--results", str(res), "--signalp-version", "6"]) == 2
    err = capsys.readouterr().err
    assert "sp.txt" in err and "1 of 2" in err and "not in the FASTA" in err and "'STRAY'" in err
    tab = tmp_path / "r.tsv"
    tab.write_text(
        "protein\trep_period\trep_n_copies\trep_coverage\nXP_1\t5\t3.0\t0.5\nSTRAY\t0\t0\t0\n"
    )
    assert main(["repeat14", *args, "--table", str(tab)]) == 2
    assert "r.tsv" in capsys.readouterr().err
    tm = tmp_path / "t.tsv"
    tm.write_text("protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nSTRAY\t11\t0\t0\t0\to\n")
    assert main(["tm", *args, "--table", str(tm)]) == 2
    err = capsys.readouterr().err
    assert "t.tsv" in err and "not in the FASTA" in err
    assert not (tmp_path / "wd" / "modules").exists()


def test_the_foreign_id_message_names_the_domain_and_blast_files(tmp_path, fasta, capsys):
    dom = tmp_path / "d.domtbl"
    write_domtbl(
        dom, "OTHER1 - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    )
    main(_pfam_args(tmp_path, fasta, _family_table(tmp_path, active=True), dom, tmp_path / "wd"))
    assert "d.domtbl" in capsys.readouterr().err


@pytest.mark.parametrize("taxon", ["0", "1", "-5"])
def test_taxon_zero_and_one_are_refused(tmp_path, fasta, capsys, taxon):
    cand = tmp_path / "c.tsv"
    cand.write_text("protein_id\ttier\tcys_frac\nXP_1\tcys_rich_sp_unassigned\t0.1\n")
    base = ["cys", "--fasta", str(fasta), "--workdir", str(tmp_path / "wd")]
    assert main([*base, "--taxon", taxon, "--candidates", str(cand)]) == 2
    assert "not a taxon of an organism" in capsys.readouterr().err
    tmap = tmp_path / "m.tsv"
    tmap.write_text("XP_1\t246410\nXP_2\t1\nBAD\t246410\n")
    assert main([*base, "--taxon-map", str(tmap), "--candidates", str(cand)]) == 2
    assert "not a taxon of an organism" in capsys.readouterr().err


def test_cys_and_expression_apply_only_to_the_applicable_taxa(tmp_path, fasta):
    cand = tmp_path / "c.tsv"
    cand.write_text("protein_id\ttier\tcys_frac\nXP_1\tcys_rich_sp_unassigned\t0.1\n")
    expr = tmp_path / "e.tsv"
    expr.write_text(
        "gene_id\tlog2fc_48h\tpadj_48h\tlog2fc_8d\nCIMG_1\t1.5\t\t\nCIMG_2\t-0.5\t0.01\t2.0\n"
    )
    pmap = tmp_path / "pm.tsv"
    pmap.write_text("protein_id\tgene_id\nXP_1\tCIMG_1\nXP_2\tCIMG_2\n")
    base = ["--fasta", str(fasta), "--workdir", str(tmp_path / "wd")]
    for taxon, state in (("746128", "not_applicable"), ("246410", "ok")):
        assert main(["cys", *base, "--taxon", taxon, "--candidates", str(cand)]) == 0
        assert read(tmp_path / "wd", "cys_rich")["XP_1"]["state"] == state
        args = ["expression", *base, "--taxon", taxon, "--table", str(expr)]
        assert main([*args, "--protein-map", str(pmap)]) == 0
        assert read(tmp_path / "wd", "expression")["XP_1"]["state"] == state
    # exact match: a species the user lists is applicable, its neighbour is not
    assert (
        main(
            ["cys", *base, "--taxon", "246411", "--applicable-taxa", "246411"]
            + ["--candidates", str(cand)]
        )
        == 0
    )
    assert read(tmp_path / "wd", "cys_rich")["XP_1"]["state"] == "ok"
    row = read(tmp_path / "wd", "expression")["XP_1"]
    assert (row["log2fc"], row["padj"], row["log2fc_8d"]) == ("1.5", "", "")  # blanks pass


def test_the_expression_table_refuses_a_bad_number_and_a_missing_column(tmp_path, fasta, capsys):
    pmap = tmp_path / "pm.tsv"
    pmap.write_text("protein_id\tgene_id\nXP_1\tCIMG_1\n")
    base = ["expression", "--fasta", str(fasta), "--workdir", str(tmp_path / "wd")]
    tail = ["--taxon", "246410", "--protein-map", str(pmap)]
    bad = tmp_path / "e.tsv"
    bad.write_text("gene_id\tlog2fc_48h\tpadj_48h\tlog2fc_8d\nCIMG_1\tnan\t0.1\t1\n")
    assert main([*base, "--table", str(bad), *tail]) == 2
    assert "not finite" in capsys.readouterr().err
    bad.write_text("gene_id\tlog2fc_48h\tpadj_48h\tlog2fc_8d\nCIMG_1\t1\t-0.1\t1\n")
    assert main([*base, "--table", str(bad), *tail]) == 2
    assert "negative" in capsys.readouterr().err
    bad.write_text("gene_id\tlog2fc_48h\nCIMG_1\t1\n")
    assert main([*base, "--table", str(bad), *tail]) == 2
    assert "missing column" in capsys.readouterr().err


def test_repeat_and_tm_identity_do_not_depend_on_the_proteome(tmp_path, fasta):
    other = tmp_path / "o.faa"
    other.write_text(">Z1\n" + "M" * 100 + "\n")
    rep = tmp_path / "r.tsv"
    rep.write_text("protein\trep_period\trep_n_copies\trep_coverage\nZ1\t0\t0\t0\nXP_1\t0\t0\t0\n")
    tm = tmp_path / "t.tsv"
    tm.write_text("protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nXP_1\t11\t0\t0\t0\to\n")
    tm2 = tmp_path / "t2.tsv"
    tm2.write_text("protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nZ1\t11\t0\t0\t0\to\n")
    ids = {}
    for key, f, rtab, ttab in (("a", fasta, rep, tm), ("b", other, rep, tm2)):
        wd = tmp_path / key
        if key == "a":
            rep.write_text("protein\trep_period\trep_n_copies\trep_coverage\nXP_1\t0\t0\t0\n")
        else:
            rep.write_text("protein\trep_period\trep_n_copies\trep_coverage\nZ1\t0\t0\t0\n")
        assert (
            main(["repeat02", "--fasta", str(f), "--workdir", str(wd), "--table", str(rtab)]) == 0
        )
        assert main(["tm", "--fasta", str(f), "--workdir", str(wd), "--table", str(ttab)]) == 0
        ids[key] = [
            json.loads((wd / "modules" / f"{m}.json").read_text()) for m in ("repeat02", "tm")
        ]
    for x, y in zip(ids["a"], ids["b"], strict=True):
        assert (x["version"], x["params_hash"], x["artefact_hash"]) == (
            y["version"],
            y["params_hash"],
            y["artefact_hash"],
        )


def test_a_protein_shorter_than_the_detector_minimum_is_not_called_not_an_error(tmp_path):
    f = tmp_path / "s.faa"
    f.write_text(">S1\nMKTAYIAKQRQ\n>L1\n" + "M" * 100 + "\n")
    tab = tmp_path / "r.tsv"
    tab.write_text("protein\trep_period\trep_n_copies\trep_coverage\nL1\t0\t0\t0\n")
    wd = tmp_path / "wd"
    assert main(["repeat14", "--fasta", str(f), "--workdir", str(wd), "--table", str(tab)]) == 0
    assert read(wd, "repeat14")["S1"]["call"] == "not_called"
    assert json.loads((wd / "modules" / "repeat14.json").read_text())["run_state"] == "ok"
