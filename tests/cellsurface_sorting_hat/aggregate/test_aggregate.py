"""The cross-proteome aggregator on SYNTHETIC engine outputs (run directories written by the helper below)."""

import csv
import gzip
import hashlib
import json
from pathlib import Path

import pytest

from cellsurface_sorting_hat.aggregate import cli, groups, inputs, tables

# a tiny taxonomy: 1 root; 100 class A (rank class) > 110 species s1, 111 species s2; 200 class B > 210 species s3; 300 no-class branch > 310 species s4
PARENT = {1: 1, 100: 1, 110: 100, 111: 100, 200: 1, 210: 200, 300: 1, 310: 300}
RANK = {
    1: "no rank",
    100: "class",
    110: "species",
    111: "species",
    200: "class",
    210: "species",
    300: "phylum",
    310: "species",
}
NAMES = {100: "ClassA", 200: "ClassB"}


def write_run(
    base,
    label,
    taxon,
    calls,
    config="cfg1",
    n_invalid=0,
    absent=(),
    gate="step1_rule@R0",
    variants=("R0",),
    extra_rows=(),
):
    """`calls` = {call: [(value, status), ...]} (ungated); a call named `g:name` is per_variant and gets one row per variant."""
    d = Path(base) / label / "out"
    d.mkdir(parents=True)
    rows = []
    n = max(len(v) for v in calls.values())
    for call, recs in calls.items():
        for i, (value, status) in enumerate(recs):
            if call.startswith("g:"):
                for v in variants:
                    rows.append((f"P{i}", call[2:], v, value, status))
            else:
                rows.append((f"P{i}", call, "", value, status))
    rows += list(extra_rows)
    with gzip.open(d / "calls.long.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            [
                "protein",
                "call",
                "variant",
                "value",
                "status",
                "status_basis",
                "other_basis",
                "derived_from",
            ]
        )
        for r in rows:
            w.writerow([*r, "", "", ""])
    with gzip.open(d / "proteins.tsv.gz", "wt", newline="") as fh:
        fh.write(
            "id\tsha256\ttaxon\tstate\n" + "".join(f"P{i}\tx\t{taxon}\tok\n" for i in range(n))
        )
    info = {
        "config_sha256": config,
        "default_gate": gate,
        "n_proteins": n,
        "n_invalid": n_invalid,
        "taxa": {str(taxon): n},
        "absent_modules": list(absent),
        "version": "t",
    }
    (d / "run.json").write_text(json.dumps(info))
    return d


def calls_for(called, notcalled, na=0, status="unvalidated", extra_na=0):
    recs = (
        [("called", status)] * called
        + [("not_called", status)] * notcalled
        + [("not_assessable", status)] * na
    )
    return recs


@pytest.fixture
def three(tmp_path):
    """s1 and s2 in class A (call X measured in s1 only), s3 in class B (no assessable protein for call Y)."""
    d1 = write_run(
        tmp_path, "s1", 110, {"X": calls_for(2, 8, status="smoke"), "Y": calls_for(1, 9)}
    )
    d2 = write_run(
        tmp_path, "s2", 111, {"X": calls_for(4, 16), "Y": calls_for(0, 20)}, absent=("modA",)
    )
    d3 = write_run(
        tmp_path,
        "s3",
        210,
        {"X": calls_for(1, 9), "Y": calls_for(0, 0, na=10)},
        absent=("modA", "modB"),
    )
    return [(str(d1), "s1"), (str(d2), "s2"), (str(d3), "s3")]


# ---- inputs
def test_read_run_reads_one_variant_only_and_ungated_calls(tmp_path):
    d = write_run(
        tmp_path, "s1", 110, {"U": calls_for(1, 1), "g:G": calls_for(1, 3)}, variants=("R0", "R1")
    )
    run = inputs.read_run(d, "s1")
    assert (
        run.variant == "R0" and len(run.calls["G"]) == 4 and len(run.calls["U"]) == 2
    )  # the R1 rows are not counted
    assert len(inputs.read_run(d, "s1", variant="R1").calls["G"]) == 4


def test_read_run_failures(tmp_path):
    d = write_run(tmp_path, "s1", 110, {"U": calls_for(1, 1)})
    (d / "run.json").write_text("{}")
    with pytest.raises(inputs.AggregateError, match="key 'config_sha256'"):
        inputs.read_run(d, "s1")
    (d / "calls.long.tsv.gz").unlink()
    with pytest.raises(inputs.AggregateError, match="not found"):
        inputs.read_run(d, "s1")
    multi = write_run(tmp_path, "m", 110, {"U": calls_for(1, 1)})
    info = json.loads((multi / "run.json").read_text())
    info["taxa"] = {"110": 1, "111": 1}
    (multi / "run.json").write_text(json.dumps(info))
    with pytest.raises(inputs.AggregateError, match="2 taxa"):
        inputs.read_run(multi, "m")
    bad = write_run(tmp_path, "b", 110, {"U": [("maybe", "unvalidated")]})
    with pytest.raises(inputs.AggregateError, match="unknown call value"):
        inputs.read_run(bad, "b")
    short = write_run(tmp_path, "c", 110, {"U": calls_for(1, 1)})
    info = json.loads((short / "run.json").read_text())
    info["n_proteins"] = 9
    (short / "run.json").write_text(json.dumps(info))
    with pytest.raises(inputs.AggregateError, match="proteins.tsv.gz has 2 rows"):
        inputs.read_run(short, "c")


def test_duplicate_labels_and_mixed_configs(tmp_path):
    a = write_run(tmp_path, "a", 110, {"U": calls_for(1, 1)}, config="c1")
    b = write_run(tmp_path, "b", 111, {"U": calls_for(1, 1)}, config="c2")
    with pytest.raises(inputs.AggregateError, match="duplicate proteome label"):
        inputs.read_runs([(str(a), "x"), (str(b), "x")])
    with pytest.raises(inputs.AggregateError, match="2 different categories.yaml"):
        inputs.read_runs([(str(a), "a"), (str(b), "b")])
    assert len(inputs.read_runs([(str(a), "a"), (str(b), "b")], allow_mixed_config=True)) == 2


def test_manifest_and_default_labels(tmp_path):
    m = tmp_path / "m.tsv"
    m.write_text("run_dir\tlabel\n/x/Scer/out\t\n/x/Calb\tCalb1\n")
    assert inputs.read_manifest(m) == [("/x/Scer/out", "Scer"), ("/x/Calb", "Calb1")]
    m.write_text("dir\n/x\n")
    with pytest.raises(inputs.AggregateError, match="run_dir"):
        inputs.read_manifest(m)


# ---- groups
def runs_of(three):
    return inputs.read_runs(three)


def test_lineage_groups_names_and_no_rank_above(three, tmp_path):
    d4 = write_run(tmp_path, "s4", 310, {"X": calls_for(1, 1)})
    runs = inputs.read_runs([*three, (str(d4), "s4")])
    g = groups.assign_groups(runs, PARENT, RANK, "class", NAMES)
    assert (
        g["s1"] == g["s2"] == ("ClassA", "lineage")
        and g["s3"] == ("ClassB", "lineage")
        and g["s4"] == ("no class above", "lineage")
    )
    assert (
        groups.assign_groups(runs, PARENT, RANK, "class")["s1"][0] == "taxon:100"
    )  # without names.dmp


def test_group_table_overrides_the_lineage_group_for_descendants(three):
    g = groups.assign_groups(
        runs_of(three), PARENT, RANK, "class", NAMES, {100: "Mine", 111: "Special"}
    )
    assert (
        g["s1"] == ("Mine", "table")
        and g["s2"] == ("Special", "table")
        and g["s3"] == ("ClassB", "lineage")
    )


def test_group_errors(three, tmp_path):
    runs = runs_of(three)
    runs[0].taxon = 424242
    with pytest.raises(inputs.AggregateError, match="not in nodes.dmp"):
        groups.assign_groups(runs, PARENT, RANK, "class")
    t = tmp_path / "g.tsv"
    for text, msg in (
        ("a\tb\n", "needs the columns"),
        ("taxon_id\tgroup\nx\tG\n", "not an integer"),
        ("taxon_id\tgroup\n1\t\n", "empty group"),
        ("taxon_id\tgroup\n1\tA\n1\tB\n", "twice"),
        ("taxon_id\tgroup\n", "no rows"),
    ):
        t.write_text(text)
        with pytest.raises(inputs.AggregateError, match=msg):
            groups.read_group_table(t)


# ---- tables
def test_counts_by_hand_na_denominator_and_status(three):
    runs = runs_of(three)
    g = groups.assign_groups(runs, PARENT, RANK, "class", NAMES)
    rows = {(r["proteome"], r["call"]): r for r in tables.count_rows(runs, g)}
    x1 = rows[("s1", "X")]
    assert (x1["called"], x1["not_called"], x1["assessable"], x1["proteins"]) == (
        2,
        8,
        10,
        10,
    ) and x1["called_per_1000_assessable"] == 200.0
    assert x1["status"] == "smoke" and x1["measured"] == "yes" and x1["unvalidated_marker"] == ""
    assert (
        rows[("s2", "X")]["unvalidated_marker"] == tables.MARK
        and rows[("s2", "X")]["measured"] == "no"
    )
    y3 = rows[("s3", "Y")]
    assert y3["assessable"] == 0 and y3["called_per_1000_assessable"] is None  # NA, never 0


def test_status_is_the_weakest_among_the_called_records_and_mixed_is_flagged(tmp_path):
    d = write_run(
        tmp_path,
        "s1",
        110,
        {
            "X": [("called", "smoke"), ("called", "unvalidated"), ("not_called", "smoke")],
            "Z": [("not_called", "estimated"), ("not_called", "smoke")],
        },
    )
    run = inputs.read_run(d, "s1")
    rows = {r["call"]: r for r in tables.count_rows([run], {"s1": ("G", "lineage")})}
    assert (
        rows["X"]["status"] == "unvalidated" and rows["X"]["status_mixed"] == "yes"
    )  # the called records decide the status
    assert (
        rows["Z"]["status"] == "smoke" and rows["Z"]["status_mixed"] == "yes"
    )  # none called: the not_called records


def test_group_rows_summaries_and_measured_only(three):
    runs = runs_of(three)
    g = groups.assign_groups(runs, PARENT, RANK, "class", NAMES)
    crows = tables.count_rows(runs, g)
    grows = {(r["group"], r["call"]): r for r in tables.group_rows(crows)}
    a = grows[("ClassA", "X")]
    assert (
        a["n_proteomes"] == 2
        and a["n_measured"] == 1
        and a["total_called"] == 6
        and a["total_assessable"] == 30
    )
    assert (
        a["mean_per_1000"] == pytest.approx(200.0)
        and a["min_per_1000"] == 200.0
        and a["max_per_1000"] == 200.0
        and a["pooled_per_1000"] == pytest.approx(200.0)
    )
    assert grows[("ClassB", "Y")]["n_rate"] == 0 and grows[("ClassB", "Y")]["mean_per_1000"] is None
    only = {(r["group"], r["call"]): r for r in tables.group_rows(crows, measured_only=True)}
    assert only[("ClassA", "X")]["n_proteomes"] == 1 and only[("ClassA", "X")]["total_called"] == 2
    assert (
        only[("ClassB", "X")]["n_proteomes"] == 0
        and only[("ClassB", "X")]["pooled_per_1000"] is None
    )


def test_warnings(three):
    runs = runs_of(three)
    g = groups.assign_groups(runs, PARENT, RANK, "class", NAMES)
    crows = tables.count_rows(runs, g)
    w = tables.warnings(runs, crows, tables.group_rows(crows), mixed_config=False)
    text = "\n".join(w)
    assert (
        "s3: no protein is assessable for 1 call(s) (rate NA): `Y`; modules without a table in that run: modA, modB"
        in text
    )
    assert "Group `ClassB` has 1 proteome(s)" in text and "Group `ClassA` has 2 proteome(s)" in text
    assert "categories.yaml files" not in text
    assert "different categories.yaml" in "\n".join(
        tables.warnings(runs, crows, tables.group_rows(crows), mixed_config=True)
    )


# ---- figures and the command
def nodes_and_names(tmp_path):
    n = tmp_path / "nodes.dmp"
    n.write_text("".join(f"{t}\t|\t{p}\t|\t{RANK[t]}\t|\n" for t, p in PARENT.items()))
    nm = tmp_path / "names.dmp"
    nm.write_text(
        "".join(f"{t}\t|\t{name}\t|\t\t|\tscientific name\t|\n" for t, name in NAMES.items())
    )
    return n, nm


def test_command_end_to_end_and_determinism(three, tmp_path):
    pytest.importorskip("matplotlib")
    n, nm = nodes_and_names(tmp_path)
    man = tmp_path / "m.tsv"
    man.write_text("run_dir\tlabel\n" + "".join(f"{d}\t{lab}\n" for d, lab in three))
    hashes = []
    for k in ("a", "b"):
        out = tmp_path / k
        assert (
            cli.main(
                [
                    "--manifest",
                    str(man),
                    "--taxdump",
                    str(n),
                    "--names",
                    str(nm),
                    "--rank",
                    "class",
                    "--out",
                    str(out),
                ]
            )
            == 0
        )
        assert {p.name for p in out.iterdir()} == {
            "proteomes.tsv",
            "counts_by_proteome.tsv",
            "counts_by_group.tsv",
            "groups.png",
            "heatmap.png",
            "status_coverage.png",
            "report.md",
        }
        hashes.append(
            {
                p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in out.iterdir()
                if p.suffix == ".png"
            }
        )
        assert all(p.stat().st_size > 1000 for p in out.iterdir() if p.suffix == ".png")
    assert hashes[0] == hashes[1]  # the figures are reproducible byte for byte
    report = (tmp_path / "a" / "report.md").read_text()
    for section in (
        "## Proteomes",
        "## Group summaries",
        "## Status coverage",
        "## Warnings",
        "## Known limits",
        "![status of each call",
        "ClassA",
        "ClassB",
    ):
        assert section in report
    with (tmp_path / "a" / "counts_by_proteome.tsv").open(newline="") as fh:
        rows = {(r["proteome"], r["call"]): r for r in csv.DictReader(fh, delimiter="\t")}
    assert (
        rows[("s3", "Y")]["called_per_1000_assessable"] == "NA"
        and rows[("s1", "X")]["called_per_1000_assessable"] == "200.0000"
    )


def test_command_stops_cleanly(three, tmp_path, capsys):
    n, nm = nodes_and_names(tmp_path)
    out = tmp_path / "exists"
    out.mkdir()
    base = ["--runs", three[0][0], "--taxdump", str(n)]
    assert cli.main([*base, "--out", str(out)]) == 2 and "exists" in capsys.readouterr().err
    assert (
        cli.main([*base, "--calls", "NOPE", "--out", str(tmp_path / "o2")]) == 2
        and "not in the runs" in capsys.readouterr().err
    )
    assert (
        cli.main(
            [
                "--runs",
                str(tmp_path / "nowhere"),
                "--taxdump",
                str(n),
                "--out",
                str(tmp_path / "o3"),
            ]
        )
        == 2
    )
    assert not (tmp_path / "o3").exists()  # nothing is written before the inputs are accepted


def test_low_assessable_cells_are_flagged_and_left_out_of_group_rates():
    def row(proteome, called, assessable, n):
        frac = assessable / n
        return {
            "proteome": proteome,
            "group": "G",
            "call": "C",
            "called": called,
            "assessable": assessable,
            "called_per_1000_assessable": 1000.0 * called / assessable,
            "assessable_fraction": frac,
            "low_assessable": "yes" if frac < tables.LOW_ASSESSABLE else "no",
            "measured": "no",
        }

    rows = [row("p1", 10, 100, 100), row("p2", 3, 3, 100)]  # p2: 3% assessable, rate 1,000
    g = tables.group_rows(rows)[0]
    assert g["n_low_assessable"] == 1 and g["n_rate"] == 1
    assert g["mean_per_1000"] == pytest.approx(100.0) and g["total_assessable"] == 100
