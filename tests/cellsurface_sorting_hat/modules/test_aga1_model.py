"""The AGA1 thresholds table and the ``aga1_model`` wrapper (SYNTHETIC score tables; the committed table is checked against its record)."""

import csv
import gzip
import json
import sys
from pathlib import Path

import pytest

from cellsurface_sorting_hat.modules import aga1_model, aga1_thresholds
from cellsurface_sorting_hat.modules.cli import main
from cellsurface_sorting_hat.taxonomy import Lineage

ROOT = Path(__file__).resolve().parents[3]
TABLE = ROOT / "data/sorting_hat/aga1_model_thresholds.tsv"
RECORD = ROOT / "docs/reports/data/aga1_proteome_test_2026-10-25/thresholds.tsv"

# a tiny taxonomy: 1 root > 4751 fungi > 147537 Saccharomycotina (the real NCBI id) > 4932 S. cerevisiae > 559292 strain; 4751 > 5052 Aspergillus > 746128 A. fumigatus
NODES = {1: 1, 4751: 1, 147537: 4751, 4932: 147537, 559292: 4932, 5052: 4751, 746128: 5052}
FASTA = ">P_all x\nMKSTA\n>P_two y\nMKSTA\n>P_one z\nMKSTA\n>P_none w\nMKSTA\n>P_missing v\nMKSTA\n>BAD\nMK*T\n"
HEADER = aga1_thresholds.HEADER


def table_text(rows):
    return "\n".join("\t".join(map(str, r)) for r in [HEADER, *rows]) + "\n"


def good_rows(extra=()):
    rows = []
    for taxon, name, cal, vals in (
        (147537, "Saccharomycotina", "yes", (0.6, 0.3, 0.1)),
        (1, "all other taxa", "no", (0.6, 0.3, 0.1)),
        *extra,
    ):
        for m, v in zip(aga1_thresholds.MODELS, vals, strict=True):
            rows.append([taxon, name, m, v, "basis", cal, "Owner", "2026-10-26", "evidence"])
    return rows


@pytest.fixture
def nodes(tmp_path):
    p = tmp_path / "nodes.dmp"
    p.write_text("".join(f"{t}\t|\t{par}\t|\tno rank\t|\n" for t, par in NODES.items()))
    return p


def write_table(tmp_path, rows):
    p = tmp_path / "thr.tsv"
    p.write_text(table_text(rows))
    return p


# ---- thresholds table
def test_committed_table_loads_and_matches_the_yeast_record():
    rows = aga1_thresholds.load_thresholds(TABLE)
    assert {r.clade_taxon for r in rows} == {1, 147537}
    with RECORD.open(newline="") as fh:
        rec = {r["model"]: float(r["score_threshold"]) for r in csv.DictReader(fh, delimiter="\t")}
    for r in rows:
        assert r.threshold == pytest.approx(
            rec[r.model], abs=1e-12
        )  # copied verbatim from the committed record
        assert (r.calibrated == "yes") == (
            r.clade_taxon == 147537
        )  # only the yeast row claims a calibration


@pytest.mark.parametrize(
    "mutate, msg",
    [
        (lambda rows: rows[:-1], "lacks rows"),
        (lambda rows: rows + [rows[0]], "duplicate row"),
        (
            lambda rows: [[*r[:3], 1.5, *r[4:]] if i == 0 else r for i, r in enumerate(rows)],
            r"in \(0, 1\)",
        ),
        (
            lambda rows: [[*r[:5], "maybe", *r[6:]] if i == 0 else r for i, r in enumerate(rows)],
            "yes or no",
        ),
        (
            lambda rows: [[*r[:6], "", *r[7:]] if i == 0 else r for i, r in enumerate(rows)],
            "signed_by is empty",
        ),
        (
            lambda rows: [
                [*r[:7], "26/10/2026", *r[8:]] if i == 0 else r for i, r in enumerate(rows)
            ],
            "YYYY-MM-DD",
        ),
        (
            lambda rows: [[*r[:2], "svm", *r[3:]] if i == 0 else r for i, r in enumerate(rows)],
            "not one of",
        ),
        (lambda rows: [[0, *r[1:]] if i == 0 else r for i, r in enumerate(rows)], "at least 1"),
        (
            lambda rows: [[*r[:4], "", *r[5:]] if i == 0 else r for i, r in enumerate(rows)],
            "basis and evidence",
        ),
    ],
)
def test_table_validation(tmp_path, mutate, msg):
    with pytest.raises(aga1_thresholds.ThresholdError, match=msg):
        aga1_thresholds.load_thresholds(write_table(tmp_path, mutate(good_rows())))


def test_table_header_and_empty(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("a\tb\n")
    with pytest.raises(aga1_thresholds.ThresholdError, match="header must be"):
        aga1_thresholds.load_thresholds(p)
    p.write_text("\t".join(HEADER) + "\n")
    with pytest.raises(aga1_thresholds.ThresholdError, match="no rows"):
        aga1_thresholds.load_thresholds(p)


def test_selection_by_lineage_deepest_clade_and_fallback(tmp_path, nodes):
    lin = Lineage.from_nodes_dmp(nodes)
    rows = aga1_thresholds.load_thresholds(
        write_table(tmp_path, good_rows(extra=[(4932, "S. cerevisiae", "yes", (0.7, 0.4, 0.2))]))
    )
    deep = aga1_thresholds.select_thresholds(
        rows, lin.ancestors(559292)
    )  # strain -> species row, deeper than Saccharomycotina
    assert {r.clade_taxon for r in deep.values()} == {4932} and deep["rf"].threshold == 0.4
    mid = aga1_thresholds.select_thresholds(rows, lin.ancestors(147537))
    assert {r.clade_taxon for r in mid.values()} == {147537}
    other = aga1_thresholds.select_thresholds(rows, lin.ancestors(746128))
    assert {r.clade_taxon for r in other.values()} == {1} and {
        r.calibrated for r in other.values()
    } == {"no"}


def test_selection_stops_without_a_fallback(tmp_path, nodes):
    lin = Lineage.from_nodes_dmp(nodes)
    rows = aga1_thresholds.load_thresholds(
        write_table(tmp_path, [r for r in good_rows() if r[0] != 1])
    )
    with pytest.raises(aga1_thresholds.ThresholdError, match="no threshold row"):
        aga1_thresholds.select_thresholds(rows, lin.ancestors(746128))


# ---- wrapper
def write_scores(tmp_path, values):
    paths = {}
    for m in aga1_thresholds.MODELS:
        p = tmp_path / f"scores_{m}.tsv"
        p.write_text(
            "accession\t"
            + f"score_{m}\n"
            + "".join(f"{acc}\t{v[m]}\n" for acc, v in values.items() if m in v)
        )
        paths[m] = p
    return paths


SCORES = {  # thresholds 0.6 / 0.3 / 0.1 for the yeast row
    "P_all": {"logreg": 0.9, "rf": 0.9, "gboost": 0.9},
    "P_two": {"logreg": 0.9, "rf": 0.3, "gboost": 0.05},  # rf exactly at its threshold counts
    "P_one": {"logreg": 0.9, "rf": 0.1, "gboost": 0.05},
    "P_none": {"logreg": 0.1, "rf": 0.1, "gboost": 0.05},
    "P_missing": {"logreg": 0.9, "rf": 0.9},  # no gboost score
}


def run(tmp_path, nodes, taxon, wd="wd", scores=SCORES, thr=None):
    fa = tmp_path / "p.faa"
    fa.write_text(FASTA)
    sc = write_scores(tmp_path, scores)
    table = thr or write_table(tmp_path, good_rows())
    code = main(
        ["aga1_model", "--fasta", str(fa), "--workdir", str(tmp_path / wd), "--taxon", str(taxon), "--taxdump", str(nodes), "--thresholds", str(table),
         "--scores", *[f"{m}={p}" for m, p in sc.items()]]
    )  # fmt: skip
    return code, tmp_path / wd


def read_module(wd):
    with gzip.open(wd / "modules" / "aga1_model.tsv.gz", "rt") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_two_of_three_rule_and_module_fields(tmp_path, nodes):
    code, wd = run(tmp_path, nodes, 559292)
    assert code == 0
    m = read_module(wd)
    assert [m[p]["call"] for p in ("P_all", "P_two", "P_one", "P_none")] == [
        "called",
        "called",
        "not_called",
        "not_called",
    ]
    assert [m[p]["n_models_call"] for p in ("P_all", "P_two", "P_one", "P_none")] == [
        "3",
        "2",
        "1",
        "0",
    ]
    assert (
        m["P_all"]["threshold_clade"] == "147537"
        and m["P_all"]["calibrated"] == "yes"
        and m["P_all"]["thr_rf"] == "0.3"
    )
    assert (
        m["P_missing"] == {**m["P_missing"], "state": "error"} and m["BAD"]["state"] == "na_invalid"
    )
    record = json.loads((wd / "modules" / "aga1_model.json").read_text())
    assert (
        record["run_state"] == "partial"
        and "2 of 3" in record["params"]["rule"]
        or "at least 2" in record["params"]["rule"]
    )


def test_fallback_row_outside_the_calibrated_clade(tmp_path, nodes):
    code, wd = run(tmp_path, nodes, 746128)
    m = read_module(wd)
    assert code == 0 and m["P_all"]["threshold_clade"] == "1" and m["P_all"]["calibrated"] == "no"


def test_unknown_taxon_stops(tmp_path, nodes, capsys):
    code, _ = run(tmp_path, nodes, 999999)
    assert code == 2 and "999999" in capsys.readouterr().err


def test_identity_changes_with_a_threshold_or_a_score_table(tmp_path, nodes):
    _, a = run(tmp_path, nodes, 559292, wd="a")
    base = json.loads((a / "modules" / "aga1_model.json").read_text())["artefact_hash"]
    rows = good_rows()
    rows[0][3] = 0.61  # one threshold of the clade that is used
    _, b = run(tmp_path, nodes, 559292, wd="b", thr=write_table(tmp_path, rows))
    assert json.loads((b / "modules" / "aga1_model.json").read_text())["artefact_hash"] != base
    changed = {**SCORES, "P_none": {"logreg": 0.11, "rf": 0.1, "gboost": 0.05}}
    _, c = run(tmp_path, nodes, 559292, wd="c", scores=changed)
    assert json.loads((c / "modules" / "aga1_model.json").read_text())["artefact_hash"] != base
    unrelated = good_rows(
        extra=[(5052, "Aspergillus", "yes", (0.7, 0.4, 0.2))]
    )  # a clade that is not used for this taxon
    _, d = run(tmp_path, nodes, 559292, wd="d", thr=write_table(tmp_path, unrelated))
    assert json.loads((d / "modules" / "aga1_model.json").read_text())["artefact_hash"] == base


@pytest.mark.parametrize(
    "text, msg",
    [
        ("accession\tscore_logreg\nP_all\t1.5\n", "outside"),
        ("accession\tscore_logreg\nP_all\tx\n", "not a number"),
        ("accession\tscore_logreg\nP_all\t0.5\nP_all\t0.6\n", "duplicate"),
        ("id\tscore_logreg\nP_all\t0.5\n", "needs the columns"),
        ("accession\tscore_logreg\n", "no rows"),
    ],
)
def test_score_table_stops(tmp_path, text, msg):
    p = tmp_path / "s.tsv"
    p.write_text(text)
    with pytest.raises(aga1_model.ScoreFormatError, match=msg):
        aga1_model.parse_scores(p, "logreg")


def test_scores_option_needs_all_three_models(tmp_path, nodes, capsys):
    fa = tmp_path / "p.faa"
    fa.write_text(FASTA)
    sc = write_scores(tmp_path, SCORES)
    code = main(["aga1_model", "--fasta", str(fa), "--workdir", str(tmp_path / "w"), "--taxon", "559292", "--taxdump", str(nodes), "--thresholds", str(TABLE),
                 "--scores", f"logreg={sc['logreg']}", f"logreg={sc['logreg']}", f"rf={sc['rf']}"])  # fmt: skip
    assert code == 2 and "exactly" in capsys.readouterr().err


def test_wrapper_imports_no_scikit_learn_or_matplotlib():
    import subprocess

    code = (
        "import sys; import cellsurface_sorting_hat.modules.aga1_model, cellsurface_sorting_hat.modules.cli; "
        "assert 'sklearn' not in sys.modules and 'matplotlib' not in sys.modules"
    )
    env = {"PYTHONPATH": str(ROOT / "src"), "PATH": "/usr/bin:/bin"}
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
