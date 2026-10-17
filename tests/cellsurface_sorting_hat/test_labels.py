"""Task 3.7 (D1, D14): ``primary_label``, one derived label per protein.

Precedence, first match wins: cell_wall_adhesion_candidate, surface_attachment_candidate,
other_surface_no_mechanism, other_not_surface, mechanism_without_class, not_assessable.
"""

import csv
import gzip
import itertools

from engine_helpers import R0, base_modules, config_with, ok, run_eval, table
from test_cli import run_cli

CWAC = "cell_wall_adhesion_candidate"
SAC = "surface_attachment_candidate"
OSN = "other_surface_no_mechanism"
ONS = "other_not_surface"
MWC = "mechanism_without_class"
NAS = "not_assessable"
ORDER = (CWAC, SAC, OSN, ONS, MWC, NAS)


def label_of(modules, cfg=None):
    from cellsurface_sorting_hat.engine import load_config
    from cellsurface_sorting_hat.labels import primary_labels

    cfg = cfg or load_config()
    records = list(run_eval(modules, cfg=cfg).values())
    return primary_labels(records, cfg)["P"]


def test_the_label_vocabulary_and_order():
    from cellsurface_sorting_hat.labels import LABELS

    assert tuple(LABELS) == ORDER


def test_no_signal_peptide_and_a_repeat_call_is_a_mechanism_without_class():
    mods = base_modules(step1="not_called", repeat02=table("repeat02", {"P": ok(call="called")}))
    assert label_of(mods) == (MWC, "tandem_repeat_protein")


def test_a_secreted_protein_with_only_a_hydrophobin_hit_is_a_surface_attachment_candidate():
    mods = base_modules(
        step1="called", pfam_hydrophobin=table("pfam_hydrophobin", {"P": ok(hit="1")})
    )
    assert label_of(mods) == (SAC, "hydrophobin_domain")


def test_a_secreted_protein_with_only_an_hsba_hit_is_a_surface_attachment_candidate():
    mods = base_modules(step1="called", pfam_hsba=table("pfam_hsba", {"P": ok(hit="1")}))
    assert label_of(mods) == (SAC, "hsba_domain")


def test_probe_p2_an_unknown_mechanism_input_gives_not_assessable():
    mods = base_modules(step1="called", pfam_hsba=table("pfam_hsba", {"P": {"state": "error"}}))
    assert label_of(mods) == (NAS, "hsba_domain")


def test_a_secreted_protein_with_only_the_antigen_rank_is_a_mechanism_without_class():
    mods = base_modules(
        step1="called", antigen_lookup=table("antigen_lookup", {"P": ok(percentile="7")})
    )
    assert label_of(mods) == (MWC, "cocci_specificity_rank_top15")


def test_a_secreted_protein_with_a_mechanism_that_no_composite_reads_is_a_mechanism_without_class(
    tmp_path,
):
    cfg = config_with(
        tmp_path,
        extra_calls=[{"name": "sowgp_ortholog", "expr": {"flag": "sowgp.hit"}}],
        mechanism_extra=["sowgp_ortholog"],
    )
    mods = base_modules(step1="called", sowgp=table("sowgp", {"P": ok(hit="1")}))
    assert label_of(mods, cfg) == (MWC, "sowgp_ortholog")


def test_the_ordinary_labels():
    repeat = {"repeat02": table("repeat02", {"P": ok(call="called")})}
    assert label_of(base_modules(step1="called", **repeat))[0] == CWAC
    assert label_of(base_modules(step1="called"))[0] == OSN
    assert label_of(base_modules(step1="not_called"))[0] == ONS
    broken = {R0: table(R0, {"P": {"state": "error"}})}
    assert label_of(base_modules(**broken))[0] == NAS
    assert label_of(base_modules(**broken, **repeat))[0] == MWC  # a mechanism holds, SP unknown


def test_cell_wall_adhesion_candidate_wins_over_surface_attachment_candidate():
    mods = base_modules(
        step1="called",
        repeat02=table("repeat02", {"P": ok(call="called")}),
        pfam_hsba=table("pfam_hsba", {"P": ok(hit="1")}),
    )
    assert label_of(mods) == (CWAC, "")


def test_a_run_without_the_default_gate_labels_every_protein_not_assessable():
    mods = base_modules()
    del mods[R0]
    assert label_of(mods)[0] == NAS


def test_every_combination_of_inputs_gets_exactly_one_label():
    """FD-A5: the label is a total function of the call values, for every input state."""
    from cellsurface_sorting_hat.engine import load_config
    from cellsurface_sorting_hat.labels import LABELS, primary_labels

    states = {
        "r0": ("called", "not_called", None),
        "rep": ("called", "not_called", None),
        "adh": (("0", "0"), ("1", "0"), ("1", "1"), None),
        "hyd": ("1", "0", None),
        "hsba": ("1", "0", None),
        "ag": ("7", "50", "na", None),
    }
    combos = list(itertools.product(*states.values()))
    ids = [f"P{i}" for i in range(len(combos))]
    tables = {n: {} for n in ("r0", "rep", "r14", "adh", "hyd", "hsba", "ag", "pa", "ah")}
    for pid, (r0, rep, adh, hyd, hsba, ag) in zip(ids, combos, strict=True):
        err = {"state": "error"}
        tables["r0"][pid] = ok(call=r0) if r0 else err
        tables["rep"][pid] = ok(call=rep) if rep else err
        tables["r14"][pid] = ok(call="not_called")
        tables["adh"][pid] = ok(hit=adh[0], hit_narrow=adh[1]) if adh else err
        tables["hyd"][pid] = ok(hit=hyd) if hyd else err
        tables["hsba"][pid] = ok(hit=hsba) if hsba else err
        tables["ag"][pid] = (
            {"state": "not_applicable"} if ag == "na" else ok(percentile=ag) if ag else err
        )
        tables["pa"][pid] = ok(hit="0")
        tables["ah"][pid] = ok(identity="0", coverage="0")
    mods = {
        R0: table(R0, tables["r0"]),
        "repeat02": table("repeat02", tables["rep"]),
        "repeat14": table("repeat14", tables["r14"]),
        "pfam_adhesion": table("pfam_adhesion", tables["adh"]),
        "pfam_hydrophobin": table("pfam_hydrophobin", tables["hyd"]),
        "pfam_hsba": table("pfam_hsba", tables["hsba"]),
        "pfam_allergen": table("pfam_allergen", tables["pa"]),
        "antigen_lookup": table("antigen_lookup", tables["ag"]),
        "allergen_homology": table("allergen_homology", tables["ah"]),
    }
    cfg = load_config()
    records = list(run_eval(mods, ids=ids, cfg=cfg).values())
    labels = primary_labels(records, cfg)
    assert sorted(labels) == sorted(ids)
    assert {label for label, _ in labels.values()} <= set(LABELS)
    value = {(r.protein, r.call, r.variant): r.value for r in records}
    for pid in ids:
        label = labels[pid][0]
        cwac, sac = value[(pid, CWAC, "R0")], value[(pid, SAC, "R0")]
        if cwac == "called":
            assert label == CWAC
        elif sac == "called":
            assert label == SAC
        else:
            assert label not in (CWAC, SAC)
        if label == OSN:
            assert value[(pid, OSN, "R0")] == "called"
        if label == ONS:
            assert value[(pid, ONS, "R0")] == "called"


# ---- outputs ------------------------------------------------------------------------------------


def test_calls_wide_ends_with_primary_label_and_primary_basis(tmp_path, write_module, nodes_dmp):
    code, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    assert code == 0
    with gzip.open(out / "calls.wide.tsv.gz", "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert list(rows[0])[-2:] == ["primary_label", "primary_basis"]
    assert len(rows) == 5 and all(r["primary_label"] in ORDER for r in rows)
    by = {r["protein"]: r for r in rows}
    assert by["SOW1"]["primary_label"] == CWAC  # repeat02 called and a signal peptide
    assert by["BAD1"]["primary_label"] == NAS  # invalid protein: no module has a value


def test_calls_long_is_not_given_a_label_column(tmp_path, write_module, nodes_dmp):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    with gzip.open(out / "calls.long.tsv.gz", "rt") as fh:
        assert "primary_label" not in fh.readline()


def test_the_report_counts_one_label_per_protein(tmp_path, write_module, nodes_dmp):
    _, out, _ = run_cli(tmp_path, write_module, nodes_dmp)
    text = (out / "report.md").read_text()
    section = text.split("## Primary label")[1].split("\n## ")[0]
    counts = {}
    for line in section.splitlines():
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) == 2 and cells[1].isdigit():
            counts[cells[0]] = int(cells[1])
    assert set(counts) == set(ORDER)  # every label has a row, also when its count is 0
    assert sum(counts.values()) == 5  # the five proteins of the toy FASTA, invalid one included
