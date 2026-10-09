"""cys8_pattern: the eight-cysteine spacing match plus the R0 signal peptide condition.

Sequences are synthetic, built from the frozen spacing in data/sorting_hat/cys8_spacing.yaml.
No truth protein and no protein seen in the discovery search is used here.
"""

from pathlib import Path

import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules import cys8

SPACING = Path(__file__).resolve().parents[3] / "data/sorting_hat/cys8_spacing.yaml"


def make_seq(gaps, lead="MKLLVAAG", tail="GSTNPAD"):
    """gaps = residues between C1..C8 (seven values, doublets are 0)."""
    s = lead + "C"
    for g in gaps:
        s += "A" * g + "C"
    return s + tail


CLASS_II = [10, 0, 11, 16, 8, 0, 10]  # Seidl-Seiboth 2011 Table 3 row II
CLASS_I = [6, 0, 35, 22, 5, 0, 16]  # inside the Seidl-Seiboth 2011 class I ranges


def prot(pid, seq, state="ok"):
    return Protein(pid, seq, "sha", state, "", 0.0)


@pytest.fixture(scope="module")
def sets():
    return cys8.load_sets(SPACING)


def test_class_ii_sequence_matches_class_ii_set(sets):
    r = cys8.classify(make_seq(CLASS_II), sets)
    assert r["pattern_match"] and "II" in r["spacing_class"].split(",")
    assert "seidl_seiboth_2011_table3" in r["sets"].split(",")


def test_class_i_sequence_matches_class_i_set(sets):
    r = cys8.classify(make_seq(CLASS_I), sets)
    assert r["pattern_match"] and "I" in r["spacing_class"].split(",")


def test_gap_out_of_every_range_does_not_match(sets):
    assert not cys8.classify(make_seq([200, 0, 11, 16, 8, 0, 10]), sets)["pattern_match"]


def test_seven_cysteines_do_not_match(sets):
    s = make_seq(CLASS_II)
    last = s.rindex("C")
    assert not cys8.classify(s[:last] + "A" + s[last + 1 :], sets)["pattern_match"]


def test_n_cys_and_length_are_reported(sets):
    s = make_seq(CLASS_II)
    r = cys8.classify(s, sets)
    assert r["n_cys"] == 8 and r["length"] == len(s)


def test_sets_missing_a_core_gap_are_skipped_and_listed(sets):
    assert "jensen_2010_c5_c6_statement" in sets.skipped
    assert all(s.name != "jensen_2010_c5_c6_statement" for s in sets.sets)


def rows(proteins, sp_calls, sets):
    return {r["id"]: r for r in cys8.cys8_rows(proteins, sets, sp_calls)}


def test_hit_needs_pattern_and_r0_called(sets):
    seq = make_seq(CLASS_II)
    p = [prot("a", seq), prot("b", seq), prot("c", seq), prot("d", "M" + "A" * 100)]
    r = rows(p, {"a": "called", "b": "not_called", "d": "called"}, sets)
    assert (r["a"]["state"], r["a"]["hit"], r["a"]["pattern_match"]) == ("ok", "1", "1")
    assert (r["b"]["state"], r["b"]["hit"], r["b"]["pattern_match"]) == ("ok", "0", "1")
    assert (r["c"]["state"], r["c"]["hit"], r["c"]["pattern_match"]) == ("ok", "", "1")
    assert (r["d"]["state"], r["d"]["hit"], r["d"]["pattern_match"]) == ("ok", "0", "0")


def test_invalid_protein_gets_the_standard_invalid_row(sets):
    r = rows([prot("x", "", state="na_invalid")], {}, sets)
    assert r["x"]["state"] == "na_invalid"


def test_params_change_with_spacing_file_and_r0_identity(tmp_path, sets):
    a = cys8.module_params(
        sets, {"module": "step1_rule@R0", "params_hash": "1", "artefact_hash": "2"}
    )
    b = cys8.module_params(
        sets, {"module": "step1_rule@R0", "params_hash": "9", "artefact_hash": "2"}
    )
    assert a != b
    other = tmp_path / "s.yaml"
    other.write_text(SPACING.read_text() + "\n# changed\n")
    c = cys8.module_params(
        cys8.load_sets(other), {"module": "step1_rule@R0", "params_hash": "1", "artefact_hash": "2"}
    )
    assert a["spacing_sha256"] != c["spacing_sha256"]
