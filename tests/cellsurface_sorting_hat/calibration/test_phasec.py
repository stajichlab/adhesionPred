import json
from pathlib import Path

import pytest

from cellsurface_sorting_hat.calibration.measure import write_status_source
from cellsurface_sorting_hat.calibration.phasec import (
    SETS,
    entries_from_phasec,
    read_names,
    species_taxid,
)
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import load_status_source
from cellsurface_sorting_hat.taxonomy import TaxonError


def cell(r, f):
    return {
        "recall": {"value": r[0], "lo": r[1], "hi": r[2]},
        "fpr": {"value": f[0], "lo": f[1], "hi": f[2]} if f else {"value": None},
    }


def make_set(label, pos, neg, c, strata=None):
    metrics = {"all": {"V-go": {"R0": c}}}
    for name, f in (strata or {}).items():
        metrics[name] = {"V-go": {"R0": {"recall": {"value": None}, "fpr": f}}}
    return {
        "label": label,
        "n_direct_positives": pos,
        "truth": {"direct": {"n": {"all": {"pos": pos, "neg": neg}}, "metrics": metrics}},
    }


def _metrics(tmp_path, scer_label="estimate"):
    data = {
        "test_sets": {
            "S1:Scer_SGD": make_set(
                scer_label,
                79,
                2000,
                cell((0.848, 0.77, 0.92), (0.03, 0.02, 0.04)),
                {
                    "N-sec": {"value": 0.084, "lo": 0.06, "hi": 0.11},
                    "N-int": {"value": 0.002, "lo": 0.001, "hi": 0.004},
                },
            ),
            "S1:Calb_CGD": make_set(
                "smoke test", 153, 2244, cell((0.477, 0.36, 0.59), (0.04, 0.03, 0.05))
            ),
            "S3-Eurotiomycetes:Afum_ASPFU": make_set(
                "smoke test", 19, 44, cell((0.947, 0.78, 1.0), (0.02, 0.0, 0.08))
            ),
        }
    }
    path = tmp_path / "metrics.json"
    path.write_text(json.dumps(data))
    return path


SET_TAXA = {"S1:Scer_SGD": [4932], "S1:Calb_CGD": [5476], "S3-Eurotiomycetes:Afum_ASPFU": [746128]}


def test_each_species_gets_its_own_entry_with_its_own_numbers(tmp_path):
    entries = entries_from_phasec(_metrics(tmp_path), SET_TAXA)
    by = {e["measure"]["calibration_set"]: e for e in entries}
    assert [e["taxa"] for e in entries] == [[4932], [5476], [746128]]
    assert by["S1:Scer_SGD"]["measure"]["sensitivity"] == {"value": 0.848, "lo": 0.77, "hi": 0.92}
    assert by["S1:Calb_CGD"]["measure"]["sensitivity"]["value"] == 0.477  # not a pooled value
    sp = by["S1:Scer_SGD"]["measure"]["specificity"]
    assert (
        sp["value"] == pytest.approx(0.97)
        and sp["lo"] == pytest.approx(0.96)
        and sp["hi"] == pytest.approx(0.98)
    )
    assert (by["S1:Scer_SGD"]["measure"]["n_pos"], by["S1:Scer_SGD"]["measure"]["n_neg"]) == (
        79,
        2000,
    )


def test_per_stratum_specificity_is_kept(tmp_path):
    scer = entries_from_phasec(_metrics(tmp_path), SET_TAXA)[0]["measure"]
    assert scer["strata"]["N-sec"]["value"] == pytest.approx(0.916)
    assert scer["strata"]["N-int"]["lo"] == pytest.approx(0.996)
    assert "strata" not in entries_from_phasec(_metrics(tmp_path), SET_TAXA)[1]["measure"]


def test_the_status_is_the_weaker_of_the_phase_c_label_and_the_rule(tmp_path):
    entries = {
        e["measure"]["calibration_set"]: e
        for e in entries_from_phasec(_metrics(tmp_path), SET_TAXA)
    }
    assert entries["S1:Scer_SGD"]["status"] == "estimated"  # label estimate and the rule agrees
    assert (
        entries["S3-Eurotiomycetes:Afum_ASPFU"]["status"] == "smoke"
    )  # label smoke test, 19 positives
    data = json.loads(_metrics(tmp_path).read_text())
    data["test_sets"]["S1:Calb_CGD"] = make_set(
        "estimate", 400, 4000, cell((0.6, 0.3, 0.9), (0.05, 0.04, 0.06))
    )
    path = tmp_path / "wide.json"
    path.write_text(json.dumps(data))
    wide = {e["measure"]["calibration_set"]: e for e in entries_from_phasec(path, SET_TAXA)}[
        "S1:Calb_CGD"
    ]
    assert wide["status"] == "smoke"  # labelled estimate, but the sensitivity interval is wide


def test_a_set_with_no_specificity_cell_is_at_most_smoke_even_if_labelled_estimate(tmp_path):
    data = json.loads(_metrics(tmp_path).read_text())
    data["test_sets"]["S1:Calb_CGD"] = make_set(
        "estimate", 400, 4000, cell((0.5, 0.48, 0.52), None)
    )
    path = tmp_path / "m2.json"
    path.write_text(json.dumps(data))
    e = {x["measure"]["calibration_set"]: x for x in entries_from_phasec(path, SET_TAXA)}[
        "S1:Calb_CGD"
    ]
    assert "specificity" not in e["measure"] and e["status"] == "smoke"


def test_a_tight_estimate_with_both_rates_is_estimated(tmp_path):
    data = json.loads(_metrics(tmp_path).read_text())
    data["test_sets"]["S1:Calb_CGD"] = make_set(
        "estimate", 400, 4000, cell((0.5, 0.48, 0.52), (0.05, 0.04, 0.06))
    )
    path = tmp_path / "m3.json"
    path.write_text(json.dumps(data))
    e = {x["measure"]["calibration_set"]: x for x in entries_from_phasec(path, SET_TAXA)}[
        "S1:Calb_CGD"
    ]
    assert e["status"] == "estimated"


def test_entries_make_a_valid_status_source(tmp_path):
    write_module(tmp_path, ModuleSpec("step1_rule@R0", "1"), [], [{"id": "A", "state": "ok"}])
    path = write_status_source(
        tmp_path, "step1_rule@R0", entries_from_phasec(_metrics(tmp_path), SET_TAXA)
    )
    assert len(load_status_source(path).entries) == 3


def test_default_sets_are_single_species_and_share_no_source():
    assert all(len(group) == 1 for group in SETS.values())
    sources = [s for group in SETS.values() for s in group]
    assert len(sources) == len(set(sources)) == 6


def test_species_taxid_reads_names_dmp_and_refuses_missing_or_ambiguous_names(tmp_path):
    names = tmp_path / "names.dmp"
    names.write_text(
        "4932\t|\tSaccharomyces cerevisiae\t|\t\t|\tscientific name\t|\n"
        "4932\t|\tbaker's yeast\t|\t\t|\tcommon name\t|\n"
        "111\t|\tTwin\t|\t\t|\tscientific name\t|\n"
        "222\t|\tTwin\t|\t\t|\tscientific name\t|\n"
    )
    table = read_names(names)
    assert species_taxid(table, "Saccharomyces cerevisiae") == 4932
    for bad in ("Twin", "Nothing here"):
        with pytest.raises(TaxonError):
            species_taxid(table, bad)


def test_a_zero_width_interval_at_a_rate_of_one_is_widened_with_wilson(tmp_path):
    data = json.loads(_metrics(tmp_path).read_text())
    data["test_sets"]["S1:Calb_CGD"] = make_set(
        "estimate", 9, 28, cell((1.0, 1.0, 1.0), (0.0, 0.0, 0.0))
    )
    path = tmp_path / "edge.json"
    path.write_text(json.dumps(data))
    m = {e["measure"]["calibration_set"]: e for e in entries_from_phasec(path, SET_TAXA)}[
        "S1:Calb_CGD"
    ]["measure"]
    assert (
        m["sensitivity"]["lo"] == pytest.approx(0.701, abs=1e-3) and m["sensitivity"]["hi"] == 1.0
    )  # 9 of 9
    assert m["specificity"]["lo"] == pytest.approx(0.879, abs=1e-3)  # 28 of 28


def _with_calb(tmp_path, name, label, pos, neg, c):
    data = json.loads(_metrics(tmp_path).read_text())
    data["test_sets"]["S1:Calb_CGD"] = make_set(label, pos, neg, c)
    path = tmp_path / name
    path.write_text(json.dumps(data))
    return path


def _calb(path):
    return {e["measure"]["calibration_set"]: e for e in entries_from_phasec(path, SET_TAXA)}[
        "S1:Calb_CGD"
    ]


def test_the_phase_c_label_caps_a_tight_measure(tmp_path):
    path = _with_calb(
        tmp_path, "l.json", "smoke test", 400, 4000, cell((0.5, 0.48, 0.52), (0.05, 0.04, 0.06))
    )
    assert (
        _calb(path)["status"] == "smoke"
    )  # an implementation that ignores the label gives estimated


def test_the_rule_boundaries_at_20_positives_and_half_width_0_10(tmp_path):
    tight = cell((0.5, 0.4, 0.6), (0.05, 0.0, 0.1))
    assert _calb(_with_calb(tmp_path, "a.json", "estimate", 20, 20, tight))["status"] == "estimated"
    assert _calb(_with_calb(tmp_path, "b.json", "estimate", 19, 20, tight))["status"] == "smoke"
    assert _calb(_with_calb(tmp_path, "c.json", "estimate", 20, 19, tight))["status"] == "smoke"
    wide = cell((0.5, 0.39, 0.61), (0.05, 0.0, 0.1))
    assert _calb(_with_calb(tmp_path, "d.json", "estimate", 20, 20, wide))["status"] == "smoke"


def test_the_notes_do_not_claim_the_cluster_floor_was_met(tmp_path):
    m = _calb(
        _with_calb(
            tmp_path, "n.json", "estimate", 400, 4000, cell((0.5, 0.48, 0.52), (0.05, 0.04, 0.06))
        )
    )["measure"]
    assert "n_clusters_pos" not in m and "cluster floor was not checked" in m["notes"]


def test_a_species_not_in_the_file_is_refused_with_path_and_key(tmp_path):
    path = _metrics(tmp_path)
    with pytest.raises(ValueError, match=r"metrics\.json.*S1:Missing"):
        entries_from_phasec(path, {**SET_TAXA, "S1:Missing": [1]})


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -0.1, 1.5, "0.5", True])
def test_a_non_finite_or_out_of_range_rate_is_refused(tmp_path, bad):
    path = _with_calb(
        tmp_path, "bad.json", "estimate", 50, 50, cell((bad, 0.4, 0.6), (0.05, 0.0, 0.1))
    )
    with pytest.raises(ValueError, match=r"bad\.json.*recall"):
        entries_from_phasec(path, SET_TAXA)


def test_a_value_without_an_interval_and_a_lo_above_the_value_are_refused(tmp_path):
    c = cell((0.5, 0.4, 0.6), (0.05, 0.0, 0.1))
    del c["recall"]["lo"]
    with pytest.raises(ValueError, match="recall"):
        entries_from_phasec(_with_calb(tmp_path, "x.json", "estimate", 50, 50, c), SET_TAXA)
    c = cell((0.5, 0.55, 0.6), (0.05, 0.0, 0.1))
    with pytest.raises(ValueError, match="lo <= value <= hi"):
        entries_from_phasec(_with_calb(tmp_path, "y.json", "estimate", 50, 50, c), SET_TAXA)


def test_an_unknown_label_and_a_bad_count_are_refused(tmp_path):
    c = cell((0.5, 0.4, 0.6), (0.05, 0.0, 0.1))
    with pytest.raises(ValueError, match="label"):
        entries_from_phasec(_with_calb(tmp_path, "u.json", "validated", 50, 50, c), SET_TAXA)
    with pytest.raises(ValueError, match=r"n\.all\.pos"):
        entries_from_phasec(_with_calb(tmp_path, "v.json", "estimate", -1, 50, c), SET_TAXA)


def test_an_unreadable_file_names_the_path(tmp_path):
    path = tmp_path / "broken.json"
    path.write_text("{not json")
    with pytest.raises(ValueError, match=r"broken\.json"):
        entries_from_phasec(path, SET_TAXA)


REAL = Path(
    "/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare/phasec/metrics.json"
)


@pytest.mark.skipif(not REAL.exists(), reason="the Phase C output is not on this machine")
def test_the_real_phase_c_file_loads_one_entry_per_species():
    taxa = {k: [i + 1] for i, k in enumerate(SETS)}
    entries = entries_from_phasec(REAL, taxa)
    by = {e["measure"]["calibration_set"]: e for e in entries}
    assert len(entries) == 6
    assert by["S3-Eurotiomycetes:Anid_EMENI"]["status"] == "estimated"
    assert all(e["status"] == "smoke" for k, e in by.items() if "Anid" not in k)
    assert by["S3-Eurotiomycetes:Afum_ASPFU"]["measure"]["n_pos"] == 19
    assert by["S1:Scer_SGD"]["measure"]["sensitivity"]["value"] == pytest.approx(0.848, abs=1e-3)
    assert by["S1:Calb_CGD"]["measure"]["sensitivity"]["value"] == pytest.approx(0.477, abs=1e-3)
