import json
from pathlib import Path

import pytest

from cellsurface_sorting_hat.calibration.intervals import wilson_p
from cellsurface_sorting_hat.calibration.measure import write_status_source
from cellsurface_sorting_hat.calibration.phasec import (
    LEAKAGE_NOTE,
    SETS,
    check_signalp,
    count_clusters,
    read_names,
    species_taxid,
)
from cellsurface_sorting_hat.calibration.phasec import entries_from_phasec as _entries_from_phasec
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import load_status_source
from cellsurface_sorting_hat.taxonomy import TaxonError

MANY = 10**6  # so many clusters that the Wilson interval does not change a file interval


def counts_from(path, c_pos=MANY, c_neg=MANY):
    """Cluster counts that agree with the protein counts of the metrics file."""
    try:
        sets = json.loads(Path(path).read_text())["test_sets"]
    except ValueError:  # a test of an unreadable file: the real function reports it
        sets = {}
    out = {s: {"pos": (0, 0), "neg": (0, 0)} for g in SETS.values() for s in g}
    for key, group in SETS.items():
        if key in sets:
            n = sets[key]["truth"]["direct"]["n"]["all"]
            out[group[0]] = {"pos": (n["pos"], c_pos), "neg": (n["neg"], c_neg)}
    return out


def entries_from_phasec(path, taxa, counts=None, **kw):
    return _entries_from_phasec(path, taxa, counts or counts_from(path), **kw)


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


def _calb_with(path, counts):
    return {
        e["measure"]["calibration_set"]: e for e in _entries_from_phasec(path, SET_TAXA, counts)
    }["S1:Calb_CGD"]


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


def test_the_notes_record_clusters_call_module_and_the_leakage_limit(tmp_path):
    path = _with_calb(
        tmp_path, "n.json", "estimate", 400, 4000, cell((0.5, 0.48, 0.52), (0.05, 0.04, 0.06))
    )
    counts = counts_from(path, c_pos=300, c_neg=900)
    m = {
        e["measure"]["calibration_set"]: e
        for e in _entries_from_phasec(path, SET_TAXA, counts, extra_notes="signalp_mode=fast")
    }["S1:Calb_CGD"]["measure"]
    assert (m["n_clusters_pos"], m["n_clusters_neg"]) == (300, 900)
    assert "cluster floor was not checked" not in m["notes"]
    assert LEAKAGE_NOTE in m["notes"]
    assert m["notes"].endswith(
        "leakage: overlap between the Phase C positives and the SignalP 6 training data "
        "was not measured"
    )
    assert "call=signal_peptide_protein; module=step1_rule@R0" in m["notes"]
    assert "signalp_mode=fast" in m["notes"]


def test_the_interval_is_the_widest_of_the_file_and_the_wilson_interval_on_clusters(tmp_path):
    path = _with_calb(
        tmp_path, "w.json", "estimate", 400, 4000, cell((0.5, 0.48, 0.52), (0.05, 0.04, 0.06))
    )
    m = _calb_with(path, counts_from(path, c_pos=40, c_neg=50))["measure"]
    _, lo, hi = wilson_p(0.5, 40)
    assert m["sensitivity"]["lo"] == pytest.approx(lo) and m["sensitivity"]["hi"] == pytest.approx(
        hi
    )
    _, lo, hi = wilson_p(0.95, 50)  # specificity = 1 - fpr = 0.95
    assert m["specificity"]["lo"] == pytest.approx(lo) and m["specificity"]["hi"] == pytest.approx(
        hi
    )
    # a file interval wider than the Wilson interval stays
    path = _with_calb(
        tmp_path, "w2.json", "estimate", 400, 4000, cell((0.5, 0.2, 0.8), (0.05, 0.0, 0.2))
    )
    m = _calb_with(path, counts_from(path, c_pos=300, c_neg=300))["measure"]
    assert (m["sensitivity"]["lo"], m["sensitivity"]["hi"]) == (0.2, 0.8)


def test_fewer_than_20_clusters_of_a_class_is_smoke_even_with_a_narrow_interval(tmp_path):
    path = _with_calb(
        tmp_path, "f.json", "estimate", 400, 4000, cell((0.99, 0.98, 1.0), (0.01, 0.0, 0.02))
    )
    assert _calb_with(path, counts_from(path, c_pos=300, c_neg=300))["status"] == "estimated"
    assert _calb_with(path, counts_from(path, c_pos=19, c_neg=300))["status"] == "smoke"
    assert _calb_with(path, counts_from(path, c_pos=300, c_neg=19))["status"] == "smoke"


def test_protein_counts_that_differ_from_metrics_are_refused_with_both_counts(tmp_path):
    path = _metrics(tmp_path)
    counts = counts_from(path)
    counts["Calb_CGD"] = {"pos": (152, 100), "neg": (2244, 100)}
    with pytest.raises(ValueError, match=r"metrics\.json.*153.*2244.*152.*2244.*Calb_CGD"):
        _entries_from_phasec(path, SET_TAXA, counts, where_counts="eval.tsv.gz joined to cl.tsv.gz")
    with pytest.raises(ValueError, match="eval.tsv.gz joined to cl.tsv.gz"):
        _entries_from_phasec(path, SET_TAXA, counts, where_counts="eval.tsv.gz joined to cl.tsv.gz")


def _write_gz(path, text):
    import gzip

    with gzip.open(path, "wt") as fh:
        fh.write(text)


EVAL_HEAD = "seq_sha256\torigin\tclass\thomology_only\tsource_ids\n"


def test_count_clusters_joins_the_two_files_and_applies_the_filters(tmp_path):
    _write_gz(
        tmp_path / "cl.tsv.gz",
        "seq_sha256\tcluster_id\n"
        + "".join(f"{h}\t{c}\n" for h, c in [("a", "1"), ("b", "1"), ("c", "2"), ("d", "3")])
        + "e\t4\nf\t5\n",
    )
    _write_gz(
        tmp_path / "ev.tsv.gz",
        EVAL_HEAD
        + "a\tgo\tpos\tno\tAnid_EMENI\n"
        + "b\tgo\tpos\tno\tAnid_EMENI\n"  # same cluster as a
        + "c\tgo\tneg\tno\tAfum_ASPFU,Anid_EMENI\n"  # in two sources
        + "d\ttc\tpos\t\tAnid_EMENI\n"  # not origin go
        + "e\tgo\tpos\tyes\tAnid_EMENI\n"  # homology only
        + "f\tgo\texcluded\tno\tAnid_EMENI\n",  # not pos or neg
    )
    got = count_clusters(tmp_path / "cl.tsv.gz", tmp_path / "ev.tsv.gz")
    assert got["Anid_EMENI"] == {"pos": (2, 1), "neg": (1, 1)}
    assert got["Afum_ASPFU"] == {"pos": (0, 0), "neg": (1, 1)}
    assert got["Scer_SGD"] == {"pos": (0, 0), "neg": (0, 0)}


def test_count_clusters_refuses_a_protein_with_no_cluster_and_a_missing_column(tmp_path):
    _write_gz(tmp_path / "cl.tsv.gz", "seq_sha256\tcluster_id\na\t1\n")
    _write_gz(tmp_path / "ev.tsv.gz", EVAL_HEAD + "z\tgo\tpos\tno\tAnid_EMENI\n")
    with pytest.raises(ValueError, match=r"ev\.tsv\.gz.*z.*cl\.tsv\.gz"):
        count_clusters(tmp_path / "cl.tsv.gz", tmp_path / "ev.tsv.gz")
    _write_gz(tmp_path / "ev2.tsv.gz", "seq_sha256\torigin\n")
    with pytest.raises(ValueError, match=r"ev2\.tsv\.gz.*missing column"):
        count_clusters(tmp_path / "cl.tsv.gz", tmp_path / "ev2.tsv.gz")


def _r0_record(tmp_path, version="6.0h-gpu", mode="fast"):
    path = tmp_path / "r0.json"
    path.write_text(json.dumps({"tools": {"signalp": version}, "params": {"mode": mode}}))
    return path


def test_check_signalp_accepts_the_matching_module_and_mode(tmp_path):
    assert check_signalp(_r0_record(tmp_path), "signalp/6-gpu", "fast") == "6.0h-gpu"


@pytest.mark.parametrize(
    "version, module, mode",
    [
        ("6.0h-gpu", "signalp/6-gpu", "slow"),  # another mode
        ("6.0h-cpu", "signalp/6-gpu", "fast"),  # CPU build, GPU module
        ("6.0h-gpu", "signalp/6", "fast"),  # GPU build, CPU module
        ("5.1b", "signalp/6-gpu", "fast"),  # another major version
    ],
)
def test_check_signalp_refuses_a_mismatch_naming_both_values(tmp_path, version, module, mode):
    with pytest.raises(ValueError, match="r0.json") as err:
        check_signalp(_r0_record(tmp_path, version), module, mode)
    assert module in str(err.value) or mode in str(err.value)


def test_check_signalp_refuses_a_module_that_is_not_signalp(tmp_path):
    with pytest.raises(ValueError, match="tmhmm/2"):
        check_signalp(_r0_record(tmp_path), "tmhmm/2", "fast")


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
    counts = count_clusters(REAL.with_name("clusters.tsv.gz"), REAL.with_name("eval_table.tsv.gz"))
    entries = _entries_from_phasec(REAL, taxa, counts)
    by = {e["measure"]["calibration_set"]: e for e in entries}
    assert len(entries) == 6
    assert by["S3-Eurotiomycetes:Anid_EMENI"]["status"] == "estimated"
    assert all(e["status"] == "smoke" for k, e in by.items() if "Anid" not in k)
    assert by["S3-Eurotiomycetes:Afum_ASPFU"]["measure"]["n_pos"] == 19
    assert by["S1:Scer_SGD"]["measure"]["sensitivity"]["value"] == pytest.approx(0.848, abs=1e-3)
    assert by["S1:Calb_CGD"]["measure"]["sensitivity"]["value"] == pytest.approx(0.477, abs=1e-3)


REAL_FORM = "module=signalp/6-gpu;signalp6_version=SignalP 6.0 (fast, eukarya) 6.0h"


def test_check_signalp_accepts_the_two_line_form_that_the_sbatch_writes(tmp_path):
    rec = _r0_record(tmp_path, REAL_FORM)
    assert check_signalp(rec, "signalp/6-gpu", "fast") == REAL_FORM


@pytest.mark.parametrize(
    "version, module",
    [
        ("module=signalp/6-cpu;signalp6_version=6.0h", "signalp/6-gpu"),
        ("module=signalp/6-gpu;signalp6_version=6.0h", "signalp/6"),
        ("module=signalp/6-gpu;signalp6_version=6.0h", "signalp/6-gpus"),
        ("module=signalp/6-gpu; call=cell_wall_adhesion_candidate", "signalp/6"),
    ],
)
def test_check_signalp_compares_a_module_token_exactly(tmp_path, version, module):
    with pytest.raises(ValueError, match="r0.json"):
        check_signalp(_r0_record(tmp_path, version), module, "fast")


def test_check_signalp_takes_the_module_token_up_to_the_semicolon(tmp_path):
    version = "module=signalp/6-gpu; call=cell_wall_adhesion_candidate"
    assert check_signalp(_r0_record(tmp_path, version), "signalp/6-gpu", "fast") == version


@pytest.mark.parametrize(
    "version, module, ok",
    [
        ("6.0h-gpu", "signalp/6-gpu", True),
        ("60-gpu", "signalp/6-gpu", False),  # 6 is not 60
        ("6.0h-cpu", "signalp/6-gpu", False),
        ("6.0h", "signalp/6", True),
        ("6.0h-gpus", "signalp/6-gpu", False),  # gpu must be a whole token
        ("6.0h-gpu", "signalp/6", False),
        ("6.0h", "signalp/6-cpu", False),  # unknown tag
    ],
)
def test_check_signalp_bare_version_fallback_is_strict(tmp_path, version, module, ok):
    rec = _r0_record(tmp_path, version)
    if ok:
        assert check_signalp(rec, module, "fast") == version
    else:
        with pytest.raises(ValueError):
            check_signalp(rec, module, "fast")


def test_a_protein_in_two_clusters_is_refused(tmp_path):
    _write_gz(tmp_path / "cl.tsv.gz", "seq_sha256\tcluster_id\na\t1\na\t2\n")
    _write_gz(tmp_path / "ev.tsv.gz", EVAL_HEAD + "a\tgo\tpos\tno\tAnid_EMENI\n")
    with pytest.raises(ValueError, match=r"cl\.tsv\.gz.*two clusters"):
        count_clusters(tmp_path / "cl.tsv.gz", tmp_path / "ev.tsv.gz")


@pytest.mark.skipif(not REAL.exists(), reason="the Phase C output is not on this machine")
def test_the_real_phase_c_cluster_counts_and_statuses():
    counts = count_clusters(REAL.with_name("clusters.tsv.gz"), REAL.with_name("eval_table.tsv.gz"))
    got = {s: (c["pos"], c["neg"]) for s, c in counts.items()}
    assert got == {
        "Anid_EMENI": ((109, 100), (164, 151)),
        "Afum_ASPFU": ((19, 17), (45, 38)),
        "Scer_SGD": ((79, 58), (3785, 3156)),
        "Calb_CGD": ((153, 113), (459, 410)),
        "Cneo_H99_GOA": ((7, 6), (32, 31)),
        "Umay_MYCMD": ((9, 9), (28, 24)),
    }
    taxa = {k: [i + 1] for i, k in enumerate(SETS)}
    status = {
        e["measure"]["calibration_set"]: e["status"]
        for e in _entries_from_phasec(REAL, taxa, counts)
    }
    assert status["S3-Eurotiomycetes:Anid_EMENI"] == "estimated"
    assert [k for k, v in status.items() if v == "estimated"] == ["S3-Eurotiomycetes:Anid_EMENI"]
