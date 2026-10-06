import json

import pytest

from cellsurface_sorting_hat.status import (
    ModuleIdentity,
    StatusEntry,
    StatusRecord,
    load_status_source,
    resolve_entry,
    resolve_status,
    weakest,
)
from cellsurface_sorting_hat.taxonomy import Lineage, TaxonError

IDENT = ModuleIdentity("step1_rule@R0", "1", "p", "a")


@pytest.fixture
def lineage(nodes_dmp):
    return Lineage.from_nodes_dmp(nodes_dmp)


def test_ancestors_and_depth(lineage):
    assert lineage.ancestors(41) == [41, 31, 20, 10, 1]
    assert lineage.depth(1) == 0
    assert lineage.depth(41) == 4


def test_unknown_taxon_raises(lineage):
    with pytest.raises(TaxonError):
        lineage.ancestors(999)


def test_descendant_or_self(lineage):
    assert lineage.is_descendant_or_self(40, 30)
    assert lineage.is_descendant_or_self(40, 40)
    assert not lineage.is_descendant_or_self(30, 40)


def test_weakest_status():
    assert weakest(["estimated", "smoke"]) == "smoke"
    assert weakest(["estimated", "unvalidated", "smoke"]) == "unvalidated"
    assert weakest([]) == "unvalidated"
    with pytest.raises(ValueError):
        weakest(["good"])


def _record(entries, ident=IDENT):
    return StatusRecord(ident, tuple(StatusEntry(taxa, status) for taxa, status in entries))


def test_status_applies_to_a_tested_taxon(lineage):
    rec = _record([((40,), "estimated")])
    assert resolve_status(rec, IDENT, 40, lineage) == ("estimated", "taxon:40")


def test_status_does_not_pass_to_a_sibling_clade_with_the_same_broad_label(lineage):
    # tested on A. fumigatus (Eurotiales); Coccidioides is Onygenales, also in Eurotiomycetes
    rec = _record([((40,), "estimated")])
    assert resolve_status(rec, IDENT, 41, lineage) == ("unvalidated", "taxon not tested")


def test_status_passes_to_descendants_of_a_tested_taxon(lineage):
    rec = _record([((30,), "smoke")])
    assert resolve_status(rec, IDENT, 42, lineage) == ("smoke", "taxon:30")


def test_most_specific_tested_taxon_wins(lineage):
    rec = _record([((20,), "smoke"), ((40,), "estimated")])
    assert resolve_status(rec, IDENT, 40, lineage)[0] == "estimated"
    assert resolve_status(rec, IDENT, 42, lineage) == ("smoke", "taxon:20")


@pytest.mark.parametrize("field", ["name", "version", "params_hash", "artefact_hash"])
def test_stale_status_source_is_refused(lineage, field):
    rec = _record([((40,), "estimated")])
    running = ModuleIdentity(**{**IDENT.__dict__, field: "changed"})
    status, basis = resolve_status(rec, running, 40, lineage)
    assert status == "unvalidated"
    assert field in basis


def test_no_status_source(lineage):
    assert resolve_status(None, IDENT, 40, lineage) == ("unvalidated", "no status_source")


def test_load_status_source(tmp_path):
    path = tmp_path / "s.json"
    path.write_text(
        json.dumps(
            {
                "module": "m",
                "version": "2",
                "params_hash": "p",
                "artefact_hash": "a",
                "entries": [{"taxa": [40, 42], "status": "estimated", "source": "x"}],
            }
        )
    )
    rec = load_status_source(path)
    assert rec.identity == ModuleIdentity("m", "2", "p", "a")
    assert rec.entries == (StatusEntry((40, 42), "estimated", "x", None),)


def test_load_status_source_rejects_unknown_status(tmp_path):
    path = tmp_path / "s.json"
    path.write_text(
        json.dumps(
            {
                "module": "m",
                "version": "2",
                "params_hash": "p",
                "artefact_hash": "a",
                "entries": [{"taxa": [40], "status": "validated"}],
            }
        )
    )
    with pytest.raises(ValueError):
        load_status_source(path)


MEASURE = {
    "calibration_set": "S1:all",
    "truth_source": "phasec/metrics.json",
    "n_pos": 232,
    "n_neg": 4244,
    "sensitivity": {"value": 0.603, "lo": 0.55, "hi": 0.65},
    "specificity": {"value": 0.963, "lo": 0.95, "hi": 0.97},
    "notes": "R0",
}


def _status_file(tmp_path, measure, extra_entry=None):
    entry = {"taxa": [40], "status": "estimated", "source": "phasec"}
    if measure is not None:
        entry["measure"] = measure
    entry.update(extra_entry or {})
    path = tmp_path / "s.json"
    path.write_text(
        json.dumps(
            {
                "module": "m",
                "version": "1",
                "params_hash": "p",
                "artefact_hash": "a",
                "entries": [entry],
            }
        )
    )
    return path


def test_measure_with_sensitivity_and_specificity_is_loaded_and_returned(tmp_path, lineage):
    rec = load_status_source(_status_file(tmp_path, MEASURE))
    ident = ModuleIdentity("m", "1", "p", "a")
    entry, tested, reason = resolve_entry(rec, ident, 40, lineage)
    assert (tested, reason) == (40, "")
    assert entry.measure["sensitivity"]["value"] == 0.603 and entry.measure["n_pos"] == 232
    assert resolve_entry(rec, ident, 41, lineage)[0] is None  # taxon 41 is not covered


def test_an_entry_without_a_measure_is_allowed(tmp_path):
    assert load_status_source(_status_file(tmp_path, None)).entries[0].measure is None


@pytest.mark.parametrize(
    "bad",
    [
        {**MEASURE, "sensitivity": {"value": 1.2, "lo": 1.0, "hi": 1.3}},
        {**MEASURE, "specificity": {"value": 0.5, "lo": 0.6, "hi": 0.7}},
        {**MEASURE, "sensitivity": {"value": 0.5}},
        {**MEASURE, "n_pos": -1},
        {**MEASURE, "surprise": 1},
    ],
)
def test_a_bad_measure_is_refused(tmp_path, bad):
    with pytest.raises(ValueError):
        load_status_source(_status_file(tmp_path, bad))


def test_an_unknown_entry_key_is_refused(tmp_path):
    with pytest.raises(ValueError, match="unknown key"):
        load_status_source(_status_file(tmp_path, None, {"surprise": 1}))


def test_a_non_integer_parent_in_nodes_dmp_names_the_line(tmp_path):
    path = tmp_path / "nodes.dmp"
    path.write_text("1\t|\t1\t|\tno rank\t|\n2\t|\tx\t|\tno rank\t|\n")
    with pytest.raises(TaxonError, match="nodes.dmp:2"):
        Lineage.from_nodes_dmp(path)


def _source(tmp_path, **over):
    data = {
        "module": "m",
        "version": "2",
        "params_hash": "p",
        "artefact_hash": "a",
        "entries": [{"taxa": [40], "status": "estimated"}],
    }
    data.update(over)
    path = tmp_path / "s.json"
    path.write_text(json.dumps(data))
    return path


@pytest.mark.parametrize(
    "over",
    [
        {"entries": [{"taxa": [True], "status": "estimated"}]},
        {"entries": [{"taxa": [40], "status": "estimated", "measure": {"n_pos": True}}]},
        {"entries": [{"taxa": [40], "status": "estimated", "measure": {"n_neg": True}}]},
        {
            "entries": [
                {
                    "taxa": [40],
                    "status": "estimated",
                    "measure": {"sensitivity": {"value": True, "lo": 0, "hi": 1}},
                }
            ]
        },
        {
            "entries": [
                {
                    "taxa": [40],
                    "status": "estimated",
                    "measure": {"specificity": {"value": 0.5, "lo": False, "hi": 1}},
                }
            ]
        },
        {"module": ""},
        {"version": ""},
        {"params_hash": ""},
        {"artefact_hash": ""},
        {
            "entries": [
                {"taxa": [40, 42], "status": "estimated"},
                {"taxa": [42], "status": "smoke"},
            ]
        },
    ],
)
def test_load_status_source_refuses(tmp_path, over):
    with pytest.raises(ValueError):
        load_status_source(_source(tmp_path, **over))


def test_per_stratum_rates_and_cluster_counts_are_accepted_and_checked(tmp_path):
    good = {
        **MEASURE,
        "n_clusters_pos": 40,
        "n_clusters_neg": 300,
        "strata": {"N-sec": {"value": 0.9, "lo": 0.8, "hi": 0.95}},
    }
    assert (
        load_status_source(_status_file(tmp_path, good)).entries[0].measure["n_clusters_pos"] == 40
    )
    bad = {**MEASURE, "strata": {"N-sec": {"value": 2, "lo": 0, "hi": 3}}}
    with pytest.raises(ValueError):
        load_status_source(_status_file(tmp_path, bad))
    with pytest.raises(ValueError):
        load_status_source(_status_file(tmp_path, {**MEASURE, "n_clusters_pos": -1}))
