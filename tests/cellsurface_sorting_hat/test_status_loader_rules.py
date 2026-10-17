"""Task 3.1 (D7, D16, D17) and the strict cluster rule of D3: the status loader checks what it reads.

A status is never higher than the measure allows. Only ``curated`` negatives count for the status
numbers. A measure-less entry is ``unvalidated``. Tests that fail on ``main`` carry a strict xfail.
"""

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
CALL_STATUS_DIR = ROOT / "docs" / "reports" / "data" / "sorting_hat" / "call_status"
XFAIL_31 = pytest.mark.xfail(strict=True, reason="Task 3.1: status loader re-validates")


def rate(value, half):
    return {"value": value, "lo": value - half, "hi": min(1.0, value + half)}


def measure(
    n_pos=25, n_neg=25, cp=22, cn=22, half=0.05, notes="leakage: none", basis="curated", **over
):
    """A measure that derives ``estimated``: 25/25 proteins, 22/22 clusters, narrow intervals."""
    m = {
        "calibration_set": "toy",
        "truth_source": "toy.tsv",
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_clusters_pos": cp,
        "n_clusters_neg": cn,
        "sensitivity": rate(0.9, half),
        "specificity": rate(0.9, half),
        "notes": notes,
    }
    if cp is None:
        del m["n_clusters_pos"]
    if cn is None:
        del m["n_clusters_neg"]
    if basis is not None:
        m["negative_basis"] = {"all": {"basis": basis, "n": n_neg}}
    m.update(over)
    return m


def write_status(tmp_path, entries, name="status.json"):
    path = tmp_path / name
    path.write_text(
        json.dumps(
            {
                "module": "step1_rule@R0",
                "version": "1",
                "params_hash": "p",
                "artefact_hash": "a",
                "entries": entries,
            }
        )
    )
    return path


def entry(status, m=None, taxa=(40,), **extra):
    e = {"taxa": list(taxa), "status": status, "source": "toy"}
    if m is not None:
        e["measure"] = m
    e.update(extra)
    return e


# ---- status_from_measure -----------------------------------------------------------------------


@XFAIL_31
def test_all_curated_negatives_give_estimated():
    from cellsurface_sorting_hat.status import status_from_measure

    assert status_from_measure(measure()) == "estimated"


@XFAIL_31
def test_a_missing_negative_basis_gives_smoke():
    from cellsurface_sorting_hat.status import status_from_measure

    assert status_from_measure(measure(basis=None)) == "smoke"


@XFAIL_31
def test_assumed_negatives_do_not_count():
    from cellsurface_sorting_hat.status import status_from_measure

    assert status_from_measure(measure(basis="assumed")) == "smoke"


def _strata(n_int=13, n_sec=12, n_pmtm=40):
    return {
        "N-int": {"basis": "curated", "n": n_int},
        "N-sec": {"basis": "curated", "n": n_sec},
        "PM-TM": {"basis": "assumed", "n": n_pmtm},
    }


@XFAIL_31
def test_a_pm_tm_stratum_does_not_move_the_status_and_has_its_own_row():
    from cellsurface_sorting_hat.status import negative_rows, status_from_measure

    with_pmtm = measure(basis=None, n_neg=25, negative_basis=_strata())
    without = measure(
        basis=None,
        n_neg=25,
        negative_basis={k: v for k, v in _strata().items() if k != "PM-TM"},
    )
    assert status_from_measure(with_pmtm) == status_from_measure(without) == "estimated"
    rows = {r["stratum"]: r for r in negative_rows(with_pmtm)}
    assert rows["PM-TM"] == {"stratum": "PM-TM", "basis": "assumed", "n": 40, "counted": False}
    assert rows["N-int"]["counted"] is True and rows["N-sec"]["counted"] is True


@XFAIL_31
def test_a_headline_n_neg_that_counts_pm_tm_is_refused():
    from cellsurface_sorting_hat.status import status_from_measure

    bad = measure(basis=None, n_neg=65, negative_basis=_strata())
    with pytest.raises(ValueError, match="n_neg"):
        status_from_measure(bad)


@XFAIL_31
def test_a_basis_other_than_curated_or_assumed_is_refused():
    from cellsurface_sorting_hat.status import status_from_measure

    with pytest.raises(ValueError, match="curated or assumed"):
        status_from_measure(measure(basis="plausible"))


@pytest.mark.parametrize("missing", ["n_clusters_pos", "n_clusters_neg"])
@XFAIL_31
def test_a_measure_without_cluster_counts_cannot_be_estimated(missing):
    from cellsurface_sorting_hat.status import status_from_measure

    kwargs = {"cp": None} if missing == "n_clusters_pos" else {"cn": None}
    assert status_from_measure(measure(**kwargs)) == "smoke"


@XFAIL_31
def test_a_measure_with_too_few_clusters_is_smoke():
    from cellsurface_sorting_hat.status import status_from_measure

    assert status_from_measure(measure(cp=19)) == "smoke"
    assert status_from_measure(measure(cn=19)) == "smoke"
    assert status_from_measure(measure(cp=20, cn=20)) == "estimated"


# ---- load_status_source ------------------------------------------------------------------------


@XFAIL_31
def test_taxon_1_with_estimated_is_refused_on_load(tmp_path):
    from cellsurface_sorting_hat.status import load_status_source

    path = write_status(tmp_path, [entry("estimated", measure(), taxa=(1,))])
    with pytest.raises(ValueError, match="taxon 1"):
        load_status_source(path)


@XFAIL_31
def test_the_minimum_taxon_is_two(tmp_path):
    from cellsurface_sorting_hat.status import MIN_TAXON, load_status_source

    assert MIN_TAXON == 2
    load_status_source(write_status(tmp_path, [entry("estimated", measure(), taxa=(2,))]))
    with pytest.raises(ValueError, match="taxon 0"):
        load_status_source(write_status(tmp_path, [entry("smoke", measure(), taxa=(0,))], "b.json"))


@XFAIL_31
def test_a_status_stronger_than_its_measure_is_refused(tmp_path):
    from cellsurface_sorting_hat.status import load_status_source

    path = write_status(tmp_path, [entry("estimated", measure(n_pos=5))])  # derives smoke
    with pytest.raises(ValueError, match="stronger than the measure allows"):
        load_status_source(path)


@XFAIL_31
def test_a_status_lower_than_its_measure_loads(tmp_path):
    from cellsurface_sorting_hat.status import load_status_source

    record = load_status_source(write_status(tmp_path, [entry("smoke", measure())]))
    assert record.entries[0].status == "smoke"
    record = load_status_source(write_status(tmp_path, [entry("unvalidated", measure())], "b.json"))
    assert record.entries[0].status == "unvalidated"


@XFAIL_31
def test_an_estimated_entry_whose_negatives_are_not_curated_is_refused(tmp_path):
    from cellsurface_sorting_hat.status import load_status_source

    path = write_status(tmp_path, [entry("estimated", measure(basis=None))])
    with pytest.raises(ValueError, match="stronger than the measure allows"):
        load_status_source(path)


@XFAIL_31
@pytest.mark.parametrize("leak", ["partial", "tuned_on_truth", "in_reference", "unknown"])
def test_leakage_above_none_allows_at_most_smoke(tmp_path, leak):
    from cellsurface_sorting_hat.status import load_status_source

    bad = write_status(tmp_path, [entry("estimated", measure(notes=f"x; leakage: {leak}"))])
    with pytest.raises(ValueError, match="leakage"):
        load_status_source(bad)
    good = write_status(tmp_path, [entry("smoke", measure(notes=f"x; leakage: {leak}"))], "b.json")
    assert load_status_source(good).entries[0].status == "smoke"


@XFAIL_31
def test_an_estimated_entry_without_a_leakage_value_is_refused(tmp_path):
    from cellsurface_sorting_hat.status import load_status_source

    path = write_status(tmp_path, [entry("estimated", measure(notes="no leakage word"))])
    with pytest.raises(ValueError, match="leakage"):
        load_status_source(path)


@XFAIL_31
def test_a_leakage_exception_lifts_the_leakage_cap(tmp_path):
    """The key is the owner's exception (D13). Without it the same entry is refused."""
    from cellsurface_sorting_hat.status import load_status_source

    exception = {"reason": "leakage unmeasured", "signer": "owner", "date": "2026-10-06"}
    notes = "x; leakage: unknown (not measured)"
    with_key = measure(notes=notes, leakage_exception=exception)
    loaded = load_status_source(write_status(tmp_path, [entry("estimated", with_key)]))
    assert loaded.entries[0].status == "estimated"
    path = write_status(tmp_path, [entry("estimated", measure(notes=notes))], "b.json")
    with pytest.raises(ValueError, match="leakage"):
        load_status_source(path)


@XFAIL_31
def test_a_leakage_exception_needs_reason_signer_and_date(tmp_path):
    from cellsurface_sorting_hat.status import load_status_source

    for bad in ({"reason": "x", "signer": "o"}, {"reason": "", "signer": "o", "date": "d"}, "yes"):
        m = measure(notes="leakage: unknown", leakage_exception=bad)
        with pytest.raises(ValueError, match="leakage_exception: needs"):
            load_status_source(write_status(tmp_path, [entry("smoke", m)]))


@XFAIL_31
def test_a_measure_less_entry_above_unvalidated_is_refused(tmp_path):
    from cellsurface_sorting_hat.status import load_status_source

    for status in ("smoke", "estimated"):
        with pytest.raises(ValueError, match="needs a measure"):
            load_status_source(write_status(tmp_path, [entry(status)]))


def test_a_measure_less_entry_with_unvalidated_loads(tmp_path):
    from cellsurface_sorting_hat.status import load_status_source

    record = load_status_source(write_status(tmp_path, [entry("unvalidated")]))
    assert record.entries[0].status == "unvalidated"


@XFAIL_31
def test_write_status_source_refuses_a_measure_less_smoke(tmp_path):
    from cellsurface_sorting_hat.calibration.measure import write_status_source
    from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module

    write_module(tmp_path, ModuleSpec("step1_rule@R0", "1"), [], [{"id": "A", "state": "ok"}])
    with pytest.raises(ValueError, match="needs a measure"):
        write_status_source(tmp_path, "step1_rule@R0", [entry("smoke")])


# ---- the committed files still load -------------------------------------------------------------


def _committed_files():
    return sorted(CALL_STATUS_DIR.rglob("*.json"))


def test_there_are_twenty_committed_status_files():
    assert len(_committed_files()) == 20


@pytest.mark.parametrize("path", _committed_files(), ids=lambda p: p.name)
def test_every_committed_status_file_loads(path):
    """All carry a measure on every entry and every entry is smoke (FD-A14)."""
    from cellsurface_sorting_hat.call_status import build_source, call_file_name
    from cellsurface_sorting_hat.engine import load_config
    from cellsurface_sorting_hat.status import load_status_source

    data = json.loads(path.read_text(encoding="utf-8-sig"))
    if "call" in data:
        source = build_source(
            data, load_config(), call_file_name(data["call"], data["variant"]), path
        )
        entries = source.entries
    else:
        entries = load_status_source(path).entries
    assert entries
    assert all(e.measure is not None and e.status == "smoke" for e in entries)
