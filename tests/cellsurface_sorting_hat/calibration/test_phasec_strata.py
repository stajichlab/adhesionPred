"""Phase C entries carry the basis of each negative stratum (D16) and the R0 leakage exception (D13).

Fixture: one species, 120 positives in 110 clusters (recall 0.90), curated negatives N-int 20 and
N-sec 20 (38 clusters, no false positive), and 10 PM-TM negatives (4 false positives). The pooled
specificity (46 of 50) is not the curated one (40 of 40).
"""

import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.calibration.phasec import SETS, count_clusters, entries_from_phasec

KEY = "S3-Eurotiomycetes:Anid_EMENI"
SOURCE = SETS[KEY][0]
XFAIL_31 = pytest.mark.xfail(strict=True, reason="Task 3.1: per-stratum basis in phasec")
XFAIL_35 = pytest.mark.xfail(strict=True, reason="Task 3.5: leakage_exception and note wording")


def _rows():
    """(sha, class, stratum, cluster) for every protein of the fixture."""
    rows = [(f"pos{i}", "pos", "", f"cp{i % 110}") for i in range(120)]
    rows += [(f"int{i}", "neg", "N-int", f"cn{i}") for i in range(20)]
    rows += [(f"sec{i}", "neg", "N-sec", f"cn{20 + i % 18}") for i in range(20)]
    rows += [(f"tm{i}", "neg", "PM-TM", f"cm{i}") for i in range(10)]
    return rows


@pytest.fixture
def phasec_files(tmp_path):
    rows = _rows()
    with gzip.open(tmp_path / "clusters.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seq_sha256", "cluster_id"])
        w.writerows((sha, cluster) for sha, _, _, cluster in rows)
    with gzip.open(tmp_path / "eval_table.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["seq_sha256", "origin", "class", "homology_only", "source_ids", "stratum"])
        w.writerows((sha, "go", cls, "no", SOURCE, st) for sha, cls, st, _ in rows)

    def cell(r, f):
        return {
            "recall": {"value": r[0], "lo": r[1], "hi": r[2]},
            "fpr": {"value": f[0], "lo": f[1], "hi": f[2]} if f else {"value": None},
        }

    def stratum(f, n):
        return {"V-go": {"R0": {"recall": {"value": None}, "fpr": f}}}, n

    metrics = {"all": {"V-go": {"R0": cell((0.9, 0.85, 0.95), (4 / 50, 0.03, 0.15))}}}
    n = {"all": {"pos": 120, "neg": 50}}
    for name, count, f in (
        ("N-int", 20, {"value": 0.0, "lo": 0.0, "hi": 0.0}),
        ("N-sec", 20, {"value": 0.0, "lo": 0.0, "hi": 0.0}),
        ("PM-TM", 10, {"value": 0.4, "lo": 0.15, "hi": 0.7}),
    ):
        metrics[name], n[name] = stratum(f, count)[0], {"neg": count}
    data = {
        "test_sets": {
            KEY: {
                "label": "estimate",
                "truth": {"direct": {"n": n, "metrics": metrics}},
            }
        }
    }
    path = tmp_path / "metrics.json"
    path.write_text(json.dumps(data))
    counts = count_clusters(tmp_path / "clusters.tsv.gz", tmp_path / "eval_table.tsv.gz")
    return path, counts


@XFAIL_31
def test_count_clusters_counts_the_curated_negatives_apart(phasec_files):
    _, counts = phasec_files
    assert counts[SOURCE]["neg"] == (50, 48)
    assert counts[SOURCE]["neg_curated"] == (40, 38)


@XFAIL_31
def test_the_headline_negatives_are_the_curated_strata_only(phasec_files):
    path, counts = phasec_files
    (entry,) = entries_from_phasec(path, {KEY: [162425]}, counts)
    m = entry["measure"]
    assert (m["n_neg"], m["n_clusters_neg"]) == (40, 38)
    assert m["negative_basis"] == {
        "N-int": {"basis": "curated", "n": 20},
        "N-sec": {"basis": "curated", "n": 20},
        "PM-TM": {"basis": "assumed", "n": 10},
    }
    assert m["specificity"]["value"] == pytest.approx(1.0)  # 40 of 40; the pooled value is 0.92
    assert m["strata"]["PM-TM"]["value"] == pytest.approx(0.6)  # kept in its own row
    assert entry["status"] == "estimated"


@XFAIL_31
def test_without_stratum_counts_the_entry_has_no_basis_and_is_at_most_smoke(phasec_files, tmp_path):
    path, _ = phasec_files
    pooled = {SOURCE: {"pos": (120, 110), "neg": (50, 48)}}
    (entry,) = entries_from_phasec(path, {KEY: [162425]}, pooled)
    assert "negative_basis" not in entry["measure"]
    assert entry["status"] == "smoke"


@XFAIL_35
def test_the_r0_entry_carries_the_leakage_exception_and_a_legal_leakage_value(phasec_files):
    from cellsurface_sorting_hat.call_status import LEAKAGE_VALUES, leakage_of

    path, counts = phasec_files
    (entry,) = entries_from_phasec(path, {KEY: [162425]}, counts)
    m = entry["measure"]
    exc = m["leakage_exception"]
    assert set(exc) == {"reason", "signer", "date"} and all(exc.values())
    assert "uncapped by owner decision" in exc["reason"]
    assert leakage_of(m) in LEAKAGE_VALUES
    assert entry["status"] == "estimated"  # the exception lifts the cap


@XFAIL_35
def test_the_reworded_note_parses_as_unknown():
    from cellsurface_sorting_hat.calibration.phasec import LEAKAGE_NOTE
    from cellsurface_sorting_hat.call_status import leakage_of

    assert LEAKAGE_NOTE.startswith("leakage: unknown (")
    assert leakage_of({"notes": "call=signal_peptide_protein; " + LEAKAGE_NOTE}) == "unknown"
