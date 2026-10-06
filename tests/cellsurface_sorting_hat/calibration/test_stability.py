import csv
import gzip

from cellsurface_sorting_hat.calibration.stability import compare_runs


def write_run(path, proteins, calls):
    path.mkdir()
    with gzip.open(path / "proteins.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["id", "sha256", "taxon", "state", "note", "trailing_stop", "ambiguous_fraction"]
        )
        for pid, sha in proteins:
            w.writerow([pid, sha, 1, "ok", "", 0, 0])
    with gzip.open(path / "calls.long.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein", "call", "variant", "value", "status", "status_basis", "other_basis"])
        for pid, call, variant, value in calls:
            w.writerow([pid, call, variant, value, "unvalidated", "", ""])


def test_only_identical_sequences_are_compared_by_sha256_not_by_id(tmp_path):
    write_run(
        tmp_path / "a",
        [("a1", "s1"), ("a2", "s2"), ("a3", "s3")],
        [("a1", "x", "", "called"), ("a2", "x", "", "called"), ("a3", "x", "", "called")],
    )
    write_run(
        tmp_path / "b",
        [("b9", "s1"), ("b8", "s2"), ("b7", "s4")],
        [("b9", "x", "", "called"), ("b8", "x", "", "not_called"), ("b7", "x", "", "called")],
    )
    r = compare_runs(tmp_path / "a", tmp_path / "b")
    assert (r["n_a"], r["n_b"], r["n_identical_sequences"]) == (3, 3, 2)
    assert dict(r["per_call"][("x", "")]) == {"agree": 1, "differ": 1}
