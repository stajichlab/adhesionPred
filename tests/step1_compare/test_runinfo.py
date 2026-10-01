import json
import platform

import pytest
import runinfo


def test_git_commit_and_python_version():
    commit = runinfo.git_commit()
    assert commit == "unknown" or len(commit) == 40
    assert runinfo.python_version() == platform.python_version()


def test_git_commit_is_unknown_outside_a_repository(monkeypatch, tmp_path):
    monkeypatch.setenv("PROJ_ROOT", str(tmp_path))
    assert runinfo.git_commit() == "unknown"


def test_require_full_accepts_true_and_skips_when_partial_allowed(tmp_path):
    log = tmp_path / "extract_log.json"
    log.write_text(json.dumps({"all_sources": True}))
    runinfo.require_full(log, "truth_set.tsv.gz", allow_partial=False)
    # With allow_partial the log is not read at all, so a missing log is accepted.
    runinfo.require_full(tmp_path / "missing.json", "x", allow_partial=True)


@pytest.mark.parametrize("value", [False, None, "true", 1])
def test_require_full_refuses_a_partial_log_with_a_clear_message(tmp_path, value):
    log = tmp_path / "extract_log.json"
    log.write_text(json.dumps({"all_sources": value}))
    with pytest.raises(runinfo.PartialInputError) as info:
        runinfo.require_full(log, "truth_set.tsv.gz", allow_partial=False)
    msg = str(info.value)
    assert "extract_log.json does not say all_sources: true" in msg
    assert f"(value: {value!r})" in msg
    assert "truth_set.tsv.gz may hold only some sources" in msg
    assert "rerun 01_extract_go_truth.py without --sources" in msg
    assert "--allow-partial-truth-set" in msg
    assert isinstance(info.value, ValueError)  # every script's main catches ValueError


def test_require_full_names_02_for_sequence_run(tmp_path):
    log = tmp_path / "sequence_run.json"
    log.write_text(json.dumps({"all_sources": False}))
    with pytest.raises(runinfo.PartialInputError, match="rerun 02_attach_sequences.py"):
        runinfo.require_full(log, "truth_sequences.tsv.gz", allow_partial=False)


@pytest.mark.parametrize("content", [None, "{not json", "[1, 2]"])
def test_require_full_refuses_a_missing_or_unreadable_log(tmp_path, content):
    log = tmp_path / "extract_log.json"
    if content is not None:
        log.write_text(content)
    with pytest.raises(runinfo.PartialInputError, match="cannot read"):
        runinfo.require_full(log, "truth_set.tsv.gz", allow_partial=False)


def test_says_all_sources(tmp_path):
    log = tmp_path / "x.json"
    assert runinfo.says_all_sources(log) is False
    log.write_text(json.dumps({"all_sources": True}))
    assert runinfo.says_all_sources(log) is True
    log.write_text("[]")
    assert runinfo.says_all_sources(log) is False


def test_atomic_write_all_writes_every_file_and_no_temp(tmp_path):
    out = tmp_path / "out"
    runinfo.atomic_write_all(
        out,
        {
            "a.txt": lambda p: p.write_text("A"),
            "b.json": lambda p: runinfo.write_json(p, {"z": 1, "a": 2}),
        },
    )
    assert (out / "a.txt").read_text() == "A"
    assert (out / "b.json").read_text() == '{\n  "a": 2,\n  "z": 1\n}\n'
    assert sorted(p.name for p in out.iterdir()) == ["a.txt", "b.json"]


def test_atomic_write_all_failure_keeps_old_outputs_and_removes_temps(tmp_path):
    (tmp_path / "a.txt").write_text("old A")
    (tmp_path / "b.txt").write_text("old B")

    def fail(path):
        path.write_text("half")
        raise OSError("disk full")

    with pytest.raises(OSError, match="disk full"):
        runinfo.atomic_write_all(
            tmp_path, {"a.txt": lambda p: p.write_text("new A"), "b.txt": fail}
        )
    assert (tmp_path / "a.txt").read_text() == "old A"
    assert (tmp_path / "b.txt").read_text() == "old B"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.txt", "b.txt"]
