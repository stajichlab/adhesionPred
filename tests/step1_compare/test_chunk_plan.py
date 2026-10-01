import json

import chunk_plan
import pytest
import seqhash
import seqsets
import truth_table
from conftest import load_script


def _unique(lengths):
    seqs = {}
    for i, n in enumerate(lengths):
        seq = ("MKV" + "ACDEFGHIKLMNPQRSTVWY"[i % 20] * n)[:n]
        seqs[seqhash.seq_sha256(seq)] = seq
    return seqsets.unique_rows(seqs)


def test_chunks_cover_every_sequence_once_and_respect_the_budget():
    unique = _unique([50, 60, 70, 300, 900, 1100, 1500, 40, 41, 2000])
    chunks = chunk_plan.build_chunks(unique, "nterm", 1200)
    rows = [r for c in chunks for r in c.rows]
    assert sorted(rows) == list(range(len(unique)))
    for c in chunks:
        assert c.residues <= 1200 or len(c.rows) == 1
    ids = [c.chunk_id for c in chunks]
    assert ids == [f"nterm_{i:04d}" for i in range(len(chunks))]


def test_cterm_chunks_hold_only_long_sequences():
    unique = _unique([50, 1100, 1500, 2000, 1022])
    chunks = chunk_plan.build_chunks(unique, "cterm", 5000)
    rows = sorted(r for c in chunks for r in c.rows)
    long_rows = sorted(int(u["row"]) for u in unique if int(u["length"]) > 1022)
    assert rows == long_rows
    assert all(c.residues == 1022 * len(c.rows) for c in chunks)


def test_chunk_membership_is_deterministic_and_length_sorted():
    unique = _unique([500, 20, 300, 800, 45, 1000, 60])
    a = chunk_plan.build_chunks(unique, "nterm", 900)
    b = chunk_plan.build_chunks(list(reversed(unique)), "nterm", 900)
    assert [c.members_sha256 for c in a] == [c.members_sha256 for c in b]
    lengths = [int(unique[r]["length"]) for c in a for r in c.rows]
    assert lengths == sorted(lengths)


def test_chunk_residues_from_rate_formula():
    assert chunk_plan.chunk_residues_from_rate(29_500) == 8_850_000
    assert chunk_plan.chunk_residues_from_rate(1.0) == chunk_plan.MIN_CHUNK
    with pytest.raises(ValueError):
        chunk_plan.chunk_residues_from_rate(0)


def test_plan_jobs_worked_example():
    # Assumed rate (review 9.1, ESM-2 150M) until J0 measures 8M and 35M.
    plan = chunk_plan.plan_jobs(
        {"esm2_t6_8M_UR50D": 29_500.0, "esm2_t12_35M_UR50D": 29_500.0},
        {"esm2_t6_8M_UR50D": 10.0, "esm2_t12_35M_UR50D": 10.0},
        37_860_229,
    )
    assert plan["total_seconds"] == pytest.approx(2586.8, abs=0.1)
    assert plan["n_jobs"] == 1
    assert plan["time_minutes"] == 75


def test_plan_jobs_splits_slow_rates_into_jobs_under_the_target():
    plan = chunk_plan.plan_jobs({"m": 2_000.0}, {"m": 0.0}, 37_860_229)
    assert plan["n_jobs"] == 6  # 18,930 s * 1.25 / 4,500 s = 5.26 -> 6
    assert plan["seconds_per_job"] * chunk_plan.SAFETY <= chunk_plan.TARGET_SECONDS


def test_chunks_for_job_round_robin():
    ids = [f"c{i}" for i in range(7)]
    jobs = [chunk_plan.chunks_for_job(ids, i, 3) for i in range(3)]
    assert jobs == [["c0", "c3", "c6"], ["c1", "c4"], ["c2", "c5"]]
    assert sorted(sum(jobs, [])) == sorted(ids)
    with pytest.raises(ValueError):
        chunk_plan.chunks_for_job(ids, 3, 3)


def _plan_inputs(tmp_path, lengths):
    out = tmp_path / "phaseb"
    out.mkdir(parents=True)
    unique = _unique(lengths)
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique)
    return out, unique


def test_plan_script_reads_j0_and_writes_a_consistent_plan(tmp_path):
    out, unique = _plan_inputs(tmp_path, [50, 400, 1100, 1300, 900])
    j0 = {
        "best": {
            "esm2_t6_8M_UR50D": {"batch_size": 32, "residues_per_s": 5.0},
            "esm2_t12_35M_UR50D": {"batch_size": 16, "residues_per_s": 4.0},
        },
        "model_load_s": {"esm2_t6_8M_UR50D": 1.0, "esm2_t12_35M_UR50D": 2.0},
    }
    (out / "j0").mkdir()
    (out / "j0" / "throughput.json").write_text(json.dumps(j0))
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path), "--chunk-residues", "1000"]) == 0
    plan = json.loads((out / "job_plan.json").read_text())
    assert plan["rate_source"] == "J0" and plan["batch_size"]["esm2_t12_35M_UR50D"] == 16
    assert plan["residues_per_model"] == 50 + 400 + 1022 + 1022 + 900 + 1022 + 1022
    chunks = chunk_plan.read_plan(out)
    assert sum(len(c.rows) for c in chunks if c.window == "nterm") == len(unique)
    assert sum(len(c.rows) for c in chunks if c.window == "cterm") == 2


def test_plan_script_stops_without_j0(tmp_path, capsys):
    _plan_inputs(tmp_path, [50, 60])
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path)]) == 2
    assert "run J0 first" in capsys.readouterr().err
    assert not (tmp_path / "phaseb" / "chunk_plan.tsv").exists()
    assert plan_script.main(["--work-dir", str(tmp_path), "--rate", "100"]) == 0
    plan = json.loads((tmp_path / "phaseb" / "job_plan.json").read_text())
    assert plan["rate_source"] == "assumed"


def test_read_plan_detects_edited_members(tmp_path):
    out, _ = _plan_inputs(tmp_path, [50, 60, 70])
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path), "--rate", "100"]) == 0
    rows = truth_table.read_tsv(out / "chunk_members.tsv.gz")
    rows[0]["seq_sha256"] = "f" * 64
    truth_table.write_tsv(out / "chunk_members.tsv.gz", chunk_plan.MEMBER_COLUMNS, rows)
    with pytest.raises(ValueError, match="members differ"):
        chunk_plan.read_plan(out)


def _planned(tmp_path):
    out, unique = _plan_inputs(tmp_path, [50, 60, 1100])
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path), "--rate", "100"]) == 0
    seq_by_row = {int(r["row"]): r["sequence"] for r in unique}
    return out, unique, seq_by_row


def test_plan_check_passes_on_the_planned_set(tmp_path):
    out, _, seq_by_row = _planned(tmp_path)
    chunk_plan.check_plan_is_current(out, chunk_plan.read_plan(out), seq_by_row)


def test_plan_check_detects_another_unique_file(tmp_path):
    out, unique, seq_by_row = _planned(tmp_path)
    truth_table.write_tsv(out / "unique_sequences.tsv.gz", seqsets.UNIQUE_COLUMNS, unique[:2])
    with pytest.raises(chunk_plan.StalePlanError, match="another unique_sequences"):
        chunk_plan.check_plan_is_current(out, chunk_plan.read_plan(out), seq_by_row)


def test_plan_check_detects_a_changed_sequence_even_with_a_matching_file_hash(tmp_path):
    import manifest

    out, _, seq_by_row = _planned(tmp_path)
    plan = json.loads((out / "job_plan.json").read_text())
    plan["unique_sequences_sha256"] = manifest.sha256_file(out / "unique_sequences.tsv.gz")
    (out / "job_plan.json").write_text(json.dumps(plan))
    seq_by_row[0] = "MKWWWW"
    with pytest.raises(chunk_plan.StalePlanError, match="does not hold the planned sequence"):
        chunk_plan.check_plan_is_current(out, chunk_plan.read_plan(out), seq_by_row)
