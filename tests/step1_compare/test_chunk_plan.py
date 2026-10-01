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


def _j0(r8=2000.0, r35=1500.0, **over):
    """A J0 throughput JSON in the Task 6 schema (device cuda, runs with status ok)."""
    j0 = {
        "schema": "step1-phaseb-j0-throughput/2",
        "device": "cuda",
        "gpu_name": "NVIDIA test GPU",
        "sample_sha256": "a" * 64,
        "best": {
            "esm2_t6_8M_UR50D": {"batch_size": 32, "residues_per_s": r8},
            "esm2_t12_35M_UR50D": {"batch_size": 16, "residues_per_s": r35},
        },
        "runs": [
            {"model": "esm2_t6_8M_UR50D", "batch_size": 32, "status": "ok"},
            {"model": "esm2_t12_35M_UR50D", "batch_size": 16, "status": "ok"},
        ],
        "model_load_s": {"esm2_t6_8M_UR50D": 1.0, "esm2_t12_35M_UR50D": 2.0},
    }
    j0.update(over)
    return j0


def _write_j0(out, j0):
    (out / "j0").mkdir(exist_ok=True)
    (out / "j0" / "throughput.json").write_text(json.dumps(j0))


def test_plan_script_reads_j0_and_writes_a_consistent_plan(tmp_path):
    out, unique = _plan_inputs(tmp_path, [50, 400, 1100, 1300, 900])
    j0 = _j0(r8=2000.0, r35=1500.0)
    _write_j0(out, j0)
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path), "--chunk-residues", "1000"]) == 0
    plan = json.loads((out / "job_plan.json").read_text())
    assert plan["rate_source"] == "J0" and plan["batch_size"]["esm2_t12_35M_UR50D"] == 16
    assert plan["j0_gpu_name"] == "NVIDIA test GPU" and plan["j0_device"] == "cuda"
    assert plan["j0_sample_sha256"] == "a" * 64
    import manifest

    assert plan["throughput_sha256"] == manifest.sha256_file(out / "j0" / "throughput.json")
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


# ---- fix round 1: job time from the largest job, rate floor, J0 checks --------------------


def test_time_limit_follows_the_largest_job_not_the_average():
    # 2 chunks of very different size, 1 job each would be needed; rates 2,000 residues/s x2.
    rates, load = {"a": 2000.0, "b": 2000.0}, {"a": 0.0, "b": 0.0}
    chunks = [40_000_000, 1_000_000]
    plan = chunk_plan.plan_jobs(rates, load, sum(chunks), chunks)
    # total = 2 * 41e6 / 2000 = 41,000 s -> ceil(41000 * 1.25 / 4500) = 12 jobs, capped at 2.
    assert plan["n_jobs"] == 2
    assert plan["longest_job_seconds"] == pytest.approx(40_000_000 * (1 / 2000 + 1 / 2000))
    assert plan["longest_job_seconds"] == 40_000.0
    assert plan["time_minutes"] * 60 >= plan["longest_job_seconds"] * 1.5
    # the old rule (average job) would give a limit shorter than the work
    assert (plan["total_seconds"] / plan["n_jobs"]) < plan["longest_job_seconds"]


def test_n_jobs_never_exceeds_the_chunk_count():
    plan = chunk_plan.plan_jobs({"m": 3.0}, {"m": 0.0}, 37_860_229, [10_000, 10_000, 10_000])
    assert plan["n_jobs"] == 3


def test_plan_jobs_with_different_model_rates():
    rates = {"fast": 40_000.0, "slow": 10_000.0}
    plan = chunk_plan.plan_jobs(rates, {"fast": 5.0, "slow": 10.0}, 1_000_000, [1_000_000])
    # 1e6/4e4 + 1e6/1e4 + 15 = 25 + 100 + 15
    assert plan["total_seconds"] == 140.0
    assert plan["longest_job_seconds"] == 140.0
    assert plan["time_minutes"] == 14  # ceil((140 * 1.5 + 600) / 60) = ceil(13.5)


def test_chunk_size_follows_the_slowest_model():
    assert chunk_plan.chunk_residues_from_rate(min(40_000.0, 20_000.0)) == 6_000_000
    plan_script = load_script("06_plan_embedding")
    rates, _, _ = plan_script.rates_from_j0(_j0(r8=40_000.0, r35=20_000.0), plan_script.MODELS)
    assert chunk_plan.chunk_residues_from_rate(min(rates.values())) == 6_000_000


def test_script_chunk_size_uses_min_rate_when_8m_is_faster(tmp_path):
    out, _ = _plan_inputs(tmp_path, [50, 400, 1100])
    _write_j0(out, _j0(r8=40_000.0, r35=20_000.0))
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path)]) == 0
    assert json.loads((out / "job_plan.json").read_text())["chunk_residues"] == 6_000_000


def test_chunk_residues_rounds_down_to_the_unit():
    # 29,555 * 300 = 8,866,500 -> 8,860,000 (floor), not 8,870,000 (round)
    assert chunk_plan.chunk_residues_from_rate(29_555) == 8_860_000


def test_sort_ties_break_by_hash_so_shuffled_input_gives_the_same_chunks():
    import random

    seqs = {}
    for i in range(30):
        seq = "M" + "ACDEFGHIKLMNPQRSTVWY"[i % 20] * 9 + "ACDEFGHIKLMNPQRSTVWY"[i // 20]
        seqs[seqhash.seq_sha256(seq)] = seq
    unique = seqsets.unique_rows(seqs)
    assert len({len(u["sequence"]) for u in unique}) == 1
    a = chunk_plan.build_chunks(unique, "nterm", 55)
    shuffled = list(unique)
    random.Random(1).shuffle(shuffled)
    b = chunk_plan.build_chunks(shuffled, "nterm", 55)
    assert len(a) > 3
    assert [c.members_sha256 for c in a] == [c.members_sha256 for c in b]
    assert [c.rows for c in a] == [c.rows for c in b]


def test_script_stops_when_chunk_residues_exceeds_the_rate_derived_size(tmp_path, capsys):
    out, _ = _plan_inputs(tmp_path, [50, 60])
    _write_j0(out, _j0(r8=2000.0, r35=2000.0))
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path), "--chunk-residues", "40000000"]) == 2
    err = capsys.readouterr().err
    assert "40000000" in err and "600000" in err
    assert not (out / "chunk_plan.tsv").exists()


@pytest.mark.parametrize(
    "j0_over, text",
    [
        ({"device": "cpu"}, "not cuda"),
        ({"schema": "step1-phaseb-j0-throughput/1"}, "schema"),
        (
            {
                "runs": [
                    {"model": "esm2_t6_8M_UR50D", "batch_size": 32, "status": "batch_failures"},
                    {"model": "esm2_t12_35M_UR50D", "batch_size": 16, "status": "ok"},
                ]
            },
            "esm2_t6_8M_UR50D at batch size 32 is not ok",
        ),
    ],
)
def test_script_stops_on_an_unusable_j0(tmp_path, capsys, j0_over, text):
    out, _ = _plan_inputs(tmp_path, [50, 60])
    _write_j0(out, _j0(**j0_over))
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path)]) == 2
    assert text in capsys.readouterr().err
    assert not (out / "chunk_plan.tsv").exists()


def test_script_stops_below_the_minimum_rate(tmp_path, capsys):
    out, _ = _plan_inputs(tmp_path, [50, 60])
    _write_j0(out, _j0(r8=364.3, r35=133.4))  # the CPU pilot numbers of the plan
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path)]) == 2
    assert "CPU or failed run" in capsys.readouterr().err
    assert not (out / "job_plan.json").exists()


@pytest.mark.parametrize("bad", ["fast", None, True, -5])
def test_script_stops_on_a_non_numeric_rate(tmp_path, capsys, bad):
    out, _ = _plan_inputs(tmp_path, [50, 60])
    _write_j0(out, _j0(r8=bad))
    plan_script = load_script("06_plan_embedding")
    assert plan_script.main(["--work-dir", str(tmp_path)]) == 2
    err = capsys.readouterr().err
    assert "STOP:" in err and "no successful run" in err and "Traceback" not in err
    assert not (out / "job_plan.json").exists()
