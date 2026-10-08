import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/map_to_proteomes.py"
spec = importlib.util.spec_from_file_location("map_to_proteomes", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def hit(q, s, pident, qcov, bits=100.0):
    return {"qseqid": q, "sseqid": s, "pident": pident, "qcovs": qcov, "bitscore": bits}


def test_threshold_keeps_96_percent_and_drops_80_percent():
    hits = [hit("A", "p1", 96.0, 95.0), hit("B", "p2", 80.0, 99.0), hit("C", "p3", 99.0, 80.0)]
    mapped, dropped = m.filter_hits(hits, min_pident=95.0, min_qcov=90.0)
    assert [r["qseqid"] for r in mapped] == ["A"]
    assert {r["qseqid"] for r in dropped} == {"B", "C"}


def test_best_hit_per_query_by_bitscore():
    hits = [hit("A", "p1", 99.0, 99.0, 150.0), hit("A", "p2", 97.0, 99.0, 120.0)]
    mapped, _ = m.filter_hits(hits, 95.0, 90.0)
    assert [(r["qseqid"], r["sseqid"]) for r in mapped] == [("A", "p1")]


def test_one_proteome_protein_is_mapped_once_and_conflict_is_listed():
    hits = [hit("A", "p1", 99.0, 99.0, 200.0), hit("B", "p1", 96.0, 99.0, 150.0)]
    mapped, dropped = m.filter_hits(hits, 95.0, 90.0)
    assert [r["qseqid"] for r in mapped] == ["A"]
    assert [(r["qseqid"], r.get("reason")) for r in dropped] == [("B", "subject taken by A")]
