"""PR 6 drift check: classes, the rule on supported hits, the shuffled null, the extraction and the job script."""

import csv
import gzip
import importlib.util
import random
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "analysis/wall_fusion_truth"
SBATCH = WF / "pr6_drift_check.sbatch"
RECORD = ROOT / "docs/reports/data/pr6_gap_uniref50_2026-10-22"


def load(name):
    spec = importlib.util.spec_from_file_location(name, WF / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


dc = load("pr6_drift_check")


def dom_line(target, tlen, qlen, full_score, ie, q_iv, t_iv, query="Q"):
    f = [
        target,
        "-",
        str(tlen),
        query,
        "-",
        str(qlen),
        "1e-30",
        str(full_score),
        "0.1",
        "1",
        "1",
        str(ie),
        str(ie),
        str(full_score),
    ]
    f += [
        "0.1",
        str(q_iv[0]),
        str(q_iv[1]),
        str(t_iv[0]),
        str(t_iv[1]),
        str(t_iv[0]),
        str(t_iv[1]),
        "0.9",
        "desc",
    ]
    return " ".join(f)


def write_dom(path, lines, trailer=True):
    body = ["# target name accession tlen query name"] + lines
    if trailer:
        body += ["# Program:         phmmer", "# [ok]"]
    path.write_text("\n".join(body) + "\n")


# --- helpers -------------------------------------------------------------------------------------------------


def test_union_length_merges_overlapping_and_adjacent_intervals():
    assert dc.union_length([(1, 10), (5, 20)]) == 20
    assert dc.union_length([(1, 10), (11, 20)]) == 20
    assert dc.union_length([(1, 5), (10, 12)]) == 8
    assert dc.union_length([]) == 0


def test_composition_measures():
    assert dc.stnqgp_fraction("SSSTTTNNNQQQ") == 1.0
    assert dc.stnqgp_fraction("ACDEFHIKLMRVWY") == 0.0
    poly = "S" * 40
    varied = "ACDEFGHIKLMNPQRSTVWY" * 2
    assert dc.low_complexity_fraction(poly) == 1.0
    assert dc.low_complexity_fraction(varied) == 0.0
    assert dc.low_complexity_fraction("ACD") == 0.0  # shorter than the window


def test_parse_domtbl_needs_the_trailer_and_unions_the_domains(tmp_path):
    p = tmp_path / "a.dom"
    write_dom(
        p,
        [
            dom_line("H1", 200, 100, 80.0, 1e-20, (1, 40), (10, 50)),
            dom_line("H1", 200, 100, 80.0, 1e-3, (30, 60), (60, 90)),
        ],
    )
    d = dc.parse_domtbl(p)["H1"]
    assert d["best_ie"] == 1e-20 and d["qlen"] == 100 and d["tlen"] == 200
    assert dc.union_length(d["q"]) == 60 and dc.union_length(d["t"]) == 72
    write_dom(
        tmp_path / "b.dom",
        [dom_line("H1", 200, 100, 80.0, 1e-20, (1, 40), (10, 50))],
        trailer=False,
    )
    with pytest.raises(ValueError, match="trailer"):
        dc.parse_domtbl(tmp_path / "b.dom")


# --- classes -------------------------------------------------------------------------------------------------


def info(score=100.0, ie=1e-20, q=((1, 80),), t=((1, 90),), qlen=100, tlen=100):
    return {"best_ie": ie, "q": list(q), "t": list(t), "score": score, "qlen": qlen, "tlen": tlen}


def test_the_three_classes_hand_made():
    sup = dc.classify({}, info(), info(score=5.0))
    assert sup["class"] == "supported_homolog" and sup["null_ratio"] == "0.050"
    bias = dc.classify(
        {}, info(), info(score=80.0)
    )  # a shuffled hit scores 80 percent of the real score
    assert bias["class"] == "bias_driven"
    edge = dc.classify(
        {}, info(), info(score=50.0)
    )  # exactly NULL_RATIO_MAX is bias_driven (the rule says LESS than)
    assert edge["class"] == "bias_driven"
    weak = dc.classify({}, info(ie=1e-3), info(score=1.0))
    assert weak["class"] == "unsupported"
    short = dc.classify(
        {}, info(q=((1, 20),), t=((1, 20),)), info(score=1.0)
    )  # covers 20 percent of both
    assert short["class"] == "unsupported" and short["query_cover"] == "0.200"
    absent = dc.classify({}, None, None)
    assert absent["class"] == "unsupported" and absent["real_score"] == ""
    # coverage of EITHER the query or the hit is enough (a short query fully inside a long hit)
    inside = dc.classify({}, info(q=((1, 90),), t=((1, 30),), tlen=300), None)
    assert inside["class"] == "supported_homolog" and inside["null_score"] == "0.0"


def row(clade, cls):
    return {"clade_fine": clade, "class": cls}


def test_rule_on_supported_hits_only_counts_named_clades():
    sup = [row("A", "supported_homolog")] * 3 + [row("B", "supported_homolog")] * 2
    v, reason, n_named, n_un = dc.verdict(sup)
    assert v == "supported_rule_met" and reason == "" and n_named == 2 and n_un == 0
    # bias-driven and unsupported hits do not count
    mixed = sup[:4] + [row("B", "bias_driven")] * 10 + [row("C", "unsupported")] * 10
    assert dc.verdict(mixed)[0] == "not_supported" and "only 4 supported" in dc.verdict(mixed)[1]
    # the (no clade) bucket counts in the total, not as a clade
    bucket = [row("A", "supported_homolog")] * 3 + [row(dc.NO_CLADE, "supported_homolog")] * 3
    v2, reason2, n2, un2 = dc.verdict(bucket)
    assert v2 == "not_supported" and n2 == 1 and un2 == 3 and "1 named clade" in reason2
    assert dc.verdict([])[0] == "not_supported"


# --- extraction ------------------------------------------------------------------------------------------------


def write_record(path):
    path.mkdir()
    rows = [
        ("Q1", "phmmer", "UniRef50_P1", "Sordariomycetes", "outside_clade"),
        ("Q1", "jackhmmer", "UniRef50_P1", "Sordariomycetes", "outside_clade"),
        ("Q1", "jackhmmer", "UniRef50_P2", "", "outside_clade"),
        ("Q1", "jackhmmer", "UniRef50_P3", "Eurotiomycetes", "inside_clade"),
        ("Q1", "jackhmmer", "UniRef50_P4", "", "non_fungal"),
        ("Q2", "jackhmmer", "UniRef50_P5", "Leotiomycetes", "outside_clade"),
        ("Q9", "jackhmmer", "UniRef50_P6", "Leotiomycetes", "outside_clade"),
    ]
    with gzip.open(path / "gap_homologs.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            ["query", "method", "hit_accession", "evalue", "score", "taxid", "clade_fine", "where"]
        )
        for q, m, h, c, wh in rows:
            w.writerow([q, m, h, "1e-9", "40.0", "5141", c, wh])


def test_extract_ids_takes_jackhmmer_outside_fungal_hits_of_the_chosen_queries(tmp_path):
    write_record(tmp_path / "rec")
    rows, ids = dc.extract_ids(tmp_path / "rec", tmp_path, ("Q1", "Q2"))
    assert [(r[0], r[1], r[2], r[6]) for r in rows] == [
        ("Q1", "UniRef50_P1", "Sordariomycetes", "yes"),
        ("Q1", "UniRef50_P2", "(no clade)", "no"),
        ("Q2", "UniRef50_P5", "Leotiomycetes", "no"),
    ]
    assert ids == ["UniRef50_P1", "UniRef50_P2", "UniRef50_P5"]


def test_extract_ids_on_the_committed_round_2_record_matches_its_counts():
    rows, ids = dc.extract_ids(RECORD, ROOT, dc.DEFAULT_QUERIES)
    per = {}
    for r in rows:
        per[r[0]] = per.get(r[0], 0) + 1
    by_gene = {
        x["genes"]: x for x in csv.DictReader(open(RECORD / "gap_per_query.tsv"), delimiter="\t")
    }
    assert per["P32323"] == int(by_gene["aga1"]["jackhmmer_outside_seq"]) == 64
    assert per["Q07730"] == int(by_gene["ece1"]["jackhmmer_outside_seq"]) == 126
    assert per["A0ACE4KFX7"] == int(by_gene["hyd1_like"]["jackhmmer_outside_seq"]) == 43
    assert len(ids) == len(set(ids))


# --- split and the shuffled null ---------------------------------------------------------------------------------


def make_fastas(tmp_path):
    hits = {
        "UniRef50_P1": "ACDEFGHIKLMNPQRSTVWY" * 3,
        "UniRef50_P2": "STSTSTSTSTSTNQNQNQ" * 3,
        "UniRef50_P5": "MKVLAAGIVG" * 5,
    }
    queries = {"Q1": "ACDEFGHIKL" * 8, "Q2": "MKVLAAGIVG" * 8}
    dc.write_fasta(tmp_path / "hits.faa", hits)
    dc.write_fasta(tmp_path / "queries.faa", queries)
    man = [
        ["Q1", "UniRef50_P1", "A", "1", "1e-9", "40", "yes"],
        ["Q1", "UniRef50_P2", dc.NO_CLADE, "1", "1e-9", "40", "no"],
        ["Q2", "UniRef50_P5", "B", "1", "1e-9", "40", "no"],
    ]
    dc.write_tsv(tmp_path / "manifest.tsv", dc.MANIFEST_HEADER, man)
    return hits, queries


def test_split_shuffles_residues_keeps_composition_and_is_deterministic(tmp_path):
    hits, _ = make_fastas(tmp_path)
    c1 = dc.split_inputs(
        tmp_path / "manifest.tsv",
        tmp_path / "hits.faa",
        tmp_path / "queries.faa",
        tmp_path / "o1",
        1,
    )
    dc.split_inputs(
        tmp_path / "manifest.tsv",
        tmp_path / "hits.faa",
        tmp_path / "queries.faa",
        tmp_path / "o2",
        1,
    )
    assert c1 == {"Q1": 2, "Q2": 1}
    n1 = dc.read_fasta(tmp_path / "o1/Q1.null.faa")
    assert (tmp_path / "o1/Q1.null.faa").read_text() == (tmp_path / "o2/Q1.null.faa").read_text()
    for h, s in n1.items():
        assert sorted(s) == sorted(hits[h]) and s != hits[h]
    # a different seed gives a different shuffle; one query's shuffle does not depend on the others
    dc.split_inputs(
        tmp_path / "manifest.tsv",
        tmp_path / "hits.faa",
        tmp_path / "queries.faa",
        tmp_path / "o3",
        2,
    )
    assert dc.read_fasta(tmp_path / "o3/Q1.null.faa") != n1
    assert dc.read_fasta(tmp_path / "o1/Q1.hits.faa") == {
        h: hits[h] for h in ("UniRef50_P1", "UniRef50_P2")
    }


def test_split_stops_on_a_missing_hit_or_query(tmp_path):
    make_fastas(tmp_path)
    dc.write_fasta(tmp_path / "hits2.faa", {"UniRef50_P1": "ACDE"})
    with pytest.raises(ValueError, match="not in"):
        dc.split_inputs(
            tmp_path / "manifest.tsv",
            tmp_path / "hits2.faa",
            tmp_path / "queries.faa",
            tmp_path / "o",
            1,
        )
    dc.write_fasta(tmp_path / "q2.faa", {"Q1": "ACDE"})
    with pytest.raises(ValueError, match="query Q2"):
        dc.split_inputs(
            tmp_path / "manifest.tsv", tmp_path / "hits.faa", tmp_path / "q2.faa", tmp_path / "o", 1
        )


def test_the_shuffle_is_a_permutation_with_a_seeded_stream():
    s = "ACDEFGHIKLMNPQRSTVWY" * 2
    a = dc.shuffled(s, random.Random("1:Q"))
    b = dc.shuffled(s, random.Random("1:Q"))
    assert a == b and sorted(a) == sorted(s) and a != s


# --- summarize end to end ------------------------------------------------------------------------------------------


def test_summarize_end_to_end_classes_rule_and_no_sequences(tmp_path):
    hits, queries = make_fastas(tmp_path)
    real, null = tmp_path / "real", tmp_path / "null"
    real.mkdir()
    null.mkdir()
    # Q1: P1 supported (real strong, null weak); P2 bias_driven (null as strong as the real one)
    write_dom(
        real / "Q1.dom",
        [
            dom_line("UniRef50_P1", 60, 80, 120.0, 1e-30, (1, 70), (1, 55), "Q1"),
            dom_line("UniRef50_P2", 54, 80, 60.0, 1e-8, (1, 60), (1, 50), "Q1"),
        ],
    )
    write_dom(
        null / "Q1.dom", [dom_line("UniRef50_P2", 54, 80, 55.0, 1e-7, (1, 40), (1, 40), "Q1")]
    )
    # Q2: P5 has no real pairwise hit: unsupported
    write_dom(real / "Q2.dom", [])
    write_dom(null / "Q2.dom", [])
    out = tmp_path / "out"
    s = dc.summarize(
        tmp_path / "manifest.tsv", tmp_path / "queries.faa", tmp_path / "hits.faa", real, null, out
    )
    assert [x[0] for x in s] == ["Q1", "Q2"]
    q1 = next(x for x in s if x[0] == "Q1")
    assert q1[1:5] == [2, 1, 1, 0] and q1[7] == "not_supported"
    q2 = next(x for x in s if x[0] == "Q2")
    assert q2[2:5] == [0, 0, 1] and q2[7] == "not_supported"
    rows = list(csv.DictReader(open(out / "drift_hits.tsv"), delimiter="\t"))
    assert {r["hit"]: r["class"] for r in rows} == {
        "UniRef50_P1": "supported_homolog",
        "UniRef50_P2": "bias_driven",
        "UniRef50_P5": "unsupported",
    }
    for p in out.iterdir():
        for line in p.read_text().splitlines():
            assert not re.fullmatch(r"[A-Z]{40,}", line)  # no sequence line in any output
    assert "PROPOSALS" in (out / "summary.txt").read_text()
    # deterministic, and a filled output directory is refused by main
    s2 = dc.summarize(
        tmp_path / "manifest.tsv",
        tmp_path / "queries.faa",
        tmp_path / "hits.faa",
        real,
        null,
        tmp_path / "out2",
    )
    assert (
        s2 == s
        and (out / "drift_hits.tsv").read_text() == (tmp_path / "out2/drift_hits.tsv").read_text()
    )
    assert (
        dc.main(
            [
                "summarize",
                "--manifest",
                str(tmp_path / "manifest.tsv"),
                "--queries-faa",
                str(tmp_path / "queries.faa"),
                "--hits-faa",
                str(tmp_path / "hits.faa"),
                "--real-dir",
                str(real),
                "--null-dir",
                str(null),
                "--out",
                str(out),
            ]
        )
        == 1
    )


def test_summarize_stops_when_a_search_did_not_finish(tmp_path):
    make_fastas(tmp_path)
    real, null = tmp_path / "real", tmp_path / "null"
    real.mkdir()
    null.mkdir()
    for q in ("Q1", "Q2"):
        write_dom(real / f"{q}.dom", [], trailer=False)
        write_dom(null / f"{q}.dom", [])
    with pytest.raises(ValueError, match="trailer"):
        dc.summarize(
            tmp_path / "manifest.tsv",
            tmp_path / "queries.faa",
            tmp_path / "hits.faa",
            real,
            null,
            tmp_path / "o",
        )


# --- the job script ---------------------------------------------------------------------------------------------------


def test_sbatch_syntax_header_and_rules():
    text = SBATCH.read_text()
    assert subprocess.run(["bash", "-n", str(SBATCH)]).returncode == 0
    assert "DRAFT, NOT SUBMITTED" in text.splitlines()[1]
    code = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("#"))
    assert "/tmp" not in code  # the scratch directory, never /tmp (the header comment may say so)
    assert "${SCRATCH:?" in text and "set -euo pipefail" in text
    assert (
        "phmmer -E 10 --domE 10" in text
        and "seqkit grep" in text
        and "module load seqkit/2.4.0 hmmer/3.4" in text
    )
    assert "jackhmmer" not in text.split("set -euo pipefail", 1)[1].replace(
        "no jackhmmer", ""
    )  # no jackhmmer run in the body
    assert '[ ! -e "$OUT" ]' in text


def test_sbatch_pins_the_round_2_file_and_the_three_queries():
    text = SBATCH.read_text()
    prov = (RECORD / "provenance.txt").read_text()
    assert re.search(r"uniprot_fasta_bytes: (\d+)", prov).group(1) == re.search(
        r"EXPECTED_UNIPROT_BYTES=(\d+)", text
    ).group(1)
    sha = re.search(r"uniprot_fasta_sha256: (\w+)", prov).group(1)
    assert sha in text
    queries = re.search(r'QUERIES="([^"]+)"', text).group(1).split()
    assert queries == list(dc.DEFAULT_QUERIES) == ["P32323", "Q07730", "A0ACE4KFX7"]


def test_freeze_values_are_the_documented_ones():
    doc = dc.__doc__
    assert dc.E_SUPPORT == 1e-5 and dc.COVER_MIN == 0.5 and dc.NULL_RATIO_MAX == 0.5
    assert (dc.MIN_SEQUENCES, dc.MIN_CLADES) == (5, 2)
    assert "PROPOSAL" in doc and "FREEZE" in doc
