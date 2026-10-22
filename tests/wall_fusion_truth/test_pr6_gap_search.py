"""PR 6 gap search drafts: pr6_gap_summary.py (fixtures, hand counts) and pr6_gap_search.sbatch (static checks)."""

import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "analysis/wall_fusion_truth"
SBATCH = WF / "pr6_gap_search.sbatch"
SCAN = ROOT / "docs/reports/data/pr6_pfam_scan_2026-10-22"


def load(name):
    spec = importlib.util.spec_from_file_location(name, WF / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


gap = load("pr6_gap_summary")

# (taxid, parent, rank, scientific name)
TAXA = [
    (1, 1, "no rank", "root"),
    (131567, 1, "no rank", "cellular organisms"),
    (2759, 131567, "domain", "Eukaryota"),
    (4751, 2759, "kingdom", "Fungi"),
    (99999, 4751, "species", "Fungus sp. unplaced"),
    (4890, 4751, "phylum", "Ascomycota"),
    (147537, 4890, "subphylum", "Saccharomycotina"),
    (4891, 147537, "class", "Saccharomycetes"),
    (4932, 4891, "species", "Saccharomyces cerevisiae"),
    (559292, 4932, "strain", "Saccharomyces cerevisiae S288C"),
    (5476, 4891, "species", "Candida albicans"),
    (147538, 4890, "subphylum", "Pezizomycotina"),
    (147550, 147538, "class", "Sordariomycetes"),
    (5141, 147550, "species", "Neurospora crassa"),
    (5140, 147550, "species", "Neurospora tetrasperma"),
    (5144, 147550, "species", "Fusarium graminearum"),
    (147545, 147538, "class", "Eurotiomycetes"),
    (162425, 147545, "species", "Aspergillus nidulans"),
    (330879, 147545, "species", "Aspergillus fumigatus"),
    (2, 131567, "domain", "Bacteria"),
    (1224, 2, "phylum", "Pseudomonadota"),
    (562, 1224, "species", "Escherichia coli"),
]


def write_dmp(tmp_path):
    nodes = tmp_path / "nodes.dmp"
    names = tmp_path / "names.dmp"
    nodes.write_text("".join(f"{t}\t|\t{p}\t|\t{r}\t|\t\n" for t, p, r, _ in TAXA))
    names.write_text(
        "".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, _, _, n in TAXA)
        + "5141\t|\tOrange bread mold\t|\t\t|\tcommon name\t|\n"
    )
    return nodes, names


def tbl_line(target, evalue, score, desc):
    cols = [
        target,
        "-",
        "query",
        "-",
        str(evalue),
        str(score),
        "0.1",
        str(evalue),
        str(score),
        "0.1",
    ]
    cols += ["1.0", "1", "0", "0", "1", "1", "1", "1", desc]
    return " ".join(cols)


def write_tbl(path, rows, trailer=True):
    lines = ["#" + " target name accession query name"]
    lines += [tbl_line(*r) for r in rows]
    if trailer:
        lines += ["# Program:         jackhmmer", "# [ok]"]
    path.write_text("\n".join(lines) + "\n")


def sp(acc, name="X_Y"):
    return f"sp|{acc}|{name}"


def desc(taxid):
    return f"Some protein OS=Some species OX={taxid} GN=g PE=1 SV=1"


@pytest.fixture
def tax(tmp_path):
    nodes, names = write_dmp(tmp_path)
    return gap.Taxonomy(
        nodes, names, {559292, 5476, 5141, 5140, 5144, 162425, 330879, 562, 4932, 99999}
    )


# --- clades -------------------------------------------------------------------------------------------------


def test_clade_rule_equals_the_one_of_map_kiss():
    map_kiss = load("map_kiss")
    cases = [
        {
            "phylum": ("Ascomycota", 1),
            "subphylum": ("Saccharomycotina", 2),
            "class": ("Saccharomycetes", 3),
        },
        {
            "phylum": ("Ascomycota", 1),
            "subphylum": ("Pezizomycotina", 2),
            "class": ("Sordariomycetes", 3),
        },
        {"phylum": ("Basidiomycota", 1), "class": ("Agaricomycetes", 3)},
        {"phylum": ("Chytridiomycota", 1)},
        {},
    ]
    for lin in cases:
        assert gap.clades_from_lineage(lin) == map_kiss.clades_from_lineage(lin)


def test_taxonomy_info_gives_clade_fungus_and_species(tax):
    s288c = tax.info(559292)
    assert (
        s288c["fungus"]
        and s288c["clade_fine"] == "Saccharomycotina"
        and s288c["species_taxid"] == 4932
    )
    assert tax.info(5141)["clade_fine"] == "Sordariomycetes"
    assert tax.info(162425)["clade_fine"] == "Eurotiomycetes"
    ecoli = tax.info(562)
    assert not ecoli["fungus"] and ecoli["clade_fine"] == "Pseudomonadota"
    assert tax.info(999999)["clade_fine"] == ""  # unknown taxon: no clade, not a fungus


# --- hits and the training rule -----------------------------------------------------------------------------


def hits_fixture():
    return [
        {
            "target": sp("Q07730"),
            "accession": "Q07730",
            "evalue": 1e-50,
            "score": 90.0,
            "taxid": 5476,
        },  # self
        {
            "target": sp("A1"),
            "accession": "A1",
            "evalue": 1e-20,
            "score": 50.0,
            "taxid": 559292,
        },  # inside
        {"target": sp("B1"), "accession": "B1", "evalue": 1e-9, "score": 30.0, "taxid": 5141},
        {"target": sp("B2"), "accession": "B2", "evalue": 1e-8, "score": 29.0, "taxid": 5140},
        {"target": sp("B3"), "accession": "B3", "evalue": 1e-7, "score": 28.0, "taxid": 5144},
        {"target": sp("C1"), "accession": "C1", "evalue": 1e-7, "score": 28.0, "taxid": 162425},
        {"target": sp("C2"), "accession": "C2", "evalue": 1e-6, "score": 27.0, "taxid": 330879},
        {
            "target": sp("E1"),
            "accession": "E1",
            "evalue": 1e-6,
            "score": 27.0,
            "taxid": 562,
        },  # bacterial
        {
            "target": sp("N1"),
            "accession": "N1",
            "evalue": 1e-6,
            "score": 27.0,
            "taxid": None,
        },  # no OX
        {
            "target": sp("B1"),
            "accession": "B1",
            "evalue": 1e-3,
            "score": 20.0,
            "taxid": 5141,
        },  # duplicate of B1
    ]


def test_classify_hits_and_the_rule_hand_counted(tax):
    rows = gap.classify_hits("Q07730", "Saccharomycotina", hits_fixture(), tax)
    where = {r["accession"]: r["where"] for r in rows}
    assert where == {
        "A1": "inside_clade",
        "B1": "outside_clade",
        "B2": "outside_clade",
        "B3": "outside_clade",
        "C1": "outside_clade",
        "C2": "outside_clade",
        "E1": "non_fungal",
        "N1": "no_taxid",
        "Q07730": "self",
    }
    # the duplicate keeps its best E-value
    assert next(r for r in rows if r["accession"] == "B1")["evalue"] == 1e-9
    rule = gap.training_rule(rows)
    # 5 outside sequences in 2 clades (Sordariomycetes 3, Eurotiomycetes 2): the rule is met; 3 species in one clade
    assert (
        rule["n_outside_seq"] == 5
        and rule["n_outside_clades"] == 2
        and rule["n_outside_species"] == 5
    )
    assert rule["rule_met"] and rule["reason"] == ""
    assert rule["per_clade"] == {"Eurotiomycetes": 2, "Sordariomycetes": 3}
    assert not rule["strict_per_clade_met"]  # 5 in each of 2 clades is not met


def test_rule_not_met_cases(tax):
    base = [{"target": sp("Q"), "accession": "Q", "evalue": 1e-50, "score": 90.0, "taxid": 5476}]

    def mk(accs_taxids):
        return base + [
            {"target": sp(a), "accession": a, "evalue": 1e-6, "score": 27.0, "taxid": t}
            for a, t in accs_taxids
        ]

    none = gap.training_rule(gap.classify_hits("Q", "Saccharomycotina", base, tax))
    assert not none["rule_met"] and none["reason"] == "no outside fungal homolog"
    four = gap.training_rule(
        gap.classify_hits(
            "Q",
            "Saccharomycotina",
            mk([("a", 5141), ("b", 5140), ("c", 162425), ("d", 330879)]),
            tax,
        )
    )
    assert not four["rule_met"] and "only 4 outside sequences" in four["reason"]
    one_clade = gap.training_rule(
        gap.classify_hits(
            "Q",
            "Saccharomycotina",
            mk([(f"s{i}", t) for i, t in enumerate([5141, 5140, 5144, 5141, 5140])]),
            tax,
        )
    )
    assert one_clade["n_outside_seq"] == 5 and not one_clade["rule_met"]
    assert "1 named clade" in one_clade["reason"]
    strict = gap.training_rule(
        gap.classify_hits(
            "Q",
            "Saccharomycotina",
            mk([(f"a{i}", 5141) for i in range(5)] + [(f"b{i}", 162425) for i in range(5)]),
            tax,
        )
    )
    assert strict["rule_met"] and strict["strict_per_clade_met"]


def test_target_status_flags_iteration_only():
    met = {"rule_met": True, "reason": "", "n_outside_seq": 6}
    no = {"rule_met": False, "reason": "only 2 outside sequences (minimum 5)", "n_outside_seq": 2}
    assert gap.target_status({"phmmer": met, "jackhmmer": met})[0] == "buildable"
    status, why = gap.target_status({"phmmer": no, "jackhmmer": met})
    assert status == "buildable_iterated_only" and "drift" in why
    assert gap.target_status({"phmmer": no, "jackhmmer": no}) == ("not_buildable", no["reason"])


def test_tblout_needs_the_ok_trailer_and_reads_the_taxon(tmp_path):
    ok = tmp_path / "ok.tbl"
    write_tbl(ok, [(sp("P1"), 1e-9, 30.0, desc(5141)), ("A0A000", 1e-5, 20.0, "no taxon here")])
    rows = gap.parse_tblout(ok)
    assert [r["accession"] for r in rows] == ["P1", "A0A000"]
    assert rows[0]["taxid"] == 5141 and rows[1]["taxid"] is None
    bad = tmp_path / "bad.tbl"
    write_tbl(bad, [(sp("P1"), 1e-9, 30.0, desc(5141))], trailer=False)
    with pytest.raises(ValueError, match="trailer"):
        gap.parse_tblout(bad)


# --- the Ecm33 check ----------------------------------------------------------------------------------------

MODEL = "NAME  Ecm33\nACC   PF12454.12\nDESC  x\nLENG  200\nGA    25.00 20.00;\nTC    25.1 20.1;\nHMM   A C D\n"


def dom_line(seq, tlen, model, score, ev, hf, ht, af, at):
    f = [
        seq,
        "-",
        str(tlen),
        model,
        "PF12454.12",
        "200",
        str(ev),
        str(score),
        "0.0",
        "1",
        "1",
        str(ev),
        str(ev),
        str(score),
        "0.0",
        str(hf),
        str(ht),
        str(af),
        str(at),
        str(af),
        str(at),
        "0.9",
        "desc",
    ]
    return " ".join(f)


def test_ecm33_check_scores_against_ga_and_reports_missing_hits(tmp_path):
    model = tmp_path / "Ecm33.hmm"
    model.write_text(MODEL)
    dom = tmp_path / "dom.txt"
    dom.write_text(
        "# hmmsearch :: search profile(s) against a sequence database\n"
        + dom_line("Q5AGC4", 420, "Ecm33", 22.5, 0.002, 10, 150, 30, 180)
        + "\n"
        + dom_line("Q5AGC4", 420, "Ecm33", 8.0, 0.9, 1, 20, 300, 320)
        + "\n# [ok]\n"
    )
    hm = gap.parse_model_header(model)
    assert hm == {"name": "Ecm33", "acc": "PF12454", "leng": 200, "ga_seq": 25.0, "ga_dom": 20.0}
    hits = gap.parse_hmmsearch_domtbl(dom, "Ecm33")
    assert hits["Q5AGC4"]["best_dom_score"] == 22.5 and hits["Q5AGC4"]["hmm_to"] == 150
    rows = gap.ecm33_rows(hm, hits, {"Q5AGC4": 420, "P38248": 406})
    by = {r[0]: r for r in rows}
    assert by["Q5AGC4"][11] == "below_GA"  # sequence score 22.5 < GA 25
    assert by["Q5AGC4"][2] == "1.034" and by["Q5AGC4"][10] == "0.705"  # 420/406; (150-10+1)/200
    assert by["P38248"][11] == "no_hit_at_E10" and by["P38248"][2] == "1.000"
    # a score over both thresholds passes
    hits["Q5AGC4"]["score_seq"], hits["Q5AGC4"]["best_dom_score"] = 30.0, 28.0
    assert {r[0]: r for r in gap.ecm33_rows(hm, hits, {"Q5AGC4": 420})}["Q5AGC4"][11] == "pass_GA"


def test_ecm33_parsers_refuse_a_model_without_ga_and_the_hmmscan_orientation(tmp_path):
    nomodel = tmp_path / "m.hmm"
    nomodel.write_text("NAME  Ecm33\nACC   PF12454.12\nLENG  200\nHMM   A\n")
    with pytest.raises(ValueError, match="GA"):
        gap.parse_model_header(nomodel)
    dom = tmp_path / "scan.txt"
    dom.write_text(dom_line("Ecm33", 200, "Q5AGC4", 30.0, 1e-5, 1, 100, 1, 100) + "\n# [ok]\n")
    with pytest.raises(ValueError, match="orientation"):
        gap.parse_hmmsearch_domtbl(dom, "Ecm33")


# --- end to end on fixtures ----------------------------------------------------------------------------------

TARGETS_HEADER = (
    "slice\tcluster_id\tgroup\tgenes\tclades_fine\tstrata_fine\taccessions\tweak_families\n"
)


def make_run(tmp_path):
    nodes, names = write_dmp(tmp_path)
    tdir = tmp_path / "in"
    tdir.mkdir()
    targets = tdir / "custom_targets.tsv"
    targets.write_text(
        TARGETS_HEADER
        + "adhesion\tKS0107\tA2\tece1\tSaccharomycotina\tclade_only\tQ07730\t\n"
        + "wall\tKS0086\tW5\tecm33\tSaccharomycotina\tclade_only\tQ5AGC4\tPF01030:Recep_L_domain\n"
    )
    ph, jk = tmp_path / "ph", tmp_path / "jk"
    ph.mkdir()
    jk.mkdir()
    one_pass = [(sp("Q07730"), 1e-50, 90.0, desc(5476))]
    iterated = one_pass + [
        (sp(f"B{i}"), 1e-6, 27.0, desc(t)) for i, t in enumerate([5141, 5140, 5144, 162425, 330879])
    ]
    write_tbl(ph / "Q07730.tbl", one_pass)
    write_tbl(jk / "Q07730.tbl", iterated)
    write_tbl(ph / "Q5AGC4.tbl", [(sp("Q5AGC4"), 1e-80, 150.0, desc(5476))] + iterated[1:])
    write_tbl(jk / "Q5AGC4.tbl", [(sp("Q5AGC4"), 1e-80, 150.0, desc(5476))] + iterated[1:])
    return targets, ph, jk, nodes, names


def run_main(tmp_path, out, extra=()):
    targets, ph, jk, nodes, names = (
        make_run(tmp_path)
        if not (tmp_path / "in").exists()
        else (
            tmp_path / "in" / "custom_targets.tsv",
            tmp_path / "ph",
            tmp_path / "jk",
            tmp_path / "nodes.dmp",
            tmp_path / "names.dmp",
        )
    )
    return gap.main(
        [
            "--targets",
            str(targets),
            "--phmmer-dir",
            str(ph),
            "--jackhmmer-dir",
            str(jk),
            "--nodes",
            str(nodes),
            "--names",
            str(names),
            "--out",
            str(out),
            *extra,
        ]
    )


def test_main_end_to_end_statuses_and_no_sequences(tmp_path, capsys):
    out = tmp_path / "out"
    assert run_main(tmp_path, out) == 0
    assert sorted(p.name for p in out.iterdir()) == sorted(gap.OUTPUTS)
    rows = [r.split("\t") for r in (out / "gap_per_query.tsv").read_text().splitlines()]
    head = rows[0]
    body = {r[3]: dict(zip(head, r, strict=True)) for r in rows[1:]}
    # ece1: one pass finds nothing outside, jackhmmer finds 5 in 2 clades: iteration only
    assert body["Q07730"]["status"] == "buildable_iterated_only"
    assert (
        body["Q07730"]["phmmer_outside_seq"] == "0"
        and body["Q07730"]["jackhmmer_outside_seq"] == "5"
    )
    assert body["Q07730"]["mode"] == "leave-cluster-out (partial)"
    # ecm33 (wall): both searches find the 5 outside sequences
    assert body["Q5AGC4"]["status"] == "buildable" and body["Q5AGC4"]["mode"] == "leave-clade-out"
    # no sequence lines anywhere in the outputs
    for p in out.iterdir():
        assert not re.search(r"^[ACDEFGHIKLMNPQRSTVWY]{30,}$", p.read_text(), re.M), p.name
    capsys.readouterr()


def test_main_is_deterministic_and_refuses_a_filled_directory(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    run_main(tmp_path, a)
    run_main(tmp_path, b)
    for name in gap.OUTPUTS:
        assert (a / name).read_bytes() == (b / name).read_bytes(), name
    with pytest.raises(SystemExit, match="refusing"):
        run_main(tmp_path, a)


def test_main_stops_when_a_query_search_is_missing(tmp_path):
    make_run(tmp_path)
    (tmp_path / "jk" / "Q5AGC4.tbl").unlink()
    with pytest.raises(SystemExit, match="missing"):
        run_main(tmp_path, tmp_path / "out")


# --- UniRef50 headers (round 2) ---------------------------------------------------------------------------


def uniref(acc, taxid, n=1, rep=None, species="Some species"):
    """(target, description) of a UniRef50 hit: the cluster id and the rest of the header line."""
    return f"UniRef50_{acc}", f"some protein n={n} Tax={species} TaxID={taxid} RepID={rep or acc}"


def test_tblout_reads_taxid_from_uniref_and_still_from_uniprot(tmp_path):
    t1, d1 = uniref("B1", 5141, n=7)
    path = tmp_path / "x.tbl"
    write_tbl(path, [(t1, 1e-9, 30.0, d1), (sp("U1"), 1e-9, 30.0, desc(5140))])
    rows = {r["target"]: r for r in gap.parse_tblout(path)}
    assert (
        rows[t1]["taxid"] == 5141 and rows[t1]["rep"] == "B1" and rows[t1]["taxid_conflict"] is None
    )
    assert rows[sp("U1")]["taxid"] == 5140 and rows[sp("U1")]["rep"] is None


def test_uniref_self_is_the_cluster_of_the_query_or_its_repid(tmp_path, tax):
    def classify(target, d):
        path = tmp_path / "t.tbl"
        write_tbl(path, [(target, 1e-9, 30.0, d)])
        rows = gap.classify_hits("Q07730", "Saccharomycotina", gap.parse_tblout(path), tax)
        return [(r["accession"], r["where"]) for r in rows]

    # the cluster of the query, and a cluster whose RepID names the query
    assert classify(*uniref("Q07730", 5476)) == [("UniRef50_Q07730", "self")]
    assert classify("UniRef50_B9WG30", "x n=2 Tax=C TaxID=5476 RepID=Q07730") == [
        ("UniRef50_B9WG30", "self")
    ]
    # a representative that is another protein is a normal hit, classified by its own taxon
    assert classify(*uniref("B9WG30", 5476)) == [("UniRef50_B9WG30", "inside_clade")]


def test_uniref_member_counts_are_not_used_by_the_rule(tmp_path, tax):
    path = tmp_path / "t.tbl"
    taxa = [5141, 5140, 5144, 162425, 330879]
    ns = [1, 1, 1, 1, 5000]  # a huge cluster counts once
    hits = [uniref(f"R{i}", t, n=n) for i, (t, n) in enumerate(zip(taxa, ns, strict=True))]
    write_tbl(path, [(target, 1e-8, 30.0, d) for target, d in hits])
    rows = gap.classify_hits("Q0", "Saccharomycotina", gap.parse_tblout(path), tax)
    rule = gap.training_rule(rows)
    assert rule["n_outside_seq"] == 5 and rule["n_outside_clades"] == 2 and rule["rule_met"]


def test_ox_and_taxid_that_disagree_give_no_taxon_and_are_reported_loudly(tmp_path, tax, capsys):
    path = tmp_path / "t.tbl"
    write_tbl(
        path,
        [
            ("UniRef50_X1", 1e-9, 30.0, "p n=1 Tax=S OX=5141 TaxID=5140 RepID=X1"),
            ("UniRef50_X2", 1e-9, 30.0, "p n=1 Tax=S OX=5141 TaxID=5141 RepID=X2"),
        ],
    )
    rows = {r["target"]: r for r in gap.parse_tblout(path)}
    assert rows["UniRef50_X1"]["taxid"] is None and rows["UniRef50_X1"]["taxid_conflict"] == (
        5141,
        5140,
    )
    assert rows["UniRef50_X2"]["taxid"] == 5141 and rows["UniRef50_X2"]["taxid_conflict"] is None
    classified = gap.classify_hits("Q0", "Saccharomycotina", list(rows.values()), tax)
    assert {r["accession"]: r["where"] for r in classified}["UniRef50_X1"] == "no_taxid"


def test_main_with_uniref_headers_and_a_taxid_conflict(tmp_path, capsys):
    make_run(tmp_path)
    ph, jk = tmp_path / "ph", tmp_path / "jk"
    own = uniref("Q5AGC4", 5476)
    outside = [
        uniref(f"R{i}", t, n=n)
        for i, (t, n) in enumerate(
            zip([5141, 5140, 5144, 162425, 330879], [1, 1, 3, 1, 900], strict=True)
        )
    ]
    rows = [(own[0], 1e-80, 150.0, own[1])] + [(t, 1e-8, 30.0, d) for t, d in outside]
    write_tbl(ph / "Q5AGC4.tbl", rows)
    write_tbl(
        jk / "Q5AGC4.tbl",
        rows + [("UniRef50_BAD", 1e-8, 30.0, "p n=1 Tax=S OX=5141 TaxID=5140 RepID=BAD")],
    )
    out = tmp_path / "out"
    assert run_main(tmp_path, out) == 0
    cap = capsys.readouterr().out
    assert "REPORT LOUD: 1 hit(s) carry OX= and TaxID= that disagree" in cap
    conflicts = (out / "taxid_conflicts.tsv").read_text().splitlines()
    assert conflicts[0] == "query\tmethod\thit\tOX\tTaxID" and conflicts[1].split("\t") == [
        "Q5AGC4",
        "jackhmmer",
        "UniRef50_BAD",
        "5141",
        "5140",
    ]
    body = [r.split("\t") for r in (out / "gap_per_query.tsv").read_text().splitlines()]
    row = dict(zip(body[0], next(r for r in body[1:] if r[3] == "Q5AGC4"), strict=True))
    assert (
        row["status"] == "buildable"
        and row["phmmer_outside_seq"] == "5"
        and row["jackhmmer_outside_seq"] == "5"
    )


def test_sbatch_is_ready_for_uniref50_and_a_longer_run():
    t = text()
    body = "\n".join(ln for ln in t.splitlines() if not ln.lstrip().startswith("#"))
    # partition, time, cores and memory are only the #SBATCH defaults; the body fixes none of them
    assert not re.search(r"--time|timeout|-p short|--mem", body)
    assert 'THREADS="${SLURM_CPUS_PER_TASK:-8}"' in t
    assert "overrides them on the sbatch command line" in t
    # the target set (28.8 GB for UniRef50) is counted once and the count is reused
    assert (
        len(re.findall(r"grep -c '\^>' \"\$UNIPROT_FASTA\"", body)) == 1 and "$N_TARGET_SEQ" in body
    )
    # raw tables are kept if the summary step fails after hours of searching
    assert "raw_after_summary_failure" in t
    # the run settings are unchanged
    assert "INC_E=1e-5" in t and "JACK_ITER=3" in t


# --- the sbatch (static checks) ------------------------------------------------------------------------------


def text():
    return SBATCH.read_text()


def test_sbatch_syntax_and_header():
    assert subprocess.run(["bash", "-n", str(SBATCH)]).returncode == 0
    t = text()
    assert "DRAFT, NOT SUBMITTED" in t.splitlines()[1]
    assert "#SBATCH -p short" in t
    assert "set -euo pipefail" in t


def test_sbatch_rules():
    t = text()
    body = "\n".join(ln for ln in t.splitlines() if not ln.lstrip().startswith("#"))
    assert "/tmp" not in body and "${SCRATCH:?" in t
    assert "module load hmmer/3.4" in t
    # no UniProt path is guessed: the set is a required argument
    assert ': "${UNIPROT_FASTA:?' in t
    assert not re.search(r"UNIPROT_FASTA=[^$]", body)
    # the thresholds are fixed in the script, with iterations capped
    assert "INC_E=1e-5" in t and "JACK_ITER=3" in t and "--incE" in t and "--incdomE" in t
    # nothing is overwritten, release is checked, trailers are checked
    assert "nothing is overwritten" in t and "Pfam${PFAM_RELEASE}/" in t and "'# \\[ok\\]'" in t
    # the Ecm33 check uses the fetched model and fails if it is not PF12454
    assert "hmmsearch --max -E 10" in t and "PF12454" in t
    # the job never submits anything itself
    assert not re.search(r"^\s*sbatch ", body, re.M)


def test_sbatch_reads_exactly_the_four_gap_queries_from_the_committed_record():
    rows = [ln.split("\t") for ln in (SCAN / "custom_targets.tsv").read_text().splitlines()]
    accs = sorted(a for r in rows[1:] for a in r[6].split(","))
    assert accs == ["A0ACE4KFX7", "P32323", "Q07730", "Q5AGC4"]
    assert "cut -f7" in text()  # the accessions column of custom_targets.tsv
    assert rows[0][6] == "accessions"


def test_sbatch_extraction_helpers_match_the_exact_header_token(tmp_path):
    t = text()
    one_seq = re.search(r"^one_seq\(\) \{.*\}$", t, re.M).group(0)
    lengths = re.search(r"awk '/\^>/ \{ if \(n\).*?' \\\n\s+\"\$TMP/ecm33.faa\"", t, re.S)
    assert lengths, "the length table awk is part of the script"
    fa = tmp_path / "x.faa"
    fa.write_text(">P38248\nMQFK\nNALT\n>P382480 other\nGGGG\n>Q5AGC4 desc\nAAAAAA\n")
    script = f'{one_seq}\none_seq "{fa}" P38248\n'
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True).stdout
    assert out == ">P38248\nMQFK\nNALT\n"  # not the P382480 record
    assert (
        subprocess.run(
            ["bash", "-c", f'{one_seq}\none_seq "{fa}" Z99999'], capture_output=True, text=True
        ).stdout
        == ""
    )
    awk_len = (
        'awk \'/^>/ { if (n) print n_id "\\t" n; n_id = substr($1, 2); n = 0; next } { n += length($0) } '
        f'END {{ print n_id "\\t" n }}\' {fa}'
    )
    assert subprocess.run(
        ["bash", "-c", awk_len], capture_output=True, text=True, check=True
    ).stdout == ("P38248\t8\nP382480\t4\nQ5AGC4\t6\n")
    assert "n_id = substr($1, 2)" in t  # the helper tested above is the one the script uses


def test_hits_without_a_fine_clade_are_not_a_clade_for_the_rule(tax):
    """The hyd1_like pattern of round 2: 3 hits in one named clade plus 2 fungal hits with no assigned fine clade."""
    assert tax.info(99999)["fungus"] and tax.info(99999)["clade_fine"] == ""
    base = [{"target": sp("Q"), "accession": "Q", "evalue": 1e-50, "score": 90.0, "taxid": 5476}]
    hits = base + [
        {"target": sp(a), "accession": a, "evalue": 1e-6, "score": 27.0, "taxid": t}
        for a, t in [("a", 5141), ("b", 5140), ("c", 5144), ("d", 99999), ("e", 99999)]
    ]
    r = gap.training_rule(gap.classify_hits("Q", "Eurotiomycetes", hits, tax))
    assert r["n_outside_seq"] == 5  # the unassigned fungal hits count in the total
    assert r["n_outside_clades"] == 1  # named clades only
    assert r["per_clade"] == {"(no clade)": 2, "Sordariomycetes": 3}
    assert not r["rule_met"] and not r["strict_per_clade_met"]
    assert "1 named clade" in r["reason"] and "2 fungal hits have no assigned clade" in r["reason"]


def test_a_second_named_clade_still_meets_the_rule_and_unassigned_hits_add_nothing_to_it(tax):
    base = [{"target": sp("Q"), "accession": "Q", "evalue": 1e-50, "score": 90.0, "taxid": 5476}]
    hits = base + [
        {"target": sp(a), "accession": a, "evalue": 1e-6, "score": 27.0, "taxid": t}
        for a, t in [("a", 5141), ("b", 5140), ("c", 162425), ("d", 99999), ("e", 99999)]
    ]
    r = gap.training_rule(gap.classify_hits("Q", "Saccharomycotina", hits, tax))
    assert r["n_outside_clades"] == 2 and r["rule_met"] and not r["strict_per_clade_met"]
    only_unassigned = base + [
        {"target": sp(f"u{i}"), "accession": f"u{i}", "evalue": 1e-6, "score": 27.0, "taxid": 99999}
        for i in range(6)
    ]
    r2 = gap.training_rule(gap.classify_hits("Q", "Saccharomycotina", only_unassigned, tax))
    assert r2["n_outside_seq"] == 6 and r2["n_outside_clades"] == 0 and not r2["rule_met"]
