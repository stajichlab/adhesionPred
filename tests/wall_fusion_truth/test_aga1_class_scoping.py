"""Tests of the AGA1-class scoping study tables (docs/reports/data/aga1_class_scoping/).

Checks: the tables exist; no sequence is committed; counts recompute from the committed truth and control files with an
independent re-implementation of the pre-declared rules; a few baseline numbers recompute from features.tsv (stdlib).
"""

import csv
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "docs/reports/data/aga1_class_scoping"
RCT = ROOT / "analysis/calibration_truth/repeat_call_truth"
CTRL = ROOT / "data/controls"

TABLES = [
    "proteins.tsv", "features.tsv", "clusters.tsv", "source_counts.tsv", "sequence_availability.tsv",
    "class_clade_counts.tsv", "lectin_overlap.tsv", "kiss_adhesion_rows.tsv", "baseline_lco.tsv", "baseline_lcl.tsv",
    "single_feature_auc.tsv", "repeat_call_vs_features.tsv", "matched_operating_point.tsv", "oof_scores.tsv",
    "repeat_rule_reproduction.tsv", "aga1_vs_named.tsv", "narrow_slices.tsv",
]  # fmt: skip


def tsv(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]  # noqa: B905


def test_tables_scripts_readme_report_exist():
    for t in TABLES:
        assert (D / t).is_file(), t
    assert (D / "README.md").is_file()
    assert (ROOT / "docs/reports/2026-10-24-aga1-class-scoping.md").is_file()
    assert len(list((D / "scripts").glob("*.py"))) >= 9


def test_no_sequence_committed():
    # a protein sequence would be a run of 25 or more amino-acid capitals in one field or line
    rx = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{25,}")
    files = list(D.glob("*.tsv")) + [
        D / "README.md",
        ROOT / "docs/reports/2026-10-24-aga1-class-scoping.md",
    ]
    for f in files:
        for i, ln in enumerate(f.read_text().splitlines(), 1):
            assert not ln.startswith(">"), f"{f.name}:{i} looks like FASTA"
            assert not rx.search(ln), f"{f.name}:{i} contains a sequence-like run"
    assert not list(D.glob("*.fa*"))
    for h in tsv(D / "features.tsv")[:1]:
        assert "sequence" not in h and "rep_unit" not in h


def test_row_counts_agree():
    P = tsv(D / "proteins.tsv")
    assert (
        len(P)
        == len({r["accession"] for r in P})
        == len(tsv(D / "features.tsv"))
        == len(tsv(D / "clusters.tsv"))
    )


def _truth():
    out = {}
    for prot in ("Scer_S288C", "Calb_SC5314"):
        for r in tsv(RCT / f"truth_v3.{prot}.annotated.tsv"):
            out[r["accession"]] = r
    return out


def test_roles_recompute_from_committed_files():
    P = {r["accession"]: r for r in tsv(D / "proteins.tsv")}
    truth = _truth()
    cons = {
        r["accession"]: r
        for r in tsv(CTRL / "repeat-mechanism/adjudication/curation_table.consensus.tsv")
    }
    cur = {r["accession"]: r for r in tsv(CTRL / "repeat-mechanism/curation_table.tsv")}
    other = {
        r["accession"]: r
        for r in tsv(CTRL / "repeat-adhesins-other-species/curation_table.other_species.tsv")
    }
    # class positives: curated adhesin AND (mechanism 2a OR truth repeat label 1)
    pos = set()
    for a, r in cur.items():
        if r["cls"] == "adhesin" and (
            cons.get(a, {}).get("final_label") == "2a" or truth.get(a, {}).get("label") == "1"
        ):
            pos.add(a)
    for a, r in truth.items():
        if r["curated_class"] == "adhesin" and (
            r["label"] == "1" or cons.get(a, {}).get("final_label") == "2a"
        ):
            pos.add(a)
    for a, r in other.items():
        if r["cls"] == "adhesin" and r["mechanism_label"] == "2a":
            pos.add(a)
    assert {a for a, r in P.items() if r["role"] == "class_positive"} == pos
    assert len(pos) == 67
    tierA = {a for a in pos if P[a]["pos_tier"] == "A"}
    assert len(tierA) == 35
    # hard negatives: union of the four sources
    hn = {r["accession"] for r in tsv(CTRL / "hard-negatives-shared/controls.tsv")}
    hn |= {
        r["accession"]
        for r in tsv(CTRL / "repeat-adhesins-other-species/hard_negatives_other_species.tsv")
    }
    hn |= {a for a, r in cur.items() if r["cls"] == "hard_negative"}
    hn |= {a for a, r in truth.items() if r["curated_class"] == "hard_negative"}
    assert {a for a, r in P.items() if r["role"] == "hard_negative"} == hn
    bg = {a for a, r in truth.items() if r["curated_class"] == "uniprot" and r["label"] == "0"}
    assert {a for a, r in P.items() if r["role"] == "background_assumed_negative"} == bg


def test_clade_counts_match_proteins():
    P = tsv(D / "proteins.tsv")
    rows = tsv(D / "class_clade_counts.tsv")
    got = {(r["subset"], r["clade_fine"]): int(r["n_proteins"]) for r in rows}
    c = Counter(
        r["clade_fine"] for r in P if r["role"] == "class_positive" and r["pos_tier"] == "A"
    )
    for clade, n in c.items():
        assert got[("pos_tierA", clade)] == n
    assert got[("pos_tierA", "ALL")] == sum(c.values()) == 35
    assert got[("pos_all", "ALL")] == 67
    assert got[("hn_all", "ALL")] == sum(1 for r in P if r["role"] == "hard_negative")


def _auc(pos, neg):
    s = 0.0
    for p in pos:
        for n in neg:
            s += 1.0 if p > n else 0.5 if p == n else 0.0
    return s / (len(pos) * len(neg))


def test_single_feature_auc_and_repeat_call_recompute():
    P = {r["accession"]: r for r in tsv(D / "proteins.tsv")}
    F = {r["accession"]: r for r in tsv(D / "features.tsv")}
    pos = [a for a, r in P.items() if r["role"] == "class_positive" and r["pos_tier"] == "A"]
    hn = [a for a, r in P.items() if r["role"] == "hard_negative"]
    got = {
        (r["feature"], r["pos_set"], r["neg_set"]): float(r["auc_raw"])
        for r in tsv(D / "single_feature_auc.tsv")
    }
    a = _auc([float(F[x]["frac_STP"]) for x in pos], [float(F[x]["frac_STP"]) for x in hn])
    assert abs(a - got[("frac_STP", "tierA", "hn_all")]) < 1e-3
    called = sum(int(F[x]["call_repeat"]) for x in pos)
    row = [
        r
        for r in tsv(D / "repeat_call_vs_features.tsv")
        if r["rule"] == "call_repeat" and r["set"] == "tierA"
    ][0]
    assert int(row["called"]) == called and int(row["n"]) == len(pos)


def test_repeat_call_reproduces_committed_status():
    # status numbers of docs/class_descriptions.yaml: S288C 0.435 / 0.987, C. albicans 0.462 / 0.976
    rows = {(r["truth"], r["z_seq_min"]): r for r in tsv(D / "repeat_rule_reproduction.tsv")}
    s = rows[("Scer_S288C", "none")]
    c = rows[("Calb_SC5314", "none")]
    assert (float(s["sensitivity"]), float(s["specificity"])) == (0.435, 0.987)
    assert (float(c["sensitivity"]), float(c["specificity"])) == (0.462, 0.976)
    # the number of truth_v3 label-1 proteins recomputes from the committed truth table
    t = tsv(RCT / "truth_v3.Scer_S288C.annotated.tsv")
    assert int(s["positives"]) == sum(r["label"] == "1" for r in t)
    assert "AGA1" in rows[("Scer_S288C", "3.5")]["gained_positives"].split(";")


def test_narrow_slices_rule_columns_are_consistent():
    for r in tsv(D / "narrow_slices.tsv"):
        clades = [x.split(":") for x in r["clades_all"].split("; ")]
        strict = sum(1 for _, n in clades if int(n) >= 5) >= 2
        weak = sum(int(n) for _, n in clades) >= 5 and len(clades) >= 2
        assert (r["rule_strict_all"] == "yes") == strict
        assert (r["rule_weak_all"] == "yes") == weak
