"""Tests of candidate_tiers.py and fetch_uniprot_annotations.py (SYNTHETIC rows) and of the committed tier record."""

import csv
import hashlib
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "analysis/aga1_class"


def load(name):
    spec = importlib.util.spec_from_file_location(name, D / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ct = load("candidate_tiers")
fu = load("fetch_uniprot_annotations")
REC = ROOT / "docs/reports/data/aga1_candidate_tiers_2026-10-26"


def row(acc, species="S1", sp=1, broad=0, narrow=0, models=0, known="", score=0.5):
    return {
        "accession": acc, "species": species, "gene": "", "entry_name": "", "protein_name": "", "length": 100, "sp_yes": sp,
        "gpi_yes": 1, "disorder_longest_run": 10, "rep_coverage": 0.0, "rep_n_copies": 0.0, "broad_flag": broad,
        "narrow_flag": narrow, "n_models_call": models, "score_logreg": score, "score_rf": score, "score_gboost": score,
        "in_850_role": known, "in_850_tier": "",
    }  # fmt: skip


def tiers_of(rows, **kw):
    out = ct.assign(pd.DataFrame(rows), **kw)
    return dict(zip(out.accession, out.tier, strict=True))


def test_tier_rules():
    rows = [
        row("N1", narrow=1), row("B2", broad=1, models=2), row("B1", broad=1, models=1), row("M2", models=2),
        row("M2nosp", models=2, sp=0), row("none", models=1), row("narrow_no_models", narrow=1, models=0),
    ]  # fmt: skip
    t = tiers_of(rows)
    assert (
        t["N1"] == "1" and t["B2"] == "1" and t["narrow_no_models"] == "1"
    )  # narrow, or broad with 2 of 3 models
    assert t["B1"] == "2" and t["M2"] == "2"  # broad alone; 2 of 3 models with a signal peptide
    assert "M2nosp" not in t or t["M2nosp"] == "3"  # no signal peptide: never tier 2
    assert t.get("none") in (None, "3")


def test_known_850_proteins_keep_a_tier_but_are_not_review_eligible_or_sampled():
    out = ct.assign(pd.DataFrame([row("K", narrow=1, known="class_positive"), row("N", narrow=1)]))
    k, n = out.set_index("accession").loc["K"], out.set_index("accession").loc["N"]
    assert k.tier == "1" and k.review_eligible == 0 and k["sample"] == ""
    assert n.review_eligible == 1 and n["sample"] == "all tier 1"


def test_cap_cuts_tier_2_by_rank_and_never_tier_1():
    rows = [
        row(f"A{i}", broad=1, score=i / 100) for i in range(1, 36)
    ]  # 35 tier-2 candidates in one species
    rows += [row("T1", narrow=1)]
    t = tiers_of(rows)
    assert t["T1"] == "1"
    assert sum(v == "2" for v in t.values()) == 29  # cap 30 minus the one tier 1 protein
    assert (
        t["A1"] == "3" and t["A35"] == "2"
    )  # the lowest-ranked tier 2 proteins overflow into tier 3
    out = ct.assign(pd.DataFrame(rows)).set_index("accession")
    assert "overflow" in out.loc["A1", "tier_reason"]


def test_tier_3_is_the_top_25_with_a_signal_peptide_and_no_flag_per_species():
    rows = [row(f"X{i}", score=i / 100) for i in range(1, 41)] + [row("nosp", sp=0, score=0.99)]
    out = ct.assign(pd.DataFrame(rows))
    t3 = set(out[out.tier == "3"].accession)
    assert (
        len(t3) == 25 and "nosp" not in t3 and "X40" in t3 and "X15" not in t3
    )  # ranks within the species


def test_rank_is_within_species_not_across_the_run():
    rows = [row(f"A{i}", species="S1", score=0.9) for i in range(30)] + [
        row(f"B{i}", species="S2", score=0.1 + i / 100) for i in range(30)
    ]
    out = ct.assign(pd.DataFrame(rows))
    assert set(out[out.species == "S2"].accession) >= {
        "B29",
        "B28",
    }  # S2's best are tier 3 although all its scores are low


def test_sample_is_seeded_spread_over_species_and_all_when_small():
    rows = [row(f"A{i}", species="S1", broad=1) for i in range(15)] + [
        row(f"B{i}", species="S2", broad=1) for i in range(15)
    ]
    one = ct.assign(pd.DataFrame(rows))
    again = ct.assign(pd.DataFrame(rows))
    pd.testing.assert_frame_equal(one, again)  # deterministic
    s = one[(one.tier == "2") & (one["sample"] != "")]
    assert len(s) == 25
    assert abs(int((s.species == "S1").sum()) - int((s.species == "S2").sum())) <= 1
    other = ct.assign(pd.DataFrame(rows), seed=1)
    assert set(other[other["sample"] != ""].accession) != set(s.accession)  # the seed matters
    small = ct.assign(pd.DataFrame([row(f"C{i}", broad=1) for i in range(5)]))
    assert (small["sample"] != "").all()  # fewer than the quota: all of them


def test_assign_failures():
    with pytest.raises(ct.TierError, match="lack columns"):
        ct.assign(pd.DataFrame([{"accession": "A"}]))
    with pytest.raises(ct.TierError, match="more than one"):
        ct.assign(pd.DataFrame([row("A"), row("A")]))


def test_annotate_requires_every_candidate():
    tiers = ct.assign(pd.DataFrame([row("A", narrow=1), row("B", narrow=1)]))
    ann = pd.DataFrame({"accession": ["A"], "pfam": ["PF1"], "interpro": [""], "domains": [""], "keywords": [""],
                        "subcellular_location": [""], "reviewed": ["unreviewed"], "uniprot_release": ["2026_03"]})  # fmt: skip
    with pytest.raises(ct.TierError, match="lack 1 candidates"):
        ct.annotate(tiers, ann)


UNIPROT_TEXT = (
    "Entry\tGene Names\tProtein names\tLength\tPfam\tInterPro\tDomain [FT]\tKeywords\tSubcellular location [CC]\tReviewed\n"
    'P1\tG1 X\tPA14 domain protein\t958\tPF10528;\tIPR018871;IPR050972;\tDOMAIN 87..249; /note="PA14"; /evidence="E"; DOMAIN 300..400; /note="GH18"\tSignal;Cell wall\t'
    "SUBCELLULAR LOCATION: Secreted, cell wall {ECO:0000256|ARBA:ARBA00004191}.\tunreviewed\n"
)


def test_uniprot_parsing_and_missing_accession(monkeypatch):
    r = fu.parse_tsv(UNIPROT_TEXT, "2026_03")["P1"]
    assert (
        r["pfam"] == "PF10528"
        and r["interpro"] == "IPR018871;IPR050972"
        and r["domains"] == "PA14;GH18"
    )
    assert r["subcellular_location"] == "Secreted, cell wall." and r["uniprot_release"] == "2026_03"
    monkeypatch.setattr(fu, "fetch_batch", lambda batch: (UNIPROT_TEXT, "2026_03"))
    assert fu.fetch(["P1"])[0]["accession"] == "P1"
    with pytest.raises(fu.FetchError, match="did not return 1"):
        fu.fetch(["P1", "P2"])
    flip = iter(["2026_03", "2026_04"])
    monkeypatch.setattr(
        fu, "fetch_batch", lambda batch: (UNIPROT_TEXT.replace("P1", batch[0]), next(flip))
    )
    with pytest.raises(fu.FetchError, match="several UniProt releases"):
        fu.fetch(["P1", "P3"], batch_size=1)


def test_cli_refuses_existing_output(tmp_path):
    f = tmp_path / "t.tsv"
    f.write_text("x\n")
    acc = tmp_path / "a.txt"
    acc.write_text("P1\n")
    r = subprocess.run(
        [
            sys.executable,
            "-I",
            str(D / "fetch_uniprot_annotations.py"),
            "--accessions",
            str(acc),
            "--out",
            str(f),
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 1 and r.stderr.startswith("FATAL:")


def test_committed_record_hashes_and_counts():
    text = (REC / "README.md").read_text()
    found = {p: s for s, p in re.findall(r"^([0-9a-f]{64})  (\S+)$", text, flags=re.M)}
    committed = {p.name for p in REC.iterdir() if p.is_file() and p.name != "README.md"}
    assert committed == set(found)
    for p in committed:
        assert hashlib.sha256((REC / p).read_bytes()).hexdigest() == found[p], p
    with (REC / "candidates_tiered.tsv").open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    n = lambda t, k=None: sum(r["tier"] == t and (k is None or r[k] == "1") for r in rows)  # noqa: E731
    assert (n("1"), n("1", "review_eligible"), n("2"), n("2", "review_eligible"), n("3")) == (
        101,
        57,
        109,
        83,
        257,
    )
    assert sum(r["sample"] != "" for r in rows) == 97
    assert {r["uniprot_release"] for r in rows} == {"2026_03"}
    with (REC / "hand_check_sample.tsv").open(newline="") as fh:
        assert [r["accession"] for r in csv.DictReader(fh, delimiter="\t")] == [
            r["accession"] for r in rows if r["sample"]
        ]
