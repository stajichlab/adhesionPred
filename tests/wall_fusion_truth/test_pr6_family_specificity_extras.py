"""PR 6 second specificity run (Big1, KNR4/SMI1, DFG10): the extras list, the sibling script and its sbatch (static checks)."""

import hashlib
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "analysis/wall_fusion_truth"
SBATCH = WF / "pr6_family_specificity_extras.sbatch"
LIST = WF / "pr6_family_specificity_extras_list.tsv"
COMMITTED = ROOT / "docs/reports/data/pr6_family_extras_candidates_2026-10-24"
FIRST_CAND = ROOT / "docs/reports/data/pr6_family_candidates_2026-10-22"
SCAN = ROOT / "docs/reports/data/pr6_pfam_scan_2026-10-22/scan_per_cluster.tsv"
MEASURE = ROOT / "data/sorting_hat/family_measure.tsv"
TABLE = ROOT / "data/sorting_hat/family_table.tsv"
ACCS = ["PF27678", "PF09346", "PF02544"]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, WF / filename)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


fs = load("pr6_family_specificity", "pr6_family_specificity.py")
ex = load("pr6_family_specificity_extras", "pr6_family_specificity_extras.py")


def tsv(path):
    return fs.read_tsv(path)


def test_the_list_has_the_three_owner_approved_rows_in_order():
    rows = ex.read_extras(LIST)
    assert [(r["order"], r["pfam_acc"], r["pfam_name"], r["cluster_id"]) for r in rows] == [
        ("1", "PF27678", "BIG1_C", "KS0087"),
        ("2", "PF09346", "SMI1_KNR4", "KS0006"),
        ("3", "PF02544", "Steroid_dh", "KS0090"),
    ]
    assert [r["panther_family"] for r in rows] == ["PTHR28285", "PTHR47432", "PTHR14624"]
    assert [r["expected"] for r in rows] == ["specific", "specific", "broad_expected"]
    assert all(r["approved_by"] == "owner decision 2026-10-24" for r in rows)


def test_the_list_agrees_with_the_committed_scan_and_misses_the_pfam_rule():
    rows = ex.read_extras(LIST)
    found = ex.check_against_scan(rows, tsv(SCAN), tsv(fs.MISSED))
    assert {a: (f["hit"]["ievalue"], f["hit"]["frac"]) for a, f in found.items()} == {
        "PF27678": (8.2e-30, 0.170),
        "PF09346": (5.2e-19, 0.323),
        "PF02544": (7.4e-9, 0.419),
    }
    assert not any(
        any(fs.rule_met(f["hit"])) for f in found.values()
    )  # why the first run left them out
    first = {r["pfam_acc"] for r in tsv(FIRST_CAND / "candidate_families.tsv")}
    assert not first & set(ACCS)


def test_a_list_that_disagrees_with_the_scan_is_refused(tmp_path):
    rows = ex.read_extras(LIST)
    for field, value in (
        ("scan_ievalue", "1e-31"),
        ("pfam_name", "Other"),
        ("cluster_id", "KS0001"),
    ):
        bad = [dict(r) for r in rows]
        bad[0][field] = value
        with pytest.raises(ValueError):
            ex.check_against_scan(bad, tsv(SCAN), tsv(fs.MISSED))
    text = LIST.read_text().splitlines()
    swapped = tmp_path / "list.tsv"
    swapped.write_text("\n".join([text[0], text[2], text[1], text[3]]) + "\n")
    with pytest.raises(ValueError):
        ex.read_extras(swapped)  # the order is part of the decision
    dropped = tmp_path / "two.tsv"
    dropped.write_text("\n".join(text[:3]) + "\n")
    with pytest.raises(ValueError):
        ex.read_extras(dropped)


def test_the_committed_extras_tables_are_the_script_output_and_hold_own_member_counts(tmp_path):
    out = tmp_path / "cand"
    rows = ex.write_extras_candidates(out)
    for name in ex.EXTRAS_OUTPUTS:
        assert (out / name).read_bytes() == (COMMITTED / name).read_bytes(), name
    assert {r["pfam_acc"]: len(r["members"]) for r in rows} == {
        "PF27678": 1,
        "PF09346": 2,
        "PF02544": 1,
    }
    own = tsv(out / "family_members.tsv")
    assert [(r["pfam_acc"], r["accession"]) for r in own] == [
        ("PF27678", "P38813"),
        ("PF09346", "A0A1D8PE53"),
        ("PF09346", "P38678"),
        ("PF02544", "P40526"),
    ]
    split = {r["accession"]: r for r in tsv(fs.MEMBERS_SPLIT)}
    assert all(split[r["accession"]]["class"] == "positive" for r in own)
    assert (out / "families.txt").read_text() == "BIG1_C\nSMI1_KNR4\nSteroid_dh\n"


def test_panther_columns_are_descriptive_metadata_only(tmp_path):
    out = tmp_path / "cand"
    ex.write_extras_candidates(out)
    meta = tsv(out / "extras_metadata.tsv")
    assert [(m["panther_family"], int(m["panther_n_fungal"]) >= 1000) for m in meta] == [
        ("PTHR28285", True),
        ("PTHR47432", True),
        ("PTHR14624", True),
    ]
    # the candidate tables read by `summarize` carry no PANTHER field: it cannot filter or label anything
    for name in ("candidate_families.tsv", "candidate_clusters.tsv", "family_members.tsv"):
        assert "PTHR" not in (out / name).read_text()


def test_the_first_run_is_unchanged(tmp_path):
    # the derivation of the 30 families is byte-identical to the committed tables, and the extras never enter it
    out = tmp_path / "first"
    sel = fs.write_candidates(fs.SCAN_RECORD, fs.MISSED, out, fs.FAMILY_TABLE)
    assert len(sel["families"]) == 30 and not set(ACCS) & set(sel["families"])
    for name in fs.CANDIDATE_OUTPUTS:
        assert (out / name).read_bytes() == (FIRST_CAND / name).read_bytes(), name
    assert fs.CANDIDATE_OUTPUTS == (
        "candidate_clusters.tsv",
        "candidate_families.tsv",
        "family_members.tsv",
        "families.txt",
    )
    assert "extras" not in (WF / "pr6_family_specificity.py").read_text().lower()


def test_a_name_missing_from_the_fetched_models_stops_the_job(tmp_path):
    cand = tmp_path / "cand"
    ex.write_extras_candidates(cand)

    def models(entries):
        return "".join(f"HMMER3/f\nNAME  {n}\nACC   {a}.12\nGA    1 1\n//\n" for a, n in entries)

    good = tmp_path / "good.hmm"
    good.write_text(models(zip(ACCS, ["BIG1_C", "SMI1_KNR4", "Steroid_dh"], strict=True)))
    assert ex.check_models(good, cand) == 3
    for label, entries in {
        "missing": zip(ACCS[:2], ["BIG1_C", "SMI1_KNR4"], strict=True),
        "renamed": zip(ACCS, ["BIG1_C", "SMI1_KNR4", "Steroid_dh_2"], strict=True),
        "extra": [*zip(ACCS, ["BIG1_C", "SMI1_KNR4", "Steroid_dh"], strict=True), ("PF00001", "X")],
    }.items():
        bad = tmp_path / f"{label}.hmm"
        bad.write_text(models(entries))
        with pytest.raises(ValueError):
            ex.check_models(bad, cand)
    rc = ex.main(
        ["check-models", "--models", str(tmp_path / "missing.hmm"), "--candidates", str(cand)]
    )
    assert rc == 2  # an unreadable file is an OSError: exit 2 with a message, never a traceback


# ---- summarize on fixtures: own members, tier cap, no committed file written ------------------------


def dom(target, acc, name="Fam"):
    return f"{target} - 300 {name} {acc}.5 100 1e-30 100.0 0.0 1 1 1e-30 1e-30 100.0 0.0 1 100 10 110 10 110 0.99 description\n"


def domtbl(rows):
    return f"# hmmsearch\n# Option settings: hmmsearch --cut_ga x\n#\n{''.join(rows)}# [ok]\n"


@pytest.fixture()
def fx(tmp_path):
    cand = tmp_path / "cand"
    ex.write_extras_candidates(
        cand
    )  # the real committed list: members P38813, A0A1D8PE53, P38678, P40526
    split = tmp_path / "members_split.tsv"
    split.write_text(
        "accession\tcluster_id\tclass\tslice\n"
        "P38813\tKS0087\tpositive\twall\nA0A1D8PE53\tKS0006\tpositive\twall\nP38678\tKS0006\tpositive\twall\n"
        "P40526\tKS0090\tpositive\twall\nX1\tKS9\tpositive\tadhesion\n"
    )
    neg = tmp_path / "neg"
    neg.mkdir()
    (neg / "controls.tsv").write_text("accession\tevidence_level\nN1\tN1\nN2\tN1\nN3\tN2\n")
    (neg / "clusters.tsv").write_text(
        "accession\tcluster_id\tset\nN1\tH1\thard_negative\nN2\tH2\thard_negative\nN3\tH3\thard_negative\n"
    )
    pos = tmp_path / "pos.dom"
    pos.write_text(
        domtbl(
            [
                dom("P38813", "PF27678", "BIG1_C"),
                dom("P38678", "PF09346", "SMI1_KNR4"),
                dom("P40526", "PF02544", "Steroid_dh"),
                dom("X1", "PF02544", "Steroid_dh"),
            ]
        )
    )
    negd = tmp_path / "neg.dom"
    negd.write_text(domtbl([dom("N1", "PF02544"), dom("N2", "PF02544"), dom("N3", "PF02544")]))
    prot = {}
    for name, n in (("P_a", 4), ("P_b", 6)):
        fa = tmp_path / f"{name}.faa"
        fa.write_text("".join(f">{name[2]}{i} annotation\nACDEFG\n" for i in range(n)))
        d = tmp_path / f"{name}.dom"
        d.write_text(domtbl([dom(f"{name[2]}0", "PF02544")]))
        prot[name] = (d, fa)
    return {
        "cand": cand,
        "split": split,
        "neg": neg,
        "pos": pos,
        "negd": negd,
        "prot": prot,
        "tmp": tmp_path,
    }


def run(fx, out):
    return ex.summarize(fx["cand"], fx["split"], fx["neg"], fx["pos"], fx["negd"], fx["prot"], out)


def test_own_member_counts_and_denominators(fx):
    rows = {r["pfam_acc"]: r for r in run(fx, fx["tmp"] / "out")}
    assert [rows[a]["n_own_members"] for a in ACCS] == [1, 2, 1]
    assert [rows[a]["own_hits"] for a in ACCS] == [1, 1, 1]  # A0A1D8PE53 has no hit in the fixture
    assert rows["PF09346"]["own_recall_circular"] == "0.500"
    assert all(rows[a]["leakage"] == "tuned_on_truth" for a in ACCS)
    d = rows[
        "PF02544"
    ]  # a broad family: hits on every negative and on the other positive; reported, not hidden
    assert (d["neg_hits_n1"], d["neg_hits_n2"], d["other_positive_hits_adhesion"]) == (2, 1, 1)
    assert d["kept_hits_total"] == 2


def test_provenance_names_the_list_not_the_rule(fx):
    out = fx["tmp"] / "out"
    run(fx, out)
    prov = {r["item"]: r["value"] for r in tsv(out / "provenance.tsv")}
    assert prov["script"].endswith("pr6_family_specificity_extras.py")
    assert prov["selection"] == ex.SELECTION and "coverage_min" not in prov
    assert prov["extras_list_sha256"] == hashlib.sha256(LIST.read_bytes()).hexdigest()
    assert prov["n_families"] == "3" and prov["tier_changes"] == "none"


def test_summarize_refuses_candidates_that_are_not_the_three_extras(fx):
    other = fx["tmp"] / "first"
    fs.write_candidates(fs.SCAN_RECORD, fs.MISSED, other, fs.FAMILY_TABLE)
    with pytest.raises(ValueError):
        ex.summarize(
            other, fx["split"], fx["neg"], fx["pos"], fx["negd"], fx["prot"], fx["tmp"] / "o2"
        )


def test_tier_is_capped_at_smoke_and_convert_is_a_dry_run(fx):
    out = fx["tmp"] / "out"
    run(fx, out)
    before = (digest(MEASURE), digest(TABLE))
    rc = ex.main(
        [
            "convert",
            "--specificity",
            str(out / "family_specificity.tsv"),
            "--source",
            "s",
            "--date",
            "2026-10-24",
        ]
    )
    assert rc == 0
    rows = fs.convert(out / "family_specificity.tsv", None, "s", "2026-10-24", "38.2", "HMMER 3.4")
    from cellsurface_sorting_hat import family_tier

    assert [r["pfam_acc"] for r in rows] == ACCS
    for r in rows:
        assert r["measure_basis"] == "descriptive_shared_sets"
        assert r["_tier"] in (family_tier.UNTESTED, family_tier.SMOKE)
        assert r["n_pos"] == r["n_pos_unknown"] and r["sensitivity"] == r["specificity"] == ""
    assert (digest(MEASURE), digest(TABLE)) == before
    with pytest.raises(SystemExit):
        ex.main(
            ["convert", "--specificity", "x", "--source", "s", "--date", "d", "--out", "y"]
        )  # no --out exists


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_outputs_hold_no_sequences_and_never_name_the_committed_files(fx):
    out = fx["tmp"] / "out"
    run(fx, out)
    for name in [*fs.SUMMARY_OUTPUTS]:
        for line in (out / name).read_text().splitlines():
            assert not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{20,}", line)
    for line in (COMMITTED / "family_members.tsv").read_text().splitlines():
        assert not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{20,}", line)
    code = [
        ln
        for ln in (WF / "pr6_family_specificity_extras.py").read_text().splitlines()
        if not ln.lstrip().startswith(("#", '"""'))
    ]
    writes = [
        ln for ln in code if re.search(r"write_text|write_tsv|write_measure|open\(.*[\"']w", ln)
    ]
    assert not any(
        "FAMILY_MEASURE" in ln or "FAMILY_TABLE" in ln or "write_measure" in ln for ln in writes
    )


def test_sbatch_follows_the_job_conventions():
    assert subprocess.run(["bash", "-n", str(SBATCH)], capture_output=True).returncode == 0
    text = SBATCH.read_text()
    code = "\n".join(
        ln
        for ln in text.splitlines()
        if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH")
    )
    assert text.startswith("#!/bin/bash -l") and "DRAFT, NOT SUBMITTED" in text.splitlines()[1]
    assert "set -euo pipefail" in code and "#SBATCH -p short" in code and '"${SCRATCH:?' in code
    assert "/tmp" not in "\n".join(ln for ln in code.splitlines() if not ln.startswith("#SBATCH"))
    assert "BASH_SOURCE" not in code and "module load hmmer/3.4" in code
    assert '[ ! -e "$OUT" ]' in code and "check_pin" in code and "sidecar_ok" in code
    assert (
        'Pfam${PFAM_RELEASE}/"' in code
        and "4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e" in code
    )
    assert "PFAM_HMM:-/srv/projects/db/pfam/current/Pfam-A.hmm" in code
    assert 'KEPT="Afum_Af293_UniProt Calb_SC5314 Cimm_RS Scer_S288C"' in code
    searches = [ln.strip() for ln in code.splitlines() if ln.strip().startswith("hmmsearch ")]
    assert len(searches) == 1 and "--cut_ga" in searches[0] and "-E " not in searches[0]
    pin = next(ln for ln in code.splitlines() if "git -C" in ln and "diff --quiet" in ln)
    for part in ('"$WF"', '"$SCAN"', '"$CAND"', '"$NEG"', '"$INPUTS"'):
        assert (
            part in pin
        )  # $WF holds the extras script and list, $CAND the committed extras tables
    assert (
        "hmmfetch -f" in code
        and 'fail "hmmfetch could not fetch' in code
        and "check-models" in code
    )
    assert "pr6_family_specificity_extras.py" in code and "pr6_family_specificity.py" not in code
    executable = "\n".join(
        ln for ln in code.splitlines() if not ln.lstrip().startswith(("echo", "}"))
    )
    assert "family_measure" not in executable
    assert "family_table.tsv" not in executable.replace(
        "FAMILY=data/sorting_hat/family_table.tsv", ""
    )
    assert (
        re.search(r"^N_FAM_EXP=3\b", code, re.M)
        and "N_POS_EXP=132" in code
        and "N_NEG_EXP=161" in code
    )
    assert "CAND=docs/reports/data/pr6_family_extras_candidates_2026-10-24" in code


def test_the_pinned_counts_match_the_committed_tables():
    n = {
        k: int(v)
        for k, v in re.findall(r"^(N_POS_EXP|N_NEG_EXP|N_FAM_EXP)=(\d+)", SBATCH.read_text(), re.M)
    }
    assert n["N_FAM_EXP"] == len(tsv(COMMITTED / "candidate_families.tsv")) == 3
    assert (
        n["N_NEG_EXP"] == len(tsv(ROOT / "data/controls/hard-negatives-shared/controls.tsv")) == 161
    )
    split = tsv(fs.MEMBERS_SPLIT)
    assert n["N_POS_EXP"] == len([r for r in split if r["class"] == "positive"]) == 132
