"""PR 6 family specificity drafts: pr6_family_specificity.py (fixtures, hand counts) and its sbatch (static checks)."""

import csv
import hashlib
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WF = ROOT / "analysis/wall_fusion_truth"
SBATCH = WF / "pr6_family_specificity.sbatch"
CAND = ROOT / "docs/reports/data/pr6_family_candidates_2026-10-22"
SCAN = ROOT / "docs/reports/data/pr6_pfam_scan_2026-10-22/scan_per_cluster.tsv"
MISSED = ROOT / "docs/reports/data/kiss_baseline_2026-10-22/clusters_missed.tsv"
PROPOSALS = ROOT / "docs/reports/data/pr6_pfam_scan_2026-10-22/panel_proposals.tsv"
MEASURE = ROOT / "data/sorting_hat/family_measure.tsv"
TABLE = ROOT / "data/sorting_hat/family_table.tsv"
INPUTS = ROOT / "scripts/sorting_hat/refresh_3_10/inputs.tsv"


def load():
    spec = importlib.util.spec_from_file_location(
        "pr6_family_specificity", WF / "pr6_family_specificity.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pr6_family_specificity"] = mod
    spec.loader.exec_module(mod)
    return mod


fs = load()


def tsv(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hit(acc="PF00001", name="Fam", frac=0.3, score=100.0, ievalue=1e-40, seq="s1"):
    return f"{acc}:{name}:frac={frac}:score={score}:ievalue={ievalue}:seq={seq}"


SCAN_HEADER = "slice\tcluster_id\tgroup\tgenes\tclades_fine\tstrata_fine\tjgi_accept\tn_members\tclass\tproposal\tpanel_models\tother_pfam_families\tweak_families\n"


def scan_row(cluster, cls, cell="", slice_="wall"):
    return f"{slice_}\t{cluster}\tG\tg\tClade\tclade_only\trow_verified\t1\t{cls}\t\t\t{cell}\t\n"


MISSED_HEADER = "slice\tcluster_id\tn_members\taccessions\tgenes\tspecies\tclades_fine\tstrata_fine\tseed_class_slice\tseed_class_all39\tbest_identity_all39\tjgi_accept\n"


def missed_row(cluster, accs, slice_="wall"):
    return f"{slice_}\t{cluster}\t{len(accs)}\t{','.join(accs)}\tg\tSp\tClade\tclade_only\tseed-independent\tseed-independent\t\trow_verified\n"


# ---- the candidate rule ----------------------------------------------------------------------------


def test_parse_other_pfam_reads_entries_including_a_jgi_local_key():
    cell = (
        hit("PF00001", "A", 0.5, 70.0, 1e-30, "P1")
        + ";"
        + hit("PF00002", "B", 0.1, 50.0, 0.0, "jgi:Neucr2:3723")
    )
    hits = fs.parse_other_pfam(cell)
    assert [h["acc"] for h in hits] == ["PF00001", "PF00002"]
    assert hits[1]["seq"] == "jgi:Neucr2:3723" and hits[1]["ievalue"] == 0.0
    assert fs.parse_other_pfam("") == []
    with pytest.raises(ValueError):
        fs.parse_other_pfam("not a hit")


@pytest.mark.parametrize(
    "frac,ievalue,cov,ev",
    [
        (0.5, 1.0, True, False),  # coverage limit is inclusive
        (0.499, 1.0, False, False),
        (0.1, 1e-30, False, True),  # E-value limit is inclusive
        (0.1, 1.1e-30, False, False),
        (0.1, 0.0, False, True),  # an E-value of 0 meets the rule
        (0.9, 1e-50, True, True),
    ],
)
def test_the_rule_at_its_boundaries(frac, ievalue, cov, ev):
    h = {"frac": frac, "ievalue": ievalue}
    assert fs.rule_met(h) == (cov, ev)


def test_a_family_is_a_candidate_only_through_its_own_hit(tmp_path):
    scan = [
        {
            "slice": "wall",
            "cluster_id": "KS1",
            "group": "G",
            "genes": "g",
            "clades_fine": "C",
            "strata_fine": "clade_only",
            "jgi_accept": "row_verified",
            "class": "hit_other_pfam",
            "other_pfam_families": hit("PF00001", "Strong", 0.1, 90.0, 1e-35, "m1")
            + ";"
            + hit("PF00002", "Weak", 0.2, 20.0, 1e-5, "m1"),
        },
        {
            "slice": "wall",
            "cluster_id": "KS2",
            "group": "G",
            "genes": "g",
            "clades_fine": "C",
            "strata_fine": "clade_only",
            "jgi_accept": "row_verified",
            "class": "hit_other_pfam",
            "other_pfam_families": hit("PF00003", "Weak3", 0.49, 20.0, 1e-5, "m2"),
        },
        {
            "slice": "wall",
            "cluster_id": "KS3",
            "group": "G",
            "genes": "g",
            "clades_fine": "C",
            "strata_fine": "clade_only",
            "jgi_accept": "row_verified",
            "class": "hit_panel_model",
            "other_pfam_families": hit("PF00004", "Panel", 0.9, 90.0, 1e-50, "m3"),
        },
        {
            "slice": "wall",
            "cluster_id": "KS4",
            "group": "G",
            "genes": "g",
            "clades_fine": "C",
            "strata_fine": "clade_only",
            "jgi_accept": "row_verified",
            "class": "no_pfam_hit",
            "other_pfam_families": "",
        },
    ]
    missed = [{"cluster_id": f"KS{i}", "accessions": f"M{i},N{i}"} for i in (1, 2, 3, 4)]
    sel = fs.select_candidates(scan, missed)
    assert [c["cluster_id"] for c in sel["clusters"]] == [
        "KS1"
    ]  # KS2 below both limits; KS3 and KS4 are other classes
    assert sorted(sel["families"]) == [
        "PF00001"
    ]  # the weak hit of the same cluster does not make PF00002 a candidate
    assert sel["families"]["PF00001"]["clusters"] == {"KS1"}


def test_the_committed_candidate_tables_are_the_rule_applied_to_the_committed_scan(tmp_path):
    out = tmp_path / "cand"
    sel = fs.write_candidates(SCAN, MISSED, out)
    for name in fs.CANDIDATE_OUTPUTS:
        assert (out / name).read_bytes() == (CAND / name).read_bytes(), name
    assert len(sel["clusters"]) == 26 and len(sel["families"]) == 30
    rows = tsv(SCAN)
    other = [r for r in rows if r["class"] == "hit_other_pfam"]
    assert len(other) == 38
    proposals = {r["cluster_id"] for r in other if r["proposal"] == "add_to_panel_proposal"}
    assert len(proposals) == 14
    assert proposals <= {
        c["cluster_id"] for c in sel["clusters"]
    }  # all 14 proposals pass (coverage >= 0.5)
    # the other 12 are partial-domain hits with an E-value of 1e-30 or better
    extra = {c["cluster_id"] for c in sel["clusters"]} - proposals
    assert len(extra) == 12
    for r in other:
        if r["cluster_id"] in extra:
            assert any(h["ievalue"] <= 1e-30 for h in fs.parse_other_pfam(r["other_pfam_families"]))
        elif r["cluster_id"] not in proposals:  # the 12 clusters left out: no hit meets either rule
            assert not any(
                any(fs.rule_met(h)) for h in fs.parse_other_pfam(r["other_pfam_families"])
            )
    assert {r["pfam_acc"] for r in tsv(PROPOSALS)} <= set(
        sel["families"]
    )  # every proposed family is tested
    assert not any(
        c["cluster_id"] in {"KS0082", "KS0116"} for c in sel["clusters"]
    )  # the two panel hits are not candidates


def test_candidates_refuses_a_non_empty_output_directory(tmp_path):
    (tmp_path / "x").write_text("keep")
    with pytest.raises(ValueError):
        fs.write_candidates(SCAN, MISSED, tmp_path)


# ---- the specificity table -------------------------------------------------------------------------


def dom(target, acc, name="Fam"):
    """One hmmsearch --domtblout line (target = sequence, query = model)."""
    return f"{target} - 300 {name} {acc}.5 100 1e-30 100.0 0.0 1 1 1e-30 1e-30 100.0 0.0 1 100 10 110 10 110 0.99 description\n"


def domtbl(rows, option="--cut_ga"):
    return f"# hmmsearch :: search profile(s) against a sequence database\n# Option settings: hmmsearch {option} x\n#\n{''.join(rows)}# [ok]\n"


@pytest.fixture()
def fx(tmp_path):
    cand = tmp_path / "cand"
    cand.mkdir()
    (cand / "candidate_families.tsv").write_text(
        "pfam_acc\tname\tn_clusters\tclusters\tslices\tn_hits\tmax_frac\tmin_ievalue\trule_coverage\trule_evalue\tin_family_table\n"
        "PF00001\tFamA\t2\tKS1;KS2\twall\t2\t0.6\t1e-40\tyes\tyes\tno\n"
        "PF00002\tFamB\t1\tKS3\tadhesion\t1\t0.2\t1e-35\tno\tyes\tno\n"
    )
    (cand / "family_members.tsv").write_text(
        "pfam_acc\tcluster_id\taccession\n"
        "PF00001\tKS1\tM1\nPF00001\tKS1\tM2\nPF00001\tKS2\tM3\nPF00002\tKS3\tM4\n"
    )
    split = tmp_path / "members_split.tsv"
    split.write_text(
        "accession\tcluster_id\tclass\tslice\n"
        "M1\tKS1\tpositive\twall\nM2\tKS1\tpositive\twall\nM3\tKS2\tpositive\twall\nM4\tKS3\tpositive\tadhesion\n"
        "M5\tKS9\tpositive\twall\nM6\tKS8\tpositive\tadhesion\nM7\tKS7\tpositive\tfusion\n"
        "N1\tKS20\tnegative\t\n"
    )
    neg = tmp_path / "neg"
    neg.mkdir()
    (neg / "controls.tsv").write_text(
        "accession\tevidence_level\nN1\tN1\nN2\tN1\nN3\tN1\nN4\tN1\nN5\tN2\nN6\tN2\n"
    )
    (neg / "clusters.tsv").write_text(
        "accession\tcluster_id\tset\n"
        "N1\tHN1\thard_negative\nN2\tHN1\thard_negative\nN3\tHN2\thard_negative\nN4\tHN3\thard_negative\n"
        "N5\tHN4\thard_negative\nN6\tHN5\thard_negative\nP9\tHN9\tadhesin_positive\n"
    )
    pos_dom = tmp_path / "pos.dom"
    # PF00001 hits M1, M3 (own) and M5 (other, wall), M6 (other, adhesion); PF00002 hits M4 (own)
    pos_dom.write_text(
        domtbl(
            [
                dom("M1", "PF00001"),
                dom("M3", "PF00001"),
                dom("M5", "PF00001"),
                dom("M6", "PF00001"),
                dom("M4", "PF00002", "FamB"),
            ]
        )
    )
    neg_dom = tmp_path / "neg.dom"
    # PF00001 hits N1, N2 (one N1 cluster), N5 (N2); PF00002 hits nothing
    neg_dom.write_text(domtbl([dom("N1", "PF00001"), dom("N2", "PF00001"), dom("N5", "PF00001")]))
    proteomes = {}
    for name, n, hits_ in (("P_a", 4, ["a1", "a2"]), ("P_b", 6, ["b1"])):
        fa = tmp_path / f"{name}.faa"
        fa.write_text(
            "".join(f">{name[2]}{i + 1} annotation of {name[2]}{i + 1}\nACDEFG\n" for i in range(n))
        )
        d = tmp_path / f"{name}.dom"
        d.write_text(domtbl([dom(h, "PF00001") for h in hits_]))
        proteomes[name] = (d, fa)
    return {
        "cand": cand,
        "split": split,
        "neg": neg,
        "pos": pos_dom,
        "negdom": neg_dom,
        "prot": proteomes,
        "tmp": tmp_path,
    }


def run_summary(fx, out, expected=None):
    return fs.summarize(
        fx["cand"], fx["split"], fx["neg"], fx["pos"], fx["negdom"], fx["prot"], out, expected
    )


def test_hand_counted_rates_and_the_family_own_denominator(fx):
    rows = {r["pfam_acc"]: r for r in run_summary(fx, fx["tmp"] / "out")}
    a, b = rows["PF00001"], rows["PF00002"]
    # the denominator is the family's own members (M1, M2, M3), not every positive
    assert (a["n_own_members"], a["own_hits"], a["own_recall_circular"]) == (3, 2, "0.667")
    assert a["n_candidate_clusters"] == 2 and a["leakage"] == "tuned_on_truth"
    # M5 and M6 hit the family but are not its own members: counted apart, by slice, outside the denominator
    assert (
        a["other_positive_hits_wall"],
        a["other_positive_hits_adhesion"],
        a["other_positive_hits_fusion"],
    ) == (1, 1, 0)
    assert (b["n_own_members"], b["own_hits"], b["other_positive_hits_wall"]) == (1, 1, 0)
    # negatives: 4 N1 rows in 3 clusters, 2 N2 rows; PF00001 hits N1, N2 (cluster HN1 once) and the N2 row N5
    assert (a["n_neg_n1"], a["neg_hits_n1"], a["n_clusters_neg_n1"], a["neg_clusters_hit_n1"]) == (
        4,
        2,
        3,
        1,
    )
    assert (a["n_neg_n2"], a["neg_hits_n2"]) == (2, 1)
    assert (b["neg_hits_n1"], b["neg_hits_n2"]) == (0, 0)
    # proteomes: 2 hits of 4 and 1 of 6 proteins: 5000.00 and 1666.67 per 10,000; kept total 3 of 10 = 3000.00
    assert (a["P_a_hits"], a["P_a_per10k"], a["P_b_hits"], a["P_b_per10k"]) == (
        2,
        "5000.00",
        1,
        "1666.67",
    )
    assert (a["kept_hits_total"], a["kept_n_proteins"], a["kept_per10k"]) == (3, 10, "3000.00")
    assert b["kept_hits_total"] == 0 and b["kept_per10k"] == "0.00"


def test_hit_tables_carry_the_header_text_and_the_own_flag(fx):
    out = fx["tmp"] / "out"
    run_summary(fx, out)
    pos = {
        (r["pfam_acc"], r["accession"]): r["own_or_other"] for r in tsv(out / "positives_hits.tsv")
    }
    assert pos[("PF00001", "M1")] == "own" and pos[("PF00001", "M5")] == "other"
    neg = {(r["pfam_acc"], r["accession"]): r for r in tsv(out / "negatives_hits.tsv")}
    assert (
        neg[("PF00001", "N5")]["evidence_level"] == "N2"
        and neg[("PF00001", "N1")]["neg_cluster"] == "HN1"
    )
    # the proteome hits keep the annotation (the FASTA header text), never the sequence
    ph = tsv(out / "proteome_hits.tsv")
    assert {r["protein"] for r in ph} == {"a1", "a2", "b1"}
    assert all(r["annotation"].startswith("annotation of") for r in ph)
    assert all(not re.search(r"[ACDEFGHIKLMNPQRSTVWY]{6,}", ";".join(r.values())) for r in ph)


def test_outputs_hold_no_sequences_and_are_deterministic(fx):
    a, b = fx["tmp"] / "o1", fx["tmp"] / "o2"
    run_summary(fx, a)
    run_summary(fx, b)
    for name in fs.SUMMARY_OUTPUTS:
        assert (a / name).read_bytes() == (b / name).read_bytes(), name
        for line in (a / name).read_text().splitlines():
            assert not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{20,}", line)


@pytest.mark.parametrize(
    "what", ["trailer", "cut_ga", "orientation", "stray_pos", "stray_neg", "count", "exists"]
)
def test_refusals(fx, what):
    out = fx["tmp"] / "out"
    expected = None
    if what == "trailer":
        fx["pos"].write_text(fx["pos"].read_text().replace("# [ok]\n", ""))
    elif what == "cut_ga":
        fx["negdom"].write_text(domtbl([dom("N1", "PF00001")], "-E 1e-5"))
    elif what == "orientation":  # a model accession in the sequence column (hmmscan layout)
        fx["pos"].write_text(
            domtbl(
                [
                    "FamA PF00001.5 300 M1 - 100 1e-30 100.0 0.0 1 1 1e-30 1e-30 100.0 0.0 1 100 10 110 10 110 0.99 d\n"
                ]
            )
        )
    elif what == "stray_pos":
        fx["pos"].write_text(domtbl([dom("N1", "PF00001")]))  # a negative among the positives
    elif what == "stray_neg":
        fx["negdom"].write_text(domtbl([dom("M1", "PF00001")]))
    elif what == "count":
        expected = {"P_a": 5, "P_b": 6}
    elif what == "exists":
        out.mkdir()
        (out / "keep").write_text("x")
    with pytest.raises(ValueError):
        run_summary(fx, out, expected)


# ---- the dry-run converter -------------------------------------------------------------------------


def test_converter_rows_are_shared_set_rows_capped_at_smoke(fx):
    out = fx["tmp"] / "out"
    run_summary(fx, out)
    rows = fs.convert(
        out / "family_specificity.tsv", None, "src", "2026-10-22", "38.2", "HMMER 3.4"
    )
    from cellsurface_sorting_hat import family_tier

    assert [r["pfam_acc"] for r in rows] == ["PF00001", "PF00002"]
    for r in rows:
        assert r["measure_basis"] == "descriptive_shared_sets"
        assert r["denominator"] == "family_own_kiss_members"
        assert (
            r["n_pos"] == r["n_pos_unknown"] == r["pos_denominator"]
        )  # nobody checked the seeds: all unknown
        assert r["sensitivity"] == r["specificity"] == ""
        assert r["_tier"] in (family_tier.UNTESTED, family_tier.SMOKE)
        assert "circular" in r["notes"] and "not owner-reviewed" in r["notes"]
    a = rows[0]
    assert (
        a["n_neg"],
        a["neg_hits"],
        a["n_clusters_neg"],
        a["neg_clusters_hit"],
        a["neg_hits_n2"],
    ) == ("4", "2", "3", "1", "1")


def test_even_a_large_measurement_is_never_estimated(tmp_path):
    from cellsurface_sorting_hat import family_tier

    spec = tmp_path / "spec.tsv"
    spec.write_text(
        "pfam_acc\tname\tn_own_members\town_hits\tn_candidate_clusters\tn_neg_n1\tneg_hits_n1\tn_clusters_neg_n1\t"
        "neg_clusters_hit_n1\tn_neg_n2\tneg_hits_n2\tleakage\tkept_hits_total\tkept_n_proteins\tkept_per10k\tP_a_hits\n"
        "PF00001\tFamA\t60\t60\t50\t137\t0\t76\t0\t24\t0\ttuned_on_truth\t0\t32491\t0.00\t0\n"
    )
    rows = fs.convert(spec, None, "src", "2026-10-22", "38.2", "HMMER 3.4")
    assert (
        rows[0]["_tier"] != family_tier.ESTIMATED
    )  # descriptive_shared_sets is capped at smoke by the derivation


def test_converter_refuses_the_committed_files_and_changes_nothing(fx):
    out = fx["tmp"] / "out"
    run_summary(fx, out)
    before = (digest(MEASURE), digest(TABLE))
    for forbidden in (MEASURE, TABLE):
        with pytest.raises(ValueError):
            fs.convert(
                out / "family_specificity.tsv", forbidden, "src", "2026-10-22", "38.2", "HMMER 3.4"
            )
    fs.convert(
        out / "family_specificity.tsv", None, "src", "2026-10-22", "38.2", "HMMER 3.4"
    )  # dry run
    proposed = fx["tmp"] / "proposed.tsv"
    fs.convert(out / "family_specificity.tsv", proposed, "src", "2026-10-22", "38.2", "HMMER 3.4")
    assert tsv(proposed)[0]["pfam_acc"] == "PF00001"
    assert (
        digest(MEASURE),
        digest(TABLE),
    ) == before  # no committed measurement or family row changed


def test_the_script_never_names_the_committed_files_as_outputs():
    text = (WF / "pr6_family_specificity.py").read_text()
    code = [ln for ln in text.splitlines() if not ln.lstrip().startswith(("#", '"""'))]
    writes = [
        ln for ln in code if re.search(r"write_text|write_tsv|write_measure|open\(.*[\"']w", ln)
    ]
    assert not any("FAMILY_MEASURE" in ln or "FAMILY_TABLE" in ln for ln in writes)


def test_a_candidate_family_that_is_already_measured_is_flagged(tmp_path):
    spec = tmp_path / "spec.tsv"
    spec.write_text(
        "pfam_acc\tname\tn_own_members\town_hits\tn_candidate_clusters\tn_neg_n1\tneg_hits_n1\tn_clusters_neg_n1\t"
        "neg_clusters_hit_n1\tn_neg_n2\tneg_hits_n2\tleakage\tkept_hits_total\tkept_n_proteins\tkept_per10k\tP_a_hits\n"
        "PF05730\tCFEM\t1\t1\t1\t137\t1\t76\t1\t24\t0\ttuned_on_truth\t0\t32491\t0.00\t0\n"
    )
    rows = fs.convert(spec, None, "src", "2026-10-22", "38.2", "HMMER 3.4")
    assert rows[0]["_already_measured"] == "yes"


# ---- the kept proteomes and the job script ---------------------------------------------------------


def test_kept_proteomes_are_one_per_species_and_match_the_inputs_table():
    counts = fs.kept_proteome_counts(INPUTS)
    assert list(counts) == ["Afum_Af293_UniProt", "Calb_SC5314", "Cimm_RS", "Scer_S288C"]
    assert sum(counts.values()) == 32491
    assert (
        "Afum_Af293_Fungi5k" not in counts
        and "Bder_ER3" not in counts
        and "Scer_S288C_nodubious" not in counts
    )
    species = {
        "Afum_Af293_UniProt": 330879,
        "Calb_SC5314": 237561,
        "Cimm_RS": 246410,
        "Scer_S288C": 559292,
    }
    taxa = {r["prot"]: int(r["taxon"]) for r in tsv(INPUTS) if r["taxon"].isdigit()}
    assert {p: taxa[p] for p in counts} == species
    # one proteome per species: A. fumigatus (330879), C. albicans, C. immitis and S. cerevisiae are four different taxa
    assert len(set(species.values())) == 4


def test_sbatch_follows_the_job_conventions():
    assert subprocess.run(["bash", "-n", str(SBATCH)], capture_output=True).returncode == 0
    text = SBATCH.read_text()
    code = "\n".join(
        ln
        for ln in text.splitlines()
        if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH")
    )
    assert text.startswith("#!/bin/bash -l")
    assert "DRAFT, NOT SUBMITTED" in text.splitlines()[1]
    assert "set -euo pipefail" in code and "#SBATCH -p short" in code
    assert '"${SCRATCH:?' in code
    assert "/tmp" not in "\n".join(ln for ln in code.splitlines() if not ln.startswith("#SBATCH"))
    assert "BASH_SOURCE" not in code
    assert "module load hmmer/3.4" in code
    searches = [ln.strip() for ln in code.splitlines() if ln.strip().startswith("hmmsearch ")]
    assert (
        len(searches) == 1 and "--cut_ga" in searches[0] and "-E " not in searches[0]
    )  # one GA search function, no E-value cut-off
    assert '[ ! -e "$OUT" ]' in code and "check_pin" in code and "sidecar_ok" in code
    assert (
        'Pfam${PFAM_RELEASE}/"' in code
        and "4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e" in code
    )
    assert 'KEPT="Afum_Af293_UniProt Calb_SC5314 Cimm_RS Scer_S288C"' in code
    for forbidden in ("Fungi5k", "Bder_ER3", "nodubious"):
        assert forbidden not in code
    # the pin covers the scripts, the committed candidate tables, the record and the negatives
    pin = next(ln for ln in code.splitlines() if "git -C" in ln and "diff --quiet" in ln)
    for part in ('"$WF"', '"$SCAN"', '"$CAND"', '"$NEG"', '"$INPUTS"'):
        assert part in pin
    executable = "\n".join(
        ln for ln in code.splitlines() if not ln.lstrip().startswith(("echo", "}"))
    )
    assert (
        "family_measure" not in executable
    )  # the job names no measurement file (only its provenance text says so)
    assert "family_table.tsv" not in executable.replace(
        "FAMILY=data/sorting_hat/family_table.tsv", ""
    )
    assert code.count("N_POS_EXP=132") == 1 and "N_FAM_EXP=30" in code and "N_NEG_EXP=161" in code


def test_the_pinned_counts_match_the_committed_tables():
    n = {
        k: int(v)
        for k, v in re.findall(r"^(N_POS_EXP|N_NEG_EXP|N_FAM_EXP)=(\d+)", SBATCH.read_text(), re.M)
    }
    split = tsv(ROOT / "docs/reports/data/kiss_cluster_split_2026-10-22/members_split.tsv")
    assert n["N_POS_EXP"] == len([r for r in split if r["class"] == "positive"]) == 132
    assert (
        n["N_NEG_EXP"] == len(tsv(ROOT / "data/controls/hard-negatives-shared/controls.tsv")) == 161
    )
    assert n["N_FAM_EXP"] == len(tsv(CAND / "candidate_families.tsv")) == 30
