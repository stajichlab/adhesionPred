import pytest

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules.pfam import (
    FAMILY_COLUMNS,
    FamilyTableError,
    NoActiveFamilyError,
    check_hit_ids,
    load_family_table,
    parse_domtblout,
    pfam_rows,
    specificity_report,
)

OPTIONS = "# Option settings:     hmmsearch --cut_ga --cpu 2 --noali fam.hmm in.fasta\n"
TRAILER = "# [ok]\n"
DOMTBL = (
    OPTIONS
    + (
        "# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias from to from to from to acc description\n"
        "P1 - 289 CFEM PF05730.17 70 1e-20 60.1 8.9 1 1 1e-21 2e-20 59.0 8.9 1 70 20 90 20 91 0.9 -\n"
        "P2 - 400 Hydrophobin PF01185.24 60 1e-12 40.0 0.0 1 1 1e-13 3e-12 39.0 0.0 1 60 5 65 5 66 0.9 -\n"
        "P3 - 300 AltA1 PF16541.11 150 1e-30 90.0 0.0 1 1 1e-31 1e-30 89.0 0.0 1 150 10 160 10 161 0.9 -\n"
        "P4 - 500 Asp PF00026.29 300 1e-40 120.0 0.0 1 1 1e-41 1e-40 119.0 0.0 1 300 10 310 10 311 0.9 -\n"
    )
    + TRAILER
)


def write_table(path, rows):
    lines = ["\t".join(FAMILY_COLUMNS)]
    for r in rows:
        lines.append("\t".join(r.get(c, "") for c in FAMILY_COLUMNS))
    path.write_text("\n".join(lines) + "\n")


def fam(acc, module="pfam_adhesion", active="yes", second="", **extra):
    base = {
        "pfam_acc": acc,
        "name": acc,
        "class": "x",
        "module": module,
        "source_pmid": "1",
        "pfam_release": "38.2",
        "specificity_note": "n",
        "second_condition": second,
        "active": active,
        "active_by": "owner" if active == "yes" else "",
        "active_date": "2026-10-05" if active == "yes" else "",
    }
    base.update(extra)
    return base


def prot(pid):
    return Protein(pid, "MKT", "x" * 64, "ok", "", 0.0)


def test_domtblout_is_parsed_without_the_accession_version(tmp_path):
    path = tmp_path / "d.domtbl"
    path.write_text(DOMTBL)
    hits = parse_domtblout(path)
    assert [(h["target"], h["acc"]) for h in hits] == [
        ("P1", "PF05730"),
        ("P2", "PF01185"),
        ("P3", "PF16541"),
        ("P4", "PF00026"),
    ]
    assert hits[0]["ievalue"] == 2e-20


def test_only_active_families_of_the_module_count(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(
        path, [fam("PF05730"), fam("PF01185", active="no"), fam("PF16541", module="pfam_allergen")]
    )
    fams = load_family_table(path)
    dpath = tmp_path / "d.domtbl"
    dpath.write_text(DOMTBL)
    hits = parse_domtblout(dpath)
    prots = [prot(p) for p in ("P1", "P2", "P3", "P4")]
    adh = {r["id"]: r["hit"] for r in pfam_rows(prots, hits, fams, "pfam_adhesion")}
    all_ = {r["id"]: r["hit"] for r in pfam_rows(prots, hits, fams, "pfam_allergen")}
    assert adh == {"P1": "1", "P2": "0", "P3": "0", "P4": "0"}  # PF01185 is inactive
    assert all_ == {"P1": "0", "P2": "0", "P3": "1", "P4": "0"}


def test_second_condition_needs_a_signal_peptide_call(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF00026", second="signal_peptide")])
    fams = load_family_table(path)
    dpath = tmp_path / "d.domtbl"
    dpath.write_text(DOMTBL)
    hits = parse_domtblout(dpath)
    p4 = [prot("P4")]
    assert pfam_rows(p4, hits, fams, "pfam_adhesion")[0]["hit"] == "0"  # no SP information
    assert pfam_rows(p4, hits, fams, "pfam_adhesion", {"P4": "not_called"})[0]["hit"] == "0"
    assert pfam_rows(p4, hits, fams, "pfam_adhesion", {"P4": "called"})[0]["hit"] == "1"


@pytest.mark.parametrize(
    "mutate,message",
    [
        (lambda r: r.update(pfam_acc="PF05730.17"), "no version"),
        (lambda r: r.update(module="other"), "module must be"),
        (lambda r: r.update(active="maybe"), "active must be"),
        (lambda r: r.update(second_condition="tm"), "second_condition"),
        (lambda r: r.update(active="yes", active_by=""), "needs active_by"),
    ],
)
def test_bad_family_rows_are_refused(tmp_path, mutate, message):
    row = fam("PF05730")
    mutate(row)
    path = tmp_path / "f.tsv"
    write_table(path, [row])
    with pytest.raises(FamilyTableError, match=message):
        load_family_table(path)


def test_duplicate_families_are_refused(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF05730"), fam("PF05730")])
    with pytest.raises(FamilyTableError, match="duplicate"):
        load_family_table(path)


def test_specificity_report_counts_and_lists_the_non_member_hits():
    rep = specificity_report(
        hit_ids={"a", "b", "x"}, member_ids={"a", "b", "c"}, universe_ids=set("abcxyz")
    )
    assert (rep["tp"], rep["fp"], rep["fn"], rep["tn"]) == (2, 1, 1, 2)
    assert rep["nonmember_hits"] == ["x"] and rep["missed_members"] == ["c"]
    assert rep["sensitivity"] == pytest.approx(2 / 3) and rep["specificity"] == pytest.approx(2 / 3)


def test_specificity_report_needs_members_inside_the_universe():
    with pytest.raises(ValueError):
        specificity_report({"a"}, {"zz"}, {"a"})


def test_no_tm_condition_drops_a_domain_in_a_protein_with_transmembrane_helices(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF05730", second="no_tm")])
    fams = load_family_table(path)
    dpath = tmp_path / "d.domtbl"
    dpath.write_text(DOMTBL)
    hits = parse_domtblout(dpath)
    p1 = [prot("P1")]
    assert pfam_rows(p1, hits, fams, "pfam_adhesion")[0]["hit"] == "0"  # no TM information
    assert (
        pfam_rows(p1, hits, fams, "pfam_adhesion", tm_counts={"P1": 7})[0]["hit"] == "0"
    )  # a receptor
    assert pfam_rows(p1, hits, fams, "pfam_adhesion", tm_counts={"P1": 0})[0]["hit"] == "1"


def test_a_domain_table_without_the_ok_trailer_or_cut_ga_is_refused(tmp_path):
    path = tmp_path / "d.domtbl"
    path.write_text(DOMTBL.replace(TRAILER, ""))
    with pytest.raises(ValueError, match="no '# \\[ok\\]' trailer"):
        parse_domtblout(path)
    path.write_text(DOMTBL.replace("--cut_ga", "-E 1e-5"))
    with pytest.raises(ValueError, match="not run with --cut_ga"):
        parse_domtblout(path)


def test_domain_table_targets_that_are_not_in_the_fasta_are_refused(tmp_path):
    path = tmp_path / "d.domtbl"
    path.write_text(DOMTBL)
    hits = parse_domtblout(path)
    with pytest.raises(ValueError, match="not in the FASTA"):
        check_hit_ids(hits, ["P1", "P2"])  # P3 and P4 are missing
    check_hit_ids(hits, ["P1", "P2", "P3", "P4", "P5"])  # no error


def test_a_module_with_no_active_family_raises_instead_of_writing_zeros(tmp_path):
    path = tmp_path / "f.tsv"
    write_table(path, [fam("PF05730", active="no")])
    with pytest.raises(NoActiveFamilyError):
        pfam_rows([prot("P1")], [], load_family_table(path), "pfam_adhesion")
