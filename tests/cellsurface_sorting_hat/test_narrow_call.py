"""Task 3.4 (D2, D4, D12): the narrow call reads a family only with evidence and an owner sign-off."""

import dataclasses
import gzip

import pytest
from test_cli import build

from cellsurface_sorting_hat.fasta import Protein
from cellsurface_sorting_hat.modules import pfam

COLS = [
    "pfam_acc",
    "name",
    "class",
    "module",
    "source_pmid",
    "pfam_release",
    "specificity_note",
    "second_condition",
    "active",
    "active_by",
    "active_date",
    "specificity_tier",
    "narrow_ok",
    "narrow_ok_by",
    "narrow_ok_date",
    "narrow_ok_evidence",
]


def fam_row(acc, module="pfam_adhesion", tier="smoke", narrow="no", **over):
    signed = narrow == "yes"
    base = {
        "pfam_acc": acc,
        "name": acc,
        "class": "x",
        "module": module,
        "source_pmid": "1",
        "pfam_release": "38.2",
        "specificity_note": "n",
        "second_condition": "",
        "active": "yes",
        "active_by": "owner",
        "active_date": "2026-10-09",
        "specificity_tier": tier,
        "narrow_ok": narrow,
        "narrow_ok_by": "owner" if signed else "",
        "narrow_ok_date": "2026-10-18" if signed else "",
        "narrow_ok_evidence": "9 of 54 E1 hits; 0 N1 hits" if signed else "",
    }
    base.update(over)
    return base


def write_table(path, rows, cols=COLS):
    lines = ["\t".join(cols)] + ["\t".join(r.get(c, "") for c in cols) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


# ---- the family table ---------------------------------------------------------------------------


def test_a_valid_table_loads_with_tier_and_narrow_flag(tmp_path):
    path = write_table(tmp_path / "t.tsv", [fam_row("PF00001", narrow="yes"), fam_row("PF00002")])
    first, second = pfam.load_family_table(path)
    assert (first.tier, first.narrow_ok) == ("smoke", True)
    assert (second.tier, second.narrow_ok) == ("smoke", False)


def test_a_table_without_specificity_tier_is_refused(tmp_path):
    old = [c for c in COLS if c not in ("specificity_tier",)]
    path = write_table(tmp_path / "t.tsv", [fam_row("PF00001")], old)
    with pytest.raises(pfam.FamilyTableError, match="specificity_tier"):
        pfam.load_family_table(path)


@pytest.mark.parametrize(
    ("over", "message"),
    [
        ({"specificity_tier": "good"}, "specificity_tier"),
        ({"narrow_ok": "maybe"}, "narrow_ok"),
        ({"specificity_tier": "untested", "narrow_ok": "yes"}, "untested"),
        ({"narrow_ok": "yes", "narrow_ok_by": ""}, "narrow_ok_by"),
        ({"narrow_ok": "yes", "narrow_ok_date": ""}, "narrow_ok_date"),
        ({"narrow_ok": "yes", "narrow_ok_evidence": ""}, "narrow_ok_evidence"),
    ],
)
def test_the_loader_refuses_these_rows(tmp_path, over, message):
    path = write_table(tmp_path / "t.tsv", [fam_row("PF00001", narrow="yes", **over)])
    with pytest.raises(pfam.FamilyTableError, match=message):
        pfam.load_family_table(path)


# ---- the Pfam rows ------------------------------------------------------------------------------


def prot(pid):
    return Protein(pid, "MKT", "x" * 64, "ok", "", 0.0)


def family(acc, tier, narrow, module="pfam_adhesion"):
    return pfam.Family(acc, acc, "x", module, "", True, tier=tier, narrow_ok=narrow)


def hit(target, acc):
    return {"target": target, "acc": acc, "ievalue": 1e-20, "score": 50.0}


def test_a_hit_in_a_family_that_is_not_narrow_sets_hit_but_not_hit_narrow():
    rows = pfam.pfam_rows(
        [prot("A"), prot("B")],
        [hit("A", "PF00001"), hit("B", "PF00002")],
        [family("PF00001", "untested", False), family("PF00002", "smoke", True)],
        "pfam_adhesion",
    )
    by = {r["id"]: r for r in rows}
    assert (by["A"]["hit"], by["A"]["hit_narrow"]) == ("1", "0")
    assert (by["B"]["hit"], by["B"]["hit_narrow"]) == ("1", "1")


def test_the_estimated_tier_is_also_narrow_when_signed():
    rows = pfam.pfam_rows(
        [prot("A")], [hit("A", "PF00001")], [family("PF00001", "estimated", True)], "pfam_adhesion"
    )
    assert rows[0]["hit_narrow"] == "1"


def test_hit_narrow_implies_hit_and_a_protein_without_a_hit_has_neither():
    families = [family("PF00001", "smoke", True), family("PF00002", "untested", False)]
    proteins = [prot(c) for c in "ABCD"]
    hits = [hit("A", "PF00001"), hit("B", "PF00002"), hit("C", "PF00001"), hit("C", "PF00002")]
    for r in pfam.pfam_rows(proteins, hits, families, "pfam_adhesion"):
        assert not (r["hit_narrow"] == "1" and r["hit"] != "1")
    by = {r["id"]: r for r in pfam.pfam_rows(proteins, hits, families, "pfam_adhesion")}
    assert (by["D"]["hit"], by["D"]["hit_narrow"]) == ("0", "0")
    assert by["C"]["families"] == "PF00001,PF00002" and by["C"]["hit_narrow"] == "1"


def test_only_pfam_adhesion_has_the_narrow_column():
    assert pfam.columns_for("pfam_adhesion") == ["hit", "hit_narrow", "families"]
    for module in ("pfam_allergen", "pfam_hydrophobin", "pfam_hsba"):
        assert pfam.columns_for(module) == ["hit", "families"]
        rows = pfam.pfam_rows(
            [prot("A")], [hit("A", "PF00001")], [family("PF00001", "smoke", True, module)], module
        )
        assert "hit_narrow" not in rows[0]


def test_the_family_digest_changes_with_tier_and_narrow_ok():
    from cellsurface_sorting_hat.modules.cli import _family_digest

    base = family("PF00001", "smoke", True)

    def digest(f):
        return _family_digest([f], "pfam_adhesion")

    assert digest(base) != digest(dataclasses.replace(base, narrow_ok=False))
    assert digest(base) != digest(dataclasses.replace(base, tier="estimated"))
    assert digest(base) == digest(dataclasses.replace(base))


OPTIONS = "# Option settings:     hmmsearch --cut_ga --cpu 2 --noali fam.hmm in.fasta\n"
DOM = "{t} - 289 {n} {a}.17 70 1e-20 60.1 8.9 1 1 1e-21 2e-20 59.0 8.9 1 70 20 90 20 91 0.9 -\n"


def test_the_pfam_command_writes_hit_narrow_for_pfam_adhesion_only(tmp_path):
    from cellsurface_sorting_hat.modules.cli import main

    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKTAYIAKQRQISFVKSHFSRQ\n>B\nMNLLPQWERTYIPASDFG\n")
    domtbl = tmp_path / "d.domtbl"
    domtbl.write_text(OPTIONS + DOM.format(t="A", n="X", a="PF00001") + "# [ok]\n")
    table = write_table(
        tmp_path / "t.tsv",
        [
            fam_row("PF00001", "pfam_adhesion", narrow="yes"),
            fam_row("PF00002", "pfam_allergen"),
            fam_row("PF00003", "pfam_hydrophobin"),
            fam_row("PF00004", "pfam_hsba"),
        ],
    )
    wd = tmp_path / "wd"
    argv = ["pfam", "--fasta", str(fasta), "--workdir", str(wd), "--domtbl", str(domtbl)]
    argv += ["--family-table", str(table), "--pfam-release", "38.2", "--pfam-sha256", "ab"]
    argv += ["--hmmer-version", "3.4"]
    assert main(argv) == 0
    for module in pfam.MODULES:
        with gzip.open(wd / "modules" / f"{module}.tsv.gz", "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
        assert ("hit_narrow" in header) == (module == "pfam_adhesion"), module
    with gzip.open(wd / "modules" / "pfam_adhesion.tsv.gz", "rt") as fh:
        lines = [ln.rstrip("\n").split("\t") for ln in fh]
    assert dict(zip(lines[0], lines[1], strict=True))["hit_narrow"] == "1"


# ---- the run ------------------------------------------------------------------------------------


def test_an_old_pfam_adhesion_table_without_hit_narrow_stops_the_run(
    tmp_path, write_module, nodes_dmp, capsys
):
    _, _, wd = build(tmp_path, write_module)
    ids = ["SOW1", "ENZ1", "DUP1", "STAR1"]
    write_module(wd, "pfam_adhesion", [{"id": i, "state": "ok", "hit": "0"} for i in ids])
    from cellsurface_sorting_hat.cli import main

    out = tmp_path / "out"
    argv = ["--fasta", str(tmp_path / "p.faa"), "--taxon-map", str(tmp_path / "taxa.tsv")]
    argv += ["--taxdump", str(nodes_dmp), "--workdir", str(wd), "--out", str(out)]
    assert main(argv) == 2
    assert "missing column 'hit_narrow'" in capsys.readouterr().err


def test_the_report_warns_when_no_family_has_narrow_ok(tmp_path, write_module, nodes_dmp):
    def report(narrow_families):
        run_dir = tmp_path / str(len(narrow_families))
        run_dir.mkdir()
        fasta, taxon_map, wd = build(run_dir, write_module)
        ids = ["SOW1", "ENZ1", "DUP1", "STAR1"]
        rows = [{"id": i, "state": "ok", "hit": "0", "hit_narrow": "0"} for i in ids]
        write_module(
            wd, "pfam_adhesion", rows, meta={"params": {"narrow_families": narrow_families}}
        )
        code, out, _ = run_cli_in(run_dir, nodes_dmp, fasta, taxon_map, wd)
        assert code == 0
        return (out / "report.md").read_text()

    assert "narrow_ok=yes" in report([]).split("## Known limits")[0]
    assert "narrow_ok=yes" not in report(["PF07691"]).split("## Known limits")[0]


def run_cli_in(tmp_path, nodes_dmp, fasta, taxon_map, wd):
    from cellsurface_sorting_hat.cli import main

    out = tmp_path / "out"
    argv = ["--fasta", str(fasta), "--taxon-map", str(taxon_map), "--taxdump", str(nodes_dmp)]
    argv += ["--workdir", str(wd), "--out", str(out)]
    return main(argv), out, wd
