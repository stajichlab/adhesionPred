"""Tests of build_control_set.py (plan Task 5.2). Fixture tables only: no network, no JGI data."""

import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "build_control_set", ROOT / "analysis/wall_fusion_truth/build_control_set.py"
)
bcs = importlib.util.module_from_spec(spec)
sys.modules["build_control_set"] = bcs
spec.loader.exec_module(bcs)

SEQ_A = "M" + "ACDEFGHIKL" * 6
SEQ_B = "M" + "LKIHGFEDCA" * 6
SEQ_C = "M" + "KLMNPQRSTV" * 6
SEQ_D = "M" + "VTSRQPNMLK" * 6
SEQ_E = "M" + "WYACDEFGHI" * 6  # a UniProt sequence that differs from the JGI one
SEQ_F = "M" + "QRSTVWYACD" * 6
SEQ_G = "M" + "EFGHIKLMNP" * 6
JGI_ONLY = "M" + "HIKLMNPQRS" * 6  # a JGI sequence: it must never reach an output


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def mrow(gene, **kw):
    base = {
        "gene_or_family_id": gene,
        "cluster_number": "1",
        "category": "wall_remodeling",
        "slice": "wall",
        "protein_id_jgi": "",
        "organism_code": "Neucr2",
        "clade_coarse": "Pezizomycotina",
        "clade_fine": "Sordariomycetes",
        "map_status": "mapped_jgi_file",
        "map_reason": "",
        "match_term": "",
        "accession": "",
        "uniprot_organism_id": "367110",
        "control_role": "kiss_row",
        "jgi_portal": "Neucr2",
        "jgi_seq_sha256": "",
        "jgi_accept": "row_verified",
    }
    base.update(kw)
    return base


ROWS = [
    # public sequence, identical to JGI
    mrow(
        "chs1",
        cluster_number="10",
        protein_id_jgi="11_Neucr2",
        accession="P00001",
        jgi_seq_sha256=sha(SEQ_A),
    ),
    # two Kiss rows for one protein: one control row
    mrow(
        "cot-2",
        cluster_number="20",
        protein_id_jgi="22_Neucr2",
        accession="P00002",
        jgi_seq_sha256=sha(SEQ_B),
    ),
    mrow(
        "gs-1",
        cluster_number="20",
        protein_id_jgi="22_Neucr2",
        accession="P00002",
        jgi_seq_sha256=sha(SEQ_B),
    ),
    # a locus-tag row of another clade (the cluster 30 is shared)
    mrow(
        "locus1",
        cluster_number="30",
        map_status="mapped_locus_tag",
        accession="P00003",
        organism_code="Aspfu1",
        clade_fine="Eurotiomycetes",
        jgi_accept="",
        jgi_portal="",
    ),
    mrow(
        "chs3",
        cluster_number="30",
        protein_id_jgi="33_Neucr2",
        accession="P00004",
        jgi_seq_sha256=sha(SEQ_D),
    ),
    # no UniProt accession: local key, the sequence stays in the JGI file
    mrow(
        "ncl1",
        cluster_number="40",
        category="adhesion",
        slice="adhesion",
        protein_id_jgi="44_Neucr2",
        jgi_seq_sha256=sha(JGI_ONLY),
        jgi_accept="portal_validated",
    ),
    # the UniProt sequence differs from the JGI one
    mrow(
        "ste2",
        cluster_number="50",
        category="adhesion",
        slice="adhesion",
        protein_id_jgi="55_Neucr2",
        accession="P00005",
        jgi_seq_sha256=sha(SEQ_C),
    ),
    # no accession by name, but exactly one cached UniProt entry has this very sequence
    mrow(
        "seqhit",
        cluster_number="90",
        category="adhesion",
        slice="adhesion",
        protein_id_jgi="99_Neucr2",
        jgi_seq_sha256=sha(SEQ_F),
        jgi_accept="portal_validated",
    ),
    # two cached UniProt entries have this sequence: the accession is ambiguous, so a local key
    mrow(
        "twin",
        cluster_number="91",
        category="adhesion",
        slice="adhesion",
        protein_id_jgi="98_Neucr2",
        jgi_seq_sha256=sha(SEQ_G),
        jgi_accept="portal_validated",
    ),
    # a name conflict: not a control
    mrow(
        "hyr1",
        cluster_number="60",
        category="adhesion",
        slice="adhesion",
        map_status="lost",
        protein_id_jgi="66_Neucr2",
        jgi_accept="",
        map_reason="jgi_header_locus_tag_entry_differs_from_name_candidate",
    ),
    # a family representative: not a control
    mrow(
        "afuX",
        cluster_number="70",
        map_status="mapped_family_rep",
        control_role="family_rep",
        organism_code="",
        clade_coarse="",
        clade_fine="",
        accession="P00006",
        match_term="AFUA_1G00010",
        uniprot_organism_id="330879",
        jgi_accept="",
        jgi_portal="",
    ),
    # a family row without a locus
    mrow(
        "NA",
        cluster_number="80",
        map_status="lost",
        control_role="unmapped_family_row",
        organism_code="",
        clade_coarse="",
        clade_fine="",
        jgi_accept="",
        jgi_portal="",
        map_reason="no_jgi_id_and_no_afu_locus",
    ),
]
ORGS = {
    "Neucr2": {
        "code": "Neucr2",
        "ncbi_taxid": "5141",
        "ncbi_name": "Neurospora crassa",
        "scope_taxid": "367110",
    },
    "Aspfu1": {
        "code": "Aspfu1",
        "ncbi_taxid": "746128",
        "ncbi_name": "Aspergillus fumigatus",
        "scope_taxid": "330879",
    },
}
UNIPROT = {
    "P00001": SEQ_A,
    "P00002": SEQ_B,
    "P00003": SEQ_C,
    "P00004": SEQ_D,
    "P00005": SEQ_E,
    "P00006": SEQ_A + "K",
    "P00007": SEQ_F,
    "P00008": SEQ_G,
    "P00009": SEQ_G,
}


def cache_text():
    head = "Entry\tReviewed\tEntry Name\tGene Names\tOrganism\tOrganism (ID)\tLength\tSequence\n"
    body = "".join(f"{a}\treviewed\t{a}_X\tg\tx\t1\t{len(s)}\t{s}\n" for a, s in UNIPROT.items())
    return {"https://rest.uniprot.org/uniprotkb/search?q=1": head + body}


@pytest.fixture
def inputs(tmp_path):
    cols = list(ROWS[0])
    with (tmp_path / "row_mapping.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in ROWS:
            w.writerow([r[c] for c in cols])
    with (tmp_path / "organisms.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["code", "ncbi_taxid", "ncbi_name", "scope_taxid"])
        for o in ORGS.values():
            w.writerow([o[k] for k in ("code", "ncbi_taxid", "ncbi_name", "scope_taxid")])
    (tmp_path / "cache.json").write_text(json.dumps(cache_text()))
    with (tmp_path / "kiss.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(
            [
                "gene_or_family_id",
                "cluster_number",
                "category",
                "source_sheet",
                "protein_id_jgi",
                "annotation",
            ]
        )
        for r in ROWS:
            w.writerow(
                [
                    r["gene_or_family_id"],
                    r["cluster_number"],
                    r["category"],
                    "Supp Data 3 (MOESM7)",
                    r["protein_id_jgi"],
                    "chitin synthase " + r["gene_or_family_id"],
                ]
            )
    return tmp_path


def run(inputs, out, *extra):
    return bcs.main(
        [
            "--out",
            str(out),
            "--row-mapping",
            str(inputs / "row_mapping.tsv"),
            "--organisms",
            str(inputs / "organisms.tsv"),
            "--cache",
            str(inputs / "cache.json"),
            "--kiss-table",
            str(inputs / "kiss.tsv"),
            *extra,
        ]
    )


def tsv(path):
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_split_into_public_local_family_and_rejected(inputs, tmp_path):
    out = tmp_path / "set"
    assert run(inputs, out) == 0
    controls = {r["accession"]: r for r in tsv(out / "controls.tsv")}
    local = {r["accession"]: r for r in tsv(out / "controls_local_sequence.tsv")}
    assert set(controls) == {"P00001", "P00002", "P00003", "P00004", "P00005", "P00007"}
    assert set(local) == {"jgi:Neucr2:44", "jgi:Neucr2:98"}  # the local key is marked jgi:
    assert local["jgi:Neucr2:44"]["sequence_in_repo"] == "no"
    assert local["jgi:Neucr2:44"]["sequence_source"] == "jgi_file_local"
    assert local["jgi:Neucr2:44"]["sequence_sha256"] == sha(JGI_ONLY)
    # the accession of seqhit was found by its sequence; the twin has two candidates and stays local
    assert (
        controls["P00007"]["gene"] == "seqhit"
        and "identical sequence" in controls["P00007"]["notes"]
    )
    assert controls["P00007"]["sequence_matches_jgi"] == "yes"
    # family representatives and unmapped rows are never controls
    assert "P00006" not in controls
    fam = tsv(out / "family_reps.tsv")
    assert [(r["accession"], r["may_be_test_control"], r["counts_toward_floor"]) for r in fam] == [
        ("P00006", "no", "no")
    ]
    rej = {r["gene_or_family_id"]: r for r in tsv(out / "rejected.tsv")}
    assert set(rej) == {"hyr1", "NA"}
    assert rej["hyr1"]["reason"].startswith("name_conflict:")
    assert "name_conflict_review" in rej["hyr1"]["see"]


def test_two_rows_of_one_protein_are_one_control(inputs, tmp_path):
    out = tmp_path / "set"
    run(inputs, out)
    row = {r["accession"]: r for r in tsv(out / "controls.tsv")}["P00002"]
    assert row["gene"] == "cot-2/gs-1"
    assert "merged Kiss rows" in row["notes"]


def test_strata_flags_and_sequence_match(inputs, tmp_path):
    out = tmp_path / "set"
    run(inputs, out)
    c = {r["accession"]: r for r in tsv(out / "controls.tsv")}
    # cluster 30 of the wall slice has a Sordariomycetes and an Eurotiomycetes row: shared
    assert c["P00003"]["cluster_stratum_fine"] == c["P00004"]["cluster_stratum_fine"] == "shared"
    assert c["P00001"]["cluster_stratum_fine"] == "clade_only"
    assert c["P00001"]["sequence_matches_jgi"] == "yes"
    assert c["P00003"]["sequence_matches_jgi"] == "not_checked"  # no JGI protein
    assert c["P00005"]["sequence_matches_jgi"] == "no"
    assert "differs from the JGI" in c["P00005"]["notes"]
    assert c["P00001"]["evidence_level"] == "row_verified"
    assert c["P00003"]["evidence_level"] == "locus_tag"
    assert c["P00001"]["leakage"] == "unknown" and c["P00001"]["tuned"] == "no"
    assert c["P00001"]["cluster_id"] == "KISS-wall-10"


def test_no_jgi_sequence_reaches_any_output_file(inputs, tmp_path):
    out = tmp_path / "set"
    run(inputs, out)
    for f in out.iterdir():
        assert JGI_ONLY not in f.read_text(), f.name
    # the public sequences are there
    assert SEQ_A in (out / "sequences.faa").read_text()


def test_guard_rejects_a_table_with_a_sequence(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("a\tb\nx\t" + JGI_ONLY + "\n")
    with pytest.raises(RuntimeError, match="looks like a protein sequence"):
        bcs.guard_tables_have_no_sequences([p])
    p.write_text("a\tb\nx\tshort text\n")
    bcs.guard_tables_have_no_sequences([p])


def test_the_set_passes_the_l1_validator(inputs, tmp_path):
    from cellsurface_sorting_hat.controls import validate_controls

    out = tmp_path / "set"
    run(inputs, out)
    rep = validate_controls(out)
    assert rep.ok, [str(d) for d in rep.defects]
    assert rep.counts["rows"] == 6 and rep.counts["sequences"] == 6


def test_two_runs_give_the_same_bytes(inputs, tmp_path):
    run(inputs, tmp_path / "a")
    run(inputs, tmp_path / "b")
    for f in sorted((tmp_path / "a").iterdir()):
        assert f.read_bytes() == (tmp_path / "b" / f.name).read_bytes(), f.name


def test_one_protein_in_two_slices_is_an_error(inputs, tmp_path):
    rows = [dict(r) for r in ROWS]
    dup = dict(rows[0])
    dup["gene_or_family_id"] = "chs1b"
    dup["slice"] = "adhesion"
    dup["category"] = "adhesion"
    rows.append(dup)
    orgs = ORGS
    seqs = {a: {"sequence": s, "entry_name": a, "genes": ""} for a, s in UNIPROT.items()}
    with pytest.raises(ValueError, match="two slices"):
        bcs.build(rows, orgs, seqs, {}, "2026-10-22")


def test_clusters_tsv_and_split_members_replace_the_proxy(inputs, tmp_path):
    clusters = tmp_path / "clusters.tsv"
    members = tmp_path / "members_split.tsv"
    accs = [
        "P00001",
        "P00002",
        "P00003",
        "P00004",
        "P00005",
        "P00007",
        "jgi:Neucr2:44",
        "jgi:Neucr2:98",
    ]
    clusters.write_text(
        "accession\tcluster_id\n" + "".join(f"{a}\tKS000{i}\n" for i, a in enumerate(accs))
    )
    members.write_text(
        "accession\tseed_class\tleakage_suggested\n"
        + "".join(f"{a}\tseed-like\tin_reference\n" for a in accs[:3])
        + "".join(f"{a}\tseed-independent\tnone\n" for a in accs[3:])
    )
    out = tmp_path / "set"
    run(inputs, out, "--clusters-tsv", str(clusters), "--split-members", str(members))
    c = {r["accession"]: r for r in tsv(out / "controls.tsv")}
    assert c["P00001"]["cluster_id"] == "KS0000" and c["P00001"]["leakage"] == "in_reference"
    assert c["P00004"]["seed_class"] == "seed-independent" and c["P00004"]["leakage"] == "none"
    assert "MMseqs2" in (out / "counts.md").read_text()
    clusters.write_text("accession\tcluster_id\nP00001\tKS0000\n")
    with pytest.raises(ValueError, match="not in the cluster file"):
        run(inputs, tmp_path / "set2", "--clusters-tsv", str(clusters))


def test_a_family_representative_without_a_sequence_is_an_error(inputs):
    seqs = {
        a: {"sequence": s, "entry_name": a, "genes": ""}
        for a, s in UNIPROT.items()
        if a != "P00006"
    }
    with pytest.raises(ValueError, match="family representative"):
        bcs.build(ROWS, ORGS, seqs, {}, "2026-10-22")


# ----------------------------------------------------------------------- the committed set


COMMITTED = ROOT / "data/controls/wall-fusion-kiss2019"
MAPPING = ROOT / "analysis/wall_fusion_truth/results/row_mapping.tsv"


def test_committed_set_passes_the_l1_validator():
    from cellsurface_sorting_hat.controls import validate_controls

    rep = validate_controls(COMMITTED)
    assert rep.ok, [str(d) for d in rep.defects]


def test_committed_set_is_what_the_script_makes(tmp_path):
    assert bcs.main(["--out", str(tmp_path)]) == 0
    for f in sorted(COMMITTED.iterdir()):
        assert f.read_bytes() == (tmp_path / f.name).read_bytes(), f.name


def test_committed_counts_are_recomputed_from_the_mapping():
    rows = tsv(MAPPING)
    strict = [r for r in rows if r["map_status"] in bcs.STRICT and r["control_role"] == "kiss_row"]
    # two Kiss rows of one protein are one control row
    keys = {
        r["accession"] or "jgi:" + r["jgi_portal"] + ":" + r["protein_id_jgi"].split("_")[0]
        for r in strict
    }
    controls = tsv(COMMITTED / "controls.tsv")
    local = tsv(COMMITTED / "controls_local_sequence.tsv")
    assert len(controls) + len(local) == len(keys)
    # a row without an accession gets a local key, unless one cached UniProt entry has its sequence
    local_keys = {r["accession"] for r in local}
    assert local_keys <= {k for k in keys if k.startswith("jgi:")}
    resolved = [r for r in controls if "identical sequence" in r["notes"]]
    assert len(local_keys) + len(resolved) == len({k for k in keys if k.startswith("jgi:")})
    assert len(tsv(COMMITTED / "family_reps.tsv")) == sum(
        r["control_role"] == "family_rep" for r in rows
    )
    rejected = tsv(COMMITTED / "rejected.tsv")
    assert len(rejected) == sum(
        r["control_role"] != "family_rep" and r["map_status"] not in bcs.STRICT for r in rows
    )
    # every strictly mapped row is in a control table: nothing is dropped, nothing is padded
    got = {r["accession"] for r in controls} | local_keys
    assert len(got) == len(keys) and {k for k in keys if not k.startswith("jgi:")} <= got
    # the owner's three name_conflict rows are not controls
    assert {"hyr1", "bhp2", "bhp3"} <= {r["gene_or_family_id"] for r in rejected}
    assert not {"hyr1", "bhp2", "bhp3"} & {
        g for r in controls + local for g in r["gene"].split("/")
    }


def test_committed_set_has_no_jgi_only_sequence():
    """A sequence in sequences.faa must be a UniProt one: its sha256 is the row's sequence_sha256."""
    seqs = {}
    cur = None
    for line in (COMMITTED / "sequences.faa").read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:]
            seqs[cur] = ""
        else:
            seqs[cur] += line
    for r in tsv(COMMITTED / "controls.tsv"):
        assert r["sequence_source"] == "uniprot"
        assert sha(seqs[r["accession"]]) == r["sequence_sha256"]
    local = {r["accession"] for r in tsv(COMMITTED / "controls_local_sequence.tsv")}
    assert not local & set(seqs)
