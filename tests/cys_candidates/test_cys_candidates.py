"""Tests for analysis/cys_candidates/cys_candidates.py.

Synthetic fixtures run everywhere. The real-data tests skip when the files are absent.
"""

import gzip
import json
import re
from pathlib import Path

import cys_candidates as cc
import pytest

SP20 = "MKLLVLSLLAAAVSAQAGSA"  # 20 residues, no Cys
CS_END = 19  # "CS pos: 19-20." -> mature starts at residue 20
SP_FILE_HEADER = (
    "# SignalP-6.0\tOrganism: Eukarya\tTimestamp: 20260930100126\n"
    "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position\n"
)
PARAMS = {"max_mature_len": 300, "min_cys_window": 8, "window": 60, "missing_models": []}


def cys_rich(n_cys, length):
    """A sequence of `length` residues with n_cys Cys packed at the start (no Cys elsewhere)."""
    assert length >= n_cys * 2
    core = "".join("CA" for _ in range(n_cys))
    return core + "G" * (length - len(core))


def make_fasta(path, records, gz=False):
    text = "".join(f">{h}\n{s}\n" for h, s in records)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        Path(path).write_text(text)
    return path


def make_signalp(path, rows, gz=False):
    """rows: (header, call, sp_prob, cs_text)"""
    lines = [SP_FILE_HEADER.rstrip("\n")]
    for h, call, prob, cs in rows:
        other = 1 - prob
        lines.append(f"{h}\t{call}\t{other:.6f}\t{prob:.6f}\t{cs}")
    text = "\n".join(lines) + "\n"
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        Path(path).write_text(text)
    return path


def make_domtbl(path, hits):
    """hits: (target, pfam_acc_with_version, name, i_evalue)"""
    lines = ["# target name accession tlen query name accession qlen ..."]
    for target, acc, name, ev in hits:
        f = [target, "-", "200", name, acc, "100", "1e-10", "50", "0.1", "1", "1", "1e-10", str(ev)]
        f += ["50", "0.1", "1", "100", "1", "100", "1", "100", "0.9", "desc"]
        lines.append(" ".join(f))
    Path(path).write_text("\n".join(lines) + "\n")
    return path


def write_provenance(path, drop=None):
    models = [
        ("PF05730.17", "CFEM"),
        ("PF04681.18", "Bys1"),
        ("PF01185.24", "Hydrophobin"),
        ("PF06766.17", "Hydrophobin_2"),
        ("PF28987.1", "DewD"),
        ("PF28404.1", "ARB_05178"),
    ]
    prov = {
        "pfam_hmm_path": "/db/current/Pfam-A.hmm",
        "pfam_hmm_resolved": "/db/2026-01-27-Pfam38.2/Pfam-A.hmm",
        "pfam_release_dir": "2026-01-27-Pfam38.2",
        "known_families_sha256": "abc123",
        "models": [{"acc": a, "name": n} for a, n in models if a != drop],
    }
    Path(path).write_text(json.dumps(prov))
    return path


def rows_for(tmp_path, recs, sp_rows, dom_hits=None, params=None):
    fa = make_fasta(tmp_path / "p.fa", recs)
    sp = cc.read_signalp(make_signalp(tmp_path / "sp.txt", sp_rows))
    dom = cc.read_domtbl(make_domtbl(tmp_path / "d.txt", dom_hits)) if dom_hits else {}
    rows = cc.build_rows("P", cc.read_fasta(fa), sp, dom, params or PARAMS)
    return {r["protein_id"]: r for r in rows}


def sp_protein(n_cys=12, mature_len=100):
    return SP20[:CS_END] + cys_rich(n_cys, mature_len)


# ------------------------------------------------------------------ SignalP parsing


def test_parse_cs_position():
    assert cc.parse_cs("CS pos: 19-20. Pr: 0.9482") == 19
    assert cc.parse_cs("CS pos: 105-106. Pr: 0.5") == 105
    assert cc.parse_cs("") is None
    assert cc.parse_cs(None) is None


def test_read_signalp_sp_and_no_cs(tmp_path):
    p = make_signalp(
        tmp_path / "s.txt",
        [
            ("A id one", "SP", 0.9996, "CS pos: 19-20. Pr: 0.9482"),
            ("B", "OTHER", 0.0004, ""),
        ],
    )
    got = cc.read_signalp(p)
    assert got["A id one"] == ("SP", pytest.approx(0.9996), 19)
    assert got["B"] == ("OTHER", pytest.approx(0.0004), None)


def test_signalp_sp_without_cs_stops(tmp_path):
    p = make_signalp(tmp_path / "s.txt", [("A", "SP", 0.9, "")])
    with pytest.raises(cc.Stop, match="CS position"):
        cc.read_signalp(p)


def test_signalp_gz(tmp_path):
    p = make_signalp(tmp_path / "s.txt.gz", [("A", "SP", 0.9, "CS pos: 5-6. Pr: 0.9")], gz=True)
    assert cc.read_signalp(p)["A"][2] == 5


def test_signalp_requires_header(tmp_path):
    p = tmp_path / "s.txt"
    p.write_text("A\tSP\t0.1\t0.9\tCS pos: 5-6. Pr: 0.9\n")
    with pytest.raises(cc.Stop, match="header"):
        cc.read_signalp(p)


# ------------------------------------------------------------------------ features


def test_cys_window_uses_mature_not_full_sequence():
    # 4 Cys inside the signal peptide (first 19 residues), 5 Cys in the mature part.
    seq = "MCCLLCLLCAAAVSAQAGS" + cys_rich(5, 60)
    assert len(seq[:CS_END]) == 19
    ft = cc.features(seq, CS_END, 60)
    assert ft["max_cys_window"] == 5
    assert ft["cys_count"] == 5
    assert ft["mature_length"] == 60
    full = cc.features(seq, None, 60)
    assert full["cys_count"] == 9


def test_mature_removal_off_by_one():
    seq = "MKLLVLSLLAAAVSAQAGS" + "C" + "AAAA"  # Cys is the first mature residue
    assert cc.features(seq, 19, 60)["cys_count"] == 1
    assert cc.features(seq, 19, 60)["mature_length"] == 5


def test_adjacent_cc_pairs_overlapping():
    assert cc.adjacent_cc("ACCCAC") == 2
    assert cc.adjacent_cc("CACAC") == 0
    assert cc.adjacent_cc("CC") == 1
    assert cc.adjacent_cc("C") == 0


def test_c_x_c():
    assert cc.c_x_c("CACAC") == 2
    assert cc.c_x_c("CCC") == 1
    assert cc.c_x_c("CAAC") == 0


def test_max_cys_window_boundary():
    w = 10
    inside = "C" + "A" * 8 + "C" + "G" * 20  # two Cys spanning exactly 10 residues
    assert cc.max_cys_window(inside, w) == (2, 1)
    outside = "C" + "A" * 9 + "C" + "G" * 20  # two Cys spanning 11 residues
    assert cc.max_cys_window(outside, w)[0] == 1


def test_max_cys_window_start_and_short_sequence():
    assert cc.max_cys_window("AAACCCAA", 3) == (3, 4)
    assert cc.max_cys_window("CAC", 60) == (2, 1)
    assert cc.max_cys_window("", 60) == (0, 1)


def test_pest_fraction():
    assert cc.pest_fraction("PSTE") == 1.0
    assert cc.pest_fraction("PSAA") == 0.5
    assert cc.pest_fraction("") == 0.0


# ---------------------------------------------------------------------------- tiers


def test_sp_cys_rich_is_unassigned(tmp_path):
    rows = rows_for(
        tmp_path,
        [("sp12", sp_protein(12, 100))],
        [("sp12", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")],
    )
    r = rows["sp12"]
    assert r["tier"] == "cys_rich_sp_unassigned"
    assert r["mature_length"] == 100
    assert r["cys_count"] == 12


def test_same_sequence_without_sp_is_no_sp(tmp_path):
    rows = rows_for(
        tmp_path,
        [("nosp", sp_protein(12, 100))],
        [("nosp", "OTHER", 0.01, "")],
    )
    assert rows["nosp"]["tier"] == "cys_rich_no_sp"
    assert rows["nosp"]["cs_end"] == "NA"
    assert rows["nosp"]["mature_length"] == 119  # full length: no SP removed


def test_cfem_hit_is_known_family(tmp_path):
    rows = rows_for(
        tmp_path,
        [("c1", sp_protein())],
        [("c1", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")],
        [("c1", "PF05730.20", "CFEM", "1e-20")],
    )
    r = rows["c1"]
    assert r["tier"] == "cys_rich_sp_known_family"
    assert r["cfem"] == 1
    assert r["cfem_evalue"] == "1e-20"


@pytest.mark.parametrize(
    "acc,name,flag",
    [
        ("PF04681.18", "Bys1", "bys1"),
        ("PF01185.24", "Hydrophobin", "hydrophobin"),
        ("PF06766.17", "Hydrophobin_2", "hydrophobin"),
    ],
)
def test_other_known_families_excluded(tmp_path, acc, name, flag):
    rows = rows_for(
        tmp_path,
        [("k", sp_protein())],
        [("k", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")],
        [("k", acc, name, "1e-9")],
    )
    assert rows["k"]["tier"] == "cys_rich_sp_known_family"
    assert rows["k"][flag] == 1


def test_pra3_like_family_stays_in_candidate_tier_and_is_flagged(tmp_path):
    rows = rows_for(
        tmp_path,
        [("p", sp_protein())],
        [("p", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")],
        [("p", "PF28404.1", "ARB_05178", "1e-30")],
    )
    assert rows["p"]["tier"] == "cys_rich_sp_unassigned"
    assert rows["p"]["pra3_like_family"] == 1
    assert rows["p"]["cfem"] == 0


def test_known_family_wins_over_pra3_like(tmp_path):
    rows = rows_for(
        tmp_path,
        [("p", sp_protein())],
        [("p", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")],
        [("p", "PF28404.1", "ARB_05178", "1e-30"), ("p", "PF05730.20", "CFEM", "1e-5")],
    )
    assert rows["p"]["tier"] == "cys_rich_sp_known_family"
    assert rows["p"]["pra3_like_family"] == 1


def test_no_sp_with_cfem_stays_no_sp_and_flagged(tmp_path):
    rows = rows_for(
        tmp_path,
        [("n", sp_protein())],
        [("n", "OTHER", 0.01, "")],
        [("n", "PF05730.20", "CFEM", "1e-5")],
    )
    assert rows["n"]["tier"] == "cys_rich_no_sp"
    assert rows["n"]["cfem"] == 1


def test_long_protein_is_other(tmp_path):
    rows = rows_for(
        tmp_path,
        [("long", sp_protein(12, 400))],
        [("long", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")],
    )
    assert rows["long"]["tier"] == "other"


def test_low_cys_is_other(tmp_path):
    rows = rows_for(
        tmp_path,
        [("low", sp_protein(3, 100))],
        [("low", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")],
    )
    assert rows["low"]["tier"] == "other"


def test_length_boundary_is_inclusive():
    f = cc.assign_tier
    assert f(True, 300, 8, set(), 300, 8) == "cys_rich_sp_unassigned"
    assert f(True, 301, 8, set(), 300, 8) == "other"
    assert f(False, 300, 8, set(), 300, 8) == "cys_rich_no_sp"
    assert f(False, 301, 8, set(), 300, 8) == "other"


def test_cys_boundary_is_inclusive():
    f = cc.assign_tier
    assert f(True, 100, 8, set(), 300, 8) == "cys_rich_sp_unassigned"
    assert f(True, 100, 7, set(), 300, 8) == "other"


def test_window_boundary_changes_tier(tmp_path):
    # 8 Cys spread over exactly 60 residues pass at W=60 and fail at W=59.
    mature = "".join("C" + "A" * 7 for _ in range(7)) + "C" + "A" * 3  # Cys at 0,8,...,56
    assert mature.count("C") == 8
    assert cc.max_cys_window(mature, 57)[0] == 8
    assert cc.max_cys_window(mature, 56)[0] == 7
    p57 = dict(PARAMS, window=57)
    p56 = dict(PARAMS, window=56)
    seq = SP20[:CS_END] + mature + "G" * 20
    sp_rows = [("w", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")]
    assert rows_for(tmp_path, [("w", seq)], sp_rows, params=p57)["w"]["tier"] != "other"
    assert rows_for(tmp_path, [("w", seq)], sp_rows, params=p56)["w"]["tier"] == "other"


def test_tier_does_not_depend_on_id(tmp_path):
    seq = sp_protein()
    sp = [("Q2TVJ9", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")]
    a = rows_for(tmp_path, [("Q2TVJ9", seq)], sp)["Q2TVJ9"]["tier"]
    sp2 = [("zzz", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")]
    b = rows_for(tmp_path, [("zzz", seq)], sp2)["zzz"]["tier"]
    assert a == b == "cys_rich_sp_unassigned"


def test_missing_models_report_na(tmp_path):
    params = dict(PARAMS, missing_models=["pra3_like_family"])
    rows = rows_for(
        tmp_path,
        [("m", sp_protein())],
        [("m", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")],
        params=params,
        dom_hits=[("m", "PF05730.20", "CFEM", "1e-5")],
    )
    assert rows["m"]["pra3_like_family"] == "NA"
    assert rows["m"]["cfem"] == 1


# ------------------------------------------------------------------------ STOP cases


def test_duplicate_fasta_id_stops(tmp_path):
    fa = make_fasta(tmp_path / "d.fa", [("a x", "MAAA"), ("a y", "MCCC")])
    with pytest.raises(cc.Stop, match="duplicate"):
        cc.read_fasta(fa)


def test_missing_fasta_stops(tmp_path):
    with pytest.raises(cc.Stop, match="FASTA not found"):
        cc.read_fasta(tmp_path / "nope.fa")


def test_signalp_id_not_in_fasta_stops(tmp_path):
    recs = [("a", "MAAAA")]
    sp = {"a": ("OTHER", 0.1, None), "ghost": ("OTHER", 0.1, None)}
    with pytest.raises(cc.Stop, match="not in the FASTA"):
        cc.build_rows("P", recs, sp, {}, PARAMS)


def test_fasta_id_without_signalp_row_stops():
    recs = [("a", "MAAAA"), ("b", "MAAAA")]
    sp = {"a": ("OTHER", 0.1, None)}
    with pytest.raises(cc.Stop, match="without a SignalP row"):
        cc.build_rows("P", recs, sp, {}, PARAMS)


def test_domtbl_target_not_in_fasta_stops():
    recs = [("a", "MAAAA")]
    sp = {"a": ("OTHER", 0.1, None)}
    with pytest.raises(cc.Stop, match="domtbl target"):
        cc.build_rows("P", recs, sp, {"ghost": {"cfem": 1e-5}}, PARAMS)


def test_unknown_model_in_domtbl_stops(tmp_path):
    p = make_domtbl(tmp_path / "d.txt", [("a", "PF00001.1", "7tm_1", "1e-5")])
    with pytest.raises(cc.Stop, match="unexpected model"):
        cc.read_domtbl(p)


def test_cs_beyond_sequence_stops():
    recs = [("a", "MAAAA")]
    sp = {"a": ("SP", 0.9, 50)}
    with pytest.raises(cc.Stop, match="beyond"):
        cc.build_rows("P", recs, sp, {}, PARAMS)


# ------------------------------------------------------------------ end to end (main)


def write_inputs(tmp_path, gz=False):
    recs = [
        ("sp12 desc one", sp_protein(12, 100)),
        ("nosp12", sp_protein(12, 100)),
        ("cfem1", sp_protein(12, 100)),
        ("long1", sp_protein(12, 400)),
    ]
    fa = make_fasta(tmp_path / ("p.fa.gz" if gz else "p.fa"), recs, gz=gz)
    sp = make_signalp(
        tmp_path / ("sp.txt.gz" if gz else "sp.txt"),
        [
            ("sp12 desc one", "SP", 0.99, "CS pos: 19-20. Pr: 0.9"),
            ("nosp12", "OTHER", 0.01, ""),
            ("cfem1", "SP", 0.99, "CS pos: 19-20. Pr: 0.9"),
            ("long1", "SP", 0.99, "CS pos: 19-20. Pr: 0.9"),
        ],
        gz=gz,
    )
    dom = make_domtbl(tmp_path / "d.txt", [("cfem1", "PF05730.20", "CFEM", "1e-9")])
    write_provenance(tmp_path / "prov.json")
    return fa, sp, dom


def read_tsv_gz(path):
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh]


@pytest.mark.parametrize("gz", [False, True])
def test_main_end_to_end(tmp_path, gz):
    fa, sp, dom = write_inputs(tmp_path, gz=gz)
    out = tmp_path / "out"
    rc = cc.main(
        [
            "--name",
            "P",
            "--fasta",
            str(fa),
            "--signalp",
            str(sp),
            "--domtbl",
            str(dom),
            "--pfam-provenance",
            str(tmp_path / "prov.json"),
            "--out-dir",
            str(out),
            "--repo",
            str(tmp_path),
        ]  # fmt: skip
    )
    assert rc == 0
    assert sorted(p.name for p in out.iterdir()) == [
        "P.tsv.gz",
        "candidates.tsv.gz",
        "run.json",
        "summary.tsv",  # fmt: skip
    ]
    tiers = {r["protein_id"]: r["tier"] for r in read_tsv_gz(out / "P.tsv.gz")}
    assert tiers == {
        "sp12": "cys_rich_sp_unassigned",
        "nosp12": "cys_rich_no_sp",
        "cfem1": "cys_rich_sp_known_family",
        "long1": "other",
    }
    cand = {r["protein_id"] for r in read_tsv_gz(out / "candidates.tsv.gz")}
    assert cand == {"sp12", "nosp12", "cfem1"}
    summ = (out / "summary.tsv").read_text().splitlines()
    assert summ[0].split("\t")[:3] == ["proteome", "n_proteins", "n_sp"]
    assert summ[1].split("\t") == ["P", "4", "3", "1", "1", "1", "1"]
    meta = json.loads((out / "run.json").read_text())
    assert meta["library"] == "none"
    assert meta["pfam"]["pfam_release_dir"] == "2026-01-27-Pfam38.2"
    assert meta["pfam"]["known_families_sha256"] == "abc123"
    assert {m["acc"] for m in meta["pfam"]["models"]} >= {"PF28404.1", "PF28987.1"}
    assert meta["parameters"]["max_mature_len"] == 300
    assert meta["inputs"]["P"]["fasta"]["sha256"] == cc.sha256_of(fa)
    assert not any("time" in k.lower() or "date" in k.lower() for k in meta)


def test_output_is_deterministic(tmp_path):
    fa, sp, dom = write_inputs(tmp_path)
    args = ["--name", "P", "--fasta", str(fa), "--signalp", str(sp), "--domtbl", str(dom),
            "--pfam-provenance", str(tmp_path / "prov.json"), "--repo", str(tmp_path)]  # fmt: skip
    assert cc.main(args + ["--out-dir", str(tmp_path / "o1")]) == 0
    assert cc.main(args + ["--out-dir", str(tmp_path / "o2")]) == 0
    for name in ("P.tsv.gz", "candidates.tsv.gz", "summary.tsv"):
        assert (tmp_path / "o1" / name).read_bytes() == (tmp_path / "o2" / name).read_bytes()


def test_thresholds_are_options(tmp_path):
    fa, sp, dom = write_inputs(tmp_path)
    out = tmp_path / "out"
    rc = cc.main(
        [
            "--name",
            "P",
            "--fasta",
            str(fa),
            "--signalp",
            str(sp),
            "--out-dir",
            str(out),
            "--max-mature-len",
            "99",
            "--repo",
            str(tmp_path),
        ]  # fmt: skip
    )
    assert rc == 0
    tiers = {r["protein_id"]: r["tier"] for r in read_tsv_gz(out / "P.tsv.gz")}
    assert tiers["sp12"] == "other"  # mature length 100 > 99
    assert tiers["nosp12"] == "other"  # full length 119 > 99


def test_stop_leaves_no_outputs(tmp_path, capsys):
    fa, sp, _ = write_inputs(tmp_path)
    out = tmp_path / "out"
    rc = cc.main(
        [
            "--name",
            "P",
            "--fasta",
            str(tmp_path / "missing.fa"),
            "--signalp",
            str(sp),
            "--out-dir",
            str(out),
        ]  # fmt: skip
    )
    assert rc == 2
    assert capsys.readouterr().err.startswith("STOP: ")
    assert not out.exists() or list(out.iterdir()) == []


def test_stop_in_second_proteome_leaves_no_outputs(tmp_path, capsys):
    fa, sp, _ = write_inputs(tmp_path)
    bad_fa = make_fasta(tmp_path / "bad.fa", [("only", "MAAA")])
    man = tmp_path / "m.tsv"
    man.write_text(f"A\t{fa}\t{sp}\nB\t{bad_fa}\t{sp}\n")
    out = tmp_path / "out"
    rc = cc.main(["--manifest", str(man), "--out-dir", str(out)])
    assert rc == 2
    assert "STOP:" in capsys.readouterr().err
    assert not out.exists() or list(out.iterdir()) == []


def test_unknown_missing_model_stops(tmp_path, capsys):
    fa, sp, _ = write_inputs(tmp_path)
    rc = cc.main(
        [
            "--name",
            "P",
            "--fasta",
            str(fa),
            "--signalp",
            str(sp),
            "--out-dir",
            str(tmp_path / "o"),
            "--missing-models",
            "nonsense",
        ]  # fmt: skip
    )
    assert rc == 2


# ------------------------------------------------------------- source-code rules


def test_source_has_no_bash_source_and_no_hardcoded_ids():
    for f in Path(cc.__file__).parent.iterdir():
        if f.suffix in (".py", ".sh", ".md"):
            assert "BASH" + "_SOURCE" not in f.read_text(), f
    src = Path(cc.__file__).read_text()
    assert not re.search(r"CIMG_\d|CPOS\d|CIB\d|QVM\d|Q2TVJ9|PRA3", src)


# ---------------------------------------------------------------- real-data tests

PRIMARY = Path("/bigdata/stajichlab/jstajich/projects/adhesionPred")
RS_FASTA = Path(
    "/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Pangenome/"
    "input_run2/CimmitisRS_FungiDB.fasta"
)
RS_SP = PRIMARY / "analysis/cocci_repeats/signalp/CimmitisRS_FungiDB/prediction_results.txt"
RS_DOM = PRIMARY / "_workdir/cys_candidates/CimmitisRS_FungiDB.domtbl.gz"
RS_PROV = PRIMARY / "_workdir/cys_candidates/pfam_provenance.json"
FASTA_DIR = PRIMARY / "_workdir/class2b_structure/fasta"
needs_rs = pytest.mark.skipif(
    not (RS_FASTA.is_file() and RS_SP.is_file() and RS_DOM.is_file() and RS_PROV.is_file()),
    reason="real C. immitis RS files absent",
)


@pytest.fixture(scope="module")
def rs_rows(tmp_path_factory):
    if not (RS_FASTA.is_file() and RS_SP.is_file() and RS_DOM.is_file() and RS_PROV.is_file()):
        pytest.skip("real C. immitis RS files absent")
    out = tmp_path_factory.mktemp("rs")
    rc = cc.main(
        [
            "--name",
            "RS",
            "--fasta",
            str(RS_FASTA),
            "--signalp",
            str(RS_SP),
            "--domtbl",
            str(RS_DOM),
            "--out-dir",
            str(out),
            "--pfam-provenance",
            str(RS_PROV),
        ]  # fmt: skip
    )
    assert rc == 0
    return {r["protein_id"]: r for r in read_tsv_gz(out / "RS.tsv.gz")}


@needs_rs
def test_real_pra3_ortholog_is_unassigned_candidate(rs_rows):
    r = rs_rows["CIMG_02492-t26_1-p1"]
    assert r["tier"] == "cys_rich_sp_unassigned"
    assert (r["length"], r["mature_length"], r["cys_count"]) == ("220", "201", "12")


@needs_rs
def test_real_cfem_controls_are_known_family(rs_rows):
    # CIMG_09696 = Ag2/PRA (Q6QJA6, fident 1.000 over 194 aa, measured with mmseqs),
    # CIMG_09560 = PRA2 (Q6K1L8, fident 0.983 over 124 aa).
    assert rs_rows["CIMG_09696-t26_1-p1"]["tier"] == "cys_rich_sp_known_family"
    assert rs_rows["CIMG_09560-t26_1-p1"]["tier"] == "cys_rich_sp_known_family"
    for r in rs_rows.values():
        if r["tier"] == "cys_rich_sp_unassigned":
            assert (r["cfem"], r["bys1"], r["hydrophobin"]) == ("0", "0", "0")


@pytest.mark.skipif(not (FASTA_DIR / "PRA3_Q2TVJ9.fasta").is_file(), reason="PRA3 fasta absent")
def test_real_uniprot_pra3_without_sp_is_no_sp(tmp_path):
    seq = cc.read_fasta(FASTA_DIR / "PRA3_Q2TVJ9.fasta")[0][1]
    assert len(seq) == 153 and seq.count("C") == 8
    rows = rows_for(
        tmp_path,
        [("Q2TVJ9", seq)],
        [("Q2TVJ9", "OTHER", 0.001, "")],
    )
    assert rows["Q2TVJ9"]["tier"] == "cys_rich_no_sp"


def test_write_failure_removes_temp_files(tmp_path, monkeypatch):
    fa, sp, _ = write_inputs(tmp_path)
    calls = {"n": 0}
    real = cc.write_gz

    def flaky(path, text):
        calls["n"] += 1
        if calls["n"] == 2:
            raise OSError("disk full")
        real(path, text)

    monkeypatch.setattr(cc, "write_gz", flaky)
    out = tmp_path / "out"
    with pytest.raises(OSError):
        cc.run([("P", str(fa), str(sp), None)], out, PARAMS, tmp_path)
    assert list(out.iterdir()) == []


# ------------------------------------------------- fix round 1: additional tests


def test_tier_uses_mature_length_not_full_length(tmp_path):
    # SP of 15 residues; full length L+5 = 305, mature length L-10 = 290.
    seq = "MKLLVLSLLAAAVSA" + cys_rich(12, 290)
    assert len(seq) == 305
    rows = rows_for(
        tmp_path,
        [("m", seq)],
        [("m", "SP", 0.99, "CS pos: 15-16. Pr: 0.9")],
    )
    assert rows["m"]["length"] == 305
    assert rows["m"]["mature_length"] == 290
    assert rows["m"]["tier"] == "cys_rich_sp_unassigned"


def test_dewd_counts_as_hydrophobin(tmp_path):
    rows = rows_for(
        tmp_path,
        [("d", sp_protein())],
        [("d", "SP", 0.99, "CS pos: 19-20. Pr: 0.9")],
        [("d", "PF28987.1", "DewD", "1e-9")],
    )
    assert rows["d"]["hydrophobin"] == 1
    assert rows["d"]["tier"] == "cys_rich_sp_known_family"


def test_internal_stop_symbol_stops(tmp_path):
    fa = make_fasta(tmp_path / "s.fa", [("p1 x", "MAAA*CCC")])
    with pytest.raises(cc.Stop, match="p1"):
        cc.read_fasta(fa)


@pytest.mark.parametrize("bad", ["MAA-CC", "MAA CC", "MAA1CC"])
def test_non_letter_stops(tmp_path, bad):
    fa = tmp_path / "s.fa"
    fa.write_text(f">p1\n{bad}\n")
    with pytest.raises(cc.Stop, match="non-letter"):
        cc.read_fasta(fa)


def test_empty_sequence_stops(tmp_path):
    fa = tmp_path / "s.fa"
    fa.write_text(">p1\n>p2\nMAAA\n")
    with pytest.raises(cc.Stop, match="empty sequence"):
        cc.read_fasta(fa)
    fa.write_text(">p1\n*\n")
    with pytest.raises(cc.Stop, match="empty sequence"):
        cc.read_fasta(fa)


def test_lowercase_and_single_trailing_stop_are_normalised(tmp_path):
    fa = tmp_path / "s.fa"
    fa.write_text(">p1\nmaacc\ncc*\n")
    assert cc.read_fasta(fa) == [("p1", "MAACCCC")]


def test_provenance_required_with_domtbl(tmp_path, capsys):
    fa, sp, dom = write_inputs(tmp_path)
    rc = cc.main(
        [
            "--name",
            "P",
            "--fasta",
            str(fa),
            "--signalp",
            str(sp),
            "--domtbl",
            str(dom),
            "--out-dir",
            str(tmp_path / "o"),
        ]  # fmt: skip
    )
    assert rc == 2
    assert "pfam_provenance" in capsys.readouterr().err
    assert not (tmp_path / "o").exists()


def test_provenance_missing_model_stops(tmp_path, capsys):
    fa, sp, dom = write_inputs(tmp_path)
    write_provenance(tmp_path / "prov.json", drop="PF28404.1")
    args = ["--name", "P", "--fasta", str(fa), "--signalp", str(sp), "--domtbl", str(dom),
            "--pfam-provenance", str(tmp_path / "prov.json"),
            "--out-dir", str(tmp_path / "o")]  # fmt: skip
    assert cc.main(args) == 2
    assert "PF28404" in capsys.readouterr().err
    # the same run is allowed when the family is declared as not searched
    assert cc.main(args + ["--missing-models", "pra3_like_family"]) == 0


def test_stale_per_proteome_files_stop(tmp_path, capsys):
    fa, sp, dom = write_inputs(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    (out / "OLD.tsv.gz").write_bytes(b"x")
    rc = cc.main(["--name", "P", "--fasta", str(fa), "--signalp", str(sp), "--out-dir", str(out)])
    err = capsys.readouterr().err
    assert rc == 2
    assert "STOP:" in err and "OLD.tsv.gz" in err
    assert sorted(p.name for p in out.iterdir()) == ["OLD.tsv.gz"]


def test_rerun_same_names_is_allowed(tmp_path):
    fa, sp, dom = write_inputs(tmp_path)
    args = ["--name", "P", "--fasta", str(fa), "--signalp", str(sp), "--out-dir",
            str(tmp_path / "o")]  # fmt: skip
    assert cc.main(args) == 0
    assert cc.main(args) == 0


def test_argparse_error_is_stop(capsys):
    with pytest.raises(SystemExit) as e:
        cc.main(["--max-mature-len", "abc", "--out-dir", "x"])
    assert e.value.code == 2
    assert capsys.readouterr().err.startswith("STOP: ")
    with pytest.raises(SystemExit) as e:
        cc.main([])
    assert e.value.code == 2


def test_missing_manifest_stops(tmp_path, capsys):
    rc = cc.main(["--manifest", str(tmp_path / "nope.tsv"), "--out-dir", str(tmp_path / "o")])
    assert rc == 2
    assert "STOP: manifest not found" in capsys.readouterr().err


def test_no_input_arguments_stops(tmp_path, capsys):
    assert cc.main(["--out-dir", str(tmp_path / "o")]) == 2
    assert capsys.readouterr().err.startswith("STOP: ")


def test_unwritable_out_dir_stops(tmp_path, capsys):
    fa, sp, _ = write_inputs(tmp_path)
    blocker = tmp_path / "afile"
    blocker.write_text("x")
    rc = cc.main(
        [
            "--name",
            "P",
            "--fasta",
            str(fa),
            "--signalp",
            str(sp),
            "--out-dir",
            str(blocker / "sub"),
        ]  # fmt: skip
    )
    assert rc == 2
    assert "STOP: cannot write" in capsys.readouterr().err


# ------------------------------------------------------- 01_known_family_hmm.sh

SCRIPT_01 = Path(cc.__file__).parent / "01_known_family_hmm.sh"


def run_01(tmp_path, env_extra, with_out_dir=True):
    import os
    import subprocess

    env = {k: v for k, v in os.environ.items() if k not in ("OUT_DIR", "PFAM_HMM", "MANIFEST")}
    if with_out_dir:
        env["OUT_DIR"] = str(tmp_path / "out")
    env.update(env_extra)
    return subprocess.run(["bash", "-l", str(SCRIPT_01)], capture_output=True, text=True, env=env)


def test_01_stops_without_out_dir(tmp_path):
    res = run_01(tmp_path, {}, with_out_dir=False)
    assert res.returncode == 2
    assert res.stderr.startswith("STOP: OUT_DIR")


def test_01_stops_on_unreadable_database_and_leaves_nothing(tmp_path):
    import shutil

    probe = shutil.which("bash")
    assert probe
    bad = tmp_path / "bad.hmm"
    bad.write_text(
        "garbage\n"
        + "".join(
            f"ACC   {a}\n"
            for a in (
                "PF05730.17",
                "PF04681.18",
                "PF01185.24",
                "PF06766.17",
                "PF28987.1",
                "PF28404.1",
            )
        )
    )
    man = tmp_path / "m.tsv"
    man.write_text("A\t/dev/null\n")
    res = run_01(tmp_path, {"PFAM_HMM": str(bad), "MANIFEST": str(man)})
    if "cannot load hmmer" in res.stderr:
        pytest.skip("hmmer module not available")
    assert res.returncode == 2
    assert "STOP:" in res.stderr
    out = tmp_path / "out"
    assert not out.exists() or [p.name for p in out.iterdir()] == []


def test_01_stops_when_model_missing(tmp_path):
    db = tmp_path / "db.hmm"
    db.write_text("ACC   PF05730.17\n")
    man = tmp_path / "m.tsv"
    man.write_text("A\t/dev/null\n")
    res = run_01(tmp_path, {"PFAM_HMM": str(db), "MANIFEST": str(man)})
    if "cannot load hmmer" in res.stderr:
        pytest.skip("hmmer module not available")
    assert res.returncode == 2
    assert "STOP:" in res.stderr and "PF04681" in res.stderr
