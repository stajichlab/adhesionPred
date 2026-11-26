"""Fixture tests for round 2 of analysis/crossclade_domain_search (offline: no network, no HPCC)."""

import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DIR = ROOT / "analysis/crossclade_domain_search"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ri = _load("round2_inputs")
rs = _load("round2_summarize")
gm = _load("gene_model_check")


def _pdb(n, plddt=80.0, low=()):
    """A minimal PDB with n residues (N and CA atoms); residues in `low` get pLDDT 30."""
    lines = []
    for i in range(1, n + 1):
        b = 30.0 if i in low else plddt
        for atom in ("N", "CA"):
            lines.append(
                f"ATOM  {len(lines) + 1:5d}  {atom:<3s} ALA A{i:4d}    {0:8.3f}{0:8.3f}{0:8.3f}  1.00{b:6.2f}           C"
            )
    return "\n".join(lines + ["END"]) + "\n"


def _write_tsv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


# round2_inputs ---------------------------------------------------------------------------------------------------------


def test_crop_trim_and_plddt():
    pdb = _pdb(10, 70.0)
    assert ri.pdb_residues(pdb) == list(range(1, 11))
    assert ri.mean_plddt(pdb) == 70.0
    c = ri.crop_pdb(pdb, 3, 5)
    assert ri.pdb_residues(c) == [3, 4, 5] and c.endswith("TER\nEND\n")
    with pytest.raises(ri.BuildError):
        ri.crop_pdb(pdb, 20, 30)
    t, kept = ri.trim_low_plddt(_pdb(10, 80.0, low={1, 2, 9}), 50)
    assert kept == 7 and ri.pdb_residues(t) == [3, 4, 5, 6, 7, 8, 10]


def test_background_order_is_seeded_whole_proteome():
    targets = [{"accession": f"A{i}", "species": "S1"} for i in range(20)] + [
        {"accession": f"B{i}", "species": "S2"} for i in range(20)
    ]
    one = ri.background_order(targets, {"A0", "B0"}, "x")
    assert one == ri.background_order(targets, {"A0", "B0"}, "x")
    assert one != ri.background_order(targets, {"A0", "B0"}, "y")
    assert len(one["S1"]) == 19 and "A0" not in one["S1"]


def test_plan_scored_roles():
    targets = [{"accession": a, "species": "S1"} for a in ("T1", "T2", "P1", "Q1")]
    hits = [{"target": "T1"}, {"target": "T2"}, {"target": "T1"}]
    truth = [{"accession": a, "in_targets": "yes"} for a in ("P1", "T2", "Q1")]
    assert ri.plan_scored(hits, targets, truth, {"Q1"}) == {
        "P1": "pf21671",
        "T1": "candidate",
        "T2": "pf21671",
    }
    with pytest.raises(ri.BuildError):
        ri.plan_scored([{"target": "ZZ"}], targets, truth, set())


def _fake_afdb(models, legacy=True):
    """models: acc -> (sequence, n_residues_in_pdb[, low-pLDDT residues]) or None."""

    def fetch(url, binary=False):
        if url.startswith(ri.AFDB):
            acc = url.rsplit("/", 1)[1]
            if models.get(acc) is None:
                raise ri.NotFoundError(url)
            seq = models[acc][0]
            e = (
                {"entryId": f"AF-{acc}-F1", "uniprotSequence": seq, "uniprotStart": 1}
                if legacy
                else {"modelEntityId": f"AF-{acc}-F1", "sequence": seq, "sequenceStart": 1}
            )
            return json.dumps([{**e, "latestVersion": 4, "pdbUrl": f"pdb://{acc}"}])
        m = models[url.split("//")[1]]
        return _pdb(m[1], 80.0, m[2] if len(m) > 2 else ()).encode()

    return fetch


def test_fetch_model_statuses():
    f = _fake_afdb({"OK": ("AAAA", 4), "SHORT": ("AAAA", 3), "NONE": None})
    assert ri.fetch_model("OK", "AAAA", f)[0] == "ok"
    assert ri.fetch_model("OK", "AAAC", f)[0] == "sequence_mismatch"
    assert ri.fetch_model("SHORT", "AAAA", f)[0] == "sequence_mismatch"
    assert ri.fetch_model("NONE", "AAAA", f)[0] == "no_model"
    assert ri.fetch_model("OK", "AAAA", _fake_afdb({"OK": ("AAAA", 4)}, legacy=False))[0] == "ok"


def test_http_get_retries_then_fails(monkeypatch):
    calls = []

    def boom(req, timeout):
        calls.append(1)
        raise ri.urllib.error.URLError("down")

    monkeypatch.setattr(ri.urllib.request, "urlopen", boom)
    with pytest.raises(ri.BuildError):
        ri.http_get("https://example.invalid/x", retries=2, backoff=0)
    assert len(calls) == 3


def _round1(tmp_path, n_bg=6):
    r1 = tmp_path / "r1"
    r1.mkdir()
    _write_tsv(
        r1 / "query_domains.tsv",
        [
            {
                "query_id": "Q1_PF1_2_5",
                "panel": "p",
                "accession": "Q1",
                "species": "S0",
                "gene": "g",
                "crop_start": 2,
                "crop_end": 5,
                "protein_length": 8,
            }
        ],
    )
    accs = ["T1", "P1"] + [f"N{i}" for i in range(n_bg)]
    _write_tsv(
        r1 / "target_proteins.tsv", [{"accession": a, "species": "S1", "length": 60} for a in accs]
    )
    _write_tsv(
        r1 / "truth_pf21671.tsv",
        [{"accession": "P1", "species": "S1", "pfam": "PF21671", "in_targets": "yes"}],
    )
    (r1 / "targets.faa").write_text("".join(f">{a}\n{'A' * 60}\n" for a in accs))
    (r1 / "queries.faa").write_text(">Q1_PF1_2_5\nCDEF\n")
    hits = tmp_path / "hits.tsv"
    _write_tsv(hits, [{"method": "foldseek", "query_id": "Q1_PF1_2_5", "target": "T1"}])
    models = {a: ("A" * 60, 60) for a in accs}
    models["Q1"] = ("MCDEFGHI", 8)
    return r1, hits, models


def test_inputs_main(tmp_path):
    r1, hits, models = _round1(tmp_path)
    models["N0"] = None  # no model: skipped, next one drawn
    models["N1"] = ("A" * 60, 60, set(range(1, 40)))  # 21 residues left after trimming: too_short
    out = tmp_path / "out"
    args = [
        "--round1-inputs",
        str(r1),
        "--round1-hits",
        str(hits),
        "--background-per-species",
        "2",
        "--out",
        str(out),
    ]
    assert ri.main(args, fetch=_fake_afdb(models), pause=0) == 0
    assert ri.pdb_residues((out / "queries/Q1_PF1_2_5.pdb").read_text()) == [2, 3, 4, 5]
    m = {r["id"]: r for r in csv.DictReader(open(out / "models.tsv"), delimiter="\t")}
    assert m["T1"]["role"] == "candidate" and m["P1"]["role"] == "pf21671"
    # the background is drawn in seeded random order until 4 usable models; N0 and N1 are drawn only if needed
    assert m.get("N0", {"status": "no_model"})["status"] == "no_model"
    assert m.get("N1", {"status": "too_short"})["status"] == "too_short"
    halves = [r["bg_half"] for r in m.values() if r["bg_half"]]
    assert sorted(halves) == ["A", "A", "B", "B"]
    assert all(m[i]["status"] == "ok" for i, r in m.items() if r["bg_half"])
    assert (
        "models.tsv" in (out / "SHA256SUMS").read_text() and not (tmp_path / "out.partial").exists()
    )
    assert ri.main(args, pause=0) == 1  # --out exists


def test_inputs_query_crop_mismatch_and_too_few_background(tmp_path):
    r1, hits, models = _round1(
        tmp_path, n_bg=2
    )  # T1 (a candidate) can be drawn too: 3 proteins for 4
    models["Q1"] = ("MCDXFGHI", 8)
    out = tmp_path / "out"
    args = [
        "--round1-inputs",
        str(r1),
        "--round1-hits",
        str(hits),
        "--background-per-species",
        "2",
        "--out",
        str(out),
    ]
    assert ri.main(args, fetch=_fake_afdb(models), pause=0) == 1
    assert not out.exists() and (tmp_path / "out.partial").exists()
    (tmp_path / "b").mkdir()
    r1b, hitsb, models_b = _round1(tmp_path / "b", n_bg=6)
    models_b["Q1"] = ("MCDXFGHI", 8)
    out_b = tmp_path / "b" / "out"
    argsb = [
        "--round1-inputs",
        str(r1b),
        "--round1-hits",
        str(hitsb),
        "--background-per-species",
        "2",
        "--out",
        str(out_b),
    ]
    assert ri.main(argsb, fetch=_fake_afdb(models_b), pause=0) == 0
    m = {r["id"]: r for r in csv.DictReader(open(out_b / "models.tsv"), delimiter="\t")}
    assert m["Q1_PF1_2_5"]["status"] == "sequence_mismatch"


# round2_summarize ------------------------------------------------------------------------------------------------------


def test_norm_label():
    known = {"M7WU28", "J9VQE5_PF21671_101_180"}
    assert rs.norm_label("/x/targets/M7WU28.pdb", known) == "M7WU28"
    assert rs.norm_label("M7WU28", known) == "M7WU28"  # Foldseek: file stem
    assert rs.norm_label("M7WU28_A", known) == "M7WU28"  # Reseek: stem + chain
    assert rs.norm_label("J9VQE5_PF21671_101_180_A", known) == "J9VQE5_PF21671_101_180"
    assert rs.norm_label("J9VQE5_PF21671_101_180.pdb", known) == "J9VQE5_PF21671_101_180"
    with pytest.raises(rs.SummaryError):
        rs.norm_label("ZZZ.pdb", known)


def test_read_scores_scale_and_best(tmp_path):
    p = tmp_path / "r.tsv"
    p.write_text("Q_A\tT_A\t1\t2\t3\t1\t2\t3\t10\t0.01\nQ_A\tT_A\t1\t2\t3\t1\t2\t3\t10\t0.0001\n")
    assert rs.read_scores(p, rs.RS_COLS, "pvalue", {"Q", "T"}, scale=10) == {
        ("Q", "T"): pytest.approx(0.001)
    }


def _m(i, role, half="", status="ok"):
    return {
        "id": i,
        "role": role,
        "bg_half": half,
        "species": "S1",
        "mean_plddt": "80",
        "residues_kept": "100",
        "status": status,
    }


def _summary_fixture():
    q = "Q1_PF_1_9"
    domains = [{"query_id": q, "panel": "basidio", "gene": "g1", "accession": "Q1"}]
    models = [
        {**_m(q, "query"), "species": "S0"},
        _m("C1", "candidate"),
        _m("C2", "candidate"),
        _m("C3", "candidate", half="A"),  # a candidate drawn into half A: left out of its own null
        _m("P1", "pf21671"),
        _m("P2", "pf21671", status="no_model"),
        _m("A1", "background", "A"),
        _m("B1", "background", "B"),
        _m("B2", "background", "B"),
    ]
    truth = [{"accession": "P1", "in_targets": "yes"}, {"accession": "P2", "in_targets": "yes"}]

    def h(method, t, e):
        return {
            "method": method,
            "query_id": q,
            "target": t,
            "evalue": e,
            "sp_yes": "1",
            "gpi_yes": "0",
            "disorder_longest_run": "0",
        }

    hits = [
        h("foldseek", "C1", "1e-5"),
        h("foldseek", "C2", "5e-4"),
        h("mmseqs", "P1", "1e-9"),
        h("foldseek", "B1", "1e-6"),
    ]
    per_query = [{"method": "foldseek", "query_id": q, "best_decoy_evalue": "1e-4"}]
    af = {
        (q, "C1"): 1e-6,
        (q, "C2"): 1e-4,
        (q, "C3"): 2e-6,
        (q, "P1"): 1e-8,
        (q, "A1"): 1e-5,
        (q, "B1"): 1e-7,
        (q, "B2"): 1e-3,
    }
    rsc = {(q, "C1"): 1e-4, (q, "P1"): 1e-6, (q, "A1"): 1e-2, (q, "B1"): 1e-3}
    return domains, models, hits, per_query, af, rsc, truth


def test_two_of_three_rule_and_half_b():
    pairs, targets, nulls, bg_hits, summary = rs.summarize(*_summary_fixture())
    p = {x["target"]: x for x in pairs}
    assert (
        p["C1"]["n_support"] == 3
    )  # prostt5 1e-5 < decoy 1e-4; af 1e-6 < A null; reseek 1e-4 < A null 1e-2
    assert p["C2"]["n_support"] == 0  # prostt5 not below decoy null; af 1e-4 not below 1e-5
    assert p["P1"]["n_support"] == 2 and p["P1"]["mmseqs_e"] == "1e-09"
    # C3 is in half A; its own score (2e-6) must not be its null (null = next best A score, A1 1e-5)
    assert p["C3"]["support_af_fs"] == "yes"
    assert "A1" not in p and "B1" not in p  # background never scored as a hit
    t = {x["target"]: x for x in targets}
    assert (
        t["C1"]["confirmed"] == "yes"
        and t["C2"]["confirmed"] == "no"
        and t["P1"]["confirmed"] == "yes"
    )
    s = {x["method"]: x for x in summary}
    assert s["2of3"]["pf21671_supported_S1"] == "1/2" and s["2of3"]["pf21671_no_model_S1"] == 1
    # half B (null from half A = min(A1 1e-5, C3 2e-6)): B1 passes prostt5 (1e-6 < 1e-4), af_fs (1e-7 < 2e-6) and reseek
    # (1e-3 <= cut, < A null 1e-2); B2 af 1e-3 is not below 2e-6
    assert s["af_fs"]["half_B_passing"] == "1/2" and s["2of3"]["half_B_passing"] == "1/2"
    n = {(x["query_id"], x["method"]): x for x in nulls}
    assert (
        n[("Q1_PF_1_9", "af_fs")]["null_n"] == "2/2"
        and n[("Q1_PF_1_9", "af_fs")]["null_best"] == "2e-06"
    )
    assert {(b["method"], b["target"]) for b in bg_hits} >= {
        ("af_fs", "A1"),
        ("af_fs", "B1"),
        ("reseek", "B1"),
    }


def test_missing_null_counts_and_unknown_target():
    d, m, h, pq, af, rsc, tr = _summary_fixture()
    rsc = {k: v for k, v in rsc.items() if k[1] not in ("A1", "C3")}  # no reseek score in half A
    _, _, nulls, _, _ = rs.summarize(d, m, h, pq, af, rsc, tr)
    n = {(x["query_id"], x["method"]): x for x in nulls}
    assert (
        n[("Q1_PF_1_9", "reseek")]["null_n"] == "0/2"
        and n[("Q1_PF_1_9", "reseek")]["null_best"] == "NA"
    )
    af[("Q1_PF_1_9", "XX")] = 1.0
    with pytest.raises(rs.SummaryError):
        rs.summarize(d, m, h, pq, af, rsc, tr)


# gene_model_check ------------------------------------------------------------------------------------------------------


def _gm_fetch(entry, members, uniref_fail_first=False):
    calls = []

    def fetch(url):
        calls.append(url)
        if "/uniprotkb/M7WU28.json" in url:
            return json.dumps(entry)
        if "/uniref/search" in url:
            if uniref_fail_first and "uniprot_id" in url:
                raise gm.CheckError("HTTP 400")
            return json.dumps({"results": [{"id": "UniRef50_A0A111"}]})
        return members

    return fetch


def test_gene_model_check():
    entry = {
        "sequence": {"value": "MKLA" + "S" * 50},
        "organism": {"scientificName": "Rhodotorula toruloides NP11"},
        "genes": [{"orfNames": [{"value": "RHTO_01234"}]}],
        "features": [
            {"type": "Signal", "location": {"start": {"value": 1}, "end": {"value": None}}}
        ],
    }
    members = (
        "Entry\tOrganism\tLength\tSignal peptide\tSequence\n"
        "M7WU28\tR. toruloides NP11\t54\t\tMKLA\n"
        "A0A111\tR. toruloides IFO0880\t80\tSIGNAL 1..20; /evidence=x\tMFFFF\n"
    )
    rows = gm.check("M7WU28", _gm_fetch(entry, members))
    assert (
        rows[0]["signal"] == "1-?"
        and rows[0]["orf"] == "RHTO_01234"
        and rows[0]["uniref50"] == "UniRef50_A0A111"
    )
    assert (
        len(rows) == 2
        and rows[1]["member"] == "A0A111"
        and rows[1]["signal"] == "1-20"
        and rows[1]["length"] == 80
    )
    assert (
        gm.check("M7WU28", _gm_fetch(entry, members, uniref_fail_first=True))[0]["uniref50"]
        == "UniRef50_A0A111"
    )
    with pytest.raises(gm.CheckError):
        gm.entry_summary({"primaryAccession": "X"})
