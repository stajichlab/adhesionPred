import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.calibration.cli import main
from cellsurface_sorting_hat.calibration.panel import panel_check
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import load_status_source

NODE_ROWS = [
    (1, 1, "no rank"),
    (5052, 1, "genus"),
    (4932, 1, "species"),
    (5476, 1, "species"),
    (746128, 5052, "species"),
    (162425, 5052, "species"),
    (5207, 1, "species"),
    (5270, 1, "species"),
    (246410, 1, "species"),
    (9999, 4932, "strain"),
    (7, 4932, "no rank"),
    (8, 1, "no rank"),
    (9, 1, "family"),
]


def write_nodes(tmp_path):
    path = tmp_path / "nodes.dmp"
    path.write_text("".join(f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\n" for t, p, r in NODE_ROWS))
    return path


OPTIONS = "# Option settings:     hmmsearch --cut_ga --cpu 2 --noali fam.hmm in.fasta\n"


def write_domtbl(path, body):
    path.write_text(OPTIONS + body + "# [ok]\n")


def write_run_json(folder):
    """``run.json`` with the identity of every module record in ``folder/modules``, as a run writes it."""
    identities = []
    for p in sorted((folder / "modules").glob("*.json")):
        rec = json.loads(p.read_text())
        identities.append(
            {k: rec[k] for k in ("module", "version", "params_hash", "artefact_hash")}
        )
    for i in identities:
        i["name"] = i.pop("module")
    (folder / "run.json").write_text(json.dumps({"module_identities": identities}))


def write_proteins(path, ids, taxon):
    """``proteins.tsv.gz`` and ``run.json`` next to the calls file, as a run writes them."""
    write_run_json(path.parent)
    with gzip.open(path.with_name("proteins.tsv.gz"), "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["id", "sha256", "taxon", "state", "note", "trailing_stop", "ambiguous_fraction"]
        )
        for pid in sorted(set(ids)):
            w.writerow([pid, "x", taxon, "ok", "", 0, "0.0000"])


def write_calls(path, rows, taxon=4932):
    write_proteins(path, [r[0] for r in rows], taxon)
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein", "call", "variant", "value", "status", "status_basis", "other_basis"])
        for p, c, v, val in rows:
            w.writerow([p, c, v, val, "unvalidated", "", ""])


def write_truth(path, rows):
    path.write_text("id\tlabel\tcluster\n" + "".join(f"{i}\t{y}\t{c}\n" for i, y, c in rows))


def test_truth_command_writes_an_entry_with_sensitivity_and_specificity(tmp_path, capsys):
    write_module(tmp_path, ModuleSpec("allergen_homology", "1"), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{i}" for i in range(40)]
    calls = [
        (i, "iuis_allergen_similarity", "", "called" if k < 12 or 20 <= k < 24 else "not_called")
        for k, i in enumerate(ids)
    ]
    calls.append(("U1", "iuis_allergen_similarity", "", "not_assessable"))
    write_calls(tmp_path / "c.tsv.gz", calls, taxon=746128)
    write_truth(
        tmp_path / "t.tsv",
        [(i, 1 if k < 20 else 0, f"c{k}") for k, i in enumerate(ids)]
        + [("U1", 1, "cu"), ("GONE", 0, "cg")],
    )
    code = main(
        [
            "truth",
            "--workdir",
            str(tmp_path),
            "--module",
            "allergen_homology",
            "--calls-long",
            str(tmp_path / "c.tsv.gz"),
            "--call",
            "iuis_allergen_similarity",
            "--truth",
            str(tmp_path / "t.tsv"),
            "--calibration-set",
            "toy",
            "--leakage",
            "none",
            "--nodes-dmp",
            str(write_nodes(tmp_path)),
            "--taxa",
            "746128",
            "--n-boot",
            "200",
        ]
    )
    assert code == 0
    entry = load_status_source(tmp_path / "status" / "allergen_homology.json").entries[0]
    m = entry.measure
    assert (m["n_pos"], m["n_neg"]) == (20, 20)
    assert m["sensitivity"]["value"] == pytest.approx(12 / 20) and m["specificity"][
        "value"
    ] == pytest.approx(16 / 20)
    assert "truth rows without a call: 1" in m["notes"] and "not assessable: 1" in m["notes"]
    assert "sensitivity if not assessable positives count as missed: 0.571" in m["notes"]
    assert entry.status == "smoke"  # 20 positives but a wide interval


def test_a_second_set_is_added_and_the_same_set_is_replaced(tmp_path):
    write_module(tmp_path, ModuleSpec("pfam_adhesion", "1"), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{i}" for i in range(10)]
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            (i, "wall_family_domain", "", "called" if k < 5 else "not_called")
            for k, i in enumerate(ids)
        ],
    )
    write_truth(tmp_path / "t.tsv", [(i, 1 if k < 5 else 0, f"c{k}") for k, i in enumerate(ids)])

    def run(name, taxon):
        write_proteins(tmp_path / "c.tsv.gz", ids, taxon)
        return main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "pfam_adhesion",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "wall_family_domain",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                name,
                "--leakage",
                "none",
                "--nodes-dmp",
                str(write_nodes(tmp_path)),
                "--taxa",
                str(taxon),
                "--n-boot",
                "50",
            ]
        )

    assert run("setA", 4932) == 0 and run("setB", 5476) == 0 and run("setA", 4932) == 0
    sets = [
        e.measure["calibration_set"]
        for e in load_status_source(tmp_path / "status" / "pfam_adhesion.json").entries
    ]
    assert sorted(sets) == ["setA", "setB"]


def test_truth_with_no_matching_protein_is_an_error(tmp_path, capsys):
    write_module(tmp_path, ModuleSpec("pfam_adhesion", "1"), [], [{"id": "A", "state": "ok"}])
    write_calls(tmp_path / "c.tsv.gz", [("X", "wall_family_domain", "", "called")])
    write_truth(tmp_path / "t.tsv", [("Y", 1, "c")])
    assert (
        main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "pfam_adhesion",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "wall_family_domain",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                "s",
                "--leakage",
                "none",
                "--nodes-dmp",
                str(write_nodes(tmp_path)),
                "--taxa",
                "4932",  # a valid species in the fixture nodes.dmp; the refusal comes from the empty match
            ]
        )
        == 2
    )
    assert "no truth protein has a call" in capsys.readouterr().err


def test_phasec_command_resolves_species_names_to_taxa(tmp_path):
    names = tmp_path / "names.dmp"
    names.write_text(
        "".join(
            f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n"
            for t, n in [
                (4932, "Saccharomyces cerevisiae"),
                (5476, "Candida albicans"),
                (746128, "Aspergillus fumigatus"),
                (162425, "Aspergillus nidulans"),
                (5207, "Cryptococcus neoformans"),
                (5270, "Ustilago maydis"),
            ]
        )
    )
    sp = tmp_path / "sets.tsv"
    sp.write_text(
        "set_key\tscientific_name\nS1:Scer_SGD\tSaccharomyces cerevisiae\nS1:Calb_CGD\tCandida albicans\n"
        "S3-Eurotiomycetes:Afum_ASPFU\tAspergillus fumigatus\nS3-Eurotiomycetes:Anid_EMENI\tAspergillus nidulans\n"
        "S3-Basidiomycota:Cneo_H99_GOA\tCryptococcus neoformans\nS3-Basidiomycota:Umay_MYCMD\tUstilago maydis\n"
    )

    def cell(r, f):
        return {
            "recall": {"value": r, "lo": r - 0.05, "hi": r + 0.05},
            "fpr": {"value": f, "lo": f / 2, "hi": f * 2},
        }

    def ts(label, pos, c):
        return {
            "label": label,
            "n_direct_positives": pos,
            "truth": {
                "direct": {
                    "n": {"all": {"pos": pos, "neg": 100}},
                    "metrics": {"all": {"V-go": {"R0": c}}},
                }
            },
        }

    metrics = tmp_path / "metrics.json"
    metrics.write_text(
        json.dumps(
            {
                "test_sets": {
                    "S1:Scer_SGD": ts("estimate", 232, cell(0.6, 0.04)),
                    "S1:Calb_CGD": ts("smoke test", 153, cell(0.5, 0.04)),
                    "S3-Eurotiomycetes:Afum_ASPFU": ts("smoke test", 19, cell(0.9, 0.01)),
                    "S3-Eurotiomycetes:Anid_EMENI": ts("estimate", 109, cell(0.7, 0.01)),
                    "S3-Basidiomycota:Cneo_H99_GOA": ts("smoke test", 7, cell(0.9, 0.08)),
                    "S3-Basidiomycota:Umay_MYCMD": ts("smoke test", 9, cell(0.9, 0.08)),
                }
            }
        )
    )
    write_r0(tmp_path / "wd")
    clusters, eval_table = write_cluster_files(
        tmp_path,
        {
            "Scer_SGD": (232, 100),
            "Calb_CGD": (153, 100),
            "Afum_ASPFU": (19, 100),
            "Anid_EMENI": (109, 100),
            "Cneo_H99_GOA": (7, 100),
            "Umay_MYCMD": (9, 100),
        },
    )
    assert (
        main(
            [
                "phasec",
                "--workdir",
                str(tmp_path / "wd"),
                "--metrics",
                str(metrics),
                "--set-species",
                str(sp),
                "--nodes-dmp",
                str(write_nodes(tmp_path)),
                "--names-dmp",
                str(names),
                *signalp_args(),
                "--clusters",
                str(clusters),
                "--eval-table",
                str(eval_table),
            ]
        )
        == 0
    )
    entries = load_status_source(tmp_path / "wd" / "status" / "step1_rule@R0.json").entries
    assert {e.taxa for e in entries} == {(4932,), (5476,), (746128,), (162425,), (5207,), (5270,)}
    assert {e.measure["calibration_set"]: e.status for e in entries}[
        "S3-Basidiomycota:Umay_MYCMD"
    ] == "smoke"


def test_a_changed_module_identity_drops_the_old_status_entries(tmp_path, capsys):
    write_module(
        tmp_path, ModuleSpec("pfam_adhesion", "1", {"a": 1}), [], [{"id": "A", "state": "ok"}]
    )
    ids = [f"P{i}" for i in range(10)]
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            (i, "wall_family_domain", "", "called" if k < 5 else "not_called")
            for k, i in enumerate(ids)
        ],
    )
    write_truth(tmp_path / "t.tsv", [(i, 1 if k < 5 else 0, f"c{k}") for k, i in enumerate(ids)])

    def run(name, taxon):
        write_proteins(tmp_path / "c.tsv.gz", ids, taxon)
        return main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "pfam_adhesion",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "wall_family_domain",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                name,
                "--nodes-dmp",
                str(write_nodes(tmp_path)),
                "--taxa",
                str(taxon),
                "--leakage",
                "none",
                "--n-boot",
                "50",
            ]
        )

    assert run("setA", 4932) == 0
    write_module(
        tmp_path, ModuleSpec("pfam_adhesion", "2", {"a": 2}), [], [{"id": "A", "state": "ok"}]
    )  # new version
    assert run("setB", 5476) == 0
    assert "module identity changed" in capsys.readouterr().err
    entries = load_status_source(tmp_path / "status" / "pfam_adhesion.json").entries
    assert [e.measure["calibration_set"] for e in entries] == [
        "setB"
    ]  # setA was measured on version 1


def test_panel_check_counts_agreement_and_leaves_known_misses_out(tmp_path):
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            ("A", "wall_family_domain", "", "called"),
            ("B", "wall_family_domain", "", "not_called"),
            ("C", "wall_family_domain", "", "called"),
        ],
    )
    panel = tmp_path / "p.tsv"
    panel.write_text(
        "protein\tcall\tvariant\texpected\tsource\nA\twall_family_domain\t\tcalled\tx\nB\twall_family_domain\t\tcalled\tx\n"
        "C\twall_family_domain\t\tknown_miss\tx\nZ\twall_family_domain\t\tcalled\tx\n"
    )
    rows, summary = panel_check(tmp_path / "c.tsv.gz", panel)
    assert summary == {"agree": 1, "disagree": 1, "known_miss": 1, "not_in_run": 1}
    assert [r["observed"] for r in rows] == ["called", "not_called", "called", "missing"]


def test_allergen_lso_command_prints_recall_per_rule(tmp_path, capsys):
    fasta = tmp_path / "a.faa"
    fasta.write_text(">Asp_f_1.0101|1\nM\n>Asp_n_1.0101|2\nM\n>Alt_a_1.0101|3\nM\n")
    blast = tmp_path / "b.tsv"
    blast.write_text(
        "Asp_f_1.0101|1\tAsp_f_1.0101|1\t100\t100\t100\t100\t500\t0\n"
        "Asp_f_1.0101|1\tAsp_n_1.0101|2\t60\t90\t100\t100\t100\t1e-20\n"
        "Asp_n_1.0101|2\tAsp_n_1.0101|2\t100\t100\t100\t100\t500\t0\n"
        "Asp_n_1.0101|2\tAsp_f_1.0101|1\t60\t90\t100\t100\t100\t1e-20\n"
        "Alt_a_1.0101|3\tAlt_a_1.0101|3\t100\t100\t100\t100\t500\t0\n"
    )
    assert main(["allergen-lso", "--blast", str(blast), "--allergen-fasta", str(fasta)]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "sequences\t3\tspecies\t3"
    assert out[1] == "iuis_allergen_similarity\t2/3" and out[2] == "iuis_allergen_homolog\t0/3"


def test_leakage_other_than_none_caps_an_estimate_at_smoke(tmp_path):
    write_module(tmp_path, ModuleSpec("antigen_lookup", "1"), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{i}" for i in range(60)]
    # 30 positives, all called; 30 negatives, none called: a tight interval that would be an estimate
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            (i, "cocci_specificity_rank_top15", "", "called" if k < 30 else "not_called")
            for k, i in enumerate(ids)
        ],
    )
    write_truth(tmp_path / "t.tsv", [(i, 1 if k < 30 else 0, f"c{k}") for k, i in enumerate(ids)])

    def run(leakage, name):
        write_proteins(tmp_path / "c.tsv.gz", ids, 246410 if name == "a" else 5476)
        return main(
            [
                "truth",
                "--workdir",
                str(tmp_path),
                "--module",
                "antigen_lookup",
                "--calls-long",
                str(tmp_path / "c.tsv.gz"),
                "--call",
                "cocci_specificity_rank_top15",
                "--truth",
                str(tmp_path / "t.tsv"),
                "--calibration-set",
                name,
                "--nodes-dmp",
                str(write_nodes(tmp_path)),
                "--taxa",
                "246410" if name == "a" else "5476",
                "--leakage",
                leakage,
                "--n-boot",
                "100",
            ]
        )

    assert run("none", "a") == 0 and run("tuned_on_truth", "b") == 0
    status = {
        e.measure["calibration_set"]: (e.status, e.measure["notes"])
        for e in load_status_source(tmp_path / "status" / "antigen_lookup.json").entries
    }
    assert status["a"][0] == "estimated" and status["b"][0] == "smoke"
    assert "leakage: tuned_on_truth" in status["b"][1]


def test_pfam_specificity_command_lists_non_member_hits_with_their_helices(tmp_path, capsys):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKT\n>B\nMKT\n>C\nMKT\n>D\nMKT\n")
    dom = tmp_path / "d.domtbl"
    row = "{} - 11 CFEM PF05730.17 70 1e-20 60 8 1 1 1e-21 2e-20 59 8 1 70 2 9 2 9 0.9 -\n"
    write_domtbl(dom, row.format("A") + row.format("B") + row.format("D"))
    members = tmp_path / "m.tsv"
    members.write_text("pfam_acc\tprotein_id\nPF05730\tA\nPF05730\tC\nPF01185\tD\n")
    tm = tmp_path / "t.tsv"
    tm.write_text(
        "protein_id\tlen\texp_aa\tfirst60\tpred_hel\ttopology\nB\t11\t150\t20\t7\to5-27i\n"
    )
    assert (
        main(
            [
                "pfam-specificity",
                "--family",
                "PF05730",
                "--domtbl",
                str(dom),
                "--members",
                str(members),
                "--universe-fasta",
                str(fasta),
                "--tm-table",
                str(tm),
            ]
        )
        == 0
    )
    out = capsys.readouterr().out.splitlines()
    assert out[0] == "family\tPF05730\thits\t3\tmembers\t2"
    assert "tp\t1" in out and "fp\t2" in out and "fn\t1" in out
    assert (
        "nonmember_hit\tB\tn_tm=7" in out
        and "nonmember_hit\tD\tn_tm=NA" in out
        and "missed_member\tC" in out
    )


def test_pfam_specificity_refuses_members_outside_the_proteome(tmp_path, capsys):
    fasta = tmp_path / "p.faa"
    fasta.write_text(">A\nMKT\n")
    dom = tmp_path / "d.domtbl"
    write_domtbl(dom, "")
    members = tmp_path / "m.tsv"
    members.write_text("pfam_acc\tprotein_id\nPF05730\tZ\n")
    assert (
        main(
            [
                "pfam-specificity",
                "--family",
                "PF05730",
                "--domtbl",
                str(dom),
                "--members",
                str(members),
                "--universe-fasta",
                str(fasta),
            ]
        )
        == 2
    )
    assert "not in the proteome" in capsys.readouterr().err


def test_panel_marks_proteins_seen_in_the_reference_as_excluded(tmp_path):
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            ("A", "iuis_allergen_homolog", "", "called"),
            ("B", "iuis_allergen_homolog", "", "called"),
        ],
    )
    panel = tmp_path / "p.tsv"
    panel.write_text(
        "protein\tcall\tvariant\texpected\tsource\ttuning\nA\tiuis_allergen_homolog\t\tcalled\tx\tin_reference\n"
        "B\tiuis_allergen_homolog\t\tcalled\tx\t\n"
    )
    _, summary = panel_check(tmp_path / "c.tsv.gz", panel)
    assert summary == {"excluded_leakage": 1, "agree": 1}


# ---- Plan 2 Global Constraints: refusals and caps added to the brief's tests ----


def truth_args(
    tmp_path,
    *extra,
    leakage="none",
    taxa="4932",
    name="s",
    module="pfam_adhesion",
    call="wall_family_domain",
):
    args = [
        "truth",
        "--workdir",
        str(tmp_path),
        "--module",
        module,
        "--calls-long",
        str(tmp_path / "c.tsv.gz"),
        "--call",
        call,
        "--truth",
        str(tmp_path / "t.tsv"),
        "--calibration-set",
        name,
        "--nodes-dmp",
        str(write_nodes(tmp_path)),
        "--taxa",
        *taxa.split(),
        "--n-boot",
        "50",
    ]
    if leakage is not None:
        args += ["--leakage", leakage]
    return args + list(extra)


def setup_truth(tmp_path, n_pos=5, n_neg=5, pos_clusters=None, record=True, taxon=4932):
    if record:
        write_module(tmp_path, ModuleSpec("pfam_adhesion", "1"), [], [{"id": "A", "state": "ok"}])
    rows = [(f"P{k}", 1, f"p{k % (pos_clusters or n_pos)}") for k in range(n_pos)]
    rows += [(f"N{k}", 0, f"n{k}") for k in range(n_neg)]
    write_calls(
        tmp_path / "c.tsv.gz",
        [(f"P{k}", "wall_family_domain", "", "called") for k in range(n_pos)]
        + [(f"N{k}", "wall_family_domain", "", "not_called") for k in range(n_neg)],
        taxon=taxon,
    )
    write_truth(tmp_path / "t.tsv", rows)


@pytest.mark.parametrize(
    "taxon", ["0", "1", "-5", "4932 1", "5052", "9", "8", "12345", "4932 5476"]
)
def test_truth_refuses_the_root_or_an_invalid_taxon_and_writes_nothing(tmp_path, capsys, taxon):
    setup_truth(tmp_path)
    assert main(truth_args(tmp_path, taxa=taxon)) == 2
    assert "taxon" in capsys.readouterr().err
    assert not (tmp_path / "status").exists()


def test_truth_requires_a_leakage_value(tmp_path):
    setup_truth(tmp_path)
    with pytest.raises(SystemExit) as err:
        main(truth_args(tmp_path, leakage=None))
    assert err.value.code == 2
    with pytest.raises(SystemExit) as err:
        main(truth_args(tmp_path, leakage="maybe"))
    assert err.value.code == 2
    assert not (tmp_path / "status").exists()


@pytest.mark.parametrize("leakage", ["partial", "tuned_on_truth", "in_reference", "unknown"])
def test_every_leakage_value_except_none_caps_at_smoke(tmp_path, leakage):
    setup_truth(tmp_path, n_pos=30, n_neg=30)
    assert main(truth_args(tmp_path, leakage=leakage)) == 0
    entry = load_status_source(tmp_path / "status" / "pfam_adhesion.json").entries[0]
    assert entry.status == "smoke" and f"leakage: {leakage}" in entry.measure["notes"]


def test_truth_with_no_negatives_is_at_most_smoke(tmp_path):
    setup_truth(tmp_path, n_pos=30, n_neg=0)
    write_truth(tmp_path / "t.tsv", [(f"P{k}", 1, f"p{k}") for k in range(30)])
    assert main(truth_args(tmp_path)) == 0
    entry = load_status_source(tmp_path / "status" / "pfam_adhesion.json").entries[0]
    assert "specificity" not in entry.measure and entry.status == "smoke"


def test_truth_with_few_clusters_is_not_an_estimate(tmp_path):
    setup_truth(tmp_path, n_pos=30, n_neg=30, pos_clusters=5)
    assert main(truth_args(tmp_path)) == 0
    entry = load_status_source(tmp_path / "status" / "pfam_adhesion.json").entries[0]
    assert entry.measure["n_clusters_pos"] == 5 and entry.status == "smoke"


def test_truth_with_too_few_negatives_is_not_an_estimate(tmp_path):
    setup_truth(tmp_path, n_pos=30, n_neg=19)
    assert main(truth_args(tmp_path)) == 0
    assert (
        load_status_source(tmp_path / "status" / "pfam_adhesion.json").entries[0].status == "smoke"
    )


def test_a_status_write_without_a_module_run_record_is_refused(tmp_path, capsys):
    setup_truth(tmp_path, record=False)
    assert main(truth_args(tmp_path)) == 2
    assert "pfam_adhesion.json" in capsys.readouterr().err
    assert not (tmp_path / "status").exists()


@pytest.mark.parametrize(
    "body, where",
    [
        ("id\tlabel\tcluster\nP0\t2\tc\n", "t.tsv:2"),
        ("id\tlabel\tcluster\nP0\tyes\tc\n", "t.tsv:2"),
        ("id\tlabel\tcluster\nP0\t1\t\n", "t.tsv:2"),
        ("id\tlabel\tcluster\nP0\t1\tc\nP0\t0\td\n", "t.tsv:3"),
        ("id\tlabel\nP0\t1\n", "t.tsv"),
    ],
)
def test_truth_parse_errors_name_the_path_and_line(tmp_path, capsys, body, where):
    setup_truth(tmp_path)
    (tmp_path / "t.tsv").write_text(body)
    assert main(truth_args(tmp_path)) == 2
    assert where in capsys.readouterr().err
    assert not (tmp_path / "status").exists()


def test_truth_refuses_an_unknown_call_value_and_a_bad_n_boot(tmp_path, capsys):
    setup_truth(tmp_path)
    write_calls(tmp_path / "c.tsv.gz", [("P0", "wall_family_domain", "", "maybe")])
    assert main(truth_args(tmp_path)) == 2
    assert "maybe" in capsys.readouterr().err
    setup_truth(tmp_path)
    assert main(truth_args(tmp_path, "--n-boot", "0")) == 2
    assert not (tmp_path / "status").exists()


def test_truth_on_calls_with_a_conflicting_duplicate_is_refused(tmp_path, capsys):
    setup_truth(tmp_path)
    write_calls(
        tmp_path / "c.tsv.gz",
        [
            ("P0", "wall_family_domain", "", "called"),
            ("P0", "wall_family_domain", "", "not_called"),
        ],
    )
    assert main(truth_args(tmp_path)) == 2
    assert "P0" in capsys.readouterr().err


def test_a_taxon_in_two_sets_is_refused_and_leaves_the_first_file_unchanged(tmp_path, capsys):
    setup_truth(tmp_path)
    assert main(truth_args(tmp_path, name="a")) == 0
    path = tmp_path / "status" / "pfam_adhesion.json"
    before = path.read_bytes()
    assert main(truth_args(tmp_path, name="b")) == 2  # same taxon 4932 under another set
    assert "4932" in capsys.readouterr().err and path.read_bytes() == before


def write_r0(workdir, version="1", signalp="6.0h-gpu", mode="fast"):
    """The R0 module record that ``cellsurface_sorting_hat_modules signalp`` writes."""
    spec = ModuleSpec(
        "step1_rule@R0",
        version,
        {"rule": "R0", "mode": mode, "organism": "eukarya"},
        (),
        {"signalp": signalp},
        artefact_digest="x" + version,
    )
    write_module(workdir, spec, [], [{"id": "A", "state": "ok"}])


def signalp_args(module="signalp/6-gpu", mode="fast"):
    return ["--phasec-signalp-module", module, "--phasec-signalp-mode", mode]


def write_cluster_files(tmp_path, sizes, n_clusters=None):
    """``clusters.tsv.gz`` and ``eval_table.tsv.gz`` with ``(n_pos, n_neg)`` proteins per source.

    Each protein of a class gets its own cluster unless ``n_clusters`` (class -> count) is given."""
    import gzip

    cl, ev = ["seq_sha256\tcluster_id"], ["seq_sha256\torigin\tclass\thomology_only\tsource_ids"]
    for source, (n_pos, n_neg) in sizes.items():
        for cls, n in (("pos", n_pos), ("neg", n_neg)):
            for k in range(n):
                h = f"{source}-{cls}-{k}"
                cluster = f"{source}-{cls}-{k % n_clusters[cls]}" if n_clusters else h
                cl.append(f"{h}\t{cluster}")
                ev.append(f"{h}\tgo\t{cls}\tno\t{source}")
    paths = []
    for name, lines in (("clusters.tsv.gz", cl), ("eval_table.tsv.gz", ev)):
        paths.append(tmp_path / name)
        with gzip.open(paths[-1], "wt") as fh:
            fh.write("\n".join(lines) + "\n")
    return paths


def phasec_inputs(
    tmp_path, names_rows=None, sets_rows=None, label_umay="smoke test", n_clusters=None
):
    names = tmp_path / "names.dmp"
    names_rows = names_rows or [
        (4932, "Saccharomyces cerevisiae"),
        (5476, "Candida albicans"),
    ]
    names.write_text("".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, n in names_rows))
    sets_rows = sets_rows or [
        ("S1:Scer_SGD", "Saccharomyces cerevisiae"),
        ("S1:Calb_CGD", "Candida albicans"),
    ]
    sp = tmp_path / "sets.tsv"
    sp.write_text("set_key\tscientific_name\n" + "".join(f"{k}\t{n}\n" for k, n in sets_rows))

    def cell(r, f):
        return {
            "recall": {"value": r, "lo": r - 0.05, "hi": r + 0.05},
            "fpr": {"value": f, "lo": f / 2, "hi": f * 2},
        }

    def ts(label, pos):
        return {
            "label": label,
            "truth": {
                "direct": {
                    "n": {"all": {"pos": pos, "neg": 100}},
                    "metrics": {"all": {"V-go": {"R0": cell(0.6, 0.04)}}},
                }
            },
        }

    metrics = tmp_path / "metrics.json"
    metrics.write_text(
        json.dumps(
            {"test_sets": {"S1:Scer_SGD": ts("estimate", 232), "S1:Calb_CGD": ts(label_umay, 153)}}
        )
    )
    clusters, eval_table = write_cluster_files(
        tmp_path, {"Scer_SGD": (232, 100), "Calb_CGD": (153, 100)}, n_clusters=n_clusters
    )
    return [
        "--metrics",
        str(metrics),
        "--set-species",
        str(sp),
        "--nodes-dmp",
        str(write_nodes(tmp_path)),
        "--names-dmp",
        str(names),
        "--clusters",
        str(clusters),
        "--eval-table",
        str(eval_table),
        *signalp_args(),
    ]


def run_phasec(tmp_path, extra):
    return main(["phasec", "--workdir", str(tmp_path), *extra])


def test_phasec_status_is_the_weaker_of_the_label_and_the_rule(tmp_path):
    write_r0(tmp_path)
    assert run_phasec(tmp_path, phasec_inputs(tmp_path)) == 0
    got = {
        e.measure["calibration_set"]: e.status
        for e in load_status_source(tmp_path / "status" / "step1_rule@R0.json").entries
    }
    # same numbers; only the Phase C label differs
    assert got == {"S1:Scer_SGD": "estimated", "S1:Calb_CGD": "smoke"}


@pytest.mark.parametrize(
    "names_rows",
    [
        [(4932, "Saccharomyces cerevisiae"), (5476, "Candida albicans"), (9, "Candida albicans")],
        [(4932, "Saccharomyces cerevisiae")],  # Candida albicans missing
        [(1, "Saccharomyces cerevisiae"), (5476, "Candida albicans")],
        [(0, "Saccharomyces cerevisiae"), (5476, "Candida albicans")],
    ],
)
def test_phasec_refuses_a_name_with_zero_or_several_ids_or_the_root(tmp_path, capsys, names_rows):
    write_r0(tmp_path)
    assert run_phasec(tmp_path, phasec_inputs(tmp_path, names_rows)) == 2
    assert capsys.readouterr().err
    assert not (tmp_path / "status").exists()


@pytest.mark.parametrize(
    "sets_rows",
    [
        [("S1:Scer_SGD", "Saccharomyces cerevisiae"), ("S1:Scer_SGD", "Candida albicans")],
        [("S9:Nothing", "Candida albicans")],
        [],
    ],
)
def test_phasec_refuses_a_set_with_several_species_an_unknown_set_or_no_set(
    tmp_path, capsys, sets_rows
):
    write_r0(tmp_path)
    args = phasec_inputs(tmp_path, sets_rows=sets_rows or [("x", "y")])
    if not sets_rows:
        (tmp_path / "sets.tsv").write_text("set_key\tscientific_name\n")
    assert run_phasec(tmp_path, args) == 2
    assert capsys.readouterr().err
    assert not (tmp_path / "status").exists()


def test_phasec_refuses_bad_set_species_rows_with_path_and_line(tmp_path, capsys):
    write_r0(tmp_path)
    args = phasec_inputs(tmp_path)
    (tmp_path / "sets.tsv").write_text(
        "set_key\tscientific_name\nS1:Scer_SGD\tSaccharomyces cerevisiae\nS1:Calb_CGD\t\n"
    )
    assert run_phasec(tmp_path, args) == 2
    assert "sets.tsv:3" in capsys.readouterr().err


def test_phasec_without_a_module_run_record_is_refused(tmp_path, capsys):
    assert run_phasec(tmp_path, phasec_inputs(tmp_path)) == 2
    assert "step1_rule@R0.json" in capsys.readouterr().err
    assert not (tmp_path / "status").exists()


def test_phasec_with_a_changed_module_identity_drops_old_entries(tmp_path, capsys):
    write_r0(tmp_path)
    assert run_phasec(tmp_path, phasec_inputs(tmp_path)) == 0
    write_r0(tmp_path, "2")
    both = phasec_inputs(tmp_path, sets_rows=[("S1:Scer_SGD", "Saccharomyces cerevisiae")])
    assert run_phasec(tmp_path, both) == 0
    assert "module identity changed" in capsys.readouterr().err
    record = load_status_source(tmp_path / "status" / "step1_rule@R0.json")
    assert [e.measure["calibration_set"] for e in record.entries] == ["S1:Scer_SGD"]
    assert record.identity.version == "2"


def test_phasec_has_no_leakage_argument(tmp_path):
    write_r0(tmp_path)
    with pytest.raises(SystemExit):
        run_phasec(tmp_path, [*phasec_inputs(tmp_path), "--leakage", "none"])


def test_panel_command_prints_verdicts_and_writes_no_status(tmp_path, capsys):
    write_calls(tmp_path / "c.tsv.gz", [("A", "wall_family_domain", "", "not_called")])
    panel = tmp_path / "p.tsv"
    panel.write_text(
        "protein\tcall\tvariant\texpected\tsource\nA\twall_family_domain\t\tcalled\tx\n"
    )
    before = sorted(p.name for p in tmp_path.iterdir())
    assert main(["panel", "--calls-long", str(tmp_path / "c.tsv.gz"), "--panel", str(panel)]) == 0
    cap = capsys.readouterr()
    assert cap.out.strip().endswith("disagree") and '"disagree": 1' in cap.err
    assert sorted(p.name for p in tmp_path.iterdir()) == before


@pytest.mark.parametrize(
    "body, where",
    [
        ("protein\tcall\tvariant\texpected\tsource\nA\tx\t\tmaybe\ts\n", "p.tsv:2"),
        ("protein\tcall\tvariant\texpected\tsource\tTuning\nA\tx\t\tcalled\ts\tsorta\n", "p.tsv:2"),
        ("protein\tcall\tvariant\tsource\nA\tx\t\ts\n", "p.tsv"),
    ],
)
def test_panel_refuses_bad_rows_with_path_and_line(tmp_path, capsys, body, where):
    write_calls(tmp_path / "c.tsv.gz", [("A", "x", "", "called")])
    panel = tmp_path / "p.tsv"
    panel.write_text(body.replace("Tuning", "tuning"))
    assert main(["panel", "--calls-long", str(tmp_path / "c.tsv.gz"), "--panel", str(panel)]) == 2
    assert where in capsys.readouterr().err


def test_panel_never_calls_a_missing_protein_an_agreement(tmp_path):
    write_calls(tmp_path / "c.tsv.gz", [("A", "x", "", "called")])
    panel = tmp_path / "p.tsv"
    panel.write_text("protein\tcall\tvariant\texpected\tsource\nZ\tx\t\tnot_called\ts\n")
    _, summary = panel_check(tmp_path / "c.tsv.gz", panel)
    assert summary == {"not_in_run": 1}


@pytest.mark.parametrize("taxon", ["4932", "9999", "7"])  # species, strain, no rank under a species
def test_truth_accepts_a_species_or_a_taxon_below_one(tmp_path, taxon):
    setup_truth(tmp_path, taxon=int(taxon))
    assert main(truth_args(tmp_path, taxa=taxon)) == 0
    assert load_status_source(tmp_path / "status" / "pfam_adhesion.json").entries[0].taxa == (
        int(taxon),
    )


def test_two_taxa_say_to_run_one_call_per_species(tmp_path, capsys):
    setup_truth(tmp_path)
    assert main(truth_args(tmp_path, taxa="4932 5476")) == 2
    assert "one call per species" in capsys.readouterr().err


def test_a_genus_is_refused_with_its_rank(tmp_path, capsys):
    setup_truth(tmp_path)
    assert main(truth_args(tmp_path, taxa="5052")) == 2
    err = capsys.readouterr().err
    assert "5052" in err and "genus" in err


def test_phasec_refuses_a_name_that_resolves_above_species(tmp_path, capsys):
    write_r0(tmp_path)
    rows = [(5052, "Saccharomyces cerevisiae"), (5476, "Candida albicans")]
    assert run_phasec(tmp_path, phasec_inputs(tmp_path, rows)) == 2
    assert "genus" in capsys.readouterr().err
    assert not (tmp_path / "status").exists()


def test_the_status_entry_records_call_module_and_variant(tmp_path):
    setup_truth(tmp_path)
    assert main(truth_args(tmp_path)) == 0
    notes = (
        load_status_source(tmp_path / "status" / "pfam_adhesion.json").entries[0].measure["notes"]
    )
    assert "call=wall_family_domain; module=pfam_adhesion; variant=" in notes
    assert "reads_modules" not in notes


def test_a_truth_set_without_negatives_says_specificity_was_not_measured(tmp_path):
    setup_truth(tmp_path, n_pos=5, n_neg=0)
    write_truth(tmp_path / "t.tsv", [(f"P{k}", 1, f"p{k}") for k in range(5)])
    assert main(truth_args(tmp_path)) == 0
    notes = (
        load_status_source(tmp_path / "status" / "pfam_adhesion.json").entries[0].measure["notes"]
    )
    assert "specificity not measured (no negatives)" in notes


def test_a_duplicate_calls_row_is_an_error_with_path_and_line(tmp_path, capsys):
    write_calls(
        tmp_path / "c.tsv.gz",
        [("A", "x", "", "called"), ("A", "x", "", "called")],
    )
    panel = tmp_path / "p.tsv"
    panel.write_text("protein\tcall\tvariant\texpected\tsource\nA\tx\t\tcalled\ts\n")
    assert main(["panel", "--calls-long", str(tmp_path / "c.tsv.gz"), "--panel", str(panel)]) == 2
    assert "c.tsv.gz:3" in capsys.readouterr().err


# ---- final review: module, taxon and reads_modules checks ----


def test_truth_refuses_a_module_that_the_call_does_not_read(tmp_path, capsys):
    setup_truth(tmp_path)
    write_module(tmp_path, ModuleSpec("antigen_lookup", "1"), [], [{"id": "A", "state": "ok"}])
    assert main(truth_args(tmp_path, module="antigen_lookup")) == 2
    err = capsys.readouterr().err
    assert "antigen_lookup" in err and "pfam_adhesion" in err
    assert not (tmp_path / "status").exists()


def test_a_call_that_reads_several_modules_is_refused_and_writes_nothing(tmp_path, capsys):
    setup_truth(tmp_path)
    write_module(tmp_path, ModuleSpec("repeat02", "1"), [], [{"id": "A", "state": "ok"}])
    write_calls(
        tmp_path / "c.tsv.gz",
        [(f"P{k}", "tandem_repeat_protein", "", "called") for k in range(5)]
        + [(f"N{k}", "tandem_repeat_protein", "", "not_called") for k in range(5)],
    )
    args = truth_args(tmp_path, module="repeat02", call="tandem_repeat_protein")
    assert main(args) == 2
    err = capsys.readouterr().err
    assert "tandem_repeat_protein" in err and "repeat02" in err and "repeat14" in err
    assert not (tmp_path / "status").exists()


def test_a_per_variant_call_counts_only_the_step_1_module_of_its_variant(tmp_path):
    write_r0(tmp_path)
    rows = [(f"P{k}", 1, f"p{k}") for k in range(5)] + [(f"N{k}", 0, f"n{k}") for k in range(5)]
    write_calls(
        tmp_path / "c.tsv.gz",
        [(r[0], "signal_peptide_protein", "R0", "called" if r[1] else "not_called") for r in rows],
    )
    write_truth(tmp_path / "t.tsv", rows)
    args = truth_args(
        tmp_path, "--variant", "R0", module="step1_rule@R0", call="signal_peptide_protein"
    )
    assert main(args) == 0
    notes = (
        load_status_source(tmp_path / "status" / "step1_rule@R0.json").entries[0].measure["notes"]
    )
    assert "call=signal_peptide_protein; module=step1_rule@R0; variant=R0" in notes


def test_a_call_that_reads_a_module_through_a_ref_names_it(tmp_path, capsys):
    # cell_wall_adhesion_candidate reads pfam_adhesion through wall_family_domain
    from cellsurface_sorting_hat.engine import load_config, modules_of_call

    reads = modules_of_call(load_config(), "cell_wall_adhesion_candidate")
    assert "pfam_adhesion" in reads and "repeat02" in reads and "step1_rule@R0" in reads


def test_truth_refuses_a_protein_of_another_taxon(tmp_path, capsys):
    setup_truth(tmp_path)  # proteins are of taxon 4932
    assert main(truth_args(tmp_path, taxa="5476")) == 2
    err = capsys.readouterr().err
    assert "10 matched truth protein(s)" in err and "'P0'" in err
    assert not (tmp_path / "status").exists()


def test_truth_refuses_when_proteins_tsv_is_missing(tmp_path, capsys):
    setup_truth(tmp_path)
    (tmp_path / "proteins.tsv.gz").unlink()
    assert main(truth_args(tmp_path)) == 2
    assert "proteins.tsv.gz" in capsys.readouterr().err
    assert not (tmp_path / "status").exists()


def test_truth_refuses_a_truth_protein_missing_from_proteins_tsv(tmp_path, capsys):
    setup_truth(tmp_path)
    write_proteins(tmp_path / "c.tsv.gz", ["P0"], 4932)
    assert main(truth_args(tmp_path)) == 2
    assert "taxon missing" in capsys.readouterr().err


def test_phasec_set_species_comment_lines_are_skipped_and_line_numbers_count_them(tmp_path, capsys):
    write_r0(tmp_path)
    args = phasec_inputs(tmp_path)
    (tmp_path / "sets.tsv").write_text(
        "# a comment\nset_key\tscientific_name\nS1:Scer_SGD\tSaccharomyces cerevisiae\n"
        "# another\nS1:Calb_CGD\t\n"
    )
    assert run_phasec(tmp_path, args) == 2
    assert "sets.tsv:5" in capsys.readouterr().err


# ---- decisions of 2026-10-06: leakage note, cluster counts, SignalP identity, run.json identity ----


def _phasec_entries(tmp_path):
    return {
        e.measure["calibration_set"]: e
        for e in load_status_source(tmp_path / "status" / "step1_rule@R0.json").entries
    }


def test_phasec_records_cluster_counts_and_applies_the_cluster_floor(tmp_path):
    write_r0(tmp_path)
    args = phasec_inputs(tmp_path, n_clusters={"pos": 10, "neg": 90})
    assert run_phasec(tmp_path, args) == 0
    scer = _phasec_entries(tmp_path)["S1:Scer_SGD"]
    assert (scer.measure["n_clusters_pos"], scer.measure["n_clusters_neg"]) == (10, 90)
    assert scer.status == "smoke"  # 10 clusters; with 200 clusters the same file is estimated


def test_phasec_with_enough_clusters_is_estimated_and_the_interval_is_the_widest(tmp_path):
    write_r0(tmp_path)
    assert run_phasec(tmp_path, phasec_inputs(tmp_path)) == 0
    scer = _phasec_entries(tmp_path)["S1:Scer_SGD"]
    assert (scer.measure["n_clusters_pos"], scer.measure["n_clusters_neg"]) == (232, 100)
    assert scer.status == "estimated"
    from cellsurface_sorting_hat.calibration.intervals import wilson_p

    _, lo, hi = wilson_p(0.6, 232)
    sens = scer.measure["sensitivity"]
    assert sens["lo"] == pytest.approx(min(0.55, lo)) and sens["hi"] == pytest.approx(max(0.65, hi))


def test_phasec_notes_carry_call_module_signalp_and_the_leakage_limit(tmp_path):
    write_r0(tmp_path)
    assert run_phasec(tmp_path, phasec_inputs(tmp_path)) == 0
    notes = _phasec_entries(tmp_path)["S1:Scer_SGD"].measure["notes"]
    assert "call=signal_peptide_protein; module=step1_rule@R0" in notes
    assert "signalp_module=signalp/6-gpu; signalp_mode=fast" in notes
    assert (
        "leakage: overlap between the Phase C positives and the SignalP 6 training data "
        "was not measured" in notes
    )
    assert "cluster floor was not checked" not in notes


@pytest.mark.parametrize("missing", ["--clusters", "--eval-table"])
def test_phasec_requires_the_cluster_files(tmp_path, missing):
    write_r0(tmp_path)
    args = phasec_inputs(tmp_path)
    i = args.index(missing)
    with pytest.raises(SystemExit):
        run_phasec(tmp_path, args[:i] + args[i + 2 :])


@pytest.mark.parametrize("missing", ["--phasec-signalp-module", "--phasec-signalp-mode"])
def test_phasec_requires_the_signalp_module_and_mode(tmp_path, missing):
    write_r0(tmp_path)
    args = phasec_inputs(tmp_path)
    i = args.index(missing)
    with pytest.raises(SystemExit):
        run_phasec(tmp_path, args[:i] + args[i + 2 :])


def test_phasec_refuses_protein_counts_that_differ_from_metrics(tmp_path, capsys):
    write_r0(tmp_path)
    args = phasec_inputs(tmp_path)
    write_cluster_files(tmp_path, {"Scer_SGD": (231, 100), "Calb_CGD": (153, 100)})
    assert run_phasec(tmp_path, args) == 2
    err = capsys.readouterr().err
    assert "metrics.json" in err and "eval_table.tsv.gz" in err and "clusters.tsv.gz" in err
    assert "232" in err and "231" in err
    assert not (tmp_path / "status").exists()


@pytest.mark.parametrize(
    "module, mode", [("signalp/6", "fast"), ("signalp/6-gpu", "slow"), ("signalp/5", "fast")]
)
def test_phasec_refuses_a_signalp_module_or_mode_that_the_record_does_not_carry(
    tmp_path, capsys, module, mode
):
    write_r0(tmp_path)
    args = phasec_inputs(tmp_path)
    i = args.index("--phasec-signalp-module")
    args[i : i + 4] = signalp_args(module, mode)
    assert run_phasec(tmp_path, args) == 2
    err = capsys.readouterr().err
    assert module in err or mode in err
    assert not (tmp_path / "status").exists()


def test_phasec_refuses_an_r0_record_run_in_another_mode(tmp_path, capsys):
    write_r0(tmp_path, mode="best")
    assert run_phasec(tmp_path, phasec_inputs(tmp_path)) == 2
    assert "best" in capsys.readouterr().err
    assert not (tmp_path / "status").exists()


def _tamper_run_json(tmp_path, module, field, value):
    path = tmp_path / "run.json"
    data = json.loads(path.read_text())
    for m in data["module_identities"]:
        if m["name"] == module:
            m[field] = value
    path.write_text(json.dumps(data))


@pytest.mark.parametrize("field", ["version", "params_hash", "artefact_hash"])
def test_truth_refuses_when_run_json_identity_differs_from_the_module_record(
    tmp_path, capsys, field
):
    setup_truth(tmp_path)
    _tamper_run_json(tmp_path, "pfam_adhesion", field, "other")
    assert main(truth_args(tmp_path)) == 2
    err = capsys.readouterr().err
    assert "pfam_adhesion" in err and field in err
    assert not (tmp_path / "status").exists()


def test_truth_refuses_when_run_json_is_missing_or_has_no_identities_or_no_module(tmp_path, capsys):
    setup_truth(tmp_path)
    (tmp_path / "run.json").unlink()
    assert main(truth_args(tmp_path)) == 2
    assert "run.json" in capsys.readouterr().err
    (tmp_path / "run.json").write_text(json.dumps({"n_proteins": 3}))
    assert main(truth_args(tmp_path)) == 2
    assert "module_identities" in capsys.readouterr().err
    (tmp_path / "run.json").write_text(json.dumps({"module_identities": []}))
    assert main(truth_args(tmp_path)) == 2
    assert "pfam_adhesion" in capsys.readouterr().err
    assert not (tmp_path / "status").exists()


def test_truth_with_a_matching_run_json_identity_is_accepted(tmp_path):
    setup_truth(tmp_path)
    assert main(truth_args(tmp_path)) == 0


def test_phasec_accepts_the_real_signalp_record_and_a_call_in_it_cannot_spoof_the_call(
    tmp_path,
):
    from cellsurface_sorting_hat.cli import _MEASURED_CALL

    write_r0(tmp_path, signalp="module=signalp/6-gpu; call=cell_wall_adhesion_candidate")
    assert run_phasec(tmp_path, phasec_inputs(tmp_path)) == 0
    notes = _phasec_entries(tmp_path)["S1:Scer_SGD"].measure["notes"]
    assert _MEASURED_CALL.findall(notes)[-1] == "signal_peptide_protein"
    assert run_phasec(tmp_path, phasec_inputs(tmp_path)) == 0


def test_phasec_refuses_a_record_with_another_module_token(tmp_path, capsys):
    write_r0(tmp_path, signalp="module=signalp/6-cpu;signalp6_version=6.0h")
    assert run_phasec(tmp_path, phasec_inputs(tmp_path)) == 2
    assert "signalp/6-cpu" in capsys.readouterr().err
    assert not (tmp_path / "status").exists()


def test_a_call_in_user_notes_does_not_spoof_the_measured_call(tmp_path):
    from cellsurface_sorting_hat.cli import StatusResolver
    from cellsurface_sorting_hat.status import ModuleIdentity
    from cellsurface_sorting_hat.taxonomy import Lineage

    setup_truth(tmp_path)
    assert main(truth_args(tmp_path, "--notes", "call=cell_wall_adhesion_candidate")) == 0
    nodes = write_nodes(tmp_path)
    rec = json.loads((tmp_path / "modules" / "pfam_adhesion.json").read_text())
    ident = {
        "pfam_adhesion": ModuleIdentity(
            "pfam_adhesion", rec["version"], rec["params_hash"], rec["artefact_hash"]
        )
    }
    resolver = StatusResolver(tmp_path, ident, Lineage.from_nodes_dmp(nodes))
    assert resolver.measured_call("pfam_adhesion", 4932) == "wall_family_domain"
