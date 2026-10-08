"""`calibrate truth --call-status`: a status for a call that reads several modules."""

import csv
import gzip
import json

import pytest

from cellsurface_sorting_hat.calibration.cli import main
from cellsurface_sorting_hat.call_status import CallStatusResolver, load_call_source
from cellsurface_sorting_hat.engine import load_config
from cellsurface_sorting_hat.modules.base import ModuleSpec, write_module
from cellsurface_sorting_hat.status import ModuleIdentity
from cellsurface_sorting_hat.taxonomy import Lineage

CALL = "tandem_repeat_protein"
READS = ("repeat02", "repeat14")
NODE_ROWS = [(1, 1, "no rank"), (5052, 1, "genus"), (4932, 1, "species"), (9999, 4932, "strain")]


def write_nodes(tmp_path):
    path = tmp_path / "nodes.dmp"
    path.write_text("".join(f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\n" for t, p, r in NODE_ROWS))
    return path


def write_run_json(folder, states=None, config_sha=None, drop=()):
    """``run.json`` as a run writes it: module identities, run states and the config hash."""
    identities, found = [], {}
    for p in sorted((folder / "modules").glob("*.json")):
        rec = json.loads(p.read_text())
        found[rec["module"]] = rec.get("run_state", "ok")
        identities.append(
            {
                "name": rec["module"],
                "version": rec["version"],
                "params_hash": rec["params_hash"],
                "artefact_hash": rec["artefact_hash"],
            }
        )
    data = {
        "module_identities": identities,
        "module_states": {**found, **(states or {})},
        "config_sha256": config_sha or load_config().sha256,
    }
    for key in drop:
        del data[key]
    (folder / "run.json").write_text(json.dumps(data))


def setup(
    tmp_path, states=None, drop=(), config_sha=None, n_pos=20, n_neg=20, called_pos=12, called_neg=1
):
    for m in READS:
        write_module(tmp_path, ModuleSpec(m, "1"), [], [{"id": "A", "state": "ok"}])
    ids = [f"P{k}" for k in range(n_pos)] + [f"N{k}" for k in range(n_neg)]
    rows = [
        (f"P{k}", CALL, "", "called" if k < called_pos else "not_called") for k in range(n_pos)
    ] + [(f"N{k}", CALL, "", "called" if k < called_neg else "not_called") for k in range(n_neg)]
    write_run_json(tmp_path, states, config_sha, drop)
    with gzip.open(tmp_path / "proteins.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["id", "sha256", "taxon", "state", "note", "trailing_stop", "ambiguous_fraction"]
        )
        for pid in ids:
            w.writerow([pid, "x", 9999, "ok", "", 0, "0.0000"])
    with gzip.open(tmp_path / "c.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein", "call", "variant", "value", "status", "status_basis", "other_basis"])
        for p, c, v, val in rows:
            w.writerow([p, c, v, val, "unvalidated", "", ""])
    (tmp_path / "t.tsv").write_text(
        "id\tlabel\tcluster\n"
        + "".join(f"P{k}\t1\tp{k}\n" for k in range(n_pos))
        + "".join(f"N{k}\t0\tn{k}\n" for k in range(n_neg))
    )


def args(tmp_path, *extra, call=CALL, taxa="4932", leakage="none", module_flag=("--call-status",)):
    return [
        "truth",
        "--workdir",
        str(tmp_path),
        *module_flag,
        "--calls-long",
        str(tmp_path / "c.tsv.gz"),
        "--call",
        call,
        "--truth",
        str(tmp_path / "t.tsv"),
        "--calibration-set",
        "toy",
        "--nodes-dmp",
        str(write_nodes(tmp_path)),
        "--taxa",
        *taxa.split(),
        "--n-boot",
        "200",
        "--leakage",
        leakage,
        *extra,
    ]


def call_file(tmp_path):
    return tmp_path / "status" / "calls" / f"{CALL}.json"


def test_a_call_status_is_written_with_the_measured_rates_and_the_leakage_cap(tmp_path):
    setup(tmp_path)
    assert main(args(tmp_path, leakage="tuned_on_truth")) == 0
    cfg = load_config()
    src = load_call_source(call_file(tmp_path), cfg)
    entry = src.entries[0]
    assert entry.taxa == (4932,)
    assert entry.status == "smoke"  # 20 positives would allow more; tuned_on_truth caps it
    m = entry.measure
    assert (m["n_pos"], m["n_neg"]) == (20, 20)
    assert m["sensitivity"]["value"] == pytest.approx(12 / 20)
    assert m["specificity"]["value"] == pytest.approx(19 / 20)
    assert "call=tandem_repeat_protein" in m["notes"] and "reads=repeat02,repeat14" in m["notes"]
    assert "leakage: tuned_on_truth" in m["notes"]
    assert [r.name for r in src.reads] == list(READS)


def test_the_written_file_is_used_for_a_strain_of_the_tested_species(tmp_path):
    setup(tmp_path)
    assert main(args(tmp_path, leakage="tuned_on_truth")) == 0
    cfg = load_config()
    ids = {m: ModuleIdentity(m, "1", *_hashes(tmp_path, m)) for m in READS}
    lineage = Lineage({1: 1, 4932: 1, 9999: 4932, 5052: 1})
    r = CallStatusResolver(tmp_path, cfg, ids, dict.fromkeys(READS, "ok"), cfg.sha256, lineage)
    assert r(CALL, "", 9999) == ("smoke", "call:tandem_repeat_protein:taxon:4932")
    assert r(CALL, "", 5052) is None


def _hashes(tmp_path, module):
    rec = json.loads((tmp_path / "modules" / f"{module}.json").read_text())
    return rec["params_hash"], rec["artefact_hash"]


@pytest.mark.parametrize("state", ["unavailable", "partial", "error", "not_run"])
def test_a_read_module_that_is_not_ok_is_refused(tmp_path, capsys, state):
    setup(tmp_path, states={"repeat14": state})
    assert main(args(tmp_path)) == 2
    err = capsys.readouterr().err
    assert "repeat14" in err and state in err
    assert not call_file(tmp_path).exists()


def test_a_run_json_without_module_states_is_refused(tmp_path, capsys):
    setup(tmp_path, drop=("module_states",))
    assert main(args(tmp_path)) == 2
    assert "module_states" in capsys.readouterr().err


def test_a_run_json_without_a_config_hash_is_refused(tmp_path, capsys):
    setup(tmp_path, drop=("config_sha256",))
    assert main(args(tmp_path)) == 2
    assert "config_sha256" in capsys.readouterr().err


def test_a_run_made_with_another_config_is_refused(tmp_path, capsys):
    setup(tmp_path, config_sha="e" * 64)
    assert main(args(tmp_path)) == 2
    assert "config" in capsys.readouterr().err
    assert not call_file(tmp_path).exists()


def test_a_module_record_that_changed_after_the_run_is_refused(tmp_path, capsys):
    setup(tmp_path)
    record = tmp_path / "modules" / "repeat02.json"
    data = json.loads(record.read_text())
    data["version"] = "2"
    record.write_text(json.dumps(data))
    assert main(args(tmp_path)) == 2
    assert "repeat02" in capsys.readouterr().err


def test_a_genus_taxon_is_refused(tmp_path, capsys):
    setup(tmp_path)
    assert main(args(tmp_path, taxa="5052")) == 2
    assert "taxon" in capsys.readouterr().err
    assert not call_file(tmp_path).exists()


@pytest.mark.parametrize(
    "call",
    [
        "wall_family_domain",
        "signal_peptide_protein",
        "cell_wall_adhesion_candidate",
        "other_not_surface",
    ],
)
def test_calls_that_cannot_have_a_call_file_are_refused(tmp_path, capsys, call):
    setup(tmp_path)
    assert main(args(tmp_path, call=call)) == 2
    assert "eligible" in capsys.readouterr().err
    assert not call_file(tmp_path).exists()


def test_exactly_one_of_module_and_call_status_is_required(tmp_path):
    setup(tmp_path)
    with pytest.raises(SystemExit) as err:
        main(args(tmp_path, module_flag=()))
    assert err.value.code == 2
    with pytest.raises(SystemExit) as err:
        main(args(tmp_path, module_flag=("--call-status", "--module", "repeat02")))
    assert err.value.code == 2


def test_the_module_path_still_refuses_a_call_that_reads_several_modules(tmp_path, capsys):
    setup(tmp_path)
    assert main(args(tmp_path, module_flag=("--module", "repeat02"))) == 2
    assert "reads 2 modules" in capsys.readouterr().err
