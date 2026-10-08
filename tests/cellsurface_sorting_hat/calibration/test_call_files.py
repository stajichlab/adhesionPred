"""Writer of call status files: validation, merge, lock, atomicity."""

import fcntl
import json

import pytest

from cellsurface_sorting_hat.calibration import call_files
from cellsurface_sorting_hat.calibration.measure import make_entry
from cellsurface_sorting_hat.call_status import load_call_source
from cellsurface_sorting_hat.engine import load_config

CALL = "tandem_repeat_protein"
CFG_SHA = "c" * 64


def measure(calibration_set="S1", n_pos=10, leakage="none"):
    return {
        "calibration_set": calibration_set,
        "truth_source": "truth.tsv",
        "n_pos": n_pos,
        "n_neg": 30,
        "n_clusters_pos": 8,
        "n_clusters_neg": 25,
        "sensitivity": {"value": 0.5, "lo": 0.3, "hi": 0.7},
        "specificity": {"value": 0.97, "lo": 0.9, "hi": 1.0},
        "notes": f"call={CALL}; leakage: {leakage}",
    }


def entry(taxon=40, calibration_set="S1", **kw):
    return make_entry([taxon], measure(calibration_set, **kw), source="truth.tsv")


@pytest.fixture
def cfg():
    return load_config()


@pytest.fixture
def workdir(tmp_path):
    folder = tmp_path / "modules"
    folder.mkdir()
    for m in ("repeat02", "repeat14"):
        # a numeric version, as an older module record may hold: it is written as a string
        (folder / f"{m}.json").write_text(
            json.dumps({"module": m, "version": 1, "params_hash": "p", "artefact_hash": "a"})
        )
    return tmp_path


def path_of(workdir):
    return workdir / "status" / "calls" / f"{CALL}.json"


def write(workdir, cfg, entries, sha=CFG_SHA, call=CALL, variant=""):
    return call_files.write_call_status(workdir, cfg, call, variant, sha, entries)


def test_writes_a_file_that_loads_and_has_string_identities(cfg, workdir):
    path = write(workdir, cfg, [entry()])
    src = load_call_source(path, cfg)
    assert [r.version for r in src.reads] == ["1", "1"]
    assert [e.status for e in src.entries] == ["smoke"]
    assert json.loads(path.read_text())["reads"][0]["version"] == "1"


def test_another_calibration_set_is_added_and_the_same_set_is_replaced(cfg, workdir):
    write(workdir, cfg, [entry(40, "S1", n_pos=10)])
    write(workdir, cfg, [entry(41, "S2", n_pos=12)])
    src = load_call_source(path_of(workdir), cfg)
    assert sorted(e.measure["calibration_set"] for e in src.entries) == ["S1", "S2"]
    write(workdir, cfg, [entry(40, "S1", n_pos=15)])
    src = load_call_source(path_of(workdir), cfg)
    sets = {e.measure["calibration_set"]: e.measure["n_pos"] for e in src.entries}
    assert sets == {"S1": 15, "S2": 12}


def test_a_changed_module_identity_drops_the_old_entries_with_a_message(cfg, workdir, capsys):
    write(workdir, cfg, [entry(40, "S1")])
    record = workdir / "modules" / "repeat02.json"
    record.write_text(
        json.dumps({"module": "repeat02", "version": "2", "params_hash": "p", "artefact_hash": "a"})
    )
    write(workdir, cfg, [entry(41, "S2")])
    assert "dropped" in capsys.readouterr().err
    src = load_call_source(path_of(workdir), cfg)
    assert [e.measure["calibration_set"] for e in src.entries] == ["S2"]


def test_a_changed_config_drops_the_old_entries(cfg, workdir, capsys):
    write(workdir, cfg, [entry(40, "S1")])
    write(workdir, cfg, [entry(41, "S2")], sha="d" * 64)
    assert "dropped" in capsys.readouterr().err
    src = load_call_source(path_of(workdir), cfg)
    assert [e.measure["calibration_set"] for e in src.entries] == ["S2"]
    assert src.config_sha256 == "d" * 64


def test_a_changed_call_hash_drops_the_old_entries(cfg, workdir, capsys, monkeypatch):
    write(workdir, cfg, [entry(40, "S1")])
    monkeypatch.setattr(call_files, "call_hash", lambda *a, **k: "f" * 64)
    write(workdir, cfg, [entry(41, "S2")])
    assert "dropped" in capsys.readouterr().err


def test_the_lock_is_held_while_the_file_is_written(cfg, workdir, monkeypatch):
    seen = {}
    original = call_files.write_atomic
    lock_path = workdir / "status" / "calls" / f"{CALL}.json.lock"

    def probe(path, data):
        with open(lock_path, "a") as other:
            try:
                fcntl.flock(other, fcntl.LOCK_EX | fcntl.LOCK_NB)
                seen["blocked"] = False
                fcntl.flock(other, fcntl.LOCK_UN)
            except BlockingIOError:
                seen["blocked"] = True
        original(path, data)

    monkeypatch.setattr(call_files, "write_atomic", probe)
    write(workdir, cfg, [entry()])
    assert seen == {"blocked": True}
    with open(lock_path, "a") as after:  # released afterwards
        fcntl.flock(after, fcntl.LOCK_EX | fcntl.LOCK_NB)


@pytest.mark.parametrize(
    "call", ["wall_family_domain", "cell_wall_adhesion_candidate", "other_not_surface", "nope"]
)
def test_ineligible_calls_are_refused(cfg, workdir, call):
    with pytest.raises(ValueError, match="eligible|unknown call"):
        write(workdir, cfg, [entry()], call=call)


def test_a_variant_for_a_plain_call_is_refused(cfg, workdir):
    with pytest.raises(ValueError, match="variant"):
        write(workdir, cfg, [entry()], variant="R0")


def test_an_entry_for_taxon_1_is_refused(cfg, workdir):
    with pytest.raises(ValueError, match="taxa"):
        write(workdir, cfg, [entry(1)])


def test_a_missing_module_record_is_refused(cfg, workdir):
    (workdir / "modules" / "repeat14.json").unlink()
    with pytest.raises(ValueError, match="repeat14"):
        write(workdir, cfg, [entry()])


def test_an_invalid_update_leaves_the_existing_file_untouched(cfg, workdir):
    path = write(workdir, cfg, [entry(40, "S1")])
    before = path.read_bytes()
    with pytest.raises(ValueError):
        write(workdir, cfg, [entry(1, "S9")])
    assert path.read_bytes() == before
