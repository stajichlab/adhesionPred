"""Call status files: loader, validity against a run, resolver (lineage, stale reasons)."""

import json

import pytest

from cellsurface_sorting_hat.call_status import (
    CallStatusError,
    CallStatusResolver,
    call_file_name,
    load_call_source,
    stale_reason,
)
from cellsurface_sorting_hat.engine import call_hash, load_config, reads_of_call
from cellsurface_sorting_hat.status import ModuleIdentity
from cellsurface_sorting_hat.taxonomy import Lineage

CALL = "tandem_repeat_protein"
# 1 root, 10 fungi, 20 class, 40 species, 400 strain below it, 41 another species, 401 strain of 41
LINEAGE = Lineage({1: 1, 10: 1, 20: 10, 40: 20, 400: 40, 41: 20, 401: 41})
CFG_SHA = "c" * 64


def measure(n_pos=10, notes="call=tandem_repeat_protein; leakage: none"):
    return {
        "calibration_set": "S1",
        "truth_source": "truth.tsv",
        "n_pos": n_pos,
        "n_neg": 30,
        "n_clusters_pos": 8,
        "n_clusters_neg": 25,
        "sensitivity": {"value": 0.5, "lo": 0.3, "hi": 0.7},
        "specificity": {"value": 0.97, "lo": 0.9, "hi": 1.0},
        "notes": notes,
    }


def entry(taxa=(40,), status="smoke", **kw):
    return {"taxa": list(taxa), "status": status, "source": "truth.tsv", "measure": measure(**kw)}


def identities(cfg, call=CALL, **changes):
    out = {}
    for m in reads_of_call(cfg, call):
        fields = {"version": "1", "params_hash": "p", "artefact_hash": "a"}
        fields.update(changes.get(m, {}))
        out[m] = ModuleIdentity(m, **fields)
    return out


def payload(cfg, call=CALL, variant="", entries=None, **over):
    data = {
        "call": call,
        "variant": variant,
        "config_sha256": CFG_SHA,
        "call_hash": call_hash(cfg, call, variant),
        "reads": [
            {"name": m, "version": "1", "params_hash": "p", "artefact_hash": "a"}
            for m in reads_of_call(cfg, call, variant)
        ],
        "entries": entries if entries is not None else [entry()],
    }
    data.update(over)
    return data


def write(tmp_path, data, name=None):
    folder = tmp_path / "status" / "calls"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / (name or call_file_name(data["call"], data["variant"]))
    path.write_text(json.dumps(data))
    return path


@pytest.fixture
def cfg():
    return load_config()


def states(cfg, call=CALL, **over):
    out = dict.fromkeys(reads_of_call(cfg, call), "ok")
    out.update(over)
    return out


def test_a_valid_file_loads(cfg, tmp_path):
    src = load_call_source(write(tmp_path, payload(cfg)), cfg)
    assert (src.call, src.variant) == (CALL, "")
    assert [e.status for e in src.entries] == ["smoke"]
    assert [r.name for r in src.reads] == ["repeat02", "repeat14"]


def renamed(cfg, call, variant=""):
    """A valid payload for the repeat call whose call and variant fields are changed (the helpers
    cannot hash an unknown, `other` or single-module call)."""
    data = payload(cfg)
    data.update(call=call, variant=variant)
    return data


def bad_entry(**changes):
    e = entry()
    e.update(changes)
    return e


@pytest.mark.parametrize(
    ("label", "build", "needle"),
    [
        (
            "no measure",
            lambda cfg: payload(cfg, entries=[{"taxa": [40], "status": "smoke"}]),
            "measure",
        ),
        (
            "estimated without room",
            lambda cfg: payload(cfg, entries=[bad_entry(status="estimated")]),
            "stronger",
        ),
        ("taxon 1", lambda cfg: payload(cfg, entries=[bad_entry(taxa=[1])]), "taxa"),
        ("taxon 0", lambda cfg: payload(cfg, entries=[bad_entry(taxa=[0])]), "taxa"),
        ("duplicate taxon", lambda cfg: payload(cfg, entries=[entry((40,)), entry((40,))]), "40"),
        ("unknown call", lambda cfg: renamed(cfg, "nope"), "unknown call"),
        ("single-module call", lambda cfg: renamed(cfg, "wall_family_domain"), "not eligible"),
        (
            "composite call",
            lambda cfg: renamed(cfg, "cell_wall_adhesion_candidate"),
            "not eligible",
        ),
        ("other call", lambda cfg: renamed(cfg, "other_not_surface"), "not eligible"),
        ("variant for a plain call", lambda cfg: renamed(cfg, CALL, "R0"), "variant"),
        ("unknown key", lambda cfg: payload(cfg, extra=1), "extra"),
    ],
)
def test_loader_refusals(cfg, tmp_path, label, build, needle):
    data = build(cfg)
    path = write(tmp_path, data, call_file_name(data["call"], data["variant"]))
    with pytest.raises(ValueError, match=needle):
        load_call_source(path, cfg)


def test_loader_refuses_a_file_whose_name_differs_from_its_content(cfg, tmp_path):
    path = write(tmp_path, payload(cfg), name="iuis_allergen_homolog.json")
    with pytest.raises(ValueError, match="file name"):
        load_call_source(path, cfg)


def test_loader_needs_the_leakage_note(cfg, tmp_path):
    e = entry()
    e["measure"]["notes"] = "call=tandem_repeat_protein"
    with pytest.raises(ValueError, match="leakage"):
        load_call_source(write(tmp_path, payload(cfg, entries=[e])), cfg)


def test_leakage_other_than_none_caps_the_status_at_smoke(cfg, tmp_path):
    big = measure()
    big.update(n_pos=40, n_neg=40, n_clusters_pos=30, n_clusters_neg=30)
    big["sensitivity"] = {"value": 0.9, "lo": 0.85, "hi": 0.95}
    big["notes"] = "call=tandem_repeat_protein; leakage: tuned_on_truth"
    e = {"taxa": [40], "status": "estimated", "source": "t", "measure": big}
    with pytest.raises(ValueError, match="leakage"):
        load_call_source(write(tmp_path, payload(cfg, entries=[e])), cfg)
    big["notes"] = "call=tandem_repeat_protein; leakage: none"
    assert (
        load_call_source(write(tmp_path, payload(cfg, entries=[e])), cfg).entries[0].status
        == "estimated"
    )


def test_a_hand_edit_that_raises_smoke_to_estimated_is_refused(cfg, tmp_path):
    data = payload(cfg)
    data["entries"][0]["status"] = "estimated"
    with pytest.raises(ValueError, match="stronger"):
        load_call_source(write(tmp_path, data), cfg)


def _stale(cfg, tmp_path, ids=None, st=None, sha=CFG_SHA, data=None):
    src = load_call_source(write(tmp_path, data or payload(cfg)), cfg)
    return stale_reason(src, cfg, ids or identities(cfg), st or states(cfg), sha)


def test_a_file_that_matches_the_run_is_not_stale(cfg, tmp_path):
    assert _stale(cfg, tmp_path) == ""


def test_a_changed_config_makes_the_file_stale(cfg, tmp_path):
    assert _stale(cfg, tmp_path, sha="d" * 64) == "config differs"


def test_reads_that_differ_from_the_modules_now_read_make_the_file_stale(cfg, tmp_path):
    data = payload(cfg)
    data["reads"] = data["reads"][:1]
    assert _stale(cfg, tmp_path, data=data) == "reads differ"


def test_a_module_without_a_run_record_makes_the_file_stale(cfg, tmp_path):
    ids = identities(cfg)
    del ids["repeat14"]
    assert _stale(cfg, tmp_path, ids=ids) == "no module run record: repeat14"


@pytest.mark.parametrize("state", ["unavailable", "error", "not_run", "partial", "missing"])
def test_a_module_that_is_not_ok_makes_the_file_stale(cfg, tmp_path, state):
    st = states(cfg)
    if state == "missing":
        del st["repeat14"]
    else:
        st["repeat14"] = state
    assert _stale(cfg, tmp_path, st=st) == f"module state not ok: repeat14 ({state})"


@pytest.mark.parametrize("field", ["version", "params_hash", "artefact_hash"])
def test_a_changed_module_identity_makes_the_file_stale(cfg, tmp_path, field):
    ids = identities(cfg, repeat02={field: "other"})
    assert _stale(cfg, tmp_path, ids=ids) == "identity differs: repeat02"


def test_a_changed_call_hash_makes_the_file_stale(cfg, tmp_path):
    data = payload(cfg, call_hash="0" * 64)
    assert _stale(cfg, tmp_path, data=data) == "call_hash differs"


def resolver(cfg, tmp_path, **kw):
    return CallStatusResolver(
        tmp_path,
        cfg,
        kw.get("ids") or identities(cfg),
        kw.get("st") or states(cfg),
        kw.get("sha", CFG_SHA),
        LINEAGE,
    )


def test_resolver_lineage(cfg, tmp_path):
    write(tmp_path, payload(cfg, entries=[entry((40,))]))
    r = resolver(cfg, tmp_path)
    assert r(CALL, "", 40) == ("smoke", "call:tandem_repeat_protein:taxon:40")
    assert r(CALL, "", 400) == ("smoke", "call:tandem_repeat_protein:taxon:40")  # strain below it
    assert r(CALL, "", 41) is None  # sibling species
    assert r(CALL, "", 401) is None  # strain of the sibling
    assert r(CALL, "", 20) is None  # an ancestor of the tested species


def test_the_most_specific_tested_taxon_wins(cfg, tmp_path):
    e_species = entry((40,), n_pos=10)
    e_strain = entry((400,), n_pos=0, notes="call=tandem_repeat_protein; leakage: none")
    e_strain["status"] = "unvalidated"
    e_strain["measure"]["sensitivity"] = None
    del e_strain["measure"]["sensitivity"]
    write(tmp_path, payload(cfg, entries=[e_species, e_strain]))
    assert resolver(cfg, tmp_path)(CALL, "", 400)[0] == "unvalidated"
    assert resolver(cfg, tmp_path)(CALL, "", 40)[0] == "smoke"


def test_no_file_gives_none(cfg, tmp_path):
    assert resolver(cfg, tmp_path)(CALL, "", 40) is None


def test_a_stale_file_gives_the_reason_and_no_status(cfg, tmp_path):
    write(tmp_path, payload(cfg))
    r = resolver(cfg, tmp_path, ids=identities(cfg, repeat02={"version": "2"}))
    assert r(CALL, "", 40) == (None, "identity differs: repeat02")


def test_a_malformed_file_stops_the_resolver_and_names_the_file(cfg, tmp_path):
    folder = tmp_path / "status" / "calls"
    folder.mkdir(parents=True)
    (folder / "tandem_repeat_protein.json").write_text("{not json")
    with pytest.raises(CallStatusError, match="tandem_repeat_protein.json"):
        resolver(cfg, tmp_path)


def test_an_orphan_file_for_a_call_that_is_not_in_the_config_stops_the_resolver(cfg, tmp_path):
    write(tmp_path, renamed(cfg, "gone"), name="gone.json")
    with pytest.raises(CallStatusError, match="gone.json"):
        resolver(cfg, tmp_path)


def test_rows_describe_valid_and_stale_files(cfg, tmp_path):
    write(tmp_path, payload(cfg, entries=[entry((40,), n_pos=10)]))
    rows = resolver(cfg, tmp_path).rows()
    assert [(r["call"], r["taxon"], r["status"], r["valid"]) for r in rows] == [
        (CALL, 40, "smoke", True)
    ]
    assert rows[0]["calibration_set"] == "S1" and rows[0]["n_pos"] == 10
    stale = resolver(cfg, tmp_path, st=states(cfg, repeat14="unavailable")).rows()
    assert stale[0]["valid"] is False and "module state not ok" in stale[0]["reason"]
