"""Helpers for engine-level tests of PR 3: module tables, a run function, config variants."""

import copy

import yaml

from cellsurface_sorting_hat.engine import ModuleTable, evaluate, load_config

R0 = "step1_rule@R0"


def table(name, rows):
    return ModuleTable(name, {k: dict(v) for k, v in rows.items()})


def ok(**fields):
    return {"state": "ok", **fields}


def base_modules(step1="not_called", **over):
    """One protein ``P``: every module present, every call not called, antigen rank not in the top."""
    mods = {
        R0: table(R0, {"P": ok(call=step1)}),
        "repeat02": table("repeat02", {"P": ok(call="not_called")}),
        "repeat14": table("repeat14", {"P": ok(call="not_called")}),
        "pfam_adhesion": table("pfam_adhesion", {"P": ok(hit="0", hit_narrow="0")}),
        "pfam_hydrophobin": table("pfam_hydrophobin", {"P": ok(hit="0")}),
        "pfam_hsba": table("pfam_hsba", {"P": ok(hit="0")}),
        "pfam_allergen": table("pfam_allergen", {"P": ok(hit="0")}),
        "antigen_lookup": table("antigen_lookup", {"P": ok(percentile="50")}),
        "allergen_homology": table("allergen_homology", {"P": ok(identity="0", coverage="0")}),
        "gpi": table("gpi", {"P": ok(call="not_called")}),
        "disorder": table("disorder", {"P": ok(longest_run="0")}),
        "aga1_model": table("aga1_model", {"P": ok(call="not_called")}),
    }
    mods.update(over)
    return mods


def run_eval(modules, status_of=None, ids=("P",), taxon=40, cfg=None):
    cfg = cfg or load_config()
    status_of = status_of or (lambda module, t: ("unvalidated", "x"))
    records = evaluate(cfg, list(ids), dict.fromkeys(ids, taxon), modules, status_of)
    return {(r.protein, r.call, r.variant): r for r in records}


def config_with(tmp_path, extra_calls=(), mechanism_extra=()):
    """A copy of the packaged config with ``extra_calls`` inserted before the ``other_*`` calls.

    ``mechanism_extra`` names calls that are added to the mechanism list of both ``other_*`` calls.
    The copy is written to ``tmp_path`` and loaded with ``load_config``."""
    from pathlib import Path

    import cellsurface_sorting_hat.engine as engine_module

    data = yaml.safe_load(
        (Path(engine_module.__file__).with_name("categories.yaml")).read_text(encoding="utf-8")
    )
    calls = copy.deepcopy(data["calls"])
    first_other = next(i for i, c in enumerate(calls) if c.get("kind") == "other")
    calls[first_other:first_other] = [copy.deepcopy(c) for c in extra_calls]
    for c in calls:
        if c.get("kind") == "other":
            c["mechanism"] = list(c["mechanism"]) + list(mechanism_extra)
    data["calls"] = calls
    path = tmp_path / "categories.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return load_config(path)
