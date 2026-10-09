"""docs/class_descriptions.yaml and docs/CLASSES.md against src/cellsurface_sorting_hat/categories.yaml (M2c of the v1 class list spec)."""

import importlib.util
import re
from pathlib import Path

import yaml

from cellsurface_sorting_hat.engine import load_config

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "make_classes_md", ROOT / "analysis/class_table/make_classes_md.py"
)
gen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gen)
DESC = yaml.safe_load((ROOT / "docs/class_descriptions.yaml").read_text())
SUFFIXES = ("_domain", "_protein", "_ortholog", "_similarity", "_homolog", "_candidate")
EXCEPTIONS = {"cocci_specificity_rank_top15"}


def test_every_call_has_a_description_and_every_description_has_a_call():
    calls = {c["name"] for c in load_config().calls}
    described = {d["call"] for d in DESC["calls"]}
    assert calls - described == set(), "calls without a description"
    assert described - calls == set(), "descriptions without a call"


def test_call_names_follow_the_naming_rule():
    for c in load_config().calls:
        name = c["name"]
        if name in EXCEPTIONS or name.startswith("other_"):
            continue
        assert name.endswith(SUFFIXES), name


def test_composites_end_in_candidate_and_other_calls_start_with_other():
    cfg = load_config()
    for c in cfg.calls:
        kind = c.get("kind")
        if kind == "other":
            assert c["name"].startswith("other_")
        elif "ref" in repr(c["expr"]):
            assert c["name"].endswith("_candidate"), c["name"]


def test_cited_status_files_exist_in_the_repository():
    for d in DESC["calls"]:
        for path in d.get("status_files", []):
            assert (ROOT / path).exists(), path


def test_each_description_has_the_required_fields():
    for d in DESC["calls"]:
        for key in ("call", "meaning", "category", "reads", "status"):
            assert d.get(key), (d["call"], key)


def test_reads_match_the_modules_the_config_uses():
    cfg = load_config()
    from cellsurface_sorting_hat.engine import modules_of_call

    for d in DESC["calls"]:
        call = cfg.call_by_name(d["call"])
        if call.get("kind") == "other":
            continue
        actual = set(modules_of_call(cfg, d["call"]))
        stated = {re.sub(r"\{step1\}", "step1_rule@R0", m) for m in d["reads"]}
        assert actual <= stated | {"step1_ml@card", "step1_rule@R1", "step1_rule@R2"}, (
            d["call"],
            actual,
            stated,
        )


def test_classes_md_is_current():
    assert (ROOT / "docs/CLASSES.md").read_text() == gen.render(DESC, load_config())
