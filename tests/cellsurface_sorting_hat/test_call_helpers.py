"""Call helpers of the engine: the modules a call reads, eligibility for a call status file, call_hash."""

import copy

import pytest
import yaml

from cellsurface_sorting_hat import engine
from cellsurface_sorting_hat.engine import (
    call_eligible,
    call_hash,
    load_config,
    reads_of_call,
)


def _custom(tmp_path, mutate):
    """The packaged config, changed by ``mutate(data)``, written to tmp_path and loaded."""
    data = yaml.safe_load(open(engine.__file__.replace("engine.py", "categories.yaml")))
    data = copy.deepcopy(data)
    mutate(data)
    path = tmp_path / "categories.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return load_config(path)


def test_reads_of_the_two_multi_module_calls():
    cfg = load_config()
    assert reads_of_call(cfg, "tandem_repeat_protein") == ["repeat02", "repeat14"]
    assert reads_of_call(cfg, "iuis_allergen_homolog") == ["allergen_homology", "pfam_allergen"]


def test_reads_of_a_per_variant_call_with_a_variant_keeps_only_that_variant():
    cfg = load_config()
    assert reads_of_call(cfg, "signal_peptide_protein", "R0") == ["step1_rule@R0"]
    assert len(reads_of_call(cfg, "signal_peptide_protein")) == 4  # no variant: every step 1 module


@pytest.mark.parametrize(
    ("call", "ok", "reason"),
    [
        ("tandem_repeat_protein", True, ""),
        ("iuis_allergen_homolog", True, ""),
        ("cell_wall_adhesion_candidate", False, "contains ref"),
        ("other_not_surface", False, "kind other"),
        ("signal_peptide_protein", False, "reads fewer than two modules"),
        ("wall_family_domain", False, "reads fewer than two modules"),
        ("no_such_call", False, "unknown call"),
    ],
)
def test_call_eligible_for_the_packaged_calls(call, ok, reason):
    assert call_eligible(load_config(), call) == (ok, reason)


def test_a_per_variant_call_with_step1_and_one_more_module_is_eligible(tmp_path):
    def mutate(d):
        d["calls"].append(
            {
                "name": "gated_pair",
                "per_variant": True,
                "expr": {"and": [{"call": "{step1}"}, {"call": "repeat02"}]},
            }
        )

    cfg = _custom(tmp_path, mutate)
    assert call_eligible(cfg, "gated_pair") == (True, "")


def test_a_call_that_names_a_step1_module_literally_is_not_eligible(tmp_path):
    def mutate(d):
        d["calls"].append(
            {
                "name": "literal_pair",
                "expr": {"or": [{"call": "step1_rule@R0"}, {"call": "repeat02"}]},
            }
        )

    cfg = _custom(tmp_path, mutate)
    assert call_eligible(cfg, "literal_pair") == (False, "literal step1 module")


def test_a_two_module_call_with_a_ref_is_not_eligible(tmp_path):
    def mutate(d):
        d["calls"].append(
            {
                "name": "ref_pair",
                "expr": {"or": [{"ref": "tandem_repeat_protein"}, {"call": "pfam_allergen"}]},
            }
        )

    cfg = _custom(tmp_path, mutate)
    assert call_eligible(cfg, "ref_pair") == (False, "contains ref")


def test_call_hash_is_stable():
    cfg = load_config()
    assert call_hash(cfg, "tandem_repeat_protein") == call_hash(cfg, "tandem_repeat_protein")
    assert len(call_hash(cfg, "tandem_repeat_protein")) == 64


def test_call_hash_changes_with_a_referenced_threshold_only(tmp_path):
    base = call_hash(load_config(), "iuis_allergen_homolog")
    changed = _custom(tmp_path, lambda d: d["thresholds"].update(allergen_identity_min=71))
    assert call_hash(changed, "iuis_allergen_homolog") != base
    other = _custom(tmp_path, lambda d: d["thresholds"].update(antigen_percentile_max=10))
    assert call_hash(other, "iuis_allergen_homolog") == base


def test_call_hash_changes_with_the_engine_semantics(monkeypatch):
    cfg = load_config()
    before = call_hash(cfg, "tandem_repeat_protein")
    monkeypatch.setattr(engine, "ENGINE_SEMANTICS", "2")
    assert call_hash(cfg, "tandem_repeat_protein") != before


def test_call_hash_changes_with_the_expression(tmp_path):
    base = call_hash(load_config(), "tandem_repeat_protein")

    def mutate(d):
        for c in d["calls"]:
            if c["name"] == "tandem_repeat_protein":
                c["expr"] = {"and": c["expr"]["or"]}

    assert call_hash(_custom(tmp_path, mutate), "tandem_repeat_protein") != base


def test_call_hash_differs_between_variants_of_a_per_variant_call(tmp_path):
    def mutate(d):
        d["calls"].append(
            {
                "name": "gated_pair",
                "per_variant": True,
                "expr": {"and": [{"call": "{step1}"}, {"call": "repeat02"}]},
            }
        )

    cfg = _custom(tmp_path, mutate)
    assert call_hash(cfg, "gated_pair", "R0") != call_hash(cfg, "gated_pair", "R2")
