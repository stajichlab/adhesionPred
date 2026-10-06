"""Category engine: evaluate the rules in categories.yaml for each protein.

Inputs are module tables (``ModuleTable``) and a function ``status_of(module, taxon)``. The engine
does no I/O except reading the config.
"""

import hashlib
import math
import operator
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from cellsurface_sorting_hat.logic import (
    CALLED,
    NOT_ASSESSABLE,
    NOT_CALLED,
    k_and,
    k_not,
    k_or,
)
from cellsurface_sorting_hat.status import UNVALIDATED, weakest

OK_STATE = "ok"
_OPS = {">=": operator.ge, "<=": operator.le, ">": operator.gt, "<": operator.lt, "==": operator.eq}
_NODE_KEYS = {"call", "flag", "test", "ref", "and", "or", "not"}


class ConfigError(ValueError):
    """categories.yaml is not valid."""


@dataclass
class ModuleTable:
    """Per-protein rows of one module: id -> dict of strings (``state``, ``call``, other fields)."""

    name: str
    rows: dict = field(default_factory=dict)

    def get(self, protein_id):
        return self.rows.get(protein_id)


@dataclass
class Config:
    default_gate: str
    step1_variants: list
    thresholds: dict
    calls: list
    evidence: list
    sha256: str

    def call_by_name(self, name):
        for c in self.calls:
            if c["name"] == name:
                return c
        raise ConfigError(f"unknown call: {name}")


@dataclass(frozen=True)
class CallRecord:
    protein: str
    call: str
    variant: str
    value: str
    status: str
    status_basis: str
    other_basis: str


def variant_label(step1_module):
    """``step1_rule@R0`` -> ``R0``."""
    return step1_module.split("@", 1)[1] if "@" in step1_module else step1_module


def load_config(path=None):
    path = Path(path) if path else Path(__file__).with_name("categories.yaml")
    raw = path.read_bytes()
    data = yaml.safe_load(raw)
    cfg = Config(
        default_gate=data["default_gate"],
        step1_variants=list(data["step1_variants"]),
        thresholds=dict(data["thresholds"]),
        calls=list(data["calls"]),
        evidence=list(data.get("evidence", [])),
        sha256=hashlib.sha256(raw).hexdigest(),
    )
    _validate(cfg)
    return cfg


def _validate(cfg):
    if cfg.default_gate not in cfg.step1_variants:
        raise ConfigError("default_gate must be one of step1_variants")
    for item in cfg.evidence:
        if "." not in item:
            raise ConfigError(f"evidence entry needs MODULE.FIELD: {item!r}")
    seen = {}
    for call in cfg.calls:
        name = call["name"]
        if name in seen:
            raise ConfigError(f"duplicate call name: {name}")
        per_variant = bool(call.get("per_variant"))
        if call.get("kind") == "other":
            if not per_variant or not seen.get(call.get("surface")):
                raise ConfigError(
                    f"{name}: an `other` call needs an earlier per_variant surface call"
                )
            if call.get("surface_is") not in (CALLED, NOT_CALLED):
                raise ConfigError(f"{name}: surface_is must be called or not_called")
            for m in call["mechanism"]:
                if m not in seen or seen[m]:
                    raise ConfigError(f"{name}: mechanism {m} must be an earlier ungated call")
        else:
            _validate_node(call["expr"], seen, per_variant, cfg.thresholds, name)
        seen[name] = per_variant


def _check_module_name(module, per_variant, where):
    if "{step1}" in module and not per_variant:
        raise ConfigError(f"{where}: {{step1}} is only allowed in a per_variant call")


def _validate_node(node, seen, per_variant, thresholds, where):
    if not isinstance(node, dict) or len(node) != 1 or next(iter(node)) not in _NODE_KEYS:
        raise ConfigError(f"{where}: bad node {node!r}")
    kind, arg = next(iter(node.items()))
    if kind in ("and", "or"):
        if not isinstance(arg, list) or not arg:
            raise ConfigError(f"{where}: {kind} needs a non-empty list")
        for child in arg:
            _validate_node(child, seen, per_variant, thresholds, where)
    elif kind == "not":
        _validate_node(arg, seen, per_variant, thresholds, where)
    elif kind == "ref":
        if arg not in seen:
            raise ConfigError(f"{where}: ref to unknown or later call {arg}")
        if seen[arg] and not per_variant:
            raise ConfigError(f"{where}: an ungated call cannot use per_variant call {arg}")
    elif kind == "test":
        _check_module_name(arg["module"], per_variant, where)
        if arg["op"] not in _OPS:
            raise ConfigError(f"{where}: bad operator {arg['op']!r}")
        value = arg["value"]
        if isinstance(value, str) and (not value.startswith("$") or value[1:] not in thresholds):
            raise ConfigError(f"{where}: unknown threshold {value!r}")
    elif kind == "flag":
        if "." not in arg:
            raise ConfigError(f"{where}: flag needs MODULE.FIELD")
        _check_module_name(arg.split(".", 1)[0], per_variant, where)
    elif kind == "call":
        _check_module_name(arg, per_variant, where)


@dataclass
class _Result:
    value: str
    contributors: frozenset


def _leaf(value, module):
    return _Result(value, frozenset([module]) if value != NOT_ASSESSABLE else frozenset())


def _row(modules, module, protein_id):
    table = modules.get(module)
    if table is None:
        return None
    row = table.get(protein_id)
    if row is None or row.get("state", OK_STATE) != OK_STATE:
        return None
    return row


def _eval(node, ctx):
    kind, arg = next(iter(node.items()))
    if kind == "call":
        module = arg.replace("{step1}", ctx["step1"] or "")
        row = _row(ctx["modules"], module, ctx["protein"])
        value = row.get("call") if row else None
        return _leaf(value if value in (CALLED, NOT_CALLED) else NOT_ASSESSABLE, module)
    if kind == "flag":
        module, fld = arg.split(".", 1)
        row = _row(ctx["modules"], module, ctx["protein"])
        raw = row.get(fld, "") if row else ""
        value = {"1": CALLED, "0": NOT_CALLED}.get(str(raw).strip(), NOT_ASSESSABLE)
        return _leaf(value, module)
    if kind == "test":
        row = _row(ctx["modules"], arg["module"], ctx["protein"])
        limit = arg["value"]
        if isinstance(limit, str):
            limit = ctx["thresholds"][limit[1:]]
        try:
            number = float(row[arg["field"]]) if row else None
        except (KeyError, TypeError, ValueError):
            number = None
        if number is None or not math.isfinite(number):  # missing, empty, NaN or infinite
            return _leaf(NOT_ASSESSABLE, arg["module"])
        return _leaf(CALLED if _OPS[arg["op"]](number, limit) else NOT_CALLED, arg["module"])
    if kind == "ref":
        callee = ctx["cfg"].call_by_name(arg)
        label = ctx["label"] if callee.get("per_variant") else ""
        return ctx["results"][(arg, label)]
    if kind == "not":
        child = _eval(arg, ctx)
        return _Result(k_not(child.value), child.contributors)
    children = [_eval(c, ctx) for c in arg]
    values = [c.value for c in children]
    if kind == "and":
        out = k_and(*values)
        decisive = [c for c in children if out != NOT_CALLED or c.value == NOT_CALLED]
    else:
        out = k_or(*values)
        decisive = [c for c in children if out != CALLED or c.value == CALLED]
    if out == NOT_ASSESSABLE:  # an unknown result has no deciding module
        return _Result(out, frozenset())
    return _Result(out, frozenset().union(*(c.contributors for c in decisive)))


def _eval_other(call, ctx):
    surface = ctx["results"][(call["surface"], ctx["label"])]
    if surface.value == NOT_ASSESSABLE:
        return _Result(NOT_ASSESSABLE, frozenset()), ""
    if surface.value != call["surface_is"]:
        return _Result(NOT_CALLED, surface.contributors), ""
    contributors, called, left_out = set(surface.contributors), set(), []
    for name in call["mechanism"]:
        res = ctx["results"][(name, "")]
        if res.value == CALLED:
            called |= res.contributors  # a false AND takes its status from all false inputs
        elif res.value == NOT_ASSESSABLE:
            left_out.append(name)
        else:
            contributors |= res.contributors
    if called:
        return _Result(NOT_CALLED, frozenset(called)), ""
    return _Result(CALLED, frozenset(contributors)), ",".join(left_out)


def available_variants(cfg, modules):
    return [v for v in cfg.step1_variants if v in modules]


def evaluate(cfg, protein_ids, taxa, modules, status_of):
    """Evaluate every call for every protein; return a list of ``CallRecord``.

    ``taxa`` maps protein ID to taxon ID. ``modules`` maps module name to ``ModuleTable``.
    ``status_of(module, taxon)`` returns ``(status, basis)``.
    """
    variants = available_variants(cfg, modules)
    status_cache, records = {}, []

    def status_for(module, taxon):
        key = (module, taxon)
        if key not in status_cache:
            status_cache[key] = status_of(module, taxon)
        return status_cache[key]

    for pid in protein_ids:
        taxon = taxa[pid]
        results = {}
        for call in cfg.calls:
            if call.get("per_variant"):
                labels = [(v, variant_label(v)) for v in variants]
            else:
                labels = [(None, "")]
            for step1, label in labels:
                ctx = {
                    "protein": pid,
                    "modules": modules,
                    "thresholds": cfg.thresholds,
                    "cfg": cfg,
                    "results": results,
                    "step1": step1,
                    "label": label,
                }
                other_basis = ""
                if call.get("kind") == "other":
                    res, other_basis = _eval_other(call, ctx)
                else:
                    res = _eval(call["expr"], ctx)
                results[(call["name"], label)] = res
                if res.contributors:
                    names = sorted(res.contributors)
                    pairs = [status_for(m, taxon) for m in names]
                    status = weakest(s for s, _ in pairs)
                    basis = ";".join(f"{m}:{b}" for m, (_, b) in zip(names, pairs, strict=True))
                else:
                    status, basis = UNVALIDATED, ""
                records.append(
                    CallRecord(pid, call["name"], label, res.value, status, basis, other_basis)
                )
    return records


def referenced_modules(cfg):
    """Names of all modules the rules use (``{step1}`` expanded to the step 1 variants)."""
    names = set()

    def walk(node):
        kind, arg = next(iter(node.items()))
        if kind in ("and", "or"):
            for child in arg:
                walk(child)
        elif kind == "not":
            walk(arg)
        elif kind == "call":
            names.add(arg)
        elif kind == "flag":
            names.add(arg.split(".", 1)[0])
        elif kind == "test":
            names.add(arg["module"])

    for call in cfg.calls:
        if call.get("kind") != "other":
            walk(call["expr"])
    expanded = set()
    for n in names:
        if "{step1}" in n:
            expanded.update(n.replace("{step1}", v) for v in cfg.step1_variants)
        else:
            expanded.add(n)
    return sorted(expanded)


def collect_evidence(cfg, protein_ids, modules):
    """Rows (protein, module, field, value) for the configured evidence fields."""
    rows = []
    for item in cfg.evidence:
        module, fld = item.split(".", 1)
        table = modules.get(module)
        if table is None:
            continue
        for pid in protein_ids:
            row = table.get(pid)
            if row and row.get("state", OK_STATE) == OK_STATE and str(row.get(fld, "")).strip():
                rows.append((pid, module, fld, str(row[fld]).strip()))
    return rows
