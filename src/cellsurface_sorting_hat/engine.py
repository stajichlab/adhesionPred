"""Category engine: evaluate the rules in categories.yaml for each protein.

Inputs are module tables (``ModuleTable``) and a function ``status_of(module, taxon)``. The engine
does no I/O except reading the config.
"""

import copy
import hashlib
import json
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
# Bump by hand when the evaluation rules change (Kleene tables, NA handling). It is part of call_hash,
# so a call status file measured under other rules is stale.
ENGINE_SEMANTICS = "1"
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
            for b in call.get("basis_calls", []):
                if b not in seen:
                    raise ConfigError(f"{name}: basis_calls entry {b!r} is not an earlier call")
                if seen[b]:
                    raise ConfigError(f"{name}: basis_calls entry {b!r} must be an ungated call")
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
    contributors: (
        frozenset  # of (module, call): the call in which the leaf that reads the module sits
    )


def _leaf(value, module, ctx):
    if value == NOT_ASSESSABLE:
        return _Result(value, frozenset())
    return _Result(value, frozenset([(module, ctx["call_name"])]))


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
        return _leaf(value if value in (CALLED, NOT_CALLED) else NOT_ASSESSABLE, module, ctx)
    if kind == "flag":
        module, fld = arg.split(".", 1)
        row = _row(ctx["modules"], module, ctx["protein"])
        raw = row.get(fld, "") if row else ""
        value = {"1": CALLED, "0": NOT_CALLED}.get(str(raw).strip(), NOT_ASSESSABLE)
        return _leaf(value, module, ctx)
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
            return _leaf(NOT_ASSESSABLE, arg["module"], ctx)
        value = CALLED if _OPS[arg["op"]](number, limit) else NOT_CALLED
        return _leaf(value, arg["module"], ctx)
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


def evaluate(
    cfg, protein_ids, taxa, modules, status_of, measured_call_of=None, call_status_of=None
):
    """Evaluate every call for every protein; return a list of ``CallRecord``.

    ``taxa`` maps protein ID to taxon ID. ``modules`` maps module name to ``ModuleTable``.
    ``status_of(module, taxon)`` returns ``(status, basis)``.

    ``measured_call_of(module, taxon)`` returns the call on which the status of the module was
    measured, or None (an entry without ``call=`` in its notes). A status counts for a call only
    when the leaf that reads the module sits in the measured call. In every other call that reads
    the module, the module is ``unvalidated`` with the basis ``module measured on call X``. Without
    ``measured_call_of`` every status counts for every call.

    ``call_status_of(call, variant, taxon)`` (optional) serves call status files, for calls that read
    several modules. It returns None (no file for the call, or the taxon is not tested), ``(status,
    basis)`` (the call was measured as a whole), or ``(None, reason)`` (a file exists and is stale: the
    module statuses are used and the basis says why). A record whose value is not assessable has no
    deciding module and never takes a call status. Contributors are grouped by the leaf call they sit
    in; a leaf with a call status contributes one status and one basis item. The record status is the
    weakest over all items. Without the hook, or when it returns None for every leaf, the result is
    exactly the one without the hook.
    """
    variants = available_variants(cfg, modules)
    status_cache, records = {}, []

    def status_for(module, taxon, leaf_call):
        key = (module, taxon)
        if key not in status_cache:
            measured = measured_call_of(module, taxon) if measured_call_of else None
            status_cache[key] = (*status_of(module, taxon), measured)
        status, basis, measured = status_cache[key]
        if measured is not None and measured != leaf_call:
            return UNVALIDATED, f"module measured on call {measured}"
        return status, basis

    call_cache = {}

    def call_status_for(leaf, label, taxon):
        variant = label if cfg.call_by_name(leaf).get("per_variant") else ""
        key = (leaf, variant, taxon)
        if key not in call_cache:
            call_cache[key] = call_status_of(leaf, variant, taxon)
        return call_cache[key]

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
                    "call_name": call["name"],
                }
                other_basis = ""
                if call.get("kind") == "other":
                    res, other_basis = _eval_other(call, ctx)
                else:
                    res = _eval(call["expr"], ctx)
                    if call.get("basis_calls") and res.value == CALLED:
                        other_basis = ",".join(
                            b for b in call["basis_calls"] if results[(b, "")].value == CALLED
                        )
                results[(call["name"], label)] = res
                if res.contributors:
                    names = sorted(res.contributors)
                    leaf_status = {}
                    if call_status_of is not None:
                        for leaf in {c for _, c in names}:
                            leaf_status[leaf] = call_status_for(leaf, label, taxon)
                    if all(v is None for v in leaf_status.values()):
                        pairs = [status_for(m, taxon, c) for m, c in names]
                        status = weakest(s for s, _ in pairs)
                        basis = ";".join(
                            f"{m}:{b}" for (m, _), (_, b) in zip(names, pairs, strict=True)
                        )
                    else:
                        status, basis = _status_with_calls(names, leaf_status, status_for, taxon)
                else:
                    status, basis = UNVALIDATED, ""
                records.append(
                    CallRecord(pid, call["name"], label, res.value, status, basis, other_basis)
                )
    return records


def _status_with_calls(names, leaf_status, status_for, taxon):
    """Status and basis of a record when at least one leaf call has a call status or a stale file."""
    items, statuses, emitted = [], [], set()
    for module, leaf in names:
        found = leaf_status.get(leaf)
        if found is not None and found[0] is not None:
            # one item for the whole leaf, where its first module would stand
            if leaf not in emitted:
                emitted.add(leaf)
                items.append(found[1])
                statuses.append(found[0])
            continue
        status, basis = status_for(module, taxon, leaf)
        if found is not None:  # a stale file: say why the module statuses are used
            basis = f"call status stale ({found[1]}); {basis}"
        items.append(f"{module}:{basis}")
        statuses.append(status)
    return weakest(statuses), ";".join(items)


def _expand_step1(cfg, names):
    expanded = set()
    for n in names:
        if "{step1}" in n:
            expanded.update(n.replace("{step1}", v) for v in cfg.step1_variants)
        else:
            expanded.add(n)
    return expanded


def _walk_modules(node, names, refs):
    kind, arg = next(iter(node.items()))
    if kind in ("and", "or"):
        for child in arg:
            _walk_modules(child, names, refs)
    elif kind == "not":
        _walk_modules(arg, names, refs)
    elif kind == "call":
        names.add(arg)
    elif kind == "flag":
        names.add(arg.split(".", 1)[0])
    elif kind == "test":
        names.add(arg["module"])
    elif kind == "ref":
        refs.add(arg)


def referenced_modules(cfg):
    """Names of all modules the rules use (``{step1}`` expanded to the step 1 variants)."""
    names = set()
    for call in cfg.calls:
        if call.get("kind") != "other":
            _walk_modules(call["expr"], names, set())
    return sorted(_expand_step1(cfg, names))


def modules_of_call(cfg, name):
    """Names of the modules that one call reads, also through ``ref`` nodes (``{step1}`` expanded)."""
    names, todo, seen = set(), [name], set()
    while todo:
        current = todo.pop()
        if current in seen:
            continue
        seen.add(current)
        call = cfg.call_by_name(current)
        if call.get("kind") == "other":
            continue
        refs = set()
        _walk_modules(call["expr"], names, refs)
        todo.extend(refs)
    return sorted(_expand_step1(cfg, names))


def reads_of_call(cfg, name, variant=""):
    """Modules that a call reads. With ``variant``, a step 1 module of another variant is left out."""
    reads = modules_of_call(cfg, name)
    if variant:
        reads = [m for m in reads if m not in cfg.step1_variants or variant_label(m) == variant]
    return reads


def _has_ref(node):
    kind, arg = next(iter(node.items()))
    if kind == "ref":
        return True
    if kind in ("and", "or"):
        return any(_has_ref(child) for child in arg)
    if kind == "not":
        return _has_ref(arg)
    return False


def _variant_labels(cfg, call):
    return [variant_label(v) for v in cfg.step1_variants] if call.get("per_variant") else [""]


def call_eligible(cfg, name):
    """``(True, "")`` when a call may have a call status file, else ``(False, reason)``.

    Eligible: the call has an expression, is not ``kind: other``, has no ``ref`` node, names no step 1
    module literally, and reads at least two modules for every variant it has.
    """
    try:
        call = cfg.call_by_name(name)
    except ConfigError:
        return False, "unknown call"
    if call.get("kind") == "other":
        return False, "kind other"
    if _has_ref(call["expr"]):
        return False, "contains ref"
    raw = set()
    _walk_modules(call["expr"], raw, set())
    if raw & set(cfg.step1_variants):
        return False, "literal step1 module"
    # Every label gives the same count while literal step 1 names are refused above; the loop is a
    # guard in case that rule is ever relaxed.
    for label in _variant_labels(cfg, call):
        if len(reads_of_call(cfg, name, label)) < 2:
            return False, "reads fewer than two modules"
    return True, ""


def _resolve_node(node, cfg, variant_module):
    """A copy of ``node`` with threshold references replaced by their numbers and {step1} by the variant."""
    kind, arg = next(iter(node.items()))
    if kind in ("and", "or"):
        return {kind: [_resolve_node(c, cfg, variant_module) for c in arg]}
    if kind == "not":
        return {kind: _resolve_node(arg, cfg, variant_module)}
    arg = copy.deepcopy(arg)
    if kind == "test":
        if isinstance(arg["value"], str) and arg["value"].startswith("$"):
            arg["value"] = cfg.thresholds[arg["value"][1:]]
        arg["module"] = arg["module"].replace("{step1}", variant_module or "{step1}")
    elif kind in ("call", "flag"):
        arg = arg.replace("{step1}", variant_module or "{step1}")
    return {kind: arg}


def call_hash(cfg, name, variant=""):
    """SHA-256 of what defines one call: its expression with the referenced thresholds as numbers,
    the step 1 module of the variant, ``per_variant``, and ``ENGINE_SEMANTICS``."""
    call = cfg.call_by_name(name)
    variant_module = ""
    if variant:
        variant_module = next(v for v in cfg.step1_variants if variant_label(v) == variant)
    payload = {
        "name": name,
        "per_variant": bool(call.get("per_variant")),
        "variant": variant,
        "variant_module": variant_module,
        "expr": _resolve_node(call["expr"], cfg, variant_module),
        "semantics": ENGINE_SEMANTICS,
    }
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode()).hexdigest()


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
