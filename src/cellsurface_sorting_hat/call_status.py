"""Status files for calls that read two or more modules.

A module status file (``status/<module>.json``) records a calibration of one module. A call that reads
several modules (for example ``tandem_repeat_protein`` = ``repeat02`` OR ``repeat14``) is measured as
a whole, so its status lives in ``status/calls/<call>[.<variant>].json``.

A call file is valid for a run only when the config, the modules the call reads, their identities and
run states, and the call definition are the ones it was measured with. A stale file is never used.
A file that is not well formed, or that names a call that does not exist or cannot have a call file,
stops the run (the caller turns ``CallStatusError`` into an input error).
"""

import json
import re
from dataclasses import dataclass
from pathlib import Path

from cellsurface_sorting_hat.engine import (
    call_eligible,
    call_hash,
    reads_of_call,
    variant_label,
)
from cellsurface_sorting_hat.status import (
    _STRENGTH,
    ENTRY_KEYS,
    SMOKE,
    ModuleIdentity,
    StatusEntry,
    _check_measure,
    best_entry,
    status_from_measure,
    weakest,
)

FILE_KEYS = {"call", "variant", "config_sha256", "call_hash", "reads", "entries"}
LEAKAGE_VALUES = ("none", "partial", "tuned_on_truth", "in_reference", "unknown")
_LEAKAGE = re.compile(r"(?:^|[\s;])leakage:\s*([a-z_]+)")
MIN_TAXON = 2  # taxon 1 is the root: a status for it would apply to every protein


class CallStatusError(ValueError):
    """A call status file is not usable."""


@dataclass(frozen=True)
class CallSource:
    path: Path
    call: str
    variant: str
    config_sha256: str
    call_hash: str
    reads: tuple  # of ModuleIdentity
    entries: tuple  # of StatusEntry


def call_file_name(call, variant=""):
    return f"{call}.{variant}.json" if variant else f"{call}.json"


def leakage_of(measure):
    found = _LEAKAGE.findall((measure or {}).get("notes", ""))
    return found[-1] if found else None


def _entries(raw, where):
    if not isinstance(raw, list) or not raw:
        raise ValueError(f"{where}: entries must be a non-empty list")
    out, listed = [], {}
    for i, e in enumerate(raw):
        at = f"{where} entries[{i}]"
        if not isinstance(e, dict) or not {"taxa", "status", "measure"} <= set(e) <= ENTRY_KEYS:
            raise ValueError(
                f"{at}: needs taxa, status and measure, and only keys {sorted(ENTRY_KEYS)}"
            )
        if e["status"] not in _STRENGTH:
            raise ValueError(f"{at}: unknown status {e['status']!r}")
        taxa = e["taxa"]
        if (
            not isinstance(taxa, list)
            or not taxa
            or any(isinstance(t, bool) or not isinstance(t, int) or t < MIN_TAXON for t in taxa)
        ):
            raise ValueError(f"{at}: taxa must be a list of integers of at least {MIN_TAXON}")
        for t in taxa:
            if t in listed:
                raise ValueError(f"{at}: taxon {t} is already in entries[{listed[t]}]")
            listed[t] = i
        measure = e["measure"]
        _check_measure(measure, at)
        allowed = status_from_measure(measure)
        if weakest([e["status"], allowed]) != e["status"]:
            raise ValueError(
                f"{at}: status {e['status']} is stronger than the measure allows ({allowed})"
            )
        leak = leakage_of(measure)
        if leak not in LEAKAGE_VALUES:
            raise ValueError(f"{at}: measure.notes needs 'leakage: <{'|'.join(LEAKAGE_VALUES)}>'")
        if leak != "none" and weakest([e["status"], SMOKE]) != e["status"]:
            raise ValueError(f"{at}: leakage {leak} allows at most smoke, not {e['status']}")
        out.append(StatusEntry(tuple(taxa), e["status"], str(e.get("source", "")), measure))
    return tuple(out)


def check_variant(cfg, call, variant, where="call status"):
    """``ValueError`` unless ``variant`` fits ``call``: one of the step 1 labels for a per-variant
    call, empty for any other call. Call this before any file name is built from the variant."""
    labels = [variant_label(v) for v in cfg.step1_variants]
    if cfg.call_by_name(call).get("per_variant"):
        if variant not in labels:
            raise ValueError(f"{where}: variant must be one of {labels} for call {call!r}")
    elif variant != "":
        raise ValueError(f"{where}: call {call!r} has no variants, but variant is {variant!r}")


def build_source(data, cfg, file_name, path=None):
    """A ``CallSource`` from parsed JSON, or ``ValueError``. ``file_name`` must be the canonical name."""
    where = str(path or file_name)
    if not isinstance(data, dict):
        raise ValueError(f"{where}: must be a JSON object")
    missing, extra = FILE_KEYS - set(data), set(data) - FILE_KEYS
    if missing or extra:
        raise ValueError(
            f"{where}: missing key(s) {sorted(missing)}, unknown key(s) {sorted(extra)}"
        )
    call, variant = data["call"], data["variant"]
    if not isinstance(call, str) or not isinstance(variant, str):
        raise ValueError(f"{where}: call and variant must be strings")
    ok, reason = call_eligible(cfg, call)
    if not ok:
        raise ValueError(
            f"{where}: call {call!r} is not eligible for a call status file ({reason})"
            if reason != "unknown call"
            else f"{where}: unknown call {call!r}"
        )
    check_variant(cfg, call, variant, where)
    if file_name != call_file_name(call, variant):
        raise ValueError(
            f"{where}: file name does not match call and variant ({call_file_name(call, variant)})"
        )
    reads = data["reads"]
    keys = {"name", "version", "params_hash", "artefact_hash"}
    if (
        not isinstance(reads, list)
        or not reads
        or any(
            not isinstance(r, dict)
            or set(r) != keys
            or any(not isinstance(v, str) for v in r.values())
            for r in reads
        )
    ):
        raise ValueError(f"{where}: reads must be a list of objects with {sorted(keys)} (strings)")
    for key in ("config_sha256", "call_hash"):
        if not isinstance(data[key], str) or not data[key]:
            raise ValueError(f"{where}: {key} must be a non-empty string")
    return CallSource(
        Path(path) if path else Path(file_name),
        call,
        variant,
        data["config_sha256"],
        data["call_hash"],
        tuple(
            ModuleIdentity(r["name"], r["version"], r["params_hash"], r["artefact_hash"])
            for r in reads
        ),
        _entries(data["entries"], where),
    )


def load_call_source(path, cfg):
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    return build_source(data, cfg, path.name, path)


def stale_reason(source, cfg, identities, states, config_sha256):
    """``""`` when the file belongs to this run, else why not (checked in a fixed order)."""
    if source.config_sha256 != config_sha256:
        return "config differs"
    now = reads_of_call(cfg, source.call, source.variant)
    if sorted(r.name for r in source.reads) != now:
        return "reads differ"
    for m in now:
        if m not in identities:
            return f"no module run record: {m}"
    for m in now:
        state = states.get(m, "missing")
        if state != "ok":
            return f"module state not ok: {m} ({state})"
    for r in source.reads:
        if identities[r.name] != r:
            return f"identity differs: {r.name}"
    if source.call_hash != call_hash(cfg, source.call, source.variant):
        return "call_hash differs"
    return ""


def _rate_text(rate):
    return (
        "not measured" if not rate else f"{rate['value']:.3f} [{rate['lo']:.3f}, {rate['hi']:.3f}]"
    )


class CallStatusResolver:
    """Loads every ``<workdir>/status/calls/*.json`` and answers ``(call, variant, taxon)``."""

    def __init__(self, workdir, cfg, identities, states, config_sha256, lineage):
        self.folder = Path(workdir) / "status" / "calls"
        self.lineage = lineage
        self._items = {}  # (call, variant) -> (source, stale reason)
        for path in sorted(self.folder.glob("*.json")):
            try:
                source = load_call_source(path, cfg)
            except (KeyError, TypeError, ValueError) as err:  # JSONDecodeError is a ValueError
                message = str(err)
                raise CallStatusError(
                    message if message.startswith(str(path)) else f"{path}: {message}"
                ) from err
            reason = stale_reason(source, cfg, identities, states, config_sha256)
            self._items[(source.call, source.variant)] = (source, reason)

    def __call__(self, call, variant, taxon):
        """None (no file, or the taxon is not tested), ``(status, basis)``, or ``(None, reason)`` (stale)."""
        item = self._items.get((call, variant))
        if item is None:
            return None
        source, reason = item
        if reason:
            return None, reason
        best = best_entry(source.entries, taxon, self.lineage)
        if best is None:
            return None
        entry, tested = best
        return entry.status, f"call:{call}:taxon:{tested}"

    def rows(self):
        """One row per entry and tested taxon of a valid file, and one row per stale file."""
        out = []
        for (call, variant), (source, reason) in sorted(self._items.items()):
            if reason:
                out.append(
                    {
                        "file": source.path.name,
                        "call": call,
                        "variant": variant,
                        "taxon": "",
                        "status": "",
                        "calibration_set": "",
                        "n_pos": "",
                        "n_neg": "",
                        "sensitivity": "",
                        "specificity": "",
                        "valid": False,
                        "reason": reason,
                    }
                )
                continue
            for entry in source.entries:
                m = entry.measure
                for taxon in entry.taxa:
                    out.append(
                        {
                            "file": source.path.name,
                            "call": call,
                            "variant": variant,
                            "taxon": taxon,
                            "status": entry.status,
                            "calibration_set": m.get("calibration_set", ""),
                            "n_pos": m.get("n_pos", ""),
                            "n_neg": m.get("n_neg", ""),
                            "sensitivity": _rate_text(m.get("sensitivity")),
                            "specificity": _rate_text(m.get("specificity")),
                            "valid": True,
                            "reason": "",
                        }
                    )
        return out
