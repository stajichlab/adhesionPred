"""Validation status of a module, as a function of (module, version, taxon).

A status applies to a protein's taxon only if that taxon is a tested taxon or a descendant of one.
A shared broad label (for example "Eurotiomycetes") is not enough.

Each entry of a status source can carry an optional ``measure`` object with the calibration set and
the sensitivity and specificity measured on it, so reports can show them.
"""

import json
from dataclasses import dataclass
from pathlib import Path

ESTIMATED, SMOKE, UNVALIDATED = "estimated", "smoke", "unvalidated"
_STRENGTH = {ESTIMATED: 2, SMOKE: 1, UNVALIDATED: 0}
ENTRY_KEYS = {"taxa", "status", "source", "measure"}
MEASURE_KEYS = {
    "calibration_set",
    "truth_source",
    "n_pos",
    "n_neg",
    "n_clusters_pos",
    "n_clusters_neg",
    "strata",
    "sensitivity",
    "specificity",
    "notes",
}


def weakest(statuses):
    """The weakest of the given statuses; ``unvalidated`` if there are none."""
    statuses = list(statuses)
    for s in statuses:
        if s not in _STRENGTH:
            raise ValueError(f"unknown status: {s!r}")
    if not statuses:
        return UNVALIDATED
    return min(statuses, key=_STRENGTH.__getitem__)


@dataclass(frozen=True)
class ModuleIdentity:
    name: str
    version: str
    params_hash: str
    artefact_hash: str


@dataclass(frozen=True)
class StatusEntry:
    taxa: tuple
    status: str
    source: str = ""
    measure: dict | None = None


@dataclass(frozen=True)
class StatusRecord:
    identity: ModuleIdentity
    entries: tuple  # of StatusEntry


def _check_rate(rate, where):
    if (
        not isinstance(rate, dict)
        or set(rate) != {"value", "lo", "hi"}
        or not all(isinstance(rate[k], int | float) and not isinstance(rate[k], bool) for k in rate)
        or not 0 <= rate["lo"] <= rate["value"] <= rate["hi"] <= 1
    ):
        raise ValueError(f"{where}: needs value, lo, hi with 0 <= lo <= value <= hi <= 1")


def _check_measure(measure, where):
    if not isinstance(measure, dict):
        raise ValueError(f"{where}: measure must be an object")
    unknown = set(measure) - MEASURE_KEYS
    if unknown:
        raise ValueError(f"{where}: unknown measure key(s): {sorted(unknown)}")
    for key in ("sensitivity", "specificity"):
        if key in measure:
            _check_rate(measure[key], f"{where}.{key}")
    strata = measure.get("strata", {})
    if not isinstance(strata, dict):
        raise ValueError(f"{where}.strata: must be an object of rates")
    for name, rate in strata.items():
        _check_rate(rate, f"{where}.strata.{name}")
    for key in ("n_pos", "n_neg", "n_clusters_pos", "n_clusters_neg"):
        if key in measure and not (
            isinstance(measure[key], int)
            and not isinstance(measure[key], bool)
            and measure[key] >= 0
        ):
            raise ValueError(f"{where}.{key}: must be a non-negative integer")


MIN_POSITIVES = 20
MIN_NEGATIVES = 20
MIN_CLUSTERS = 20
MAX_HALF_WIDTH = 0.10
# Floating-point subtraction can give 0.10000000000000003 for an interval of exactly 0.10.
_TOLERANCE = 1e-9


def status_from_measure(
    measure,
    min_positives=MIN_POSITIVES,
    min_negatives=MIN_NEGATIVES,
    max_half_width=MAX_HALF_WIDTH,
    min_clusters=MIN_CLUSTERS,
):
    """``estimated`` needs at least 20 positives and 20 negatives (and, when recorded, at least 20
    independent clusters of each) and a 95% interval half-width of at most 0.10 for both
    sensitivity and specificity. Phase C uses the same floor for recall. A measure with any
    positive but without a specificity is at most ``smoke``: a rule that was never tested on
    negatives is not an estimate. No sensitivity or no positive gives ``unvalidated``.

    A measure with a non-finite or out-of-range rate or count is refused with ``ValueError``."""
    _check_measure(measure, "measure")
    sens, spec = measure.get("sensitivity"), measure.get("specificity")
    if not sens or measure.get("n_pos", 0) < 1:
        return UNVALIDATED
    narrow = (sens["hi"] - sens["lo"]) / 2 <= max_half_width + _TOLERANCE
    if spec:
        narrow = narrow and (spec["hi"] - spec["lo"]) / 2 <= max_half_width + _TOLERANCE
    enough = measure["n_pos"] >= min_positives and measure.get("n_neg", 0) >= min_negatives
    # independent clusters, when the measure records them (older files may not)
    for key in ("n_clusters_pos", "n_clusters_neg"):
        if key in measure and measure[key] < min_clusters:
            enough = False
    if spec and enough and narrow:
        return ESTIMATED
    return SMOKE


def load_status_source(path):
    data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    identity = ModuleIdentity(
        data["module"], data["version"], data["params_hash"], data["artefact_hash"]
    )
    for name in ("module", "version", "params_hash", "artefact_hash"):
        value = getattr(identity, name if name != "module" else "name")
        if not isinstance(value, str) or not value:
            raise ValueError(f"{path}: {name} must be a non-empty string")
    entries, listed = [], {}
    for i, e in enumerate(data["entries"]):
        where = f"{path} entries[{i}]"
        unknown = set(e) - ENTRY_KEYS
        if unknown:
            raise ValueError(f"{where}: unknown key(s): {sorted(unknown)}")
        if e["status"] not in _STRENGTH:
            raise ValueError(f"unknown status {e['status']!r} in {path}")
        measure = e.get("measure")
        if measure is not None:
            _check_measure(measure, where)
        if isinstance(e["taxa"], str | bytes) or any(isinstance(t, bool) for t in e["taxa"]):
            raise ValueError(f"{where}: taxa must be a list of integers")
        taxa = tuple(int(t) for t in e["taxa"])
        for t in taxa:
            if t in listed:
                raise ValueError(f"{where}: taxon {t} is already in entries[{listed[t]}]")
            listed[t] = i
        entries.append(StatusEntry(taxa, e["status"], str(e.get("source", "")), measure))
    return StatusRecord(identity, tuple(entries))


def best_entry(entries, taxon, lineage):
    """``(entry, tested_taxon)`` for the most specific tested taxon that ``taxon`` is, or descends
    from, or None. Used for module status files and for call status files."""
    best = None  # (depth, tested taxon, entry)
    for entry in entries:
        for tested in entry.taxa:
            if lineage.is_descendant_or_self(taxon, tested):
                depth = lineage.depth(tested)
                if best is None or depth > best[0]:
                    best = (depth, tested, entry)
    return None if best is None else (best[2], best[1])


def resolve_entry(record, running, taxon, lineage):
    """Return (entry, tested_taxon, reason). ``entry`` is None when nothing applies.

    A stale record is refused: name, version, params hash and artefact hash must match.
    """
    if record is None:
        return None, None, "no status_source"
    for field in ("name", "version", "params_hash", "artefact_hash"):
        if getattr(record.identity, field) != getattr(running, field):
            return None, None, f"status_source stale: {field} differs"
    best = best_entry(record.entries, taxon, lineage)
    if best is None:
        return None, None, "taxon not tested"
    return best[0], best[1], ""


def resolve_status(record, running, taxon, lineage):
    """Return (status, basis)."""
    entry, tested, reason = resolve_entry(record, running, taxon, lineage)
    if entry is None:
        return UNVALIDATED, reason
    return entry.status, f"taxon:{tested}"
