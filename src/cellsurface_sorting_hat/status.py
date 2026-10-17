"""Validation status of a module, as a function of (module, version, taxon).

A status applies to a protein's taxon only if that taxon is a tested taxon or a descendant of one.
A shared broad label (for example "Eurotiomycetes") is not enough.

Each entry of a status source can carry an optional ``measure`` object with the calibration set and
the sensitivity and specificity measured on it, so reports can show them.
"""

import json
import re
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
    "negative_basis",
    "leakage_exception",
    "sensitivity",
    "specificity",
    "notes",
}
CURATED, ASSUMED = "curated", "assumed"
NEGATIVE_BASES = (CURATED, ASSUMED)
LEAKAGE_VALUES = ("none", "partial", "tuned_on_truth", "in_reference", "unknown")
_LEAKAGE = re.compile(r"(?:^|[\s;])leakage:\s*([a-z_]+)")
MIN_TAXON = 2  # taxon 1 is the root: a status for it would apply to every protein


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


def leakage_of(measure):
    """The last ``leakage: <value>`` word in ``measure.notes``, or None."""
    found = _LEAKAGE.findall((measure or {}).get("notes", ""))
    return found[-1] if found else None


def _check_negative_basis(measure, where):
    """``negative_basis`` is ``{stratum: {"basis": curated|assumed, "n": count}}``.

    Only ``curated`` strata count for the status numbers. When a measure has a curated stratum,
    the headline ``n_neg`` is the sum over the curated strata (an ``assumed`` stratum such as PM-TM
    has its own row and is not in it). A measure with no curated stratum has headline ``n_neg``
    equal to the sum over all strata; it can reach ``smoke`` only."""
    basis = measure["negative_basis"]
    if not isinstance(basis, dict) or not basis:
        raise ValueError(f"{where}.negative_basis: must be an object of strata")
    for name, item in basis.items():
        at = f"{where}.negative_basis.{name}"
        if not isinstance(item, dict) or set(item) != {"basis", "n"}:
            raise ValueError(f"{at}: needs basis and n")
        if item["basis"] not in NEGATIVE_BASES:
            raise ValueError(f"{at}: basis must be curated or assumed (got {item['basis']!r})")
        if not (isinstance(item["n"], int) and not isinstance(item["n"], bool) and item["n"] >= 0):
            raise ValueError(f"{at}.n: must be a non-negative integer")
    curated = sum(v["n"] for v in basis.values() if v["basis"] == CURATED)
    total = sum(v["n"] for v in basis.values())
    expected = curated if any(v["basis"] == CURATED for v in basis.values()) else total
    if "n_neg" in measure and measure["n_neg"] != expected:
        raise ValueError(
            f"{where}.n_neg: is {measure['n_neg']} but negative_basis gives {expected} "
            "(the headline counts the curated strata only)"
        )


def _check_leakage_exception(measure, where):
    exc = measure["leakage_exception"]
    keys = {"reason", "signer", "date"}
    if (
        not isinstance(exc, dict)
        or set(exc) != keys
        or not all(isinstance(exc[k], str) and exc[k].strip() for k in keys)
    ):
        raise ValueError(
            f"{where}.leakage_exception: needs reason, signer and date (non-empty text)"
        )


def negative_rows(measure):
    """One row per negative stratum: ``stratum``, ``basis``, ``n`` and whether it ``counted`` for
    the status numbers. An assumed stratum (for example PM-TM) is listed apart and never counted."""
    basis = (measure or {}).get("negative_basis") or {}
    return [
        {"stratum": k, "basis": v["basis"], "n": v["n"], "counted": v["basis"] == CURATED}
        for k, v in basis.items()
    ]


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
    if "negative_basis" in measure:
        _check_negative_basis(measure, where)
    if "leakage_exception" in measure:
        _check_leakage_exception(measure, where)


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
    """``estimated`` needs at least 20 positives, at least 20 curated negatives and at least 20
    independent clusters of each class, and a 95% interval half-width of at most 0.10 for both
    sensitivity and specificity. Phase C uses the same floor for recall.

    Only negatives in a ``curated`` stratum of ``negative_basis`` count. A measure without
    ``negative_basis`` has no curated negatives, so it is at most ``smoke``. A cluster count that is
    missing does not meet the floor (an older rule skipped the floor for a missing count).
    A measure with any positive but without a specificity is at most ``smoke``: a rule that was
    never tested on negatives is not an estimate. No sensitivity or no positive gives
    ``unvalidated``.

    A measure with a non-finite or out-of-range rate or count is refused with ``ValueError``."""
    _check_measure(measure, "measure")
    sens, spec = measure.get("sensitivity"), measure.get("specificity")
    if not sens or measure.get("n_pos", 0) < 1:
        return UNVALIDATED
    narrow = (sens["hi"] - sens["lo"]) / 2 <= max_half_width + _TOLERANCE
    if spec:
        narrow = narrow and (spec["hi"] - spec["lo"]) / 2 <= max_half_width + _TOLERANCE
    curated = sum(r["n"] for r in negative_rows(measure) if r["counted"])
    enough = measure["n_pos"] >= min_positives and curated >= min_negatives
    for key in ("n_clusters_pos", "n_clusters_neg"):
        if measure.get(key, 0) < min_clusters:
            enough = False
    if spec and enough and narrow:
        return ESTIMATED
    return SMOKE


def check_entry(where, status, measure, taxa):
    """``ValueError`` unless ``status`` is allowed for ``measure`` and ``taxa``.

    A status is never higher than the status that ``status_from_measure`` derives. An entry without
    a measure is ``unvalidated`` only. Taxon 1 and 0 are refused. Leakage other than ``none`` (or no
    ``leakage:`` value) allows at most ``smoke``, unless the leakage is ``unknown`` (not measured)
    and the measure carries ``leakage_exception`` (the owner's recorded exception, for rule R0).
    The exception never lifts the cap of a measured leakage."""
    for t in taxa:
        if t < MIN_TAXON:
            raise ValueError(f"{where}: taxon {t} is below the minimum taxon {MIN_TAXON}")
    if measure is None:
        if status != UNVALIDATED:
            raise ValueError(f"{where}: status {status} needs a measure (without one: unvalidated)")
        return
    allowed = status_from_measure(measure)
    if weakest([status, allowed]) != status:
        hint = ""
        if "leakage_exception" not in measure and "negative_basis" not in measure:
            # an entry that an older `calibrate phasec` wrote, before the per-stratum rules
            hint = (
                f"; the entry in {where} carries no leakage_exception and no per-stratum "
                "negative_basis, so the file is out of date: re-run `calibrate phasec`"
            )
        raise ValueError(
            f"{where}: status {status} is stronger than the measure allows ({allowed}){hint}"
        )
    leak = leakage_of(measure)
    if leak not in LEAKAGE_VALUES:
        leak = "unknown"
    # The owner's exception covers the unmeasured overlap only (D13); a measured leakage
    # (partial, tuned_on_truth, in_reference) keeps its cap even when the key is present.
    excepted = leak == "unknown" and "leakage_exception" in measure
    if leak != "none" and not excepted and weakest([status, SMOKE]) != status:
        raise ValueError(f"{where}: leakage {leak} allows at most smoke, not {status}")


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
        check_entry(where, e["status"], measure, taxa)
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


def exception_text(entry):
    """`` (<reason>)`` when the entry carries the owner's leakage exception, else an empty string.
    The reason goes into the basis column so that a reader sees that the cap was lifted."""
    exc = (entry.measure or {}).get("leakage_exception")
    return f" ({exc['reason']})" if exc else ""


def resolve_status(record, running, taxon, lineage):
    """Return (status, basis)."""
    entry, tested, reason = resolve_entry(record, running, taxon, lineage)
    if entry is None:
        return UNVALIDATED, reason
    return entry.status, f"taxon:{tested}{exception_text(entry)}"
