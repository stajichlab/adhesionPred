"""Build ``measure`` objects and write status sources (the files that record calibration)."""

import json
from numbers import Integral
from pathlib import Path

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.calibration.intervals import cluster_bootstrap
from cellsurface_sorting_hat.status import (  # noqa: F401  (re-exported: older code imports these here)
    _STRENGTH,
    ENTRY_KEYS,
    ESTIMATED,
    MAX_HALF_WIDTH,
    MIN_CLUSTERS,
    MIN_NEGATIVES,
    MIN_POSITIVES,
    SMOKE,
    UNVALIDATED,
    _check_measure,
    load_status_source,
    status_from_measure,
    weakest,
)


def build_measure(calibration_set, truth_source, y, call, clusters, notes="", n_boot=2000, seed=1):
    rates = cluster_bootstrap(y, call, clusters, n_boot=n_boot, seed=seed)
    y = list(y)
    measure = {
        "calibration_set": calibration_set,
        "truth_source": truth_source,
        "n_pos": sum(1 for v in y if v == 1),
        "n_neg": sum(1 for v in y if v == 0),
        "n_clusters_pos": rates["n_clusters_pos"],
        "n_clusters_neg": rates["n_clusters_neg"],
    }
    for key in ("sensitivity", "specificity"):
        if rates[key] is not None:
            measure[key] = rates[key]
    if notes:
        measure["notes"] = notes
    return measure


def _taxa(taxa):
    if isinstance(taxa, str | bytes):
        raise ValueError("taxa must be a list of integers")
    out = []
    for t in taxa:
        if isinstance(t, bool) or not isinstance(t, Integral) or t < 0:
            raise ValueError(f"taxa must be non-negative integers (got {t!r})")
        out.append(int(t))
    return out


def make_entry(taxa, measure, source="", cap=None):
    """A status entry. ``cap`` limits the status (for example ``cap="smoke"`` when the truth set
    overlaps the data that set the rule). The entry takes the weaker of the measured status and
    the cap; a cap never raises a status."""
    if cap is not None and cap not in _STRENGTH:
        raise ValueError(f"unknown cap: {cap!r}")
    taxa = _taxa(taxa)
    status = status_from_measure(measure)
    if cap is not None:
        status = weakest([status, cap])
    return {"taxa": taxa, "status": status, "source": source, "measure": measure}


def _check_entries(entries):
    seen = {}
    for i, e in enumerate(entries):
        where = f"entries[{i}]"
        if not isinstance(e, dict) or not {"taxa", "status"} <= set(e) <= ENTRY_KEYS:
            raise ValueError(f"{where}: needs taxa and status, and only keys {sorted(ENTRY_KEYS)}")
        if e["status"] not in _STRENGTH:
            raise ValueError(f"{where}: unknown status {e['status']!r}")
        taxa = _taxa(e["taxa"])
        measure = e.get("measure")
        if measure is None:
            if e["status"] == ESTIMATED:
                raise ValueError(f"{where}: status estimated needs a measure")
        else:
            _check_measure(measure, where)
            allowed = status_from_measure(measure)
            if weakest([e["status"], allowed]) != e["status"]:
                raise ValueError(
                    f"{where}: status {e['status']} is stronger than the measure allows ({allowed})"
                )
        label = (measure or {}).get("calibration_set", where)
        for t in taxa:
            if t in seen:
                raise ValueError(f"taxon {t} appears in two entries ({seen[t]} and {label})")
            seen[t] = label


def write_status_source(workdir, module, entries):
    """Write ``<workdir>/status/<module>.json`` using the identity in ``modules/<module>.json``.

    Entries are validated before anything is written (valid measures, a status no stronger than the
    measure allows, no tested taxon in two entries because the engine would pick the first of equal
    depth, no NaN). The file is also read back with ``load_status_source``.
    """
    record = json.loads((Path(workdir) / "modules" / f"{module}.json").read_text())
    entries = list(entries)
    _check_entries(entries)
    data = {k: record[k] for k in ("module", "version", "params_hash", "artefact_hash")}
    data["entries"] = entries
    payload = (json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    path = Path(workdir) / "status" / f"{module}.json"
    write_atomic(path, payload)
    load_status_source(path)  # raises if the file is not valid
    return path
