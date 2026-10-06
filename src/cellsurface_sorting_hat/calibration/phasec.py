"""Turn the Phase C ``metrics.json`` into status entries for the parameter-free rule R0 (one per species).

R0 ("SignalP calls a signal peptide") has no fitted parameter, so its sensitivity and false-positive
rate on a test set are measurements of the rule itself. Only the single-species test sets in
``SETS`` are used, so no tested taxon appears in two entries. Specificity is ``1 - FPR`` over all
negative classes of the set, with the interval ends swapped.

The file is validated before any entry is returned: a missing key or a non-finite or out-of-range
number is refused with a ``ValueError`` that names the path and the key.
"""

import json
import math
from pathlib import Path

from cellsurface_sorting_hat.calibration.intervals import wilson
from cellsurface_sorting_hat.calibration.measure import make_entry
from cellsurface_sorting_hat.taxonomy import TaxonError

# Phase C species test set -> the source ID (row of species.tsv). One entry per species: a pooled
# estimate describes neither species (S. cerevisiae and C. albicans, A. fumigatus and A. nidulans
# differ strongly).
SETS = {
    "S1:Scer_SGD": ("Scer_SGD",),
    "S1:Calb_CGD": ("Calb_CGD",),
    "S3-Eurotiomycetes:Afum_ASPFU": ("Afum_ASPFU",),
    "S3-Eurotiomycetes:Anid_EMENI": ("Anid_EMENI",),
    "S3-Basidiomycota:Cneo_H99_GOA": ("Cneo_H99_GOA",),
    "S3-Basidiomycota:Umay_MYCMD": ("Umay_MYCMD",),
}
STRATA = ("N-int", "N-sec", "PM-TM")
LABELS = ("estimate", "smoke test")


def read_names(names_dmp):
    """``names.dmp`` -> ``{scientific name: [taxon IDs]}``."""
    out = {}
    for line in Path(names_dmp).read_text(encoding="utf-8-sig").splitlines():
        f = [x.strip() for x in line.split("|")]
        if len(f) >= 4 and f[3] == "scientific name":
            out.setdefault(f[1], []).append(int(f[0]))
    return out


def species_taxid(names, scientific_name):
    """The one taxon ID for a scientific name; refuses a missing or ambiguous name."""
    ids = names.get(scientific_name, [])
    if len(ids) != 1:
        raise TaxonError(f"{scientific_name!r} has {len(ids)} taxon IDs in names.dmp; expected 1")
    return ids[0]


def _get(obj, key, where):
    if not isinstance(obj, dict) or key not in obj:
        raise ValueError(f"{where}: missing key {key!r}")
    return obj[key]


def _number(x, where):
    if isinstance(x, bool) or not isinstance(x, int | float) or not math.isfinite(x):
        raise ValueError(f"{where}: must be a finite number (got {x!r})")
    if not 0 <= x <= 1:
        raise ValueError(f"{where}: must be between 0 and 1 (got {x!r})")
    return float(x)


def _count(x, where):
    if isinstance(x, bool) or not isinstance(x, int) or x < 0:
        raise ValueError(f"{where}: must be a non-negative integer (got {x!r})")
    return x


def rate_from_phasec(cell, where):
    """A Phase C rate cell -> ``{value, lo, hi}``, or None when the cell has no value (the class
    has no such proteins). A cell with a value but no interval is refused."""
    if cell is None or _get(cell, "value", where) is None:
        return None
    value = _number(cell["value"], f"{where}.value")
    lo = _number(_get(cell, "lo", where), f"{where}.lo")
    hi = _number(_get(cell, "hi", where), f"{where}.hi")
    if not lo <= value <= hi:
        raise ValueError(f"{where}: needs lo <= value <= hi (got {lo}, {value}, {hi})")
    return {"value": value, "lo": lo, "hi": hi}


def entries_from_phasec(metrics_json, set_taxa, candidate="R0", variant="V-go", truth="direct"):
    """Return a list of entries (see ``make_entry``), one per test set in ``set_taxa``.

    ``set_taxa`` maps a Phase C test set key to its species-level taxon IDs. The status is the
    weaker of the Phase C label (``estimate`` or ``smoke test``) and ``status_from_measure``. A
    Phase C measure records no cluster counts, so the cluster floor is not applied and is not
    claimed; the entry says so in its notes. The rates of the negative classes (N-int, N-sec,
    PM-TM) are kept in ``strata``, because the pooled specificity depends on the mix of negatives.
    """
    path = Path(metrics_json)
    try:
        m = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f"{path}: cannot read as JSON: {exc}") from exc
    sets = _get(m, "test_sets", str(path))
    entries = []
    for key, taxa in set_taxa.items():
        where = f"{path} test_sets[{key!r}]"
        ts = _get(sets, key, f"{path} test_sets")
        label = _get(ts, "label", where)
        if label not in LABELS:
            raise ValueError(f"{where}.label: must be one of {LABELS} (got {label!r})")
        t = _get(_get(ts, "truth", where), truth, f"{where}.truth")
        counts = _get(_get(t, "n", f"{where}.truth.{truth}"), "all", f"{where}.n")
        n_pos = _count(_get(counts, "pos", f"{where}.n.all"), f"{where}.n.all.pos")
        n_neg = _count(_get(counts, "neg", f"{where}.n.all"), f"{where}.n.all.neg")
        metrics = _get(t, "metrics", f"{where}.truth.{truth}")

        sens = _cell(metrics, "all", "recall", variant, candidate, where)
        fpr = _cell(metrics, "all", "fpr", variant, candidate, where)
        measure = {
            "calibration_set": key,
            "truth_source": f"{path.name}:test_sets/{key}/truth/{truth}",
            "n_pos": n_pos,
            "n_neg": n_neg,
            "notes": (
                f"Phase C {label}; rule {candidate}, variant {variant}; GO direct evidence; "
                "cluster counts not recorded, so the cluster floor was not checked"
            ),
        }
        if sens:
            measure["sensitivity"] = widen_at_boundary(sens, n_pos)
        if fpr:
            measure["specificity"] = widen_at_boundary(_spec_from_fpr(fpr), n_neg)
        strata = {}
        for name in STRATA:
            f = _cell(metrics, name, "fpr", variant, candidate, where)
            if f:
                # the per-stratum count is optional; without it the interval is not widened
                n_raw = t["n"].get(name, {}).get("neg", 0)
                n_stratum = _count(n_raw, f"{where}.n.{name}.neg")
                strata[name] = widen_at_boundary(_spec_from_fpr(f), n_stratum)
        if strata:
            measure["strata"] = strata
        entry = make_entry(taxa, measure, source=measure["truth_source"])
        if label != "estimate" and entry["status"] == "estimated":
            entry["status"] = "smoke"
        entries.append(entry)
    return entries


def _spec_from_fpr(fpr):
    return {"value": 1 - fpr["value"], "lo": 1 - fpr["hi"], "hi": 1 - fpr["lo"]}


def widen_at_boundary(rate, n):
    """A percentile bootstrap gives a zero-width interval at a rate of 0 or 1 (every resample
    agrees). Combine it with the Wilson interval on ``n`` counted proteins, which is wider."""
    if not rate or n <= 0 or rate["value"] not in (0.0, 1.0):
        return rate
    k = int(round(rate["value"] * n))
    _, lo, hi = wilson(k, n)
    return {"value": rate["value"], "lo": min(rate["lo"], lo), "hi": max(rate["hi"], hi)}


def _cell(metrics, stratum, rate, variant, candidate, where):
    """One rate of one stratum; None when the stratum has no such cell (not an error, except for
    ``all``, which every set must have)."""
    c = metrics.get(stratum, {}).get(variant, {}).get(candidate)
    if c is None:
        if stratum == "all":
            raise ValueError(f"{where}: no metrics for all/{variant}/{candidate}")
        return None
    return rate_from_phasec(_get(c, rate, f"{where}.{stratum}"), f"{where}.{stratum}.{rate}")
