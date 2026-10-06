"""Turn the Phase C ``metrics.json`` into status entries for the parameter-free rule R0 (one per species).

R0 ("SignalP calls a signal peptide") has no fitted parameter, so its sensitivity and false-positive
rate on a test set are measurements of the rule itself. Only the single-species test sets in
``SETS`` are used, so no tested taxon appears in two entries. Specificity is ``1 - FPR`` over all
negative classes of the set, with the interval ends swapped.

Cluster counts come from the Phase C cluster file (``count_clusters``). The interval of each rate is
the union of the interval in ``metrics.json`` and the Wilson interval on the number of clusters of
that class (the same widest-of rule as ``cluster_bootstrap``). A set whose protein counts differ from
``metrics.json`` is refused.

The file is validated before any entry is returned: a missing key or a non-finite or out-of-range
number is refused with a ``ValueError`` that names the path and the key.
"""

import csv
import gzip
import json
import math
from pathlib import Path

from cellsurface_sorting_hat.calibration.intervals import wilson, wilson_p
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
LEAKAGE_NOTE = "leakage: overlap between the Phase C positives and the SignalP 6 training data was not measured"
# the call and module that the R0 measurement belongs to (decision of 2026-10-06)
R0_CALL = "signal_peptide_protein"
R0_VARIANT = "R0"
R0_MODULE = f"step1_rule@{R0_VARIANT}"


def _open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, encoding="utf-8-sig", newline="")


def _read_tsv(path, columns):
    try:
        with _open_text(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            missing = [c for c in columns if c not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(f"{path}: missing column(s) {missing}")
            yield from reader
    except (OSError, EOFError) as exc:
        raise ValueError(f"{path}: cannot be read: {exc}") from exc


def count_clusters(clusters_path, eval_table_path):
    """Proteins and distinct clusters per Phase C source and class.

    Rows of ``eval_table`` with origin ``go``, ``homology_only`` ``no`` and class ``pos`` or ``neg``
    are joined to ``clusters`` on ``seq_sha256``. A row belongs to every source listed in
    ``source_ids`` (comma separated). Returns ``{source: {"pos": (n_proteins, n_clusters),
    "neg": (...)}}`` for the sources in ``SETS``. A protein with no cluster is refused."""
    cluster_of = {}
    for r in _read_tsv(clusters_path, ("seq_sha256", "cluster_id")):
        if cluster_of.setdefault(r["seq_sha256"], r["cluster_id"]) != r["cluster_id"]:
            raise ValueError(f"{clusters_path}: {r['seq_sha256']} is in two clusters")
    wanted = {s for group in SETS.values() for s in group}
    proteins = {(s, c): 0 for s in wanted for c in ("pos", "neg")}
    clusters = {key: set() for key in proteins}
    columns = ("seq_sha256", "origin", "class", "homology_only", "source_ids")
    for r in _read_tsv(eval_table_path, columns):
        if r["origin"] != "go" or r["homology_only"] != "no" or r["class"] not in ("pos", "neg"):
            continue
        for source in {x.strip() for x in r["source_ids"].split(",")} & wanted:
            cluster = cluster_of.get(r["seq_sha256"])
            if cluster is None:
                raise ValueError(
                    f"{eval_table_path}: {r['seq_sha256']} ({source}) is not in {clusters_path}"
                )
            proteins[(source, r["class"])] += 1
            clusters[(source, r["class"])].add(cluster)
    return {
        s: {c: (proteins[(s, c)], len(clusters[(s, c)])) for c in ("pos", "neg")} for s in wanted
    }


def check_signalp(record_path, module, mode):
    """Compare the SignalP module and mode of Phase C with the work directory's R0 record.

    The record holds the SignalP version (``tools.signalp``, for example ``6.0h-gpu``) and the mode
    (``params.mode``). It does not hold the name of the environment module (``signalp/6-gpu``), so
    the module name is compared through its major version and its ``gpu`` tag. Returns the version
    string of the record."""
    try:
        record = json.loads(Path(record_path).read_text())
        version = str(record["tools"]["signalp"])
        record_mode = str(record["params"]["mode"])
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise ValueError(
            f"{record_path}: no SignalP version or mode in the record: {exc!r}"
        ) from exc
    name, _, tag = module.partition("/")
    major, _, suffix = tag.partition("-")
    if name != "signalp" or not major:
        raise ValueError(f"--phasec-signalp-module {module!r}: expected a name like signalp/6-gpu")
    if not version.startswith(major) or (suffix == "gpu") != ("gpu" in version):
        raise ValueError(
            f"--phasec-signalp-module {module!r} does not match SignalP version {version!r} "
            f"in {record_path}"
        )
    if mode != record_mode:
        raise ValueError(
            f"--phasec-signalp-mode {mode!r} does not match mode {record_mode!r} in {record_path}"
        )
    return version


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


def entries_from_phasec(
    metrics_json,
    set_taxa,
    cluster_counts,
    extra_notes="",
    candidate="R0",
    variant="V-go",
    truth="direct",
    where_counts="the cluster files",
):
    """Return a list of entries (see ``make_entry``), one per test set in ``set_taxa``.

    ``set_taxa`` maps a Phase C test set key to its species-level taxon IDs. The status is the
    weaker of the Phase C label (``estimate`` or ``smoke test``) and ``status_from_measure``.
    ``cluster_counts`` is the result of ``count_clusters``; its protein counts must equal the counts
    in ``metrics.json``. The rates of the negative classes (N-int, N-sec,
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
        (p_pos, c_pos), (p_neg, c_neg) = (cluster_counts[SETS[key][0]][c] for c in ("pos", "neg"))
        if (p_pos, p_neg) != (n_pos, n_neg):
            raise ValueError(
                f"{where}: {path} has n_pos={n_pos}, n_neg={n_neg}; {where_counts} give "
                f"{p_pos} and {p_neg} proteins for {SETS[key][0]}"
            )

        sens = _cell(metrics, "all", "recall", variant, candidate, where)
        fpr = _cell(metrics, "all", "fpr", variant, candidate, where)
        measure = {
            "calibration_set": key,
            "truth_source": f"{path.name}:test_sets/{key}/truth/{truth}",
            "n_pos": n_pos,
            "n_neg": n_neg,
            "n_clusters_pos": c_pos,
            "n_clusters_neg": c_neg,
            "notes": (
                f"Phase C {label}; rule {candidate}, variant {variant}; GO direct evidence; "
                f"interval is the widest of the file's interval and the Wilson interval on the "
                f"cluster count; call={R0_CALL}; module={R0_MODULE}; variant={R0_VARIANT}; "
                f"{extra_notes + '; ' if extra_notes else ''}{LEAKAGE_NOTE}"
            ),
        }
        if sens:
            measure["sensitivity"] = widest_with_clusters(widen_at_boundary(sens, n_pos), c_pos)
        if fpr:
            spec = widen_at_boundary(_spec_from_fpr(fpr), n_neg)
            measure["specificity"] = widest_with_clusters(spec, c_neg)
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


def widest_with_clusters(rate, n_clusters):
    """Union of ``rate`` and the Wilson interval with ``p = value`` on ``n_clusters`` clusters."""
    if not rate or n_clusters < 1:
        return rate
    _, lo, hi = wilson_p(rate["value"], n_clusters)
    return {
        "value": rate["value"],
        "lo": min(rate["lo"], lo, rate["value"]),
        "hi": max(rate["hi"], hi, rate["value"]),
    }


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
