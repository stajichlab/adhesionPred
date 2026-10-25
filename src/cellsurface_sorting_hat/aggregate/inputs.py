"""Read the outputs of several engine runs (one proteome per run) for aggregation.

A run directory holds ``calls.long.tsv.gz``, ``proteins.tsv.gz`` and ``run.json``. Calls are read for ONE step 1 variant only: the file has one row per
variant for a per_variant call, and counting all variants would multiply the counts. The default is the gate variant of the run (``run.json``
``default_gate``); ungated calls have no variant.
"""

import csv
import gzip
import json
from dataclasses import dataclass, field
from pathlib import Path

from cellsurface_sorting_hat.engine import variant_label

FILES = ("calls.long.tsv.gz", "proteins.tsv.gz", "run.json")
VALUES = ("called", "not_called", "not_assessable", "not_applicable")


class AggregateError(ValueError):
    """An input of the aggregation is missing or not usable."""


@dataclass
class RunData:
    label: str
    run_dir: str
    taxon: int
    n_proteins: int
    n_invalid: int
    config_sha256: str
    version: str
    variant: str
    absent_modules: tuple
    calls: dict = field(default_factory=dict)  # call -> list of (value, status)
    call_order: list = field(default_factory=list)


def read_manifest(path):
    """TSV with ``run_dir`` and an optional ``label``; the label defaults to the name of the run directory's parent when the directory is called ``out``."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "run_dir" not in reader.fieldnames:
            raise AggregateError(f"{path}: needs a 'run_dir' column")
        out = []
        for n, r in enumerate(reader, 2):
            d = (r["run_dir"] or "").strip()
            if not d:
                raise AggregateError(f"{path}:{n}: empty run_dir")
            out.append((d, (r.get("label") or "").strip() or default_label(d)))
    if not out:
        raise AggregateError(f"{path}: no runs")
    return out


def default_label(run_dir):
    p = Path(run_dir)
    return p.parent.name if p.name in ("out", "output") else p.name


def _read_gz_tsv(path):
    try:
        with gzip.open(path, "rt", newline="", encoding="utf-8-sig") as fh:
            yield from csv.DictReader(fh, delimiter="\t")
    except (OSError, EOFError, UnicodeDecodeError, gzip.BadGzipFile) as err:
        raise AggregateError(f"{path}: cannot be read ({err.__class__.__name__}: {err})") from err


def read_run(run_dir, label, variant=None):
    d = Path(run_dir)
    for name in FILES:
        if not (d / name).is_file():
            raise AggregateError(
                f"{d}: {name} not found; give an output directory of cellsurface_sorting_hat"
            )
    try:
        info = json.loads((d / "run.json").read_text(encoding="utf-8-sig"))
    except ValueError as err:
        raise AggregateError(f"{d / 'run.json'}: not valid JSON ({err})") from err
    for key in ("config_sha256", "default_gate", "n_proteins", "taxa"):
        if key not in info:
            raise AggregateError(f"{d / 'run.json'}: key {key!r} is missing")
    taxa = {int(t): n for t, n in info["taxa"].items()}
    if len(taxa) != 1:
        raise AggregateError(
            f"{d}: the run has {len(taxa)} taxa ({sorted(taxa)}); aggregate runs of one proteome each (run the engine per species)"
        )
    (taxon,) = taxa
    gate = variant if variant is not None else variant_label(info["default_gate"])
    run = RunData(
        label, str(d), taxon, int(info["n_proteins"]), int(info.get("n_invalid", 0)), info["config_sha256"], str(info.get("version", "")),
        gate, tuple(info.get("absent_modules", [])),
    )  # fmt: skip
    for r in _read_gz_tsv(d / "calls.long.tsv.gz"):
        if r["variant"] not in ("", gate):
            continue
        if r["value"] not in VALUES:
            raise AggregateError(f"{d}: unknown call value {r['value']!r} for {r['call']!r}")
        if r["call"] not in run.calls:
            run.calls[r["call"]] = []
            run.call_order.append(r["call"])
        run.calls[r["call"]].append((r["value"], r["status"]))
    if not run.calls:
        raise AggregateError(f"{d}: no call rows for variant {gate!r}")
    n_taxon = sum(1 for r in _read_gz_tsv(d / "proteins.tsv.gz"))
    if n_taxon != run.n_proteins:
        raise AggregateError(
            f"{d}: proteins.tsv.gz has {n_taxon} rows but run.json says {run.n_proteins}"
        )
    return run


def read_runs(entries, variant=None, allow_mixed_config=False):
    runs, seen = [], set()
    for run_dir, label in entries:
        if label in seen:
            raise AggregateError(f"duplicate proteome label {label!r}")
        seen.add(label)
        runs.append(read_run(run_dir, label, variant))
    configs = {r.config_sha256 for r in runs}
    if len(configs) > 1 and not allow_mixed_config:
        raise AggregateError(
            f"the runs were made with {len(configs)} different categories.yaml files; use --allow-mixed-config to aggregate them anyway "
            "(the report lists the mixture)"
        )
    return runs
