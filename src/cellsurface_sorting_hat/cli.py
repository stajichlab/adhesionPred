"""Command line driver: read module tables from a workdir, evaluate the rules, write outputs.

Version 1 of the driver does not run modules. It reads their outputs from
``<workdir>/modules/<module>.tsv.gz`` (columns ``id``, ``state``, ``call`` and module fields) and
optional ``<module>.json`` (run record) and ``<workdir>/status/<module>.json`` (status source).
"""

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import zlib
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from cellsurface_sorting_hat import __version__
from cellsurface_sorting_hat.engine import (
    ConfigError,
    ModuleTable,
    collect_evidence,
    evaluate,
    load_config,
    referenced_modules,
)
from cellsurface_sorting_hat.fasta import NA_INVALID, FastaError, read_fasta
from cellsurface_sorting_hat.outputs import (
    RunInfo,
    render_report,
    write_atomic_text,
    write_evidence,
    write_long,
    write_proteins,
    write_run_json,
    write_wide,
)
from cellsurface_sorting_hat.status import (
    UNVALIDATED,
    ModuleIdentity,
    load_status_source,
    resolve_entry,
)
from cellsurface_sorting_hat.taxonomy import Lineage, TaxonError

UNUSABLE_RUN_STATES = {"unavailable", "not_run", "error"}
RUN_STATES = {"ok", "partial", "unavailable", "not_run", "error"}
CALL_VALUES = {"called", "not_called"}
FLAG_VALUES = {"0", "1"}


class RunError(RuntimeError):
    """The run cannot continue; the message is shown to the user."""


class InputError(RunError):
    """An input file is malformed; the message names the file."""


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="cellsurface_sorting_hat", description=__doc__)
    p.add_argument("--fasta", required=True)
    p.add_argument("--taxon", type=int, help="NCBI taxon ID for every protein not in --taxon-map")
    p.add_argument("--taxon-map", help="TSV: protein ID, taxon ID (overrides --taxon per protein)")
    p.add_argument("--taxdump", required=True, help="path to nodes.dmp from the NCBI taxonomy")
    p.add_argument("--workdir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--config", help="categories.yaml (default: the packaged file)")
    return p.parse_args(argv)


def read_taxon_map(path):
    out, first_line = {}, {}
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except UnicodeDecodeError as err:
        raise InputError(f"{path}: cannot be read as text") from err
    for n, line in enumerate(text.splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = [f.strip() for f in line.split("\t")]
        if len(fields) < 2 or not fields[0]:
            raise InputError(f"{path}:{n}: expected 'protein ID<TAB>taxon ID'")
        try:
            taxon = int(fields[1])
        except ValueError:
            raise InputError(f"{path}:{n}: taxon ID is not an integer: {fields[1]!r}") from None
        if fields[0] in out:
            raise InputError(
                f"{path}:{n}: ID {fields[0]!r} already on line {first_line[fields[0]]}"
            )
        out[fields[0]] = taxon
        first_line[fields[0]] = n
    return out


def assign_taxa(proteins, default_taxon, taxon_map):
    taxa = {}
    for p in proteins:
        taxon = taxon_map.get(p.id, default_taxon)
        if taxon is None:
            raise RunError(f"no taxon for protein {p.id}: give --taxon or list it in --taxon-map")
        taxa[p.id] = taxon
    return taxa


@dataclass
class LoadedModules:
    tables: dict = field(default_factory=dict)
    states: dict = field(default_factory=dict)
    identities: dict = field(default_factory=dict)
    notes: dict = field(default_factory=dict)
    state_counts: dict = field(default_factory=dict)
    unmatched: dict = field(default_factory=dict)


def _read_json(path):
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except (json.JSONDecodeError, UnicodeDecodeError) as err:
        raise InputError(f"{path}: not valid JSON ({err.__class__.__name__})") from err
    if not isinstance(data, dict):
        raise InputError(f"{path}: must be a JSON object")
    return data


def _read_module_table(path):
    rows = {}
    try:
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            columns = next(reader, [])
            for needed in ("id", "state"):
                if needed not in columns:
                    raise InputError(f"{path.name}: missing column {needed!r}")
            for n, values in enumerate(reader, 2):
                if len(values) != len(columns):
                    raise InputError(
                        f"{path.name}:{n}: {len(values)} fields, the header has {len(columns)}"
                    )
                row = dict(zip(columns, values, strict=True))
                if row["id"] in rows:
                    raise InputError(f"{path.name}: duplicate row for ID {row['id']!r}")
                rows[row["id"]] = row
    except (OSError, EOFError, zlib.error, UnicodeDecodeError) as err:
        raise InputError(f"{path.name}: cannot be read ({err.__class__.__name__})") from err
    return columns, rows


def required_columns(cfg):
    """Module name -> columns the rules read, from the config."""
    need = defaultdict(set)

    def add(module, column):
        names = cfg.step1_variants if "{step1}" in module else [module]
        for name in names:
            need[module.replace("{step1}", name) if "{step1}" in module else name].add(column)

    def walk(node):
        kind, arg = next(iter(node.items()))
        if kind in ("and", "or"):
            for child in arg:
                walk(child)
        elif kind == "not":
            walk(arg)
        elif kind == "call":
            add(arg, "call")
        elif kind == "flag":
            add(arg.split(".", 1)[0], arg.split(".", 1)[1])
        elif kind == "test":
            add(arg["module"], arg["field"])

    for call in cfg.calls:
        if "expr" in call:
            walk(call["expr"])
    return need


def load_modules(workdir, protein_ids, invalid_ids, required=None):
    """Read module tables. Missing IDs, invalid proteins and bad values become unknown rows."""
    loaded = LoadedModules()
    folder = Path(workdir) / "modules"
    if not folder.exists():
        return loaded
    names = sorted(
        {p.name[: -len(".tsv.gz")] for p in folder.glob("*.tsv.gz")}
        | {p.stem for p in folder.glob("*.json")}
    )
    wanted = set(protein_ids)
    for name in names:
        meta_path = folder / f"{name}.json"
        meta = _read_json(meta_path) if meta_path.exists() else {}
        table_path = folder / f"{name}.tsv.gz"
        loaded.identities[name] = ModuleIdentity(
            name,
            str(meta.get("version", "")),
            str(meta.get("params_hash", "")),
            str(meta.get("artefact_hash", "")),
        )
        state = meta.get("run_state", "ok" if table_path.exists() else "error")
        if state not in RUN_STATES:
            raise InputError(f"{meta_path}: run_state {state!r} is not one of {sorted(RUN_STATES)}")
        loaded.states[name] = state
        if not table_path.exists():
            loaded.notes[name] = "no result table"
            continue
        if state in UNUSABLE_RUN_STATES:
            continue
        columns, rows = _read_module_table(table_path)
        for column in sorted((required or {}).get(name, ())):
            if column not in columns:
                raise InputError(f"{table_path.name}: missing column {column!r}")
        matched = wanted & set(rows)
        if not matched:
            loaded.states[name] = "error"
            loaded.notes[name] = (
                f"no FASTA ID matches the {len(rows)} row(s); check the ID format of the module"
            )
            continue
        if len(rows) > len(matched):
            loaded.unmatched[name] = len(set(rows) - wanted)
        table, counts = ModuleTable(name), Counter()
        for pid in protein_ids:
            if pid in invalid_ids:
                row = {"state": NA_INVALID}
            elif pid not in rows:
                row = {"state": "error"}
            else:
                row = dict(rows[pid])
                if row["state"] == "ok" and _bad_value(columns, row):
                    row["state"] = "bad_value"
            table.rows[pid] = row
            if row["state"] != "ok":
                counts[row["state"]] += 1
        missing = sum(1 for pid in protein_ids if pid not in rows and pid not in invalid_ids)
        if missing:
            loaded.states[name] = "partial"
            loaded.notes[name] = f"{missing} protein(s) have no row"
        loaded.tables[name] = table
        loaded.state_counts[name] = dict(counts)
    return loaded


def _bad_value(columns, row):
    if "call" in columns and row.get("call") not in CALL_VALUES:
        return True
    if "hit" in columns and row.get("hit") not in FLAG_VALUES:
        return True
    return False


def check_identical_sequences(proteins, tables):
    """Per module, count groups of identical sequences whose module rows differ."""
    groups = defaultdict(list)
    for p in proteins:
        if p.state != NA_INVALID:
            groups[p.sha256].append(p.id)
    groups = [ids for ids in groups.values() if len(ids) > 1]
    out = {}
    for name, table in tables.items():
        bad = 0
        for ids in groups:
            seen = {
                tuple(sorted((k, v) for k, v in table.rows[i].items() if k != "id")) for i in ids
            }
            bad += len(seen) > 1
        if bad:
            out[name] = bad
    return out


class StatusResolver:
    """Resolve the status of a module for a taxon, and the matched calibration entry."""

    def __init__(self, workdir, identities, lineage):
        self.folder = Path(workdir) / "status"
        self.identities = identities
        self.lineage = lineage
        self._records = {}

    def _record(self, module):
        if module not in self._records:
            path = self.folder / f"{module}.json"
            try:
                self._records[module] = load_status_source(path) if path.exists() else None
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as err:
                raise InputError(f"{path}: not a valid status source ({err})") from err
        return self._records[module]

    def entry(self, module, taxon):
        """(entry, tested taxon, reason); reason is empty when an entry applies."""
        identity = self.identities.get(module)
        if identity is None:
            return None, None, "no module run record"
        return resolve_entry(self._record(module), identity, taxon, self.lineage)

    def __call__(self, module, taxon):
        entry, tested, reason = self.entry(module, taxon)
        if entry is None:
            return UNVALIDATED, reason
        return entry.status, f"taxon:{tested}"


def _rate_text(rate):
    if not rate:
        return "not measured"
    return f"{rate['value']:.3f} [{rate['lo']:.3f}, {rate['hi']:.3f}]"


def calibration_rows(resolver, modules, taxa):
    """One row per (module, taxon): status, calibration set, counts, sensitivity, specificity."""
    rows = []
    for module in sorted(modules):
        for taxon in sorted(set(taxa)):
            entry, tested, reason = resolver.entry(module, taxon)
            measure = (entry.measure if entry else None) or {}
            rows.append(
                {
                    "module": module,
                    "taxon": taxon,
                    "status": entry.status if entry else UNVALIDATED,
                    "matched_taxon": tested if entry else "",
                    "reason": reason,
                    "calibration_set": measure.get("calibration_set", ""),
                    "n_pos": measure.get("n_pos", ""),
                    "n_neg": measure.get("n_neg", ""),
                    "sensitivity": _rate_text(measure.get("sensitivity")),
                    "specificity": _rate_text(measure.get("specificity")),
                }
            )
    return rows


def run(args):
    cfg = load_config(args.config)
    lineage = Lineage.from_nodes_dmp(args.taxdump)
    try:
        proteins = read_fasta(args.fasta)
    except FastaError as err:
        message = str(err)
        if not message.startswith(str(args.fasta)):
            message = f"{args.fasta}: {message}"
        raise InputError(message) from err
    taxon_map = read_taxon_map(args.taxon_map) if args.taxon_map else {}
    if args.taxon is None and not taxon_map:
        raise RunError("give --taxon or --taxon-map")
    taxa = assign_taxa(proteins, args.taxon, taxon_map)
    for taxon in set(taxa.values()):
        lineage.ancestors(taxon)  # raises TaxonError for an unknown taxon
    ids = [p.id for p in proteins]
    invalid = {p.id for p in proteins if p.state == NA_INVALID}
    loaded = load_modules(args.workdir, ids, invalid, required_columns(cfg))
    tables = loaded.tables
    resolver = StatusResolver(args.workdir, loaded.identities, lineage)
    records = evaluate(cfg, ids, taxa, tables, resolver)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    counts = Counter(taxa.values())
    info = RunInfo(
        version=__version__,
        taxa=dict(counts),
        taxonomy_sha256=hashlib.sha256(Path(args.taxdump).read_bytes()).hexdigest(),
        config_sha256=cfg.sha256,
        default_gate=cfg.default_gate,
        thresholds=cfg.thresholds,
        n_proteins=len(proteins),
        n_invalid=len(invalid),
        invalid=[[p.id, p.note] for p in proteins if p.state == NA_INVALID],
        n_trailing_stop=sum(p.trailing_stop for p in proteins),
        module_states=loaded.states,
        module_notes=loaded.notes,
        state_counts=loaded.state_counts,
        absent_modules=[
            m
            for m in referenced_modules(cfg)
            if m not in tables and m not in cfg.step1_variants and m not in loaded.states
        ],
        unavailable_variants=[v for v in cfg.step1_variants if v not in tables],
        inconsistent=check_identical_sequences(proteins, tables),
        unmatched_module_ids=loaded.unmatched,
        map_ids_not_in_fasta=len(set(taxon_map) - set(ids)),
        calibration=calibration_rows(resolver, set(tables) | {cfg.default_gate}, taxa.values()),
    )
    report = render_report(info, records)  # render first: a failure leaves no partial output
    write_long(out / "calls.long.tsv.gz", records)
    write_wide(out / "calls.wide.tsv.gz", records, ids)
    write_evidence(out / "evidence.tsv.gz", collect_evidence(cfg, ids, tables))
    write_proteins(out / "proteins.tsv.gz", proteins, taxa)
    write_atomic_text(out / "report.md", report)
    write_run_json(out / "run.json", info)
    return records


def main(argv=None):
    try:
        run(parse_args(argv))
    except (
        RunError,
        FastaError,
        TaxonError,
        ConfigError,
        OSError,
        subprocess.CalledProcessError,
    ) as err:
        print(f"cellsurface_sorting_hat: error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
