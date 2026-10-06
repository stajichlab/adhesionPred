"""Shared helpers: hashes and the writer for a module table plus its run record."""

import csv
import gzip
import hashlib
import io
import json
from dataclasses import dataclass, field
from pathlib import Path

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.fasta import NA_INVALID


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def params_hash(params):
    """Hash of a parameter dict (canonical JSON)."""
    return hashlib.sha256(json.dumps(params, sort_keys=True, default=str).encode()).hexdigest()


def artefact_hash(paths):
    """Hash of the databases or models a module uses: sorted ``name:sha256`` lines."""
    lines = sorted(f"{Path(p).name}:{sha256_file(p)}" for p in paths)
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


@dataclass(frozen=True)
class ModuleSpec:
    name: str
    version: str
    params: dict = field(default_factory=dict)
    artefacts: tuple = ()
    tools: dict = field(default_factory=dict)
    artefact_digest: str = (
        ""  # use instead of hashing large files (for example a recorded Pfam sha256)
    )


def invalid_row(protein):
    return {"id": protein.id, "state": NA_INVALID}


def write_module(workdir, spec, columns, rows, run_state="ok", note=""):
    """Write ``<workdir>/modules/<name>.tsv.gz`` and ``<name>.json`` (see the module contract).

    ``columns`` are the fields after ``id`` and ``state``. A row that lacks a column gets an empty
    value. Every row needs ``id`` and ``state``.
    """
    folder = Path(workdir) / "modules"
    header = ["id", "state"] + list(columns)
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    seen = set()
    for r in rows:
        if r["id"] in seen:
            raise ValueError(f"{spec.name}: duplicate row for ID {r['id']!r}")
        seen.add(r["id"])
        writer.writerow([r.get(h, "") for h in header])
    write_atomic(folder / f"{spec.name}.tsv.gz", gzip.compress(buf.getvalue().encode(), mtime=0))
    record = {
        "module": spec.name,
        "version": spec.version,
        "params_hash": params_hash(spec.params),
        "artefact_hash": spec.artefact_digest or artefact_hash(spec.artefacts),
        "run_state": run_state,
        "params": spec.params,
        "tools": spec.tools,
        "n_rows": len(seen),
        "note": note,
    }
    write_atomic(
        folder / f"{spec.name}.json", (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
    )
    return record
