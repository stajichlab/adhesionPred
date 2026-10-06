"""Module result cache and atomic writes.

Results are stored per sequence sha256 in one table per module identity. The identity key covers
module name, version, params hash, artefact hash and tool versions, so a change in any of them
starts a new table.

Concurrency: ``ModuleCache`` takes an exclusive ``fcntl`` lock on ``<table>.lock`` for an update and
a shared lock for a read, so a reader never sees a data file with a sidecar that belongs to another
version. A table that fails its checksum (for example after a crash between the two renames) is
treated as a cache miss and is rebuilt by the next update. ``fcntl`` locks may not be reliable on
every network file system; use a node-local or otherwise tested directory for parallel runs.
"""

import csv
import fcntl
import gzip
import hashlib
import io
import json
import os
import uuid
from contextlib import contextmanager
from pathlib import Path


class CacheError(RuntimeError):
    """A cached file is missing its checksum or does not match it."""


def identity_key(identity, tool_versions=None):
    payload = {
        "name": identity.name,
        "version": identity.version,
        "params_hash": identity.params_hash,
        "artefact_hash": identity.artefact_hash,
        "tools": dict(sorted((tool_versions or {}).items())),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _tmp_name(path):
    return path.with_name(f"{path.name}.tmp{os.getpid()}-{uuid.uuid4().hex[:8]}")


def write_atomic(path, data):
    """Write bytes to ``path`` through a temporary file, then write ``path + '.sha256'``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = _tmp_name(path)
    tmp.write_bytes(data)
    os.replace(tmp, path)
    side = path.with_name(path.name + ".sha256")
    side_tmp = _tmp_name(side)
    side_tmp.write_text(hashlib.sha256(data).hexdigest() + "\n")
    os.replace(side_tmp, side)


def read_verified(path):
    path = Path(path)
    side = path.with_name(path.name + ".sha256")
    if not side.exists():
        raise CacheError(f"{path} has no checksum file")
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != side.read_text().strip():
        raise CacheError(f"{path} does not match its checksum")
    return data


@contextmanager
def _locked(lock_path, exclusive):
    """Hold a lock on ``lock_path``. A reader that cannot open the lock file (a read-only
    directory) reads without a lock. A steady stream of readers can delay a writer."""
    if exclusive:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fh = open(lock_path, "a")
    except OSError:
        if exclusive:
            raise
        yield
        return
    with fh:
        fcntl.flock(fh, fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


class ModuleCache:
    """A table of module rows keyed by sequence sha256, for one module identity."""

    def __init__(self, directory, key):
        self.path = Path(directory) / f"{key}.tsv.gz"
        self._lock = self.path.with_name(self.path.name + ".lock")

    def _load_unlocked(self):
        if not self.path.exists():
            return {}
        try:
            text = gzip.decompress(read_verified(self.path)).decode()
        except (CacheError, OSError, EOFError):
            return {}  # a damaged table is a cache miss
        return {row["sha256"]: row for row in csv.DictReader(io.StringIO(text), delimiter="\t")}

    def load(self):
        with _locked(self._lock, exclusive=False):
            return self._load_unlocked()

    def update(self, rows):
        """Merge ``rows`` (dicts with a ``sha256`` key) into the table and rewrite it atomically."""
        with _locked(self._lock, exclusive=True):
            table = self._load_unlocked()
            for row in rows:
                table[row["sha256"]] = row
            columns = ["sha256"] + sorted({k for r in table.values() for k in r} - {"sha256"})
            buf = io.StringIO()
            writer = csv.DictWriter(buf, columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for sha in sorted(table):
                writer.writerow(table[sha])
            write_atomic(self.path, gzip.compress(buf.getvalue().encode(), mtime=0))
