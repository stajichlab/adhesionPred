"""Write ``status/calls/<call>[.<variant>].json`` (see ``cellsurface_sorting_hat.call_status``)."""

import fcntl
import json
import sys
from contextlib import contextmanager
from pathlib import Path

from cellsurface_sorting_hat.cache import write_atomic
from cellsurface_sorting_hat.call_status import build_source, call_file_name, load_call_source
from cellsurface_sorting_hat.engine import call_eligible, call_hash, reads_of_call


@contextmanager
def _locked(path):
    """Exclusive ``fcntl`` lock on ``<file>.lock``. It is best effort on a network file system
    (see cache.py); one writer per work directory is the supported use."""
    lock = path.with_name(path.name + ".lock")
    lock.parent.mkdir(parents=True, exist_ok=True)
    with open(lock, "a") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def _reads(workdir, cfg, call, variant):
    """Identities of the modules the call reads, from the module run records (all values strings)."""
    reads = []
    for module in reads_of_call(cfg, call, variant):
        path = Path(workdir) / "modules" / f"{module}.json"
        if not path.is_file():
            raise ValueError(
                f"{path}: module run record not found; run the module {module!r} first"
            )
        record = json.loads(path.read_text(encoding="utf-8-sig"))
        reads.append(
            {
                "name": module,
                "version": str(record.get("version", "")),
                "params_hash": str(record.get("params_hash", "")),
                "artefact_hash": str(record.get("artefact_hash", "")),
            }
        )
    return reads


def write_call_status(workdir, cfg, call, variant, config_sha256, entries):
    """Add ``entries`` to the call status file, validate, and write it atomically.

    Entries of other calibration sets stay and entries of the same set are replaced. All old entries
    are dropped (with a message on stderr) when the modules read, their identities, the call
    definition or the config differ from the old file. The result is checked before anything is
    written, so an invalid update leaves the existing file as it was.
    """
    ok, reason = call_eligible(cfg, call)
    if not ok:
        raise ValueError(
            f"unknown call {call!r}"
            if reason == "unknown call"
            else f"call {call!r} is not eligible for a call status file ({reason})"
        )
    name = call_file_name(call, variant)
    path = Path(workdir) / "status" / "calls" / name
    entries = list(entries)
    with _locked(path):
        data = {
            "call": call,
            "variant": variant,
            "config_sha256": config_sha256,
            "call_hash": call_hash(cfg, call, variant),
            "reads": _reads(workdir, cfg, call, variant),
        }
        kept = []
        if path.exists():
            try:
                old = json.loads(path.read_text(encoding="utf-8-sig"))
                old_entries = old["entries"]
                for e in old_entries:
                    e["measure"]["calibration_set"]
            except (ValueError, KeyError, TypeError) as err:
                raise ValueError(
                    f"{path}: cannot read the existing call status file: {err!r}"
                ) from err
            same = all(old.get(k) == data[k] for k in ("reads", "call_hash", "config_sha256"))
            if same:
                replaced = {e["measure"]["calibration_set"] for e in entries}
                kept = [e for e in old_entries if e["measure"]["calibration_set"] not in replaced]
            else:
                print(
                    f"{path}: module, call or config changed; {len(old_entries)} old entr(ies) dropped",
                    file=sys.stderr,
                )
        data["entries"] = kept + entries
        build_source(data, cfg, name, path)  # raises before anything is written
        payload = (json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
        write_atomic(path, payload)
        load_call_source(path, cfg)  # read back (still under the lock)
    return path
