"""Shared helpers of the Phase C scripts: STOP error, run JSON, output hashes.

Every Phase C script writes its outputs with `write_outputs`: each output goes to a temp name,
the run JSON is written last and records the SHA-256 of every other output
(`outputs_sha256`), then all files are moved into place (runinfo.atomic_write_all). The next
script calls `require_current` and stops when a file differs from its recorded SHA-256.
"""

import importlib
import json
import platform
from pathlib import Path

import manifest
import runinfo

PHASEC_DIR = Path(__file__).resolve().parent
SEED = 20261001
VARIANTS = ("V-go", "V-kw")


class StopError(ValueError):
    """An input is missing, stale or inconsistent. The script prints STOP and exits 2."""


def out_dir(work) -> Path:
    """Phase C outputs go to $STEP1_WORKDIR/phasec/."""
    return Path(work) / "phasec"


def library_versions() -> dict:
    """Python and library versions for the run JSON (None when a library is missing)."""
    out = {"python": platform.python_version()}
    for name in ("numpy", "scipy", "sklearn"):
        try:
            out[name] = importlib.import_module(name).__version__
        except ImportError:
            out[name] = None
    return out


def read_json(path) -> dict:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError) as exc:
        raise StopError(f"cannot read {path}: {exc}") from exc


def write_outputs(out, writers: dict, run_name: str, log: dict) -> dict:
    """Write the outputs and then the run JSON, which records `outputs_sha256`.

    The run JSON writer reads the temp files `.tmp.<name>` that runinfo.atomic_write_all has
    written before it (writers run in dict order). No file gets its final name before all
    writers have succeeded."""
    out = Path(out)
    names = list(writers)

    def write_run(path):
        log["outputs_sha256"] = {n: manifest.sha256_file(out / f".tmp.{n}") for n in names}
        runinfo.write_json(path, log)

    runinfo.atomic_write_all(out, {**writers, run_name: write_run})
    return log


def require_current(out, run_name: str, names, producer: str) -> dict:
    """Stop unless every file in `names` has the SHA-256 that `run_name` recorded."""
    out = Path(out)
    log = read_json(out / run_name)
    recorded = log.get("outputs_sha256", {})
    for name in names:
        path = out / name
        if not path.exists():
            raise StopError(f"{path} is missing; run {producer}")
        if recorded.get(name) != manifest.sha256_file(path):
            raise StopError(
                f"{name} differs from the SHA-256 in {run_name} (stale or edited input); "
                f"re-run {producer}"
            )
    return log


def float_or_stop(text: str, what: str) -> float:
    """float(text) for a finite number; StopError for empty, NaN or inf."""
    try:
        value = float(text)
    except (TypeError, ValueError) as exc:
        raise StopError(f"{what}: {text!r} is not a number") from exc
    if value != value or value in (float("inf"), float("-inf")):
        raise StopError(f"{what}: {text!r} is not a finite number")
    return value
