"""Shared run helpers for scripts 01-04: provenance, the all_sources guard, atomic outputs.

Standard library only.
"""

import json
import os
import platform
import subprocess
from collections.abc import Callable
from pathlib import Path

import paths

# The script that writes each run log, named in the message of require_full.
PRODUCERS = {
    "extract_log.json": "01_extract_go_truth.py",
    "sequence_run.json": "02_attach_sequences.py",
}


class PartialInputError(ValueError):
    """A run log is missing, unreadable, or does not say all_sources: true."""


def git_commit() -> str:
    """HEAD commit of the repository, or 'unknown'. Does not depend on the time."""
    try:
        done = subprocess.run(
            ["git", "-C", str(paths.repo_root()), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return done.stdout.strip() or "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def python_version() -> str:
    return platform.python_version()


def _all_sources(log_path: Path):
    return json.loads(Path(log_path).read_text()).get("all_sources")


def says_all_sources(log_path: Path) -> bool:
    """True if the run log exists and says all_sources: true."""
    try:
        return _all_sources(log_path) is True
    except (OSError, ValueError, AttributeError):
        return False


def require_full(log_path: Path, what: str, allow_partial: bool) -> None:
    """Stop unless the run log says that `what` covers all sources.

    With allow_partial the log is not read. Raise PartialInputError (a ValueError)."""
    if allow_partial:
        return
    log_path = Path(log_path)
    try:
        all_sources = _all_sources(log_path)
    except (OSError, ValueError, AttributeError) as exc:
        raise PartialInputError(
            f"cannot read {log_path} ({exc}); use --allow-partial-truth-set to skip this check"
        ) from exc
    if all_sources is not True:
        producer = PRODUCERS.get(log_path.name, "the script that wrote it")
        raise PartialInputError(
            f"{log_path.name} does not say all_sources: true (value: {all_sources!r}), so "
            f"{what} may hold only some sources; rerun {producer} without --sources or use "
            "--allow-partial-truth-set"
        )


def write_json(path: Path, obj) -> None:
    """Write JSON with sorted keys, indent 2 and a final newline (byte-stable)."""
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def atomic_write_all(out_dir: Path, writers: dict[str, Callable[[Path], None]]) -> None:
    """Write each output to `.tmp.<name>` with its writer, then os.replace all of them.

    The files are moved in the order of `writers` only after every writer has succeeded. A
    failure leaves the earlier outputs untouched, and the temp files are always removed."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    temps = {name: out_dir / f".tmp.{name}" for name in writers}
    try:
        for name, writer in writers.items():
            writer(temps[name])
        for name in writers:
            os.replace(temps[name], out_dir / name)
    finally:
        for temp in temps.values():
            temp.unlink(missing_ok=True)
