"""Locate the repository root and the step 1 work directory. Standard library only.

The repository root comes from the PROJ_ROOT environment variable. If PROJ_ROOT is not set,
the root is two levels above this file. The work directory comes from STEP1_WORKDIR. If
STEP1_WORKDIR is not set, it is `<workdir from config/site.yaml>/step1_compare`.
"""

import os
from pathlib import Path

STEP1_DIR = Path(__file__).resolve().parent


def repo_root() -> Path:
    env = os.environ.get("PROJ_ROOT")
    if env:
        return Path(env).resolve()
    return STEP1_DIR.parents[1]


def site_value(key: str, site_yaml: Path | None = None) -> str:
    """Read one top-level `key: value` line from config/site.yaml without PyYAML."""
    path = site_yaml or repo_root() / "config" / "site.yaml"
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(f"{key}:"):
                return line.split(":", 1)[1].strip()
    raise KeyError(f"{key} not found in {path}")


def workdir() -> Path:
    env = os.environ.get("STEP1_WORKDIR")
    if env:
        return Path(env)
    return Path(site_value("workdir")) / "step1_compare"


def downloads_dir() -> Path:
    return workdir() / "downloads"
