"""Centralized site paths for adhesionPred analysis scripts.

Reads config/site.yaml so no script hardcodes /bigdata/stajichlab/... paths.
To point at a different site, edit config/site.yaml only.
"""

from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SITE_YAML = REPO_ROOT / "config" / "site.yaml"

with open(SITE_YAML) as _fh:
    _RAW = yaml.safe_load(_fh)

PATHS = {k: Path(v) for k, v in _RAW.items()}

FUNGI5K_INPUT = PATHS["fungi5k_input"]
FUNGI5K_DUCKDB = PATHS["fungi5k_duckdb"]
FUNGI5K_SAMPLES = PATHS["fungi5k_samples"]
COCCI_PANGENOME = PATHS["cocci_pangenome"]
COCCI_LONGREAD = PATHS["cocci_longread"]
SPHERULE_RNASEQ = PATHS["spherule_rnaseq"]
RHODOTORULA_BIOFILM = PATHS["rhodotorula_biofilm"]
WORKDIR = PATHS["workdir"]


def get(key: str) -> Path:
    return PATHS[key]
