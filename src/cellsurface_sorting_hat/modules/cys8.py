"""Eight-cysteine spacing match plus the step 1 (R0) signal peptide call -> module ``cys8_pattern``.

The spacing comes from ``data/sorting_hat/cys8_spacing.yaml`` (published sets, kept separate). A protein
matches when its sequence fits any one set as stated (decision D9 of the hydrophobin spec). A set needs the
five core gaps C1_C2, C3_C4, C4_C5, C5_C6 and C7_C8. C2_C3 and C6_C7 are 0 when a set does not list them
(the cysteine doublets). A set without a core gap is skipped and listed.

``hit`` is 1 when the sequence matches and the R0 call is ``called``. It is 0 when the sequence does not
match, or matches and R0 is ``not_called``. It is empty (not assessable) when the sequence matches and the R0
row is missing or not usable. ``pattern_match`` records the sequence match alone.
"""

import re
from dataclasses import dataclass

import yaml

from cellsurface_sorting_hat.modules.base import invalid_row, sha256_file

COLUMNS = ["hit", "pattern_match", "spacing_class", "sets", "n_cys", "length"]
CORE = ("C1_C2", "C3_C4", "C4_C5", "C5_C6", "C7_C8")
ORDER = ("C1_C2", "C2_C3", "C3_C4", "C4_C5", "C5_C6", "C6_C7", "C7_C8")
CLASSES = {"class_I": "I", "class_II": "II", "other": "other"}


@dataclass(frozen=True)
class SpacingSet:
    cls: str
    name: str
    gaps: tuple  # ((min, max), ...) for ORDER
    regex: object


@dataclass(frozen=True)
class SpacingSets:
    sets: tuple
    skipped: tuple
    sha256: str


def _regex(gaps):
    parts = ["C"]
    for lo, hi in gaps:
        parts.append(f".{{{lo},{hi}}}C")
    return re.compile("".join(parts))


def load_sets(path):
    data = yaml.safe_load(open(path))
    sets, skipped = [], []
    for key, label in CLASSES.items():
        for name, gaps in (data.get(key) or {}).items():
            if any(g not in gaps for g in CORE):
                skipped.append(name)
                continue
            full = tuple(
                (int(gaps[g]["min"]), int(gaps[g]["max"])) if g in gaps else (0, 0) for g in ORDER
            )
            sets.append(SpacingSet(label, name, full, _regex(full)))
    return SpacingSets(tuple(sets), tuple(skipped), sha256_file(path))


def classify(sequence, spacing):
    seq = sequence.upper()
    hit_sets = [s for s in spacing.sets if s.regex.search(seq)]
    classes = sorted({s.cls for s in hit_sets}, key=("I", "II", "other").index)
    return {
        "pattern_match": bool(hit_sets),
        "spacing_class": ",".join(classes),
        "sets": ",".join(s.name for s in hit_sets),
        "n_cys": seq.count("C"),
        "length": len(seq),
    }


def cys8_rows(proteins, spacing, sp_calls):
    """``sp_calls`` maps id to the R0 ``call`` value for rows whose R0 state is ``ok``."""
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
            continue
        r = classify(p.sequence, spacing)
        if not r["pattern_match"]:
            hit = "0"
        else:
            sp = sp_calls.get(p.id)
            hit = "" if sp is None else ("1" if sp == "called" else "0")
        rows.append(
            {
                "id": p.id,
                "state": "ok",
                "hit": hit,
                "pattern_match": "1" if r["pattern_match"] else "0",
                "spacing_class": r["spacing_class"],
                "sets": r["sets"],
                "n_cys": str(r["n_cys"]),
                "length": str(r["length"]),
            }
        )
    return rows


def module_params(spacing, condition):
    """Parameters that decide the result: the spacing file hash, the sets used and the R0 identity."""
    return {
        "spacing_sha256": spacing.sha256,
        "sets": [s.name for s in spacing.sets],
        "skipped_sets": list(spacing.skipped),
        "rule": "any one published set; hit = pattern match and R0 called",
        "conditions": {"sp_module": condition},
    }
