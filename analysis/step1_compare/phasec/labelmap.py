"""Phase C class mapping (Phase C spec 3.1): one function maps a truth row to its class.

| label     | d8_class        | class    | reporting stratum                   |
|-----------|-----------------|----------|-------------------------------------|
| P-ext     | "" or P-gpi     | pos      | subset (wall, extracellular-only)   |
| P-ext     | PM-TM           | neg      | PM-TM                               |
| P-ext     | pm-unresolved   | excluded | pm-unresolved (ruling C-1)          |
| N-int     | any             | neg      | N-int                               |
| N-sec     | any             | neg      | N-sec                               |
| ambiguous | any             | excluded | ambiguous                           |
| unlabelled| any             | excluded | unlabelled (never a negative)       |

Read `subset`, not `stratum`, for positives: 03_triage_pm.py writes d8_class over stratum.
"""

POS = "pos"
NEG = "neg"
EXCLUDED = "excluded"
CLASSES = (POS, NEG, EXCLUDED)
D8_CLASSES = ("", "P-gpi", "PM-TM", "pm-unresolved")
LABELS = ("P-ext", "N-int", "N-sec", "ambiguous", "unlabelled")
SUBSETS = ("wall", "extracellular-only")


def class_of(label: str, d8_class: str) -> str:
    """Return pos, neg or excluded. ValueError for a label or P-ext d8_class not in the table."""
    if label == "P-ext":
        if d8_class in ("", "P-gpi"):
            return POS
        if d8_class == "PM-TM":
            return NEG
        if d8_class == "pm-unresolved":
            return EXCLUDED
        raise ValueError(f"unknown d8_class {d8_class!r} for a P-ext gene")
    if label in ("N-int", "N-sec"):
        return NEG
    if label in ("ambiguous", "unlabelled"):
        return EXCLUDED
    raise ValueError(f"unknown label {label!r}")


def stratum_of(label: str, subset: str, d8_class: str) -> str:
    """Reporting stratum of one truth row (the table above)."""
    cls = class_of(label, d8_class)
    if label == "P-ext" and cls == POS:
        if subset not in SUBSETS:
            raise ValueError(f"P-ext gene with subset {subset!r}; expected one of {SUBSETS}")
        return subset
    if label == "P-ext":
        return d8_class
    return label
