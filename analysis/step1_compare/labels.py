"""Step 1 label rules (spec section 2.2) on GO ancestor sets. Standard library only.

`label_terms` are the ancestors of the rows that an evidence policy accepts. `any_terms` are the
ancestors of all rows of the gene, IEA included. The rules copy /tmp/glyco_spec/d1_count.py.
"""

from gaf import EXPERIMENTAL_CODES, HOMOLOGY_CODES

WALL = "GO:0005618"
EXTRACELLULAR = "GO:0005576"
PLASMA_MEMBRANE = "GO:0005886"
INTERNAL = frozenset({"GO:0005829", "GO:0005634", "GO:0005739"})  # cytosol, nucleus, mito
SECRETORY = frozenset({"GO:0012505", PLASMA_MEMBRANE, "GO:0005773"})  # endomembrane, PM, vacuole
ANY_SECRETORY = SECRETORY | {"GO:0016020", "GO:0071944", WALL, EXTRACELLULAR}
SURFACE = frozenset({WALL, EXTRACELLULAR})
HIGH_THROUGHPUT_CODES = frozenset({"HDA", "HMP", "HEP", "HGI", "HTP"})

P_EXT = "P-ext"
AMBIGUOUS = "ambiguous"
N_INT = "N-int"
N_SEC = "N-sec"
UNLABELLED = "unlabelled"
LABELS = (P_EXT, AMBIGUOUS, N_INT, N_SEC, UNLABELLED)

POLICIES = ("non_iea", "no_homology", "experimental")


def policy_accepts(policy: str, evidence: str) -> bool:
    if policy == "non_iea":
        return evidence != "IEA"
    if policy == "no_homology":
        return evidence != "IEA" and evidence not in HOMOLOGY_CODES
    if policy == "experimental":
        return evidence in EXPERIMENTAL_CODES
    raise ValueError(f"unknown evidence policy {policy!r}")


def classify(label_terms: frozenset[str] | set[str], any_terms: frozenset[str] | set[str]) -> str:
    surface = bool(label_terms & SURFACE)
    internal = bool(label_terms & INTERNAL)
    if surface:
        return AMBIGUOUS if internal else P_EXT
    if internal and not (any_terms & ANY_SECRETORY):
        return N_INT
    if label_terms & SECRETORY and not (any_terms & SURFACE):
        return N_SEC
    return UNLABELLED


def subset_of(label: str, label_terms: frozenset[str] | set[str]) -> str:
    if label != P_EXT:
        return ""
    return "wall" if WALL in label_terms else "extracellular-only"


def htp_only(internal_codes: set[str]) -> bool:
    """True if the gene has non-IEA internal evidence and all of it is high-throughput (R-A)."""
    return bool(internal_codes) and internal_codes <= HIGH_THROUGHPUT_CODES


def is_pm_candidate(label: str, label_terms: frozenset[str] | set[str]) -> bool:
    """P-ext with a non-IEA plasma-membrane term: the input set of D8 triage."""
    return label == P_EXT and PLASMA_MEMBRANE in label_terms
