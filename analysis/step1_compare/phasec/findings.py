"""Estimate label and machine-checked findings (Phase C spec 4, rulings C-8, C-11).

All functions read values that 11_evaluate.py has already put into metrics.json, so every
boolean in findings.json can be traced to a number there.
"""

ESTIMATE_HALF_WIDTH = 0.10
MIN_DIRECT_POSITIVES = 20  # owner decision 2026-10-01 (ruling C-8)
SATURATION_AUC = 0.99
LABEL_CANDIDATES = ("R2", "M8", "M35", "M8-C", "M35-C", "H")  # ruling C-8
ML = ("M8", "M35", "M8-C", "M35-C", "H")
COMPARATORS = ("B1", "R2")
N_S2 = 3  # S2-Calb_CGD, S2-Scer_SGD, S2-Spom_PomBase


def floor_met(n_direct_positives: int) -> bool:
    """The test set has at least MIN_DIRECT_POSITIVES direct-evidence positives."""
    return n_direct_positives >= MIN_DIRECT_POSITIVES


def estimate_label(recall_half_widths: dict, n_direct_positives: int) -> str:
    """`estimate` when the recall half-width is <= 0.10 for every listed candidate and the test
    set has at least 20 direct-evidence positives, else `smoke test` (also when a half-width is
    not defined, for example without positives). The count floor stops a zero-width interval
    of a small set (for example 7 of 7 positives found in every resample) from passing."""
    values = [recall_half_widths.get(c) for c in LABEL_CANDIDATES]
    widths_ok = all(v is not None and v <= ESTIMATE_HALF_WIDTH + 1e-12 for v in values)
    if widths_ok and floor_met(n_direct_positives):
        return "estimate"
    return "smoke test"


def finding_a(b1_auc: dict) -> dict:
    """(a) B1 does not saturate: ROC-AUC below 0.99 on every listed test set."""
    holds = bool(b1_auc) and all(v is not None and v < SATURATION_AUC for v in b1_auc.values())
    return {"threshold": SATURATION_AUC, "b1_roc_auc": b1_auc, "holds": holds}


def beats(diff: dict) -> bool:
    """The comparator's N-sec FPR minus the ML candidate's: the interval lies above 0."""
    return diff.get("lo") is not None and diff["lo"] > 0


def finding_b(diffs: dict) -> dict:
    """(b) for one test set: diffs[ml][comparator] = {value, lo, hi}. ML beats B1 and R2 when
    both intervals lie above 0."""
    out = {}
    for ml in ML:
        if ml not in diffs:
            continue
        flags = {f"beats_{c}": beats(diffs[ml].get(c, {})) for c in COMPARATORS}
        out[ml] = {
            **{c: diffs[ml].get(c) for c in COMPARATORS},
            **flags,
            "holds": all(flags.values()),
        }
    return {"candidates": out, "holds_for": [m for m in out if out[m]["holds"]]}


def finding_c(per_test_set: dict) -> dict:
    """(c) the statement of (b) on every S2 test set. per_test_set[ts] = finding_b(...)."""
    names = sorted(per_test_set)
    complete = len(names) == N_S2  # fewer S2 sets than expected: the statement does not hold
    holds_for = [
        m for m in ML if complete and all(m in per_test_set[ts]["holds_for"] for ts in names)
    ]
    return {"test_sets": per_test_set, "n_s2": len(names), "expected_s2": N_S2,
            "holds_for": holds_for}  # fmt: skip
