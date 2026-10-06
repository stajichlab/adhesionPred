"""Kleene three-valued logic for category calls.

Values are the strings used in every output file: ``called`` (true), ``not_called`` (false) and
``not_assessable`` (unknown).
"""

CALLED = "called"
NOT_CALLED = "not_called"
NOT_ASSESSABLE = "not_assessable"
VALUES = (CALLED, NOT_CALLED, NOT_ASSESSABLE)


def _check(*values):
    for v in values:
        if v not in VALUES:
            raise ValueError(f"not a call value: {v!r}")


def k_not(a):
    _check(a)
    if a == CALLED:
        return NOT_CALLED
    if a == NOT_CALLED:
        return CALLED
    return NOT_ASSESSABLE


def k_and(*values):
    """False if any input is false; else unknown if any is unknown; else true."""
    _check(*values)
    if NOT_CALLED in values:
        return NOT_CALLED
    if NOT_ASSESSABLE in values:
        return NOT_ASSESSABLE
    return CALLED


def k_or(*values):
    """True if any input is true; else unknown if any is unknown; else false."""
    _check(*values)
    if CALLED in values:
        return CALLED
    if NOT_ASSESSABLE in values:
        return NOT_ASSESSABLE
    return NOT_CALLED
