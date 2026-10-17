"""Kleene three-valued logic for category calls, with a fourth value for species-limited calls.

Values are the strings used in every output file: ``called`` (true), ``not_called`` (false),
``not_assessable`` (unknown) and ``not_applicable`` (the call cannot hold for this protein, for
example a rule that exists for one species only).

``not_applicable`` acts as ``not_called`` inside ``and`` and ``or``, so a call can never become true
because an input does not apply. ``and`` and ``or`` never return it. ``not`` of ``not_applicable``
is ``not_applicable``: the value must not turn into ``called``. A leaf call that reads a row that
is not applicable returns ``not_applicable`` (the engine does that).
"""

CALLED = "called"
NOT_CALLED = "not_called"
NOT_ASSESSABLE = "not_assessable"
NOT_APPLICABLE = "not_applicable"
VALUES = (CALLED, NOT_CALLED, NOT_ASSESSABLE, NOT_APPLICABLE)


def _check(*values):
    for v in values:
        if v not in VALUES:
            raise ValueError(f"not a call value: {v!r}")


def _as_inputs(values):
    """Check the values; not_applicable acts as not_called inside and/or."""
    _check(*values)
    return tuple(NOT_CALLED if v == NOT_APPLICABLE else v for v in values)


def k_not(a):
    _check(a)
    if a == CALLED:
        return NOT_CALLED
    if a == NOT_CALLED:
        return CALLED
    return a  # not_assessable stays unknown; not_applicable stays not applicable


def k_and(*values):
    """False if any input is false; else unknown if any is unknown; else true."""
    values = _as_inputs(values)
    if NOT_CALLED in values:
        return NOT_CALLED
    if NOT_ASSESSABLE in values:
        return NOT_ASSESSABLE
    return CALLED


def k_or(*values):
    """True if any input is true; else unknown if any is unknown; else false."""
    values = _as_inputs(values)
    if CALLED in values:
        return CALLED
    if NOT_ASSESSABLE in values:
        return NOT_ASSESSABLE
    return NOT_CALLED
