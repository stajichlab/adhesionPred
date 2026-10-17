"""The four-value truth table of spec 3.5 (D18) and the ``other_*`` rule (D5).

``not_applicable`` (NA) acts as ``not_called`` inside ``and`` and ``or``. ``and`` and ``or`` never
return NA. ``not(NA)`` is NA (owner decision of 2026-10-17). The table below is written by hand
from the spec; it is not read from the code.

Spec 3.5 prints the ``not`` column of the NA row as ``C``. That is the first-draft value. The text
under the table and conflict 13 of spec section 9 say ``not(NA)`` is ``not_applicable``; this file
follows the text (see the note in scratchpad pr3/conflicts.md).
"""

import pytest

C, N, U, NA = "called", "not_called", "not_assessable", "not_applicable"

SECOND = (C, N, U, NA)  # order of the second input in the rows below
NOT = {C: N, N: C, U: U, NA: NA}
AND_ROWS = {  # first input -> result for each second input
    C: (C, N, U, N),
    N: (N, N, N, N),
    U: (U, N, U, N),
    NA: (N, N, N, N),
}
OR_ROWS = {
    C: (C, C, C, C),
    N: (C, N, U, N),
    U: (C, U, U, U),
    NA: (C, N, U, N),
}
AND = {(a, b): r for a, row in AND_ROWS.items() for b, r in zip(SECOND, row, strict=True)}
OR = {(a, b): r for a, row in OR_ROWS.items() for b, r in zip(SECOND, row, strict=True)}

XFAIL_NA = pytest.mark.xfail(strict=True, reason="Task 3.2: not_applicable in logic.py")


def _marked(cases):
    for key, want in cases:
        key_t = key if isinstance(key, tuple) else (key,)
        marks = [XFAIL_NA] if NA in key_t else []
        yield pytest.param(key, want, marks=marks, id="-".join(key_t) + "->" + want)


@pytest.mark.parametrize(("a", "want"), list(_marked(NOT.items())))
def test_not(a, want):
    from cellsurface_sorting_hat.logic import k_not

    assert k_not(a) == want


@pytest.mark.parametrize(("pair", "want"), list(_marked(AND.items())))
def test_and(pair, want):
    from cellsurface_sorting_hat.logic import k_and

    assert k_and(*pair) == want


@pytest.mark.parametrize(("pair", "want"), list(_marked(OR.items())))
def test_or(pair, want):
    from cellsurface_sorting_hat.logic import k_or

    assert k_or(*pair) == want


@pytest.mark.xfail(strict=True, reason="Task 3.2: not_applicable in logic.py")
def test_and_and_or_never_return_not_applicable():
    from cellsurface_sorting_hat.logic import k_and, k_or

    for pair in AND:
        assert k_and(*pair) != NA
        assert k_or(*pair) != NA
    assert k_and(NA) == N and k_or(NA) == N  # one operand: NA acts as not_called


@pytest.mark.xfail(strict=True, reason="Task 3.2: not_applicable in logic.py")
def test_the_value_vocabulary_has_four_values():
    from cellsurface_sorting_hat import logic

    assert set(logic.VALUES) == {C, N, U, NA}
    assert logic.NOT_APPLICABLE == NA
