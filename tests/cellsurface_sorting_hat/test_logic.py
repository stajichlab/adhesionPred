"""The Kleene tables of spec section 3.3, every case."""

import itertools

import pytest

from cellsurface_sorting_hat.logic import (
    CALLED as T,
)
from cellsurface_sorting_hat.logic import (
    NOT_ASSESSABLE as U,
)
from cellsurface_sorting_hat.logic import (
    NOT_CALLED as F,
)
from cellsurface_sorting_hat.logic import (
    k_and,
    k_not,
    k_or,
)

AND = {
    (T, T): T,
    (T, F): F,
    (T, U): U,
    (F, T): F,
    (F, F): F,
    (F, U): F,
    (U, T): U,
    (U, F): F,
    (U, U): U,
}
OR = {
    (T, T): T,
    (T, F): T,
    (T, U): T,
    (F, T): T,
    (F, F): F,
    (F, U): U,
    (U, T): T,
    (U, F): U,
    (U, U): U,
}


@pytest.mark.parametrize("a,b", list(itertools.product([T, F, U], repeat=2)))
def test_and_or_tables(a, b):
    assert k_and(a, b) == AND[(a, b)]
    assert k_or(a, b) == OR[(a, b)]


@pytest.mark.parametrize("a,expected", [(T, F), (F, T), (U, U)])
def test_not(a, expected):
    assert k_not(a) == expected


def test_a_true_input_keeps_an_or_true_when_another_input_is_unknown():
    assert k_or(U, T) == T


def test_a_false_input_makes_an_and_false_when_another_input_is_unknown():
    assert k_and(U, F) == F


def test_empty_and_or():
    assert k_and() == T
    assert k_or() == F


def test_rejects_values_that_are_not_call_values():
    with pytest.raises(ValueError):
        k_and("yes", T)
