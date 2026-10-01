import pytest
import seqwindow


def _protein(n: int) -> str:
    # A sequence whose residue at position i (0-based) is unambiguous for slicing checks.
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(alphabet[i % 20] for i in range(n))


def test_cterm_window_takes_last_1022():
    seq = "M" * 500 + _protein(1022)
    got = seqwindow.cterm_window(seq)
    assert len(got) == 1022
    assert got == seq[500:]
    assert got[-1] == seq[-1] and got[0] == seq[500]


def test_nterm_window_takes_first_1022():
    seq = _protein(1500)
    assert seqwindow.nterm_window(seq) == seq[:1022]


@pytest.mark.parametrize("n", [1, 40, 1021, 1022])
def test_short_protein_windows_are_the_whole_sequence(n):
    seq = _protein(n)
    assert seqwindow.nterm_window(seq) == seq
    assert seqwindow.cterm_window(seq) == seq
    assert not seqwindow.needs_cterm(seq)


def test_needs_cterm_only_above_the_limit():
    assert not seqwindow.needs_cterm("A" * 1022)
    assert seqwindow.needs_cterm("A" * 1023)


def test_window_dispatch_and_unknown_name():
    seq = _protein(1100)
    assert seqwindow.window(seq, "nterm") == seq[:1022]
    assert seqwindow.window(seq, "cterm") == seq[78:]
    with pytest.raises(ValueError, match="unknown window"):
        seqwindow.window(seq, "middle")


def test_bad_characters_lists_what_esm_would_change():
    assert seqwindow.bad_characters("MKXBUZO") == set()
    assert seqwindow.bad_characters("MK-J.*") == {"-", "J", ".", "*"}


def test_constants_match_the_package():
    # The package imports torch; skip where it is not installed (py3.12 stdlib runs).
    pytest.importorskip("torch")
    pytest.importorskip("esm")
    from surface_glyco import card, embeddings

    assert seqwindow.MAX_RESIDUES == card.MAX_RESIDUES
    assert set(embeddings.ESM_RESIDUES) == set(seqwindow.ESM_RESIDUES)
