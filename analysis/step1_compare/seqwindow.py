"""ESM-2 input windows for Phase B (spec 5, decision Q10). Standard library only.

ESM-2 takes at most MAX_RESIDUES residues (1,024 tokens with BOS and EOS). The N-terminal
window is the first MAX_RESIDUES residues; `surface_glyco.embeddings.get_esm_embeddings` cuts
there too. The C-terminal window is the last MAX_RESIDUES residues. It exists only for proteins
longer than MAX_RESIDUES; for shorter proteins both windows are the whole sequence.

Every sequence must hold only ESM_RESIDUES. Then `sanitize_sequence` in the package changes
nothing, and a window cut here is exactly the input the model sees.
"""

MAX_RESIDUES = 1022  # same value as surface_glyco.card.MAX_RESIDUES
ESM_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWYXBUZO")  # surface_glyco.embeddings.ESM_RESIDUES
WINDOWS = ("nterm", "cterm")


def bad_characters(sequence: str) -> set[str]:
    """Characters of `sequence` that ESM-2 does not accept as they are."""
    return set(sequence) - ESM_RESIDUES


def needs_cterm(sequence: str) -> bool:
    return len(sequence) > MAX_RESIDUES


def nterm_window(sequence: str) -> str:
    return sequence[:MAX_RESIDUES]


def cterm_window(sequence: str) -> str:
    return sequence[-MAX_RESIDUES:]


def window(sequence: str, name: str) -> str:
    if name == "nterm":
        return nterm_window(sequence)
    if name == "cterm":
        return cterm_window(sequence)
    raise ValueError(f"unknown window {name!r}; valid: {', '.join(WINDOWS)}")
