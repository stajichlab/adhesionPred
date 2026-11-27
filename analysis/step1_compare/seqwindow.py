"""ESM-2 input windows for Phase B (spec 5, decision Q10). Standard library only.

ESM-2 takes at most MAX_RESIDUES residues (1,024 tokens with BOS and EOS). The N-terminal
window is the first MAX_RESIDUES residues; `esm_embed.get_esm_embeddings` cuts
there too. The C-terminal window is the last MAX_RESIDUES residues. It exists only for proteins
longer than MAX_RESIDUES; for shorter proteins both windows are the whole sequence.

Every sequence must hold only ESM_RESIDUES. Then `esm_embed.sanitize_sequence` changes
nothing, and a window cut here is exactly the input the model sees.
"""

# ESM-2 context is 1,024 tokens including BOS and EOS; esm_embed imports these three.
MAX_RESIDUES = 1022
ESM_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWYXBUZO")  # residues the ESM-2 alphabet accepts
# Layer pooled by Phase B: the final layer of the 8M model, the middle layer of the 35M model.
DEFAULT_REPR_LAYER = 6
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
