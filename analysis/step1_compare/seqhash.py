"""Sequence cleaning and hashing shared by D1 step 2 and D10. Standard library only.

`clean` matches surface_glyco.io.clean_sequence (J -> L, '*' removed) and also upper-cases
and removes whitespace, so the same protein gives the same hash from any FASTA source.
"""

import hashlib


def clean(sequence: str) -> str:
    return "".join(sequence.split()).upper().replace("J", "L").replace("*", "")


def seq_sha256(sequence: str) -> str:
    return hashlib.sha256(clean(sequence).encode("ascii")).hexdigest()
