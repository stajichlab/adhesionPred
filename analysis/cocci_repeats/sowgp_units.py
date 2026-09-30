"""Shared phase-anchored repeat-unit logic for the SOWgp visualizations (09-13).

A SOWgp repeat unit is 47 aa and starts with the motif PTDCYGDC. Units are located
by that anchor, not by dividing the repeat region into equal-length pieces (which
drifts out of phase whenever the profiler's region edges are not unit boundaries).

Unit classes use two diagnostic sites (0-based offsets within the unit):
  offset 8        : E vs K                 (PTDCYGDC[E/K]DG...)
  offsets 20-33   : DDYDG vs DYDDG block
The last unit runs into the C-terminus (..GSPPPKETK..) and is classed 'terminal'.
"""

import re

import pandas as pd

ANCHOR = "PTDCYGDC"
PERIOD = 47
CLASSES = ["E-DDYDG", "E-DYDDG", "K-DYDDG", "K-DDYDG", "terminal"]
# categorical slots 1-4 of the validated default palette; terminal = neutral
CLASS_COLOR = {
    "E-DDYDG": "#2a78d6",
    "E-DYDDG": "#eb6834",
    "K-DYDDG": "#1baf7a",
    "K-DDYDG": "#eda100",
    "terminal": "#b5b4ad",
}
OTHER_COLOR = "#ffffff"  # class not resolved: drawn as an outlined empty block
SPECIES_INK = {"immitis": "#1c5cab", "posadasii": "#a8431c"}
INK, MUTED, GRID = "#222220", "#6b6a64", "#e6e5df"


def unit_class(u):
    if "SPPP" in u[34:40] and "ETK" in u[36:47]:
        return "terminal"
    s9 = u[8] if len(u) > 8 else "?"
    win = u[20:34]
    blk = "DDYDG" if "DDYDG" in win else ("DYDDG" if "DYDDG" in win else "?")
    return f"{s9}-{blk}"


def anchor_units(seq):
    """Return [(start, inferred)] for each unit.

    A unit whose anchor motif is mutated is inferred when two anchor hits are
    2 periods (+-1 aa) apart.
    """
    pos = [m.start() for m in re.finditer(ANCHOR, seq)]
    out = []
    for a, b in zip(pos, pos[1:] + [None]):
        out.append((a, False))
        if b is not None and abs((b - a) - 2 * PERIOD) <= 1:
            out.append((a + PERIOD, True))
    return out


def unit_seqs(seq, units=None):
    """Unit sequences: from each anchor to the next anchor, capped at PERIOD aa."""
    units = anchor_units(seq) if units is None else units
    starts = [p for p, _ in units]
    res = []
    for i, p in enumerate(starts):
        nxt = starts[i + 1] if i + 1 < len(starts) else len(seq)
        res.append(seq[p : min(p + PERIOD, nxt)])
    return res


def is_regular(units):
    starts = [p for p, _ in units]
    return all(abs((b - a) - PERIOD) <= 1 for a, b in zip(starts, starts[1:]))


def ungapped_identity(a, b):
    n = min(len(a), len(b))
    return sum(x == y for x, y in zip(a[:n], b[:n])) / n if n else 0.0


def distinct_alleles(copies_tsv):
    """Collapse the per-strain copies table to distinct protein sequences."""
    d = pd.read_csv(copies_tsv, sep="\t")
    d["strain"] = d["strain"].str.replace(".proteins", "", regex=False)
    hap = (
        d.groupby("seq")
        .agg(
            n=("strain", "size"),
            species=("species", "first"),
            length=("length", "first"),
            strains=("strain", lambda s: ",".join(sorted(s))),
        )
        .reset_index()
    )
    rows = []
    for _, r in hap.iterrows():
        us = anchor_units(r.seq)
        if not us:
            continue
        seqs = unit_seqs(r.seq, us)
        rows.append(
            {
                "seq": r.seq,
                "n": r.n,
                "species": r.species,
                "length": r.length,
                "strains": r.strains,
                "units": us,
                "unit_seqs": seqs,
                "n_units": len(us),
                "regular": is_regular(us),
                "classes": [unit_class(u) for u in seqs],
            }
        )
    H = pd.DataFrame(rows)
    H = H.sort_values(
        ["species", "n_units", "length", "n"], ascending=[True, True, True, False]
    ).reset_index(drop=True)
    H["allele"] = [f"{sp[:3]}{k}u_{i:02d}" for i, (sp, k) in enumerate(zip(H.species, H.n_units))]
    return H


def unit_diff(a, b):
    """aa differences between two units: mismatches + gap columns in a global alignment."""
    from Bio.Align import PairwiseAligner

    pa = PairwiseAligner()
    pa.mode = "global"
    pa.match_score, pa.mismatch_score = 2, -1
    pa.open_gap_score, pa.extend_gap_score = -4, -1
    aln = pa.align(a, b)[0]
    ta, tb = aln[0], aln[1]
    return sum(x != y for x, y in zip(ta, tb))
