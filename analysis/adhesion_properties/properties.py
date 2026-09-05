"""Sequence property computation for adhesion-protein analysis."""

_KYTE_DOOLITTLE = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}
_AROMATIC = set("FWY")
_POSITIVE = set("KR")
_NEGATIVE = set("DE")

_NAN_RESULT = {
    "pct_ser": float("nan"),
    "pct_thr": float("nan"),
    "pct_pro": float("nan"),
    "pct_ser_thr_pro": float("nan"),
    "pct_cys": float("nan"),
    "aromaticity": float("nan"),
    "mean_hydrophobicity": float("nan"),
    "net_charge_ph7": float("nan"),
}


def sequence_properties(peptide: str) -> dict[str, float]:
    """Compute composition/property features for one peptide sequence.

    All percentage features are in [0, 100]. Unrecognized residue codes
    (X, *, etc.) are counted toward length/percentages but contribute 0 to
    hydrophobicity, not excluded — a coarse approximation, not a crash.
    net_charge_ph7 is a simple Lys+Arg minus Asp+Glu residue count, not a
    full pKa model — stated as an approximation in the report.
    """
    seq = peptide.upper()
    n = len(seq)
    if n == 0:
        return dict(_NAN_RESULT)

    counts: dict[str, int] = {}
    for aa in seq:
        counts[aa] = counts.get(aa, 0) + 1

    pct_ser = 100.0 * counts.get("S", 0) / n
    pct_thr = 100.0 * counts.get("T", 0) / n
    pct_pro = 100.0 * counts.get("P", 0) / n
    pct_cys = 100.0 * counts.get("C", 0) / n
    aromaticity = 100.0 * sum(counts.get(aa, 0) for aa in _AROMATIC) / n
    hydrophobicity_sum = sum(_KYTE_DOOLITTLE.get(aa, 0.0) * c for aa, c in counts.items())
    mean_hydrophobicity = hydrophobicity_sum / n
    net_charge = sum(counts.get(aa, 0) for aa in _POSITIVE) - sum(
        counts.get(aa, 0) for aa in _NEGATIVE
    )

    return {
        "pct_ser": pct_ser,
        "pct_thr": pct_thr,
        "pct_pro": pct_pro,
        "pct_ser_thr_pro": pct_ser + pct_thr + pct_pro,
        "pct_cys": pct_cys,
        "aromaticity": aromaticity,
        "mean_hydrophobicity": mean_hydrophobicity,
        "net_charge_ph7": float(net_charge),
    }
