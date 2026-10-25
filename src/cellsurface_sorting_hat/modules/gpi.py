"""PredGPI gff3 results -> module ``gpi`` (a GPI-anchor call and its confidence class).

PredGPI writes one gff3 row per protein: type ``Chain`` (no GPI anchor; start 1, end = protein length) or ``GPI-anchor`` (start = end = omega
site). The score column is the confidence class of the tool (1.0, 0.70 or 0.55 for an estimated false positive rate of at most 0.0015, 0.005
or 0.01), not a probability. PredGPI never scores a sequence of 40 residues or fewer: it writes a silent ``Chain`` row, which this wrapper turns
into ``not_assessable`` (never ``not_called``).
"""

import math
from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

MODULE = "gpi"
MIN_LEN = 41  # PredGPI scores a sequence only if it has more than 40 residues
CLASSES = {  # gff3 score -> (confidence class, upper bound of the estimated false positive rate)
    "1.0": ("high", "0.0015"),
    "0.7": ("probable", "0.005"),
    "0.70": ("probable", "0.005"),
    "0.55": ("weak", "0.01"),
}
COLUMNS = ["call", "omega_site", "confidence_class", "fpr_upper_bound"]


class GpiFormatError(ValueError):
    """The PredGPI file does not have the expected shape."""


def parse_gpi(path, lengths):
    """Return ``{protein ID: {"kind": "Chain"|"GPI-anchor", "omega": int|None, "cls": str, "ub": str}}``.

    ``lengths`` maps a protein ID to its length (a ``Chain`` row must span the protein). Stops on a row with the
    wrong number of columns, another source or type, a ``GPI-anchor`` with start != end or outside the protein,
    a ``Chain`` that does not span the protein, an unknown confidence score, a duplicate ID or an ID that is not in
    ``lengths``.
    """
    out = {}
    for n, raw in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        f = raw.split("\t")
        where = f"{path}:{n}"
        if len(f) != 9:
            raise GpiFormatError(f"{where}: expected 9 TAB separated gff3 columns, got {len(f)}")
        pid, source, kind, start, end, score = f[:6]
        if not pid.strip() or pid != pid.strip():
            raise GpiFormatError(f"{where}: empty or padded ID {pid!r}")
        if pid in out:
            raise GpiFormatError(f"{where}: duplicate ID {pid!r}")
        if source != "PredGPI":
            raise GpiFormatError(f"{where}: source column is {source!r}, expected PredGPI")
        if pid not in lengths:
            raise GpiFormatError(f"{where}: {pid!r} is not a protein of the FASTA")
        try:
            start, end = int(start), int(end)
        except ValueError:
            raise GpiFormatError(f"{where}: start and end must be integers") from None
        length = lengths[pid]
        if kind == "Chain":
            if start != 1 or end != length:
                raise GpiFormatError(
                    f"{where}: Chain {start}-{end} does not span {pid!r} (length {length})"
                )
            out[pid] = {"kind": kind, "omega": None, "cls": "", "ub": ""}
        elif kind == "GPI-anchor":
            if start != end or not 1 <= start <= length:
                raise GpiFormatError(
                    f"{where}: GPI-anchor {start}-{end} is not one site inside {pid!r} (length {length})"
                )
            try:
                value = float(score)
            except ValueError:
                raise GpiFormatError(
                    f"{where}: confidence score {score!r} is not a number"
                ) from None
            if not math.isfinite(value) or score not in CLASSES:
                raise GpiFormatError(
                    f"{where}: unknown confidence score {score!r} (expected 1.0, 0.70 or 0.55)"
                )
            cls, ub = CLASSES[score]
            out[pid] = {"kind": kind, "omega": start, "cls": cls, "ub": ub}
        else:
            raise GpiFormatError(
                f"{where}: unknown row type {kind!r} (expected Chain or GPI-anchor)"
            )
    if not out:
        raise GpiFormatError(f"{path}: no rows")
    return out


def gpi_rows(proteins, parsed):
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in parsed:
            rows.append({"id": p.id, "state": "error"})
        elif len(p.sequence) < MIN_LEN:
            rows.append({"id": p.id, "state": "not_assessable"})
        else:
            hit = parsed[p.id]
            anchored = hit["kind"] == "GPI-anchor"
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "call": "called" if anchored else "not_called",
                    "omega_site": hit["omega"] if anchored else "",
                    "confidence_class": hit["cls"],
                    "fpr_upper_bound": hit["ub"],
                }
            )
    return rows
