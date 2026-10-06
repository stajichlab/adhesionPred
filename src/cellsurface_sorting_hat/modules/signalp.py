"""SignalP 6 results -> module ``step1_rule@R0`` (rule R0: SignalP calls a signal peptide, ``SP``)."""

from pathlib import Path

from cellsurface_sorting_hat.modules.base import invalid_row

MODULE = "step1_rule@R0"


class SignalPFormatError(ValueError):
    """The SignalP file does not have the expected shape."""


def parse_signalp(path):
    """Return ``{protein ID: {"prediction": str, "sp_prob": float}}``.

    ``prediction_results.txt`` has two ``#`` header lines, then one TAB separated row per protein:
    ID (the FASTA header, the ID is its first token), prediction, OTHER probability, SP(Sec/SPI)
    probability, cleavage site. Prediction ``SP`` is the Sec/SPI signal peptide that rule R0 uses.
    """
    out = {}
    for n, line in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        if line.startswith("#"):
            if "Organism:" in line and "Eukarya" not in line:
                raise SignalPFormatError(f"{path}:{n}: not a Eukarya run: {line.strip()!r}")
            continue
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise SignalPFormatError(f"{path}:{n}: expected at least 4 TAB separated fields")
        pid = fields[0].split()[0]
        try:
            sp_prob = float(fields[3])
        except ValueError:
            raise SignalPFormatError(f"{path}:{n}: SP probability is not a number") from None
        if pid in out:
            raise SignalPFormatError(f"{path}:{n}: duplicate ID {pid!r}")
        out[pid] = {"prediction": fields[1], "sp_prob": sp_prob}
    if not out:
        raise SignalPFormatError(f"{path}: no prediction rows")
    return out


def signalp_rows(proteins, parsed):
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in parsed:
            rows.append({"id": p.id, "state": "error"})
        else:
            hit = parsed[p.id]
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "call": "called" if hit["prediction"] == "SP" else "not_called",
                    "sp_prob": f"{hit['sp_prob']:.6f}",
                    "prediction": hit["prediction"],
                }
            )
    return rows


COLUMNS = ["call", "sp_prob", "prediction"]
