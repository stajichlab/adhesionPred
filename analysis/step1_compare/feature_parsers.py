"""Parse SignalP 6 and PredGPI outputs, and compute Ser+Thr fraction. Standard library only.

SignalP 6 (`--organism eukarya --mode fast`) writes `prediction_results.txt`: two `#` header
lines, the second `# ID<TAB>Prediction<TAB>OTHER<TAB>SP(Sec/SPI)<TAB>CS Position`, then one TAB
separated row per protein. The CS column is `CS pos: 24-25. Pr: 0.5487` for SP rows and empty
for OTHER rows. Columns are found by header name, never by whitespace splitting (a whitespace
split mis-columns the CS field; see analysis/cocci_repeats/01_signalp.sh). `output.gff3` has
one `signal_peptide` line per SP protein (start 1, end = last residue of the signal peptide).

PredGPI rows come from jobs/predgpi_scores.py (columns id, length, gpi_call, gpi_prob, omega,
fpr, svm). Both readers accept plain or gzip files (gzip found by its magic bytes).
"""

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from gaf import open_text

SP_COLUMNS = ("ID", "Prediction", "OTHER", "SP(Sec/SPI)", "CS Position")
_CS = re.compile(r"^CS pos: (\d+)-(\d+)\. Pr: ([0-9.]+)$")
GPI_COLUMNS = ("id", "length", "gpi_call", "gpi_prob", "omega", "fpr", "svm")
GPI_CALLS = ("highly_probable", "probable", "weakly", "none", "too_short")


class OutputFormatError(ValueError):
    """A tool output does not have the expected format."""


def _num(convert, text: str, where: str, what: str):
    """convert(text), or raise OutputFormatError('<where>: <what> <text> is not a number')."""
    try:
        return convert(text)
    except ValueError as exc:
        raise OutputFormatError(f"{where}: {what} {text!r} is not a valid number") from exc


@dataclass(frozen=True, slots=True)
class SignalPCall:
    prediction: str
    other_prob: float
    sp_prob: float
    cs_end: int | None
    cs_prob: float | None


def parse_signalp(path: str | Path) -> dict[str, SignalPCall]:
    calls: dict[str, SignalPCall] = {}
    header = None
    with open_text(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if line.startswith("# ID"):
                header = line[2:].split("\t")
                missing = [c for c in SP_COLUMNS if c not in header]
                if missing:
                    raise OutputFormatError(f"{path}: header lacks {missing}")
                continue
            if line.startswith("#") or not line:
                continue
            if header is None:
                raise OutputFormatError(f"{path}:{lineno}: data row before the '# ID' header")
            fields = line.split("\t")
            if len(fields) != len(header):
                raise OutputFormatError(f"{path}:{lineno}: {len(fields)} fields, not {len(header)}")
            row = dict(zip(header, fields, strict=True))
            pred, cs = row["Prediction"], row["CS Position"].strip()
            if pred not in ("SP", "OTHER"):
                raise OutputFormatError(f"{path}:{lineno}: unexpected prediction {pred!r}")
            cs_end = cs_prob = None
            if pred == "SP":
                match = _CS.match(cs)
                if not match:
                    raise OutputFormatError(f"{path}:{lineno}: SP row with CS field {cs!r}")
                where = f"{path}:{lineno}"
                cs_end = _num(int, match.group(1), where, "CS position")
                cs_prob = _num(float, match.group(3), where, "CS probability")
            elif cs:
                raise OutputFormatError(f"{path}:{lineno}: OTHER row with a CS field")
            if row["ID"] in calls:
                raise OutputFormatError(f"{path}:{lineno}: duplicate id {row['ID']}")
            where = f"{path}:{lineno}"
            calls[row["ID"]] = SignalPCall(
                pred,
                _num(float, row["OTHER"], where, "OTHER probability"),
                _num(float, row["SP(Sec/SPI)"], where, "SP(Sec/SPI) probability"),
                cs_end,
                cs_prob,
            )
    if header is None:
        raise OutputFormatError(f"{path}: no '# ID' header line")
    return calls


def parse_signalp_gff(path: str | Path) -> dict[str, int]:
    """id -> end of the signal_peptide feature."""
    ends: dict[str, int] = {}
    with open_text(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            if raw.startswith("#") or not raw.strip():
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "signal_peptide":
                raise OutputFormatError(f"{path}:{lineno}: unexpected GFF line {raw[:60]!r}")
            if f[0] in ends:
                raise OutputFormatError(f"{path}:{lineno}: duplicate id {f[0]}")
            ends[f[0]] = _num(int, f[4], f"{path}:{lineno}", "GFF end")
    return ends


def check_signalp_consistency(calls: dict[str, SignalPCall], gff_ends: dict[str, int]) -> None:
    sp = {k: c.cs_end for k, c in calls.items() if c.prediction == "SP"}
    if sp != gff_ends:
        diff = sorted(set(sp.items()) ^ set(gff_ends.items()))[:5]
        raise OutputFormatError(f"prediction_results.txt and output.gff3 disagree: {diff}")


@dataclass(frozen=True, slots=True)
class GpiCall:
    call: str
    prob: float
    omega: int | None
    fpr: float | None
    svm: float | None


# gpi_prob and the allowed estimated false positive rate (fpr) of each call. The bounds are
# inclusive on both sides because jobs/predgpi_scores.py writes fpr with 6 significant digits.
GPI_PROB = {"highly_probable": 1.0, "probable": 0.70, "weakly": 0.55, "none": 0.0, "too_short": 0.0}
GPI_FPR_RANGE = {
    "highly_probable": (0.0, 0.0015),
    "probable": (0.0015, 0.005),
    "weakly": (0.005, 0.01),
    "none": (0.01, 1.0),
}


def _check_gpi_row(where: str, call: GpiCall) -> None:
    if call.prob != GPI_PROB[call.call]:
        raise OutputFormatError(
            f"{where}: gpi_call {call.call} has gpi_prob {call.prob}, expected {GPI_PROB[call.call]}"
        )
    if call.call == "too_short":
        if call.omega is not None or call.fpr is not None or call.svm is not None:
            raise OutputFormatError(f"{where}: too_short row has omega, fpr or svm")
        return
    if call.fpr is None or call.svm is None:
        raise OutputFormatError(f"{where}: {call.call} row lacks fpr or svm")
    low, high = GPI_FPR_RANGE[call.call]
    if not low <= call.fpr <= high:
        raise OutputFormatError(
            f"{where}: gpi_call {call.call} with fpr {call.fpr}, expected {low} to {high}"
        )
    if call.call == "none":
        if call.omega is not None:
            raise OutputFormatError(f"{where}: none row has omega {call.omega}")
    elif call.omega is None or call.omega < 0:
        raise OutputFormatError(f"{where}: {call.call} row has omega {call.omega}")


def parse_predgpi_scores(path: str | Path) -> dict[str, GpiCall]:
    calls: dict[str, GpiCall] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != GPI_COLUMNS:
            raise OutputFormatError(f"{path}: columns {reader.fieldnames} != {GPI_COLUMNS}")
        for r in reader:
            where = f"{path}:{reader.line_num}"
            if r["gpi_call"] not in GPI_CALLS:
                raise OutputFormatError(f"{where}: unknown gpi_call {r['gpi_call']!r}")
            if r["id"] in calls:
                raise OutputFormatError(f"{where}: duplicate id {r['id']}")
            call = GpiCall(
                r["gpi_call"],
                _num(float, r["gpi_prob"], where, "gpi_prob"),
                _num(int, r["omega"], where, "omega") if r["omega"] else None,
                _num(float, r["fpr"], where, "fpr") if r["fpr"] else None,
                _num(float, r["svm"], where, "svm") if r["svm"] else None,
            )
            _check_gpi_row(where, call)
            calls[r["id"]] = call
    return calls


def ser_thr_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.count("S") + sequence.count("T")) / len(sequence)
