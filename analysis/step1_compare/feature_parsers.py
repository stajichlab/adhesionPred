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
                cs_end, cs_prob = int(match.group(1)), float(match.group(3))
            elif cs:
                raise OutputFormatError(f"{path}:{lineno}: OTHER row with a CS field")
            if row["ID"] in calls:
                raise OutputFormatError(f"{path}:{lineno}: duplicate id {row['ID']}")
            calls[row["ID"]] = SignalPCall(
                pred, float(row["OTHER"]), float(row["SP(Sec/SPI)"]), cs_end, cs_prob
            )
    if header is None:
        raise OutputFormatError(f"{path}: no '# ID' header line")
    return calls


def parse_signalp_gff(path: str | Path) -> dict[str, int]:
    """id -> end of the signal_peptide feature."""
    ends: dict[str, int] = {}
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            f = raw.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "signal_peptide":
                raise OutputFormatError(f"{path}: unexpected GFF line {raw[:60]!r}")
            ends[f[0]] = int(f[4])
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


def parse_predgpi_scores(path: str | Path) -> dict[str, GpiCall]:
    calls: dict[str, GpiCall] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != GPI_COLUMNS:
            raise OutputFormatError(f"{path}: columns {reader.fieldnames} != {GPI_COLUMNS}")
        for r in reader:
            if r["gpi_call"] not in GPI_CALLS:
                raise OutputFormatError(f"{path}: unknown gpi_call {r['gpi_call']!r}")
            if r["id"] in calls:
                raise OutputFormatError(f"{path}: duplicate id {r['id']}")
            calls[r["id"]] = GpiCall(
                r["gpi_call"],
                float(r["gpi_prob"]),
                int(r["omega"]) if r["omega"] else None,
                float(r["fpr"]) if r["fpr"] else None,
                float(r["svm"]) if r["svm"] else None,
            )
    return calls


def ser_thr_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.count("S") + sequence.count("T")) / len(sequence)
