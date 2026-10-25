"""Per-residue disorder scores (AIUPred 2.x, IUPred3) -> module ``disorder`` (longest disordered run, fraction, mean).

The input is the raw per-residue text of the tool: record lines ``#>ID`` (AIUPred 2.x) or ``>ID`` (IUPred3 through the job script), score rows
``pos aa score`` (separated by white space) and other ``#`` lines or blank lines, which are ignored. The file is READ AS A STREAM: for 27,811 proteins it has
13.76 million rows, and only one protein's scores are in memory at a time. A residue counts as disordered at a score of at least ``RESIDUE_THRESHOLD``.
The module has no ``call`` column; calls read ``longest_run`` through ``{test: ...}`` in ``categories.yaml``. The file carries the residue letters of the
proteins, so it is compared with the FASTA and never copied into the work directory.
"""

import math

from cellsurface_sorting_hat.modules.base import invalid_row

MODULE = "disorder"
RESIDUE_THRESHOLD = 0.5
SCORE_SLACK = (
    0.05  # the tools overshoot [0, 1] slightly (AIUPred up to 1.023); further out is an error
)
COLUMNS = ["longest_run", "frac", "mean", "n_res"]


class DisorderFormatError(ValueError):
    """The per-residue file does not have the expected shape."""


def _finish(name, residues, scores, expected, where, threshold, out):
    if name in out:
        raise DisorderFormatError(f"{where}: duplicate record {name!r}")
    if not scores:
        raise DisorderFormatError(f"{where}: record {name!r} has no residues")
    if name in expected:
        seq = expected[name]
        if len(scores) != len(seq):
            raise DisorderFormatError(
                f"{where}: record {name!r} has {len(scores)} residues, the FASTA sequence has {len(seq)}"
            )
        if "".join(residues) != seq:
            raise DisorderFormatError(
                f"{where}: the residues of record {name!r} differ from the FASTA sequence"
            )
    longest = run = n_dis = 0
    for s in scores:
        if s >= threshold:
            run += 1
            n_dis += 1
            longest = max(longest, run)
        else:
            run = 0
    out[name] = {
        "longest_run": longest,
        "frac": n_dis / len(scores),
        "mean": sum(scores) / len(scores),
        "n_res": len(scores),
    }


def parse_disorder_stream(lines, expected=None, threshold=RESIDUE_THRESHOLD, label="disorder"):
    """Return ``{protein ID: {"longest_run", "frac", "mean", "n_res"}}`` from an iterable of lines.

    ``expected`` maps a protein ID to its FASTA sequence (checked for length and residues; a record that is not in it is
    returned and rejected by the caller). Stops on scores more than ``SCORE_SLACK`` outside [0, 1], positions out of order,
    a malformed row, a row before any record, a duplicate record, a record without residues or a residue-count or residue
    mismatch with ``expected``.
    """
    expected = expected or {}
    out, name, residues, scores = {}, None, [], []
    for n, line in enumerate(lines, 1):
        text = line.strip()
        where = f"{label}:{n}"
        if not text:
            continue
        if text.startswith("#>") or text.startswith(">"):
            if name is not None:
                _finish(name, residues, scores, expected, where, threshold, out)
            head = text[2:] if text.startswith("#>") else text[1:]
            parts = head.split()
            if not parts:
                raise DisorderFormatError(f"{where}: empty record name")
            name, residues, scores = parts[0], [], []
            continue
        if text.startswith("#"):
            continue
        if name is None:
            raise DisorderFormatError(f"{where}: a score row before any record line")
        f = text.split()
        if len(f) != 3 or len(f[1]) != 1:
            raise DisorderFormatError(f"{where}: expected 'pos aa score', got {text!r}")
        try:
            pos, score = int(f[0]), float(f[2])
        except ValueError:
            raise DisorderFormatError(
                f"{where}: position or score is not a number: {text!r}"
            ) from None
        if pos != len(scores) + 1:
            raise DisorderFormatError(
                f"{where}: position {pos} out of order (expected {len(scores) + 1})"
            )
        if not math.isfinite(score) or not -SCORE_SLACK <= score <= 1.0 + SCORE_SLACK:
            raise DisorderFormatError(f"{where}: score {f[2]} is outside [0, 1]")
        residues.append(f[1])
        scores.append(score)
    if name is not None:
        _finish(name, residues, scores, expected, f"{label}:end", threshold, out)
    if not out:
        raise DisorderFormatError(f"{label}: no records")
    return out


def parse_disorder(path, expected=None, threshold=RESIDUE_THRESHOLD):
    with open(path, encoding="utf-8-sig") as fh:
        return parse_disorder_stream(fh, expected, threshold, label=str(path))


def disorder_rows(proteins, parsed):
    rows = []
    for p in proteins:
        if p.state != "ok":
            rows.append(invalid_row(p))
        elif p.id not in parsed:
            rows.append({"id": p.id, "state": "error"})
        else:
            h = parsed[p.id]
            rows.append(
                {
                    "id": p.id,
                    "state": "ok",
                    "longest_run": h["longest_run"],
                    "frac": f"{h['frac']:.4f}",
                    "mean": f"{h['mean']:.4f}",
                    "n_res": h["n_res"],
                }
            )
    return rows
