"""Raw SignalP 6 / PredGPI / disorder outputs -> one table per protein for the AGA1-class prediction step (stdlib only).

SignalP and PredGPI formats were checked against the real outputs of 2026-10-24 (jobs 30451696, 30453935); the disorder format is
still UNVERIFIED (see analysis/aga1_class/README.md):
  SignalP 6 `prediction_results.txt`: 2 `#` lines (the first holds `Organism: Eukarya`), then TAB rows
      `ID<TAB>Prediction<TAB>OTHER<TAB>SP(Sec/SPI)<TAB>CS Position`, CS like `CS pos: 24-25. Pr: 0.9770` (empty for OTHER).
  PredGPI `predgpi_raw.gff3` (read directly, no adapter): 9 TAB columns `ID PredGPI <type> start end score . . attributes`, one row
      per protein. Type `Chain` = not GPI-anchored (start 1, end = protein length; the score 1.0 means nothing). Type
      `GPI-anchor` = GPI-anchored (start = end = omega site; the score is the confidence class of the tool: 1.0 = estimated false
      positive rate <= 0.0015, 0.70 = <= 0.005, 0.55 = <= 0.01; the tool calls GPI at <= 0.01). The gff3 holds no false positive
      rate and no raw score, so none is invented: the table gives `gpi_confidence_class` (high / probable / weak) and
      `gpi_fpr_upper_bound`. PredGPI never scores a sequence of 40 aa or fewer (it returns a silent "no"), so those get `NA`.
      The protein lengths (lengths.tsv of stage_inputs.py) are therefore required.
  Disorder (AIUPred or IUPred3): per-residue text, one record per protein: a line `>ID ...` followed by rows `pos aa score`
      (whitespace separated, pos 1..n); `#` lines are ignored. If the raw tool prints another layout, convert it first.
The disorder block is optional but explicit: pass either --disorder FILE --disorder-tool T, or --no-disorder (the disorder columns
are then left out of the table entirely, never written as zero or NA).
Fails loudly (non-zero exit, message on stderr) on a missing accession, a duplicate, an empty file, an unparseable row, an
unexpected accession, a length mismatch. Nothing is ever written as a silent zero: a value that does not
apply is written as NA.
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

GPI_CLASSES = {
    "1.0": ("high", 0.0015),
    "0.7": ("probable", 0.005),
    "0.70": ("probable", 0.005),
    "0.55": ("weak", 0.01),
}
PREDGPI_MIN_LEN = 41  # PredGPI scores a sequence only if len(seq) > 40
SHORT_NOTE = "length<=40: PredGPI does not score sequences this short"
CS_RE = re.compile(r"^CS pos:\s*(\d+)-(\d+)\.\s*Pr:\s*([0-9.eE+-]+)\s*$")
OUT_COLS = [
    "accession", "signal_peptide_pred", "sp_cleavage_pos", "sp_cleavage_prob", "sp_prob", "sp_prediction_raw",
    "gpi_pred", "gpi_omega_site", "gpi_confidence_class", "gpi_fpr_upper_bound", "gpi_note",
    "disorder_tool", "disorder_threshold", "disorder_n_res", "disorder_mean", "disorder_frac", "disorder_longest_run",
]  # fmt: skip


class ConvertError(Exception):
    """A tool output does not have the expected shape."""


def _num(text, where, what):
    try:
        v = float(text)
    except ValueError:
        raise ConvertError(f"{where}: {what} is not a number: {text!r}") from None
    if not math.isfinite(v):
        raise ConvertError(f"{where}: {what} is not finite: {text!r}")
    return v


def _read_lines(path, label):
    p = Path(path)
    if not p.is_file():
        raise ConvertError(f"{label}: file {path} does not exist")
    text = p.read_text(encoding="utf-8-sig")
    if not text.strip():
        raise ConvertError(f"{label}: file {path} is empty")
    return text.splitlines()


def read_manifest(path):
    ids = [ln.strip() for ln in _read_lines(path, "manifest") if ln.strip()]
    if len(set(ids)) != len(ids):
        raise ConvertError("manifest: duplicate accessions")
    return ids


def read_lengths(path):
    rows = list(csv.reader(_read_lines(path, "lengths"), delimiter="\t"))
    if rows[0][:2] != ["accession", "length"]:
        raise ConvertError("lengths: header must start with accession<TAB>length")
    out = {}
    for r in rows[1:]:
        if r[0] in out:
            raise ConvertError(f"lengths: duplicate accession {r[0]}")
        out[r[0]] = int(_num(r[1], f"lengths {r[0]}", "length"))
    return out


def parse_signalp(path):
    out, organism = {}, False
    for n, line in enumerate(_read_lines(path, "signalp"), 1):
        where = f"signalp:{path}:{n}"
        if line.startswith("#"):
            if "Organism:" in line:
                if line.split("Organism:", 1)[1].split()[:1] != ["Eukarya"]:
                    raise ConvertError(f"{where}: not a Eukarya run: {line.strip()!r}")
                organism = True
            continue
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 5:
            raise ConvertError(f"{where}: expected 5 TAB separated fields, got {len(f)}")
        pid = f[0].split()[0] if f[0].strip() else ""
        if not pid:
            raise ConvertError(f"{where}: empty ID")
        if pid in out:
            raise ConvertError(f"{where}: duplicate ID {pid}")
        pred = f[1].strip()
        p_other, p_sp = _num(f[2], where, "OTHER probability"), _num(f[3], where, "SP probability")
        cs = f[4].strip()
        if pred == "SP":
            m = CS_RE.match(cs)
            if not m:
                raise ConvertError(
                    f"{where}: prediction SP without a parsable 'CS pos: a-b. Pr: p' field: {cs!r}"
                )
            cpos, cprob = int(m.group(1)), _num(m.group(3), where, "cleavage probability")
        elif pred == "OTHER":
            if cs:
                raise ConvertError(
                    f"{where}: prediction OTHER but a cleavage field is present: {cs!r}"
                )
            cpos, cprob = "NA", "NA"
        else:
            raise ConvertError(
                f"{where}: unexpected prediction {pred!r} in a Eukarya run (expected SP or OTHER)"
            )
        out[pid] = {"pred": "yes" if pred == "SP" else "no", "cpos": cpos, "cprob": cprob,
                    "p_sp": p_sp, "p_other": p_other, "raw": pred}  # fmt: skip
    if not out:
        raise ConvertError(f"signalp: {path}: no prediction rows")
    if not organism:
        raise ConvertError(f"signalp: {path}: header lacks 'Organism: Eukarya'")
    return out


def parse_predgpi(path, lengths):
    """Read the raw PredGPI gff3 (one row per protein) against the protein lengths."""
    out = {}
    for n, line in enumerate(_read_lines(path, "predgpi"), 1):
        where = f"predgpi:{path}:{n}"
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) != 9:
            raise ConvertError(f"{where}: expected 9 TAB separated gff3 columns, got {len(f)}")
        pid, source, kind, start, end, score = f[:6]
        if not pid.strip() or pid != pid.strip():
            raise ConvertError(f"{where}: empty or padded ID {pid!r}")
        if pid in out:
            raise ConvertError(f"{where}: duplicate ID {pid}")
        if source != "PredGPI":
            raise ConvertError(f"{where}: source column is {source!r}, expected PredGPI")
        if pid not in lengths:
            raise ConvertError(f"{where}: {pid} is not in the lengths table")
        try:
            start, end = int(start), int(end)
        except ValueError:
            raise ConvertError(f"{where}: start/end are not integers: {f[3]!r} {f[4]!r}") from None
        short = lengths[pid] < PREDGPI_MIN_LEN
        if kind == "Chain":
            if start != 1 or end != lengths[pid]:
                raise ConvertError(
                    f"{where}: Chain row {start}-{end} for {pid} does not span the protein (length {lengths[pid]})"
                )
            out[pid] = {"pred": "NA" if short else "no", "omega": "NA", "cls": "NA", "ub": "NA",
                        "note": SHORT_NOTE if short else ""}  # fmt: skip
        elif kind == "GPI-anchor":
            if start != end:
                raise ConvertError(
                    f"{where}: GPI-anchor row for {pid} has start {start} != end {end}"
                )
            if not 1 <= start <= lengths[pid]:
                raise ConvertError(
                    f"{where}: omega site {start} outside the protein (length {lengths[pid]})"
                )
            if short:
                raise ConvertError(
                    f"{where}: {pid} is {lengths[pid]} aa; PredGPI should not score it, but a GPI-anchor row exists"
                )
            if score not in GPI_CLASSES:
                raise ConvertError(
                    f"{where}: unknown GPI-anchor confidence score {score!r} (known: 1.0, 0.70, 0.55)"
                )
            cls, ub = GPI_CLASSES[score]
            out[pid] = {"pred": "yes", "omega": start, "cls": cls, "ub": ub, "note": ""}
        else:
            raise ConvertError(
                f"{where}: unknown feature type {kind!r} (expected Chain or GPI-anchor)"
            )
    if not out:
        raise ConvertError(f"predgpi: {path}: no rows")
    return out


# The real tools overshoot [0, 1] slightly (850 proteins: AIUPred up to 1.0036, IUPred3 from -0.0123 to 1.0096; 27,811 proteins:
# AIUPred 273 of 13.76 million rows up to 1.023, 2026-10-24/25); more than this is still an error. Values are used as written.
SCORE_SLACK = 0.05


def parse_disorder(path, threshold):
    out, cur = {}, None
    for n, line in enumerate(_read_lines(path, "disorder"), 1):
        where = f"disorder:{path}:{n}"
        # AIUPred 2.x echoes each record header as '#>ID' (verified on the real output of 2026-10-24); its other '#' lines are
        # the banner and the column header
        if line.startswith("#>"):
            line = line[1:]
        if line.startswith("#") or not line.strip():
            continue
        if line.startswith(">"):
            toks = line[1:].split()
            if not toks:
                raise ConvertError(f"{where}: empty record name")
            if toks[0] in out:
                raise ConvertError(f"{where}: duplicate record {toks[0]}")
            cur = out[toks[0]] = []
            continue
        if cur is None:
            raise ConvertError(f"{where}: a score row before any '>' record line")
        f = line.split()
        if len(f) != 3:
            raise ConvertError(f"{where}: expected 'pos aa score', got {line.strip()!r}")
        pos = int(_num(f[0], where, "position"))
        if pos != len(cur) + 1:
            raise ConvertError(f"{where}: position {pos} out of order (expected {len(cur) + 1})")
        s = _num(f[2], where, "score")
        if not -SCORE_SLACK <= s <= 1.0 + SCORE_SLACK:
            raise ConvertError(f"{where}: score {s} outside [0, 1]")
        cur.append(s)
    if not out:
        raise ConvertError(f"disorder: {path}: no records")
    summary = {}
    for name, sc in out.items():
        if not sc:
            raise ConvertError(f"disorder: record {name} has no residues")
        longest = run = 0
        for s in sc:
            run = run + 1 if s >= threshold else 0
            longest = max(longest, run)
        summary[name] = {"n": len(sc), "mean": sum(sc) / len(sc),
                         "frac": sum(s >= threshold for s in sc) / len(sc), "longest": longest}  # fmt: skip
    return summary


def check_ids(label, have, want):
    missing = [a for a in want if a not in have]
    extra = sorted(set(have) - set(want))
    if missing:
        raise ConvertError(
            f"{label}: {len(missing)} manifest accessions missing, first: {missing[:5]}"
        )
    if extra:
        raise ConvertError(
            f"{label}: {len(extra)} accessions not in the manifest, first: {extra[:5]}"
        )


def convert(manifest, signalp, predgpi, disorder, tool, threshold, lengths):
    """`disorder` None means the disorder block is not requested (its columns are then absent)."""
    want = read_manifest(manifest)
    ln = read_lengths(lengths)
    check_ids("lengths", ln, want)
    sp, gp = parse_signalp(signalp), parse_predgpi(predgpi, ln)
    checks = [("signalp", sp), ("predgpi", gp)]
    dis = None
    if disorder is not None:
        dis = parse_disorder(disorder, threshold)
        checks.append(("disorder", dis))
    for label, have in checks:
        check_ids(label, have, want)
    if dis is not None:
        bad = [a for a in want if dis[a]["n"] != ln[a]]
        if bad:
            raise ConvertError(f"disorder: {len(bad)} records differ in length from lengths.tsv, first: "
                               f"{bad[0]} has {dis[bad[0]]['n']} residues, expected {ln[bad[0]]}")  # fmt: skip
    rows = []
    for a in want:
        s, g = sp[a], gp[a]
        row = {
            "accession": a, "signal_peptide_pred": s["pred"], "sp_cleavage_pos": s["cpos"],
            "sp_cleavage_prob": s["cprob"], "sp_prob": s["p_sp"], "sp_prediction_raw": s["raw"],
            "gpi_pred": g["pred"], "gpi_omega_site": g["omega"], "gpi_confidence_class": g["cls"],
            "gpi_fpr_upper_bound": g["ub"], "gpi_note": g["note"],
        }  # fmt: skip
        if dis is not None:
            d = dis[a]
            row.update({
                "disorder_tool": tool, "disorder_threshold": threshold, "disorder_n_res": d["n"],
                "disorder_mean": f"{d['mean']:.4f}", "disorder_frac": f"{d['frac']:.4f}",
                "disorder_longest_run": d["longest"],
            })  # fmt: skip
        rows.append(row)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--signalp", required=True)
    ap.add_argument("--predgpi", required=True)
    ap.add_argument("--disorder", help="per-residue disorder file (needs --disorder-tool)")
    ap.add_argument("--disorder-tool", choices=["aiupred", "iupred3"])
    ap.add_argument(
        "--no-disorder",
        action="store_true",
        help="explicitly leave the disorder block out (no disorder columns)",
    )
    ap.add_argument("--disorder-threshold", type=float, default=0.5)
    ap.add_argument(
        "--lengths",
        required=True,
        help="lengths.tsv of stage_inputs.py (PredGPI checks, residue counts)",
    )
    ap.add_argument("--out", required=True, help="output TSV; must not exist")
    a = ap.parse_args(argv)
    if a.no_disorder == bool(a.disorder):
        ap.error("give exactly one of --disorder FILE (with --disorder-tool) or --no-disorder")
    if a.disorder and not a.disorder_tool:
        ap.error("--disorder needs --disorder-tool")
    if a.no_disorder and a.disorder_tool:
        ap.error("--no-disorder and --disorder-tool contradict each other")
    if Path(a.out).exists():
        print(f"FATAL: {a.out} exists; nothing is overwritten", file=sys.stderr)
        return 1
    try:
        rows = convert(
            a.manifest,
            a.signalp,
            a.predgpi,
            None if a.no_disorder else a.disorder,
            a.disorder_tool,
            a.disorder_threshold,
            a.lengths,
        )
    except ConvertError as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        cols = [c for c in OUT_COLS if c in rows[0]]
        w.writerow(cols)
        w.writerows([[r[c] for c in cols] for r in rows])
    print(f"wrote {a.out}: {len(rows)} proteins")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
