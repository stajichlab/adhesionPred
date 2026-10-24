"""Raw SignalP 6 / PredGPI / disorder outputs -> one table per protein for the AGA1-class prediction step (stdlib only).

FORMATS ARE UNVERIFIED until the first real output is seen (see analysis/aga1_class/README.md):
  SignalP 6 `prediction_results.txt` (format as parsed by src/cellsurface_sorting_hat/modules/signalp.py, plus the cleavage column):
      '# SignalP-6.0 ... Organism: Eukarya ...' / '# ID<TAB>Prediction<TAB>OTHER<TAB>SP(Sec/SPI)<TAB>CS Position'
      then 'ID<TAB>SP<TAB>p_other<TAB>p_sp<TAB>CS pos: 22-23. Pr: 0.9876' (the CS field is empty for OTHER).
  PredGPI: NORMALISED intermediate TSV (the raw PredGPI output format is not known to the author; adapt the raw output to this
      once seen): header `id<TAB>prediction<TAB>omega_site<TAB>fp_rate[<TAB>score]`, prediction `GPI-anchored` or
      `Not GPI-anchored`, omega_site and score may be empty or NA.
  Disorder (AIUPred or IUPred3): per-residue text, one record per protein: a line `>ID ...` followed by rows `pos aa score`
      (whitespace separated, pos 1..n); `#` lines are ignored. If the raw tool prints another layout, convert it first.
Fails loudly (non-zero exit, message on stderr) on a missing accession, a duplicate, an empty file, an unparseable row, an
unexpected accession, a length mismatch (with --lengths). Nothing is ever written as a silent zero: a value that does not
apply is written as NA.
"""

import argparse
import csv
import math
import re
import sys
from pathlib import Path

CS_RE = re.compile(r"^CS pos:\s*(\d+)-(\d+)\.\s*Pr:\s*([0-9.eE+-]+)\s*$")
OUT_COLS = [
    "accession", "signal_peptide_pred", "sp_cleavage_pos", "sp_cleavage_prob", "sp_prob", "sp_prediction_raw",
    "gpi_pred", "gpi_omega_site", "gpi_fp_rate", "gpi_score",
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


def parse_predgpi(path):
    lines = [ln for ln in _read_lines(path, "predgpi") if ln.strip() and not ln.startswith("#")]
    if not lines:
        raise ConvertError(f"predgpi: {path}: no rows")
    head = lines[0].split("\t")
    need = ["id", "prediction", "omega_site", "fp_rate"]
    if head[:4] != need:
        raise ConvertError(f"predgpi: header must start with {need}, got {head}")
    out = {}
    for n, line in enumerate(lines[1:], 2):
        where = f"predgpi:{path}:row{n}"
        f = line.split("\t")
        if len(f) < 4:
            raise ConvertError(f"{where}: expected at least 4 TAB separated fields")
        pid = f[0].strip()
        if not pid:
            raise ConvertError(f"{where}: empty id")
        if pid in out:
            raise ConvertError(f"{where}: duplicate id {pid}")
        pred = {"GPI-anchored": "yes", "Not GPI-anchored": "no"}.get(f[1].strip())
        if pred is None:
            raise ConvertError(f"{where}: unexpected prediction {f[1]!r}")
        omega = f[2].strip()
        omega = "NA" if omega in ("", "NA") else int(_num(omega, where, "omega site"))
        fp = _num(f[3], where, "FP rate")
        if fp < 0:
            raise ConvertError(f"{where}: negative FP rate")
        score = f[4].strip() if len(f) > 4 else ""
        score = "NA" if score in ("", "NA") else _num(score, where, "score")
        out[pid] = {"pred": pred, "omega": omega, "fp": fp, "score": score}
    return out


def parse_disorder(path, threshold):
    out, cur = {}, None
    for n, line in enumerate(_read_lines(path, "disorder"), 1):
        where = f"disorder:{path}:{n}"
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
        if not 0.0 <= s <= 1.0:
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


def convert(manifest, signalp, predgpi, disorder, tool, threshold, lengths=None):
    want = read_manifest(manifest)
    sp, gp, dis = (
        parse_signalp(signalp),
        parse_predgpi(predgpi),
        parse_disorder(disorder, threshold),
    )
    for label, have in (("signalp", sp), ("predgpi", gp), ("disorder", dis)):
        check_ids(label, have, want)
    if lengths:
        ln = read_lengths(lengths)
        check_ids("lengths", ln, want)
        bad = [a for a in want if dis[a]["n"] != ln[a]]
        if bad:
            raise ConvertError(f"disorder: {len(bad)} records differ in length from lengths.tsv, first: "
                               f"{bad[0]} has {dis[bad[0]]['n']} residues, expected {ln[bad[0]]}")  # fmt: skip
    rows = []
    for a in want:
        s, g, d = sp[a], gp[a], dis[a]
        rows.append({
            "accession": a, "signal_peptide_pred": s["pred"], "sp_cleavage_pos": s["cpos"],
            "sp_cleavage_prob": s["cprob"], "sp_prob": s["p_sp"], "sp_prediction_raw": s["raw"],
            "gpi_pred": g["pred"], "gpi_omega_site": g["omega"], "gpi_fp_rate": g["fp"], "gpi_score": g["score"],
            "disorder_tool": tool, "disorder_threshold": threshold, "disorder_n_res": d["n"],
            "disorder_mean": f"{d['mean']:.4f}", "disorder_frac": f"{d['frac']:.4f}", "disorder_longest_run": d["longest"],
        })  # fmt: skip
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--signalp", required=True)
    ap.add_argument("--predgpi", required=True)
    ap.add_argument("--disorder", required=True)
    ap.add_argument("--disorder-tool", required=True, choices=["aiupred", "iupred3"])
    ap.add_argument("--disorder-threshold", type=float, default=0.5)
    ap.add_argument("--lengths", help="lengths.tsv of stage_inputs.py (checks the residue counts)")
    ap.add_argument("--out", required=True, help="output TSV; must not exist")
    a = ap.parse_args(argv)
    if Path(a.out).exists():
        print(f"FATAL: {a.out} exists; nothing is overwritten", file=sys.stderr)
        return 1
    try:
        rows = convert(
            a.manifest,
            a.signalp,
            a.predgpi,
            a.disorder,
            a.disorder_tool,
            a.disorder_threshold,
            a.lengths,
        )
    except ConvertError as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(OUT_COLS)
        w.writerows([[r[c] for c in OUT_COLS] for r in rows])
    print(f"wrote {a.out}: {len(rows)} proteins")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
