"""Ser/Thr stalk and C-terminal tail features per protein from a per-residue AIUPred file (stdlib only).

Input is the raw AIUPred output of the AGA1-class prediction step (`#>ID` record lines, `pos<TAB>aa<TAB>score` rows). That file
carries the residue letters, so it is NOT in the repository; this script writes no residues, only per-protein numbers.

Stalk rule (owner decision 2026-10-24): a disordered run is a maximal stretch of residues with an AIUPred score at or above
--threshold (0.5). A run holds a stalk if some window of --window residues (30) lying wholly inside it is at least --min-st-frac
(0.5) Ser/Thr. `stalk_run_len` is the length of the longest run that holds a stalk (0 if none); `max_st_window_frac` is the
highest Ser/Thr fraction of any such window over all runs of at least --window residues (0 if there is none).

C-terminal tail (exploratory, for proteins PredGPI does not call): `ctail_hydro` is the highest mean Kyte-Doolittle value of a
15-residue window inside the last 40 residues, `ctail_after` the number of residues after that window; both NA for a protein
shorter than 15 residues.

Usage: python3 -I stalk_features.py --raw AIUPRED_RAW.txt --manifest manifest_accessions.txt --out FEATURES.tsv
         [--window 30] [--min-st-frac 0.5] [--threshold 0.5]
Refuses an existing --out. Fails loudly (FATAL:, exit 1) on an accession set that differs from the manifest, an out-of-order
position, a bad score row, a duplicate record or a parameter out of range.
"""

import argparse
import csv
import sys
from pathlib import Path

KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}  # fmt: skip
TAIL_SPAN, TAIL_WIN = 40, 15
OUT_COLS = [
    "accession", "n_res", "stalk_window", "stalk_min_frac", "stalk_threshold", "stalk_run_len", "max_st_window_frac",
    "ctail_hydro", "ctail_after",
]  # fmt: skip


class StalkError(Exception):
    pass


def read_raw(path):
    """accession -> (residues, scores)"""
    out, cur = {}, None
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if line.startswith("#>"):
            name = line[2:].split()[0] if line[2:].split() else ""
            if not name or name in out:
                raise StalkError(f"{path}:{n}: empty or duplicate record {name!r}")
            cur = out[name] = ([], [])
            continue
        if line.startswith("#") or not line.strip():
            continue
        if cur is None:
            raise StalkError(f"{path}:{n}: a score row before any record line")
        f = line.split()
        if len(f) != 3 or len(f[1]) != 1:
            raise StalkError(f"{path}:{n}: expected 'pos aa score', got {line!r}")
        try:
            pos, score = int(f[0]), float(f[2])
        except ValueError as e:
            raise StalkError(f"{path}:{n}: {e}") from None
        if pos != len(cur[0]) + 1:
            raise StalkError(
                f"{path}:{n}: position {pos} out of order (expected {len(cur[0]) + 1})"
            )
        cur[0].append(f[1])
        cur[1].append(score)
    if not out or any(not r[0] for r in out.values()):
        raise StalkError(f"{path}: no records, or a record without residues")
    return out


def runs(scores, threshold):
    i, n = 0, len(scores)
    while i < n:
        if scores[i] >= threshold:
            j = i
            while j + 1 < n and scores[j + 1] >= threshold:
                j += 1
            yield i, j
            i = j + 1
        else:
            i += 1


def stalk(residues, scores, window, min_frac, threshold):
    """(longest run holding a stalk window, best Ser/Thr fraction of any window inside a run of >= window)"""
    st = [1 if a in "ST" else 0 for a in residues]
    longest, best = 0, 0.0
    for i, j in runs(scores, threshold):
        if j - i + 1 < window:
            continue
        c = sum(st[i : i + window])
        top = c
        for k in range(i + 1, j - window + 2):
            c += st[k + window - 1] - st[k - 1]
            top = max(top, c)
        best = max(best, top / window)
        if top / window >= min_frac:
            longest = max(longest, j - i + 1)
    return longest, best


def ctail(residues):
    t = residues[-TAIL_SPAN:]
    if len(t) < TAIL_WIN:
        return "NA", "NA"
    top, at = None, 0
    for i in range(len(t) - TAIL_WIN + 1):
        m = sum(KD.get(c, 0.0) for c in t[i : i + TAIL_WIN]) / TAIL_WIN
        if top is None or m > top:
            top, at = m, i
    return f"{top:.3f}", len(t) - (at + TAIL_WIN)


def features(raw, manifest, window, min_frac, threshold):
    if window < 1 or not 0 < min_frac <= 1 or not 0 <= threshold <= 1:
        raise StalkError(
            "parameters out of range (window >= 1, 0 < min-st-frac <= 1, 0 <= threshold <= 1)"
        )
    want = Path(manifest).read_text().split()
    if set(raw) != set(want) or len(want) != len(set(want)):
        raise StalkError(
            f"the accessions of the raw file differ from the manifest ({len(raw)} vs {len(want)})"
        )
    rows = []
    for a in want:
        res, sc = raw[a]
        run, frac = stalk(res, sc, window, min_frac, threshold)
        h, after = ctail(res)
        rows.append({
            "accession": a, "n_res": len(res), "stalk_window": window, "stalk_min_frac": min_frac,
            "stalk_threshold": threshold, "stalk_run_len": run, "max_st_window_frac": f"{frac:.3f}",
            "ctail_hydro": h, "ctail_after": after,
        })  # fmt: skip
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--raw", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=30)
    ap.add_argument("--min-st-frac", type=float, default=0.5)
    ap.add_argument("--threshold", type=float, default=0.5)
    a = ap.parse_args(argv)
    if Path(a.out).exists():
        print(f"FATAL: {a.out} exists; nothing is overwritten", file=sys.stderr)
        return 1
    try:
        rows = features(read_raw(a.raw), a.manifest, a.window, a.min_st_frac, a.threshold)
    except (StalkError, OSError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(
        f"wrote {a.out}: {len(rows)} proteins, window {a.window}, min Ser/Thr fraction {a.min_st_frac}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
