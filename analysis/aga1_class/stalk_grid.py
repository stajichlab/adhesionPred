"""Sensitivity grid of the broad flag over the stalk-rule settings (stdlib only; same raw AIUPred file as stalk_features.py).

For every window x minimum Ser/Thr fraction x minimum run, counts per role_group the proteins that have a signal peptide, a
PredGPI GPI call and a stalk run of at least the minimum (the broad flag with the stalk rule and no tail rescue). The setting
`none` is the flag without the stalk rule. Reads the AIUPred record table for the signal peptide and GPI calls.

Usage: python3 -I stalk_grid.py --raw AIUPRED_RAW.txt --record RECORD.tsv --manifest manifest_accessions.txt --out GRID.tsv
Refuses an existing --out; FATAL: and exit 1 on a mismatch of accessions.
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stalk_features as sf  # noqa: E402

WINDOWS, FRACS, MIN_RUNS, THRESHOLD = (20, 30, 50), (0.4, 0.5, 0.6), (100, 200), 0.5


def grid(raw, record):
    parts = {r["accession"]: r for r in record}
    if set(parts) != set(raw):
        raise sf.StalkError("the record table and the raw file hold different accessions")
    groups = sorted({r["role_group"] for r in record})
    n = Counter(r["role_group"] for r in record)
    settings = [("none", "none")] + [(w, f) for w in WINDOWS for f in FRACS]
    out = []
    for n_min in MIN_RUNS:
        for w, f in settings:
            yes = Counter()
            for a, r in parts.items():
                if r["signal_peptide_pred"] != "yes" or r["gpi_pred"] != "yes":
                    continue
                run = (
                    int(r["disorder_longest_run"])
                    if w == "none"
                    else sf.stalk(*raw[a], w, f, THRESHOLD)[0]
                )
                yes[r["role_group"]] += run >= n_min
            for g in groups:
                out.append(
                    {
                        "min_run": n_min,
                        "stalk_window": w,
                        "stalk_min_frac": f,
                        "role_group": g,
                        "n": n[g],
                        "yes": yes[g],
                    }
                )
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("raw", "record", "manifest", "out"):
        ap.add_argument(f"--{k}", required=True)
    a = ap.parse_args(argv)
    if Path(a.out).exists():
        print(f"FATAL: {a.out} exists; nothing is overwritten", file=sys.stderr)
        return 1
    try:
        raw = sf.read_raw(a.raw)
        if set(raw) != set(Path(a.manifest).read_text().split()):
            raise sf.StalkError("the raw file differs from the manifest")
        with open(a.record, newline="") as fh:
            rows = grid(raw, list(csv.DictReader(fh, delimiter="\t")))
    except (sf.StalkError, OSError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    cols = ["min_run", "stalk_window", "stalk_min_frac", "role_group", "n", "yes"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {a.out}: {len(rows)} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
