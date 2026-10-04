#!/usr/bin/env python3.12
"""Calibrate the significance threshold in 14_repeat_detect_general.py against its true null.

Written for the 2026-09-30 review of REPORT_2026-09-29_repeat_detector_divergence.md.

14's rule is "the best period must reach z >= Z_MIN". The per-period statistic is well
calibrated -- but the threshold is applied to the ARGMAX over every period scanned, so the
relevant null is the maximum of ~77 correlated draws, not one N(0,1) draw. This script
measures both, on random sequence where no repeat exists by construction.

Three questions:
  1. Is z flat across the period range? The exclusion rule in period_z() removes far more
     periods when p is small (25% of periods survive at p=4 vs 96% at p=47), which could
     inflate z at short periods.
  2. What is the per-period false-positive rate at Z_MIN?
  3. What is the PER-PROTEIN rate, i.e. how often does max(z) clear Z_MIN by chance?

Run:  /usr/bin/python3.12 19_zthreshold_probe.py            (from this directory)
      /usr/bin/python3.12 19_zthreshold_probe.py --trials 2000 --length 400

Takes about 90 s at the defaults. Results are printed, not written -- this is a calibration
probe, not a pipeline stage.

Caveat: residues are drawn uniformly from the 20 amino acids. Real proteome composition is
more skewed and therefore a harder null, so the rates here are floors.
"""

import argparse
import importlib.util
import random
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROBE_PERIODS = (4, 5, 6, 8, 10, 12, 14, 17, 20, 24, 34, 47, 60, 75)


def load_detector():
    """Import 14_repeat_detect_general.py, whose name is not a valid module name."""
    path = HERE / "14_repeat_detect_general.py"
    spec = importlib.util.spec_from_file_location("repeat_detect_general", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def null_set_sizes(m):
    """Fraction of periods left in the comparison set, as a function of p."""
    periods = list(range(m.MIN_PERIOD, m.MAX_PERIOD + 1))
    out = {}
    for p in periods:
        rel = [q for q in periods if not (q % p in (0, 1) or (q + 1) % p == 0 or p % q == 0)]
        out[p] = len(rel) / len(periods)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=600)
    ap.add_argument("--length", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    m = load_detector()
    random.seed(args.seed)
    np.random.seed(args.seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    print(f"detector: MIN_PERIOD={m.MIN_PERIOD} MAX_PERIOD={m.MAX_PERIOD} Z_MIN={m.Z_MIN}")
    print(f"null: {args.trials} random sequences of {args.length} aa, uniform over 20 aa\n")

    print("1. Null-set size by period (fraction of periods left after excluding")
    print("   multiples and divisors of p):")
    frac = null_set_sizes(m)
    print("   " + "  ".join(f"p={p}:{frac[p]:.2f}" for p in PROBE_PERIODS if p in frac))
    print()

    by_period, max_z = {}, []
    for _ in range(args.trials):
        s = "".join(random.choice(aa) for _ in range(args.length))
        raw = m.raw_scan(m.encode(s), m.SIM)
        zs = []
        for p in raw:
            z = m.period_z(raw, p)
            if z is None:
                continue
            by_period.setdefault(p, []).append(z)
            zs.append(z)
        if zs:
            max_z.append(max(zs))

    print("2. Per-period z on random sequence (ideal: mean 0, sd 1, flat in p):")
    print(f"   {'p':>4} {'mean':>8} {'sd':>7} {'frac>=Z_MIN':>12}")
    for p in PROBE_PERIODS:
        if p not in by_period:
            continue
        a = np.array(by_period[p])
        print(f"   {p:>4} {a.mean():>8.3f} {a.std():>7.3f} {(a >= m.Z_MIN).mean():>12.5f}")
    allz = np.concatenate([np.array(v) for v in by_period.values()])
    print(f"   {'all':>4} {allz.mean():>8.3f} {allz.std():>7.3f} {(allz >= m.Z_MIN).mean():>12.5f}")
    print()

    mz = np.array(max_z)
    n_per = len(by_period)
    print(f"3. PER-PROTEIN null: max z over the ~{n_per} periods scanned.")
    print("   This is the distribution Z_MIN is actually compared against.")
    print(
        f"   mean {mz.mean():.2f}  sd {mz.std():.2f}  median {np.median(mz):.2f}  "
        f"95th {np.percentile(mz, 95):.2f}  max {mz.max():.2f}"
    )
    for t in (3.0, 3.5, 4.0, 4.5, 5.0):
        flag = "  <- current Z_MIN" if abs(t - m.Z_MIN) < 1e-9 else ""
        print(f"     frac of random proteins with max z >= {t}: {(mz >= t).mean():.4f}{flag}")
    print()
    print("   A per-period rate near 0 does NOT mean a per-protein rate near 0: the threshold")
    print("   is applied to the best of many periods. Calibrate Z_MIN against section 3.")


if __name__ == "__main__":
    main()
