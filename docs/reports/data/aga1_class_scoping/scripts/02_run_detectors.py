#!/usr/bin/env python3
"""Run the two existing repeat detectors (analysis/cocci_repeats 02 and 14, the ones behind `tandem_repeat_protein`)
locally on every protein of proteins.tsv. Local CPU only, no HPCC. The FASTA and raw detector outputs stay in the
cache OUTSIDE the repo; 03_features.py reads them and writes numbers only.

Needs numpy + biopython (detector 14): run with the repo environment python, e.g.
  /Users/jstajich/projects/adhesionPred/.venv/bin/python 02_run_detectors.py
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import CACHE, OUT, ROOT, all_sequences, read_tsv  # noqa: E402

DETECTORS = {"repeat02": "02_repeat_profile.py", "repeat14": "14_repeat_detect_general.py"}


def main():
    seqs, _ = all_sequences()
    ids = [r["accession"] for r in read_tsv(OUT / "proteins.tsv")]
    work = CACHE / "detectors"
    work.mkdir(parents=True, exist_ok=True)
    fa = work / "all.faa"
    with open(fa, "w") as fh:
        for a in ids:
            fh.write(f">{a}\n{seqs[a]}\n")
    for name, script in DETECTORS.items():
        subprocess.run(
            [sys.executable, str(ROOT / "analysis/cocci_repeats" / script), str(fa), "--out", str(work / f"{name}.tsv"), "--min-len", "80"],
            check=True,
        )  # fmt: skip
        print("done", name)


if __name__ == "__main__":
    main()
