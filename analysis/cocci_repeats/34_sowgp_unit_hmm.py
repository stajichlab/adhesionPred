#!/usr/bin/python3.12
"""Build the SOWgp repeat-unit profile HMM and pick the modal unit.

Input   sowgp_repeat_viz.units.tsv (from 09_sowgp_repeat_viz.py; gunzip the .gz if absent).
        One row per anchored unit of every full-length SOWgp copy in the pangenome.
Output  sowgp_unit_set.fa     the distinct unit sequences (one record each, with copy count)
        sowgp_unit.msa.fa     MAFFT alignment of them
        sowgp_unit.hmm        hmmbuild profile (amino, --wpb, default entropy weighting)
        sowgp_unit_modal.fa   the modal unit (most frequent over all copies)

Units are 47 aa from each PTDCYGDC anchor (sowgp_units.py). The terminal unit (class
`terminal`) is kept: its first ~34 aa are homologous to the others and the last 13 aa are
not, which the profile scores as a weak tail. Units of < 40 aa (partial, at a sequence end)
are dropped.

Needs mafft and hmmer/3.4 on PATH. Run inside `module load mafft hmmer/3.4`, e.g.
    srun -p short -c 2 --mem 4G -t 10 bash -lc 'module load mafft hmmer/3.4; \
        /usr/bin/python3.12 34_sowgp_unit_hmm.py'
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
MIN_UNIT = 40


def main():
    f = HERE / "sowgp_repeat_viz.units.tsv"
    u = pd.read_csv(f if f.exists() else str(f) + ".gz", sep="\t")
    u = u[u.unit_seq.str.len() >= MIN_UNIT]
    cnt = u.unit_seq.value_counts()
    modal = cnt.index[0]
    print(f"{len(u)} units from {u.protein.nunique()} proteins; {len(cnt)} distinct")
    print(f"modal unit ({cnt.iloc[0]} copies): {modal}")
    with open(HERE / "sowgp_unit_set.fa", "w") as fh:
        for i, (s, n) in enumerate(cnt.items()):
            fh.write(f">unit{i:02d}_n{n}\n{s}\n")
    (HERE / "sowgp_unit_modal.fa").write_text(f">SOWgp_modal_unit\n{modal}\n")
    aln = subprocess.run(
        [
            "mafft",
            "--localpair",
            "--maxiterate",
            "1000",
            "--quiet",
            str(HERE / "sowgp_unit_set.fa"),
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    (HERE / "sowgp_unit.msa.fa").write_text(aln)
    subprocess.run(
        [
            "hmmbuild",
            "--amino",
            "-n",
            "SOWgp_unit",
            str(HERE / "sowgp_unit.hmm"),
            str(HERE / "sowgp_unit.msa.fa"),
        ],
        check=True,
    )
    print("wrote sowgp_unit_set.fa, sowgp_unit.msa.fa, sowgp_unit.hmm, sowgp_unit_modal.fa")


if __name__ == "__main__":
    sys.exit(main())
