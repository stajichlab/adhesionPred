#!/usr/bin/env python
"""Class 2b: what is PRA/Ag2's actual fold, and who else has it?

The "fold" worth searching for is NOT the whole 194 aa protein. AlphaFold DB already has a
model for PRA (UniProt Q6QJA6, AF-Q6QJA6-F1, mean pLDDT 62 -- moderate/low overall), and its
per-residue confidence lines up exactly with the domain architecture reported in Zhu et al.
1996 (Gene 181:121-125): the N- and C-terminal signal peptides and the Pro/Thr tetrapeptide-
repeat region (aa 89-141) are all low-confidence (pLDDT ~51-52, consistent with disorder), while
the N-terminal region (aa 20-84) is confidently folded (pLDDT 83.7) -- and InterPro/Pfam
independently call that same region a CFEM domain (PF05730, e=1.5e-13).

So class 2b's real question is: which other Coccidioides/Onygenales proteins carry a CFEM
domain? This is answered from the Fungi_5k functionalDB's precomputed Pfam scan (no new HMMER
or structure-prediction run needed) rather than from a Foldseek structural search, which is not
available on this HPCC node (no foldseek install; see NOTE at the bottom).

Run on HPCC (needs function.duckdb, so srun, not the login node):
    srun -p epyc -c 2 --mem 10G -t 15 python3 pra_cfem_survey.py --out cfem_onygenales.tsv
"""

import argparse
import sys
from pathlib import Path

import duckdb
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from analysis._common.paths import FUNGI5K_DUCKDB, FUNGI5K_SAMPLES  # noqa: E402

CFEM_PFAM = "PF05730"


def load_onygenales_prefixes():
    df = pd.read_csv(FUNGI5K_SAMPLES)
    ony = df[df["ORDER"] == "Onygenales"]
    return dict(zip(ony["LOCUSTAG"], ony["SPECIESIN"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="output TSV of CFEM hits")
    args = ap.parse_args()

    prefixes = load_onygenales_prefixes()
    con = duckdb.connect(str(FUNGI5K_DUCKDB), read_only=True)
    con.execute("PRAGMA memory_limit='6GB'")

    placeholders = ",".join(f"'{p}'" for p in prefixes)
    pfam = con.execute(
        f"""
        SELECT protein_id, species_prefix, pfam_acc, full_seq_e_value, domain_i_evalue,
               hmm_from, hmm_to, ali_from, ali_to
        FROM pfam
        WHERE pfam_acc ILIKE '{CFEM_PFAM}%' AND species_prefix IN ({placeholders})
        ORDER BY species_prefix, domain_i_evalue
        """
    ).fetchdf()
    pfam["species"] = pfam["species_prefix"].map(prefixes)

    ids = tuple(pfam["protein_id"])
    id_placeholders = ",".join(f"'{i}'" for i in ids)
    sp = con.execute(
        f"SELECT protein_id, probability AS signalp_prob FROM signalp WHERE protein_id IN ({id_placeholders})"
    ).fetchdf()
    pfam = pfam.merge(sp, on="protein_id", how="left")
    pfam["secreted"] = pfam["signalp_prob"].apply(
        lambda p: "yes" if pd.notna(p) and p > 0.5 else "no_signalp_call"
    )

    pfam.to_csv(args.out, sep="\t", index=False)
    print(
        f"{len(pfam)} CFEM (PF05730) hits across {pfam['protein_id'].nunique()} proteins, "
        f"{pfam['species'].nunique()} Onygenales genomes -> {args.out}"
    )
    print(pfam.groupby("species").size().sort_values(ascending=False).head(10).to_string())
    print()
    cocci = pfam[pfam["species"].str.contains("Coccidioides", na=False)]
    print(f"Coccidioides: {len(cocci)} CFEM proteins")
    print(
        cocci[
            ["protein_id", "species", "domain_i_evalue", "ali_from", "ali_to", "secreted"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
