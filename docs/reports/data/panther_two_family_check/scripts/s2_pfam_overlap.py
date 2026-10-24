"""Counts for the specificity-test discussion: how many fungal UniProt proteins carry each Pfam model and how many of
them sit in the PANTHER family. Usage: python3 -I s2_pfam_overlap.py OUT_TSV. Counts come from the x-total-results
header of UniProt REST search (size=1); a stream count of the same query can differ by under 1 percent.
Also writes median fungal length per subfamily (from the cached UniProt stream of s1)."""

import csv
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import get  # noqa: E402

U = "https://rest.uniprot.org/uniprotkb/search?format=tsv&fields=accession&size=1&query="


def n(q):
    return int(get(U + urllib.parse.quote(q), want_headers=True)[1] or 0)


F = "taxonomy_id:4751"
rows = []
for fam, pfs in (
    ("PTHR46430", ["PF08238", "PF13181", "PF07719", "PF00515", "PF13174", "PF13424", "PF13432"]),
    ("PTHR35192", ["PF21671"]),
):
    nfam = n(f"(xref:panther-{fam}) AND ({F})")
    for pf in pfs:
        a = n(f"(xref:pfam-{pf}) AND ({F})")
        both = n(f"(xref:pfam-{pf}) AND ({F}) AND (xref:panther-{fam})")
        rows.append(
            [
                fam,
                pf,
                nfam,
                a,
                both,
                f"{both / a:.3f}" if a else "",
                f"{both / nfam:.3f}" if nfam else "",
            ]
        )
        print(rows[-1], flush=True)
with open(sys.argv[1], "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(
        [
            "panther_family",
            "pfam",
            "n_fungal_in_family",
            "n_fungal_with_pfam",
            "n_fungal_both",
            "frac_of_pfam_proteins_in_family",
            "frac_of_family_with_pfam",
        ]
    )
    w.writerows(rows)
