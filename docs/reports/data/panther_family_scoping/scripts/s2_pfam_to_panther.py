"""Step 2: which PANTHER families do the proteins of each Pfam family of the panel fall in (UniProt REST).

Usage: python3 -I s2_pfam_to_panther.py REPO_ROOT OUT_TSV
Pfam families: the 30 candidates of the specificity test (docs/reports/data/pr6_family_candidates_2026-10-22/candidate_families.tsv),
the active rows of data/sorting_hat/family_table.tsv and the 21 models of analysis/wall_fusion_report/wall_fusion_models.tsv.
Per Pfam: the fungal protein count (x-total-results), then a sample = all reviewed fungal proteins (up to 2000) plus the
first 500 unreviewed ones returned by UniProt; PANTHER families (the base id PTHRxxxxx) are tallied over the sample.
"""

import collections
import csv
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import get  # noqa: E402

root, out = Path(sys.argv[1]), Path(sys.argv[2])
sets = collections.defaultdict(set)
names = {}


def rd(path, accol, setname, namecol=None, flt=None):
    with (root / path).open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if flt and not flt(r):
                continue
            sets[r[accol]].add(setname)
            if namecol:
                names[r[accol]] = r[namecol]


rd(
    "docs/reports/data/pr6_family_candidates_2026-10-22/candidate_families.tsv",
    "pfam_acc",
    "candidate30",
    "name",
)
rd(
    "data/sorting_hat/family_table.tsv",
    "pfam_acc",
    "family_table_active",
    "name",
    lambda r: r["active"] == "yes",
)
rd("analysis/wall_fusion_report/wall_fusion_models.tsv", "id", "wall_fusion_models", "name")

U = "https://rest.uniprot.org/uniprotkb/"
F = "accession,reviewed,organism_id,xref_pfam,xref_panther"


def q(pf, extra):
    return urllib.parse.quote(f"(xref:pfam-{pf}) AND (taxonomy_id:4751){extra}")


def parse(tsv):
    lines = tsv.strip().split("\n")[1:]
    res = []
    for ln in lines:
        c = ln.split("\t")
        if len(c) < 5:
            c += [""] * (5 - len(c))
        res.append((c[0], c[1], {x.split(":")[0] for x in c[4].split(";") if x}))
    return res


cols = [
    "pfam_acc",
    "pfam_name",
    "panel_sets",
    "n_fungal_pfam_uniprot",
    "n_reviewed_fungal",
    "sample_size",
    "panther_family",
    "n_in_sample",
    "frac_of_sample",
    "n_sample_without_panther",
]
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for pf in sorted(sets):
        _, tot = get(
            f"{U}search?query={q(pf, '')}&format=tsv&fields=accession&size=1", want_headers=True
        )
        _, nrev = get(
            f"{U}search?query={q(pf, ' AND (reviewed:true)')}&format=tsv&fields=accession&size=1",
            want_headers=True,
        )
        sample = []
        if int(nrev or 0):
            sample += parse(
                get(f"{U}stream?query={q(pf, ' AND (reviewed:true)')}&format=tsv&fields={F}")
            )
        sample += parse(
            get(f"{U}search?query={q(pf, ' AND (reviewed:false)')}&format=tsv&fields={F}&size=500")
        )
        cnt = collections.Counter()
        none = 0
        for _, _, fam in sample:
            if not fam:
                none += 1
            for x in fam:
                cnt[x] += 1
        n = len(sample)
        for fam, c in cnt.most_common(6):
            w.writerow(
                [
                    pf,
                    names.get(pf, ""),
                    "+".join(sorted(sets[pf])),
                    tot,
                    nrev,
                    n,
                    fam,
                    c,
                    f"{c / n:.3f}" if n else "",
                    none,
                ]
            )
        if not cnt:
            w.writerow(
                [pf, names.get(pf, ""), "+".join(sorted(sets[pf])), tot, nrev, n, "", 0, "0", none]
            )
        print(pf, tot, nrev, n, cnt.most_common(2), flush=True)
