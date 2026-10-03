#!/usr/bin/env python3.12
"""Compare SOWgp unit counts and gene models at the 5 UArizona long-read loci.

For each assembly and each SOWgp locus (scaffold):
  - genomic units: number of tblastn first-part hits (query aa 1-31) of one RS unit, from
    44_longread_unit_tblastn.sh. Each unit has two hits (aa 1-31, aa 32-47) split by a ~56 nt intron.
  - funannotate models overlapping the array (GFF3 in the UArizona For_Marc folder) and their
    anchor counts (PTDCYGDC) from sowgp_anchored_longread.tsv.
  - miniprot model (43_longread_sowgp_miniprot.sh): best seed hit by identity over the locus.
Output: longread_locus_summary.tsv
"""

import csv
import glob
import re
from collections import defaultdict

LR = "/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc"
STRAINS = ["CiB10637", "CiB10992", "VFC140", "Cpos1038", "Cpos3700"]

anch = list(csv.DictReader(open("sowgp_anchored_longread.tsv"), delimiter="\t"))

out = []
for s in STRAINS:
    # genomic units from tblastn
    firsts = defaultdict(list)
    for line in open(f"longread_unit_tblastn/{s}.tsv"):
        f = line.split("\t")
        qs, ss, se = int(f[4]), int(f[6]), int(f[7])
        if qs == 1:
            firsts[f[1]].append(min(ss, se))
    # funannotate models in the array window
    gff = glob.glob(f"{LR}/{s}/*.gff3")[0]
    genes = []
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "gene":
            genes.append((f[0], int(f[3]), int(f[4]), re.search(r"ID=([^;]+)", f[8]).group(1)))
    # miniprot mRNA records per scaffold
    mp = defaultdict(list)
    for line in open(f"longread_miniprot/{s}.sowgp.gff"):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "mRNA":
            ident = float(re.search(r"Identity=([\d.]+)", f[8]).group(1))
            tgt = re.search(r"Target=(\S+) (\d+) (\d+)", f[8])
            if re.search(r"CIB10637|SOWgp66", tgt.group(1)):
                mp[f[0]].append(
                    (
                        ident,
                        tgt.group(1),
                        int(f[3]),
                        int(f[4]),
                        int(tgt.group(3)),
                        "Frameshift" in f[8],
                    )
                )
    for sc, pos in firsts.items():
        pos.sort()
        lo, hi = pos[0], pos[-1] + 200
        models = [g for g in genes if g[0] == sc and g[2] >= lo - 500 and g[1] <= hi + 500]
        ids = [m[3] for m in models]
        n_anch = sum(
            int(a["n_anchor"])
            for a in anch
            if a["protein"].rsplit("-T", 1)[0] in ids or a["protein"].replace("-T1", "") in ids
        )
        best = max(mp.get(sc, [(0, "", 0, 0, 0, False)]))
        out.append(
            {
                "strain": s,
                "scaffold": sc,
                "array_first": pos[0],
                "genomic_units": len(pos),
                "unit_spacing": ",".join(str(b - a) for a, b in zip(pos, pos[1:], strict=False)),
                "funannotate_models": ";".join(f"{m[3]}:{m[1]}-{m[2]}" for m in models),
                "n_funannotate_models": len(models),
                "funannotate_anchors": n_anch,
                "miniprot_best": f"{best[1]}:{best[2]}-{best[3]}",
                "miniprot_ident": best[0],
                "miniprot_aa": best[4],
                "miniprot_frameshift": best[5],
            }
        )

cols = list(out[0])
with open("longread_locus_summary.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in out:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")
for r in out:
    print(
        r["strain"],
        r["scaffold"],
        "genomic_units",
        r["genomic_units"],
        "| funannotate models",
        r["n_funannotate_models"],
        "anchors",
        r["funannotate_anchors"],
        "| miniprot",
        r["miniprot_best"],
        r["miniprot_ident"],
        r["miniprot_aa"],
        "fs" if r["miniprot_frameshift"] else "",
    )
