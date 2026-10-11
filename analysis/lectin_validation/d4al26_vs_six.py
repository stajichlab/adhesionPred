#!/usr/bin/env python3.12
"""Compare Swiss-Prot D4AL26 (Arthroderma benhamiae 'Fucose-specific lectin') with the six Onygenales
PF07938 proteins of spec F27 (analysis/lectin_class_report/f27_six/six.faa).
Usage (mmseqs on PATH): d4al26_vs_six.py UNIPROT_CANDIDATES_JSON_GZ SIX_FAA OUT_TSV TMP_DIR
"""

import gzip
import json
import subprocess
import sys

js, six, out, tmp = sys.argv[1:5]
d = {e["primaryAccession"]: e for e in json.load(gzip.open(js, "rt"))}
import os  # noqa: E402

os.makedirs(tmp, exist_ok=True)
with open(tmp + "/q.faa", "w") as fh:
    fh.write(">D4AL26\n" + d["D4AL26"]["sequence"]["value"] + "\n")
subprocess.run(
    [
        "mmseqs", "easy-search", tmp + "/q.faa", six, out, tmp + "/t",
        "--min-seq-id", "0.3", "-s", "7.5", "-v", "1",
        "--format-output", "query,target,pident,alnlen,qcov,tcov,evalue",
    ],
    check=True,
)  # fmt: skip
print(open(out).read())
