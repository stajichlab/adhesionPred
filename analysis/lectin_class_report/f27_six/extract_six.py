"""Extract the six F27 PF07938 proteins (spec F27, group 1) from the source proteomes.

Usage: /usr/bin/python3.12 extract_six.py POOLS_DIR SELECTION_TSV WHOLE_PROTEOME_HITS_TSV OUT_FAA OUT_TSV

Selection: proteins that the discovery pools table lists with PF07938 in the pfam column, in a
proteome of group onygenales, with a Pfam 38.2 PF07938 sequence score of 100 or more in the
whole-proteome hit table. The script stops unless there are exactly 6, each in a different genus,
and each with the same length in the pools table and in the source FASTA.
"""

import csv
import gzip
import os
import sys

pools_dir, sel_tsv, hits_tsv, out_faa, out_tsv = sys.argv[1:6]
sel = list(csv.DictReader(open(sel_tsv), delimiter="\t"))
score = {}
for r in csv.DictReader(open(hits_tsv), delimiter="\t"):
    if r["model"] == "Fungal_lectin":
        score[(r["proteome"], r["protein"])] = float(r["seq_score"])
found = []
for s in sel:
    if s["group"] != "onygenales":
        continue
    with gzip.open(os.path.join(pools_dir, s["name"] + ".pools.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pf = (r["pfam"] or "").split(";")
            if "PF07938" in pf and score.get((s["name"], r["protein"]), 0) >= 100:
                found.append((s, r))
assert len(found) == 6, len(found)
assert len({s["genus"] for s, _ in found}) == 6
with open(out_faa, "w") as fa, open(out_tsv, "w") as tv:
    tv.write("accession\tgenus\tproteome\tlength\tpf07938_score\tpool\tsource_file\n")
    for s, r in found:
        want = r["protein"]
        seq, cur, got = [], None, None
        for line in open(s["fasta"]):
            if line.startswith(">"):
                if cur == want:
                    break
                cur = line[1:].split()[0]
                seq = []
            elif cur == want:
                seq.append(line.strip())
        if cur == want:
            got = "".join(seq).rstrip("*")
        assert got is not None, want
        assert len(got) == int(r["length"]), (want, len(got), r["length"])
        fa.write(f">{want} {s['name']}\n{got}\n")
        row = [want, s["genus"], s["name"], str(len(got)), "%.1f" % score[(s["name"], want)]]
        tv.write("\t".join(row + [r["pool"], s["fasta"]]) + "\n")
print("OK: 6 proteins")
