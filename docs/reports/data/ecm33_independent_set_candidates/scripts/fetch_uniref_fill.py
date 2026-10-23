# ruff: noqa
# Scoping script kept exactly as run on 2026-10-23 (stdlib only).
import csv, sys, urllib.request, time

have = {l.split("\t")[0] for l in open(sys.argv[1]) if not l.startswith("accession")}
ids = [l.strip() for l in open(sys.argv[2]) if l.strip() and l.strip() not in have]
n = 0
with open(sys.argv[1], "a") as f:
    for a in ids:
        u = (
            "https://rest.uniprot.org/uniref/search?query=(identity:0.5)+AND+(uniprot_id:%s)&format=tsv&fields=id"
            % a
        )
        for k in range(5):
            try:
                t = (
                    urllib.request.urlopen(
                        urllib.request.Request(
                            u, headers={"User-Agent": "adhesionPred-scoping/1.0"}
                        ),
                        timeout=40,
                    )
                    .read()
                    .decode()
                    .split("\n")
                )
                break
            except Exception as e:
                time.sleep(2)
                t = []
        r = [x for x in t[1:] if x]
        if r:
            f.write("%s\t%s\n" % (a, r[0].strip()))
            n += 1
print(len(ids), n)
