import collections
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
rep = {}
for line in open(sys.argv[2]):
    a, r = line.split()
    rep[a] = r
asco = {
    "Aleuria aurantia",
    "Arthrobotrys oligospora",
    "Aspergillus fumigatus",
    "Aspergillus oryzae",
    "Pseudallescheria apiosperma (Scedosporium apiospermum)",
    "Cladosporium fulvum",
    "Zymoseptoria tritici",
    "Sclerotinia sclerotiorum",
    "Gibberella zeae",
    "Magnaporthe oryzae",
    "Neurospora crassa",
    "Tuber borchii",
}
keep = [
    r
    for r in rows
    if r["uniprot"] in rep
    and not (
        "PA14" in r["unilectin_class"]
        or "Emp" in r["unilectin_class"]
        or "Ralstonia" in r["species"]
    )
]


def hits(r):
    return {h.split(":")[0] for h in r["pfam_ga_hits"].split(";") if h}


for name, M in [
    ("five", {"PF07938", "PF21691", "PF00337", "PF07367", "PF21595"}),
    ("PF07938", {"PF07938"}),
]:
    for sub, f in [("all", lambda r: True), ("asco", lambda r: r["species"] in asco)]:
        k = [r for r in keep if f(r)]
        h = [r for r in k if M & hits(r)]
        print(
            name,
            sub,
            "seqs",
            len(h),
            "/",
            len(k),
            "clusters",
            len({rep[r["uniprot"]] for r in h}),
            "/",
            len({rep[r["uniprot"]] for r in k}),
            [r["uniprot"] for r in h],
        )
c = collections.defaultdict(set)
for r in keep:
    c[r["species"]].add(rep[r["uniprot"]])
print(sorted(((len(v), s) for s, v in c.items()), reverse=True)[:4])
for r in rows:
    print(r["uniprot"], r["unilectin_class"], r["pfam_ga_hits"], rep.get(r["uniprot"]), sep="\t")
