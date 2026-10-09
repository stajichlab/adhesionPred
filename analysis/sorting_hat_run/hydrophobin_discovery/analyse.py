import collections
import csv
import gzip
import re

W = "_workdir/sorting_hat/calibration/hydrophobin_discovery"
S = "_workdir/sorting_hat"
PROTS = [
    "Afum_Af293_UniProt",
    "Afum_A1163",
    "Afum_W72310",
    "Scer_S288C",
    "Calb_SC5314",
    "Cimm_RS",
    "Bder_ER3",
    "Cneo_H99",
]
HYD = {"Hydrophobin", "Hydrophobin_2", "Eas", "DewD", "Hyd1F", "Hydrophobin_D", "Hydrophobin_like"}
# 8-Cys hydrophobin signature, generous spacing for class I and II: C-Xn-CC-Xn-C-Xn-C-Xn-CC-Xn-C-Xn
MOTIF = re.compile(r"C.{3,12}CC.{8,45}C.{5,30}C.{3,12}CC.{4,25}C.{0,20}")


def fasta(p):
    """Return {id: (sequence, header)}."""
    d, n, hdr, b = {}, None, "", []
    for line in open(p):
        line = line.rstrip()
        if line.startswith(">"):
            if n:
                d[n] = ("".join(b).rstrip("*"), hdr)
            n = line[1:].split()[0]
            hdr = line[1:]
            b = []
        else:
            b.append(line)
    if n:
        d[n] = ("".join(b).rstrip("*"), hdr)
    return d


summary = collections.defaultdict(lambda: collections.defaultdict(set))
for p in PROTS:
    for line in open(f"{W}/{p}.domtbl"):
        if line.startswith("#"):
            continue
        c = line.split(None, 22)
        summary[p][c[3]].add(c[0])
print("proteins with a hit, per model (cut_ga):")
models = sorted({m for p in PROTS for m in summary[p]})
print("model".ljust(18) + " ".join(p[:6].rjust(7) for p in PROTS))
for m in models:
    print(m.ljust(18) + " ".join(str(len(summary[p].get(m, ()))).rjust(7) for p in PROTS))
print()
print("Hydrophobin-class models: proteins hit by a model other than PF01185 (Hydrophobin):")
for p in PROTS:
    seqs = fasta(f"{S}/{p}.faa")
    h01185 = summary[p].get("Hydrophobin", set())
    others = (
        set().union(*[summary[p].get(m, set()) for m in HYD - {"Hydrophobin"}]) if HYD else set()
    )
    extra = others - h01185
    allhyd = h01185 | others
    sp = {}
    try:
        for r in csv.DictReader(
            gzip.open(f"{S}/{p}/modules/step1_rule@R0.tsv.gz", "rt"), delimiter="\t"
        ):
            sp[r["id"]] = r["call"] == "called"
    except OSError:
        pass
    print(
        f"== {p}: PF01185 hits {len(h01185)}; hit by another hydrophobin-class model only: {len(extra)} {[(e, [m for m in HYD if e in summary[p].get(m, ())]) for e in sorted(extra)][:6]}"
    )
    # Cys-motif scan among secreted, short proteins, not hit by any hydrophobin model
    cand = []
    for pid, (s, hdr) in seqs.items():
        if pid in allhyd or len(s) > 400 or len(s) < 60:
            continue
        if sp and not sp.get(pid, False):
            continue
        if MOTIF.search(s) and s.count("C") >= 8:
            cand.append((pid, len(s), s.count("C"), hdr[:80]))
    print(f"   Cys-motif, secreted, <=400 aa, no hydrophobin-model hit: {len(cand)}")
    for c in cand[:8]:
        print("     ", c)
