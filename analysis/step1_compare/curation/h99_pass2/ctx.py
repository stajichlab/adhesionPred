import re
import sys

f = sys.argv[1]
t = open("texts/" + f).read()
n = int(sys.argv[2]) if sys.argv[2].isdigit() else 250
pats = sys.argv[3:] if sys.argv[2].isdigit() else sys.argv[2:]
for p in pats:
    print("##", p)
    for m in list(re.finditer(p, t, re.I))[:4]:
        print("  >>", t[max(0, m.start() - n) : m.end() + n].replace("\n", " | "))
