#!/usr/bin/python3.12
"""Derived HMM file for the relaxed Pfam level: the seven hydrophobin-class Pfam models with the GA line rewritten
to one full-sequence cutoff and a domain cutoff of -1000 (so the domain test never binds).

Usage:
  make_relaxed_hmm.py template --pfam PATH --out FILE      fetch the seven models (hmmfetch) unchanged
  make_relaxed_hmm.py final --template FILE --cutoff X --out FILE
Needs hmmer on PATH for `template`.
"""

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

MODELS = [
    "Hydrophobin",
    "Hydrophobin_2",
    "Eas",
    "DewD",
    "Hyd1F",
    "Hydrophobin_D",
    "Hydrophobin_like",
]
DOMAIN_GA = -1000.0


def rewrite_ga(text, cutoff):
    if not re.search(r"^GA\s", text, re.M):
        raise ValueError("no GA line in the model file")
    return re.sub(r"^GA\s.*$", f"GA    {cutoff:.2f} {DOMAIN_GA:.2f};", text, flags=re.M)


def score_table(lines):
    best = {}
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        c = line.split()
        s = float(c[5])
        if c[0] not in best or s > best[c[0]]:
            best[c[0]] = s
    return best


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("cmd", choices=["template", "final"])
    ap.add_argument("--pfam")
    ap.add_argument("--template")
    ap.add_argument("--cutoff", type=float)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "template":
        pfam = Path(a.pfam)
        keys = Path(a.out).with_suffix(".keys")
        keys.write_text("\n".join(MODELS) + "\n")
        with open(a.out, "w") as f:
            subprocess.run(["hmmfetch", "-f", str(pfam), str(keys)], stdout=f, check=True)
        n = sum(1 for line in open(a.out) if line.startswith("NAME"))
        if n != len(MODELS):
            raise SystemExit(f"hmmfetch gave {n} models, expected {len(MODELS)}")
        prov = {
            "pfam_path": str(pfam),
            "pfam_resolved": str(pfam.resolve()),
            "pfam_sha256": sha(pfam.resolve()),
            "models": MODELS,
            "template_sha256": sha(a.out),
        }
        json.dump(prov, open(Path(a.out).with_suffix(".provenance.json"), "w"), indent=1)
        print(json.dumps(prov, indent=1))
    else:
        text = rewrite_ga(open(a.template).read(), a.cutoff)
        open(a.out, "w").write(text)
        print(f"wrote {a.out} with GA {a.cutoff:.2f} {DOMAIN_GA:.2f}; sha256 {sha(a.out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
