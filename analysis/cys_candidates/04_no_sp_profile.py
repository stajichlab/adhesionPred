#!/usr/bin/env python3.12
"""What does the cys_rich_no_sp tier hold? Two steps, read-only, standard library only.

Step 1 (--write-fasta OUT.faa): write the cys_rich_no_sp proteins of all proteomes to a FASTA
  (ids are <proteome>|<protein_id>). Then run, on a login node (about 3 minutes for 154
  proteins against the full Pfam-A of the central database):
    hmmsearch --cut_ga --noali -o /dev/null --domtblout nosp.domtbl PFAM_A_HMM OUT.faa
Step 2 (--domtbl nosp.domtbl): print counts: proteins, with a Pfam hit (GA cut-off), without,
  with a zinc-finger / Fe-S / RING-type domain (name rule below), first residue M, maximum SP
  probability, and the number of each Pfam name.

The name rule is a keyword list, not a curated annotation: see ZN_FE_RULE.
"""

import argparse
import gzip
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cys_candidates as cc  # noqa: E402

ZN_FE_RULE = re.compile(
    r"^(zf-|Fer4|ANAPC|BBOX|ANCHR-like_BBOX|BUD31|MYND|IBR|Rcat|PHF5|Nab2|Zn_clus|GFA)"
)


def load(path):
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, ln.rstrip("\n").split("\t"), strict=True)) for ln in fh]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--table-dir", required=True)
    ap.add_argument("--write-fasta")
    ap.add_argument("--domtbl")
    a = ap.parse_args()
    try:
        jobs = cc.read_manifest(a.manifest)
        rows = {}
        for name, fasta, _sp, _dom in jobs:
            tbl = Path(a.table_dir) / f"{name}.tsv.gz"
            if not tbl.is_file():
                raise cc.Stop(f"table not found: {tbl}")
            seqs = dict(cc.read_fasta(fasta))
            for r in load(tbl):
                if r["tier"] == "cys_rich_no_sp":
                    rows[f"{name}|{r['protein_id']}"] = (r, seqs[r["header"]])
        if a.write_fasta:
            tmp = Path(a.write_fasta + ".tmp")
            tmp.write_text("".join(f">{k}\n{s}\n" for k, (_r, s) in rows.items()))
            tmp.replace(a.write_fasta)
            print(f"wrote {len(rows)} proteins to {a.write_fasta}")
        if a.domtbl:
            names = {}
            with cc.open_text(a.domtbl) as fh:
                for ln in fh:
                    if ln.startswith("#") or not ln.strip():
                        continue
                    f = ln.split(None, 4)
                    if f[0] not in rows:
                        raise cc.Stop(f"domtbl target {f[0]!r} is not a cys_rich_no_sp protein")
                    names.setdefault(f[0], set()).add(f[3])
            n = len(rows)
            hit = [k for k in rows if k in names]
            zn = [k for k in hit if any(ZN_FE_RULE.match(x) for x in names[k])]
            met = sum(1 for _r, s in rows.values() if s.startswith("M"))
            maxp = max(float(r["sp_prob"]) for r, _s in rows.values())
            print(f"proteins\t{n}")
            print(f"with_pfam_hit_cut_ga\t{len(hit)}")
            print(f"without_pfam_hit\t{n - len(hit)}")
            print(f"with_zn_fe_ring_type_name\t{len(zn)}")
            print(f"hit_but_not_zn_fe_ring_type\t{len(hit) - len(zn)}")
            print(f"starts_with_M\t{met}")
            print(f"max_sp_prob\t{maxp:.4f}")
            cnt = Counter(x for k in hit for x in names[k])
            print("pfam_name\tn_proteins")
            for k, v in sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0])):
                print(f"{k}\t{v}")
    except cc.Stop as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
