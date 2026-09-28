#!/usr/bin/env python3
"""Annotate the intersected Coccidioides antigen candidates.

Adds, for each candidate: length, signal peptide, TM helices, Pfam domains, from the Fungi_5k
functional database (via the FungiDB -> Fungi_5k id map built by 01_build_inputs.sh).

A serodiagnostic antigen should be secreted or surface-exposed. Candidates with no signal
peptide and no TM helix are intracellular and are flagged, not silently kept.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import duckdb

DB = Path("/bigdata/stajichlab/shared/projects/Fungi_5k/functionalDB/function.duckdb")
RS_PREFIX = "FA2214EC"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--idmap", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.candidates), delimiter="\t"))
    # FungiDB protein -> best Fungi_5k protein
    best = {}
    for line in open(args.idmap):
        f = line.rstrip("\n").split("\t")
        q, t, pid = f[0], f[1], float(f[2])
        if q not in best or pid > best[q][1]:
            best[q] = (t, pid)
    print(f"{len(rows)} candidates; id map covers {len(best)} reference proteins")

    ids = [best[r["protein"]][0] for r in rows if r["protein"] in best]
    print(f"{len(ids)} candidates map into Fungi_5k")

    con = duckdb.connect(str(DB), read_only=True)
    con.execute("PRAGMA memory_limit='3GB'")
    con.register("_w", __import__("pandas").DataFrame({"pid": ids}))

    def fetch(table, cols="protein_id"):
        q = (
            f"SELECT {cols} FROM {table} t JOIN _w w ON "
            f"replace(t.protein_id,'.protein','') = w.pid "
            f"WHERE t.protein_id LIKE '{RS_PREFIX}%'"
        )
        return con.execute(q).fetchdf()

    sp = set(fetch("signalp").protein_id.str.replace(".protein", "", regex=False))
    tm = set(fetch("tmhmm").protein_id.str.replace(".protein", "", regex=False))
    pf = fetch("pfam", "protein_id, pfam_id")
    doms = defaultdict(set)
    for p, d in zip(pf.protein_id.str.replace(".protein", "", regex=False), pf.pfam_id):
        doms[p].add(d)
    lens = fetch("gene_proteins", "protein_id, length")
    length = dict(zip(lens.protein_id.str.replace(".protein", "", regex=False), lens.length))

    out = []
    for r in rows:
        fid = best.get(r["protein"], ("", 0))[0]
        d = sorted(doms.get(fid, []))
        has_sp, has_tm = fid in sp, fid in tm
        out.append(
            {
                **r,
                "fungi5k_id": fid,
                "length": length.get(fid, ""),
                "signal_peptide": "yes" if has_sp else "no",
                "tm_helix": "yes" if has_tm else "no",
                "pfam_domains": ";".join(d),
                "n_pfam": len(d),
                "surface_plausible": "yes" if (has_sp or has_tm) else "NO - intracellular?",
            }
        )
    out.sort(key=lambda r: -float(r.get("spherule48h_TPM", 0) or 0))
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    n_sp = sum(1 for r in out if r["signal_peptide"] == "yes")
    n_nodom = sum(1 for r in out if r["n_pfam"] == 0)
    n_surf = sum(1 for r in out if r["surface_plausible"] == "yes")
    print(f"\n{len(out)} candidates annotated")
    print(
        f"   signal peptide: {n_sp}   TM helix: {sum(1 for r in out if r['tm_helix'] == 'yes')}"
        f"   surface-plausible: {n_surf}"
    )
    print(f"   no Pfam domain at all: {n_nodom}")
    print(f"\n{'protein':<21}{'len':>5}{'SP':>4}{'TM':>4}{'sph48h':>9}{'log2FC':>8}  domains")
    for r in out:
        print(
            f"{r['protein'][:20]:<21}{str(r['length']):>5}{r['signal_peptide'][:1].upper():>4}"
            f"{r['tm_helix'][:1].upper():>4}{float(r.get('spherule48h_TPM', 0)):>9.0f}"
            f"{float(r.get('log2FC_spherule48h', 0)):>8.2f}  {r['pfam_domains'][:46] or '(none)'}"
        )


if __name__ == "__main__":
    main()
