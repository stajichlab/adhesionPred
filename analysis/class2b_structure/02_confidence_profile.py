#!/usr/bin/env python3.12
"""Per-residue pLDDT profile and confident-core extraction.

A structural search is only meaningful over the part of a model that is actually folded.
AlphaFold stores per-residue pLDDT in the B-factor column of its mmCIF. This script reads it,
reports the confidence distribution per protein, and writes a "confident core" PDB containing
only residues with pLDDT >= the cutoff (default 70, AlphaFold's own confident threshold).

Residues below 70 are treated as not-a-fold, not as a low-quality fold. Foldseek alignments
that run through such residues are not evidence of a fold relationship.

    /usr/bin/python3.12 02_confidence_profile.py --workdir <workdir>
"""

import argparse
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# 3-letter to 1-letter, only what AlphaFold emits
AA3 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def read_cif_atoms(path):
    """Return (list of atom dicts, ordered list of (resnum, resname, plddt))."""
    lines = path.read_text().splitlines()
    # find the atom_site loop header
    cols, atoms, in_loop, header = [], [], False, False
    for ln in lines:
        if ln.startswith("_atom_site."):
            cols.append(ln.strip().split(".")[1])
            header = True
            continue
        if header and (ln.startswith("ATOM") or ln.startswith("HETATM")):
            in_loop = True
        elif in_loop and (ln.startswith("#") or not ln.strip()):
            break
        if in_loop:
            f = ln.split()
            atoms.append(dict(zip(cols, f, strict=False)))
    residues = {}
    for a in atoms:
        rn = int(a["label_seq_id"])
        residues.setdefault(rn, (a["label_comp_id"], float(a["B_iso_or_equiv"])))
    ordered = [(rn, residues[rn][0], residues[rn][1]) for rn in sorted(residues)]
    return atoms, ordered


def write_pdb(atoms, out, keep_resnums=None):
    """Write a minimal PDB. AlphaFold models are single-chain, so chain A throughout."""
    with open(out, "w") as fh:
        n = 0
        for a in atoms:
            rn = int(a["label_seq_id"])
            if keep_resnums is not None and rn not in keep_resnums:
                continue
            n += 1
            # strict PDB column layout: name 13-16, altLoc 17, resName 18-20, chain 22,
            # resSeq 23-26, iCode 27, coords from 31. Foldseek silently returns nothing if
            # these columns are off by one.
            name = a["label_atom_id"]
            name = name if len(name) == 4 else (" " + name).ljust(4)
            fh.write(
                "ATOM  {:>5d} {:4s} {:>3s} A{:>4d}    "
                "{:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2s}\n".format(
                    n,
                    name,
                    a["label_comp_id"],
                    rn,
                    float(a["Cartn_x"]),
                    float(a["Cartn_y"]),
                    float(a["Cartn_z"]),
                    float(a["occupancy"]),
                    float(a["B_iso_or_equiv"]),
                    a["type_symbol"],
                )
            )
        fh.write("TER\nEND\n")
    return n


def segments(resnums):
    """Collapse a sorted residue-number list into contiguous ranges."""
    out, start, prev = [], None, None
    for r in resnums:
        if start is None:
            start = prev = r
            continue
        if r == prev + 1:
            prev = r
        else:
            out.append((start, prev))
            start = prev = r
    if start is not None:
        out.append((start, prev))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--candidates", default=str(HERE / "candidates.tsv"))
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--cutoff", type=float, default=70.0)
    ap.add_argument(
        "--min-seg",
        type=int,
        default=5,
        help="drop confident segments shorter than this many residues",
    )
    args = ap.parse_args()

    wd = Path(args.workdir)
    (wd / "core_pdb").mkdir(parents=True, exist_ok=True)
    (wd / "full_pdb").mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(open(args.candidates), delimiter="\t"))
    out_rows = []
    per_res = open(wd / "plddt_per_residue.tsv", "w")
    per_res.write("label\taccession\tresnum\taa\tplddt\n")

    for row in rows:
        acc, label = row["accession"], row["label"]
        cif = wd / "cif" / f"{label}_{acc}.cif"
        if not cif.exists():
            continue
        atoms, ordered = read_cif_atoms(cif)
        pl = np.array([p for _, _, p in ordered])
        for rn, rname, p in ordered:
            per_res.write(f"{label}\t{acc}\t{rn}\t{AA3.get(rname, 'X')}\t{p:.2f}\n")

        keep = {rn for rn, _, p in ordered if p >= args.cutoff}
        segs = [s for s in segments(sorted(keep)) if s[1] - s[0] + 1 >= args.min_seg]
        keep = {r for a, b in segs for r in range(a, b + 1)}

        write_pdb(atoms, wd / "full_pdb" / f"{label}.pdb")
        n_core = write_pdb(atoms, wd / "core_pdb" / f"{label}.pdb", keep) if keep else 0
        if not keep:
            (wd / "core_pdb" / f"{label}.pdb").unlink(missing_ok=True)

        out_rows.append(
            {
                "label": label,
                "accession": acc,
                "length": len(ordered),
                "mean_plddt": round(float(pl.mean()), 1),
                "median_plddt": round(float(np.median(pl)), 1),
                "frac_plddt_ge70": round(float((pl >= 70).mean()), 3),
                "frac_plddt_ge90": round(float((pl >= 90).mean()), 3),
                "frac_plddt_lt50": round(float((pl < 50).mean()), 3),
                "core_residues": len(keep),
                "core_segments": ";".join(f"{a}-{b}" for a, b in segs) or "none",
                "core_atoms": n_core,
                "role": row["role"],
                "pfam": row["pfam"],
            }
        )
    per_res.close()

    fields = list(out_rows[0].keys())
    with open(wd / "confidence_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    hdr = [
        "label",
        "role",
        "length",
        "mean_plddt",
        "frac_plddt_ge70",
        "core_residues",
        "core_segments",
    ]
    print("\t".join(hdr))
    for r in out_rows:
        print("\t".join(str(r[h]) for h in hdr))


if __name__ == "__main__":
    main()
