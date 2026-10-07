#!/usr/bin/env python3
"""Write DRAFT member lists for the Pfam specificity review (Plan 2, Task 13 step 3).

A member is a protein that belongs to a family by function or by its name, not because Pfam
found a domain in it. These lists are drafts. The names come from UniProt and SGD headers and from
papers checked in docs/paper/05-literature-verification.md (CfmA-C, CalA). They were not reviewed
by an expert. Do not report sensitivity or specificity from them (see the plan, Task 13 step 3).

Output per proteome: members.tsv (columns pfam_acc, protein_id) and members.rationale.tsv.
Usage: make_draft_members.py --proteome NAME --fasta FILE --out-dir DIR
"""

import argparse
import gzip
import sys

# family -> list of (key, label, source note). key = UniProt accession or systematic ORF name.
AF293 = {
    "PF05730": [
        ("Q4WLB9", "CfmA", "PMID 24361821"),
        ("Q4WMA6", "CfmB", "PMID 24361821"),
        ("Q4WNE1", "CfmC", "PMID 24361821"),
    ],
    "PF04681": [("Q4WXJ1", "CalA", "PMID 27841851; docs/reports/2026-09-29-class2b-structure.md")],
    "PF01185": [
        (a, n, "UniProt gene name rod*")
        for a, n in (
            ("P41746", "RodA"),
            ("E9QT94", "RodB"),
            ("Q4WBR8", "RodC"),
            ("Q4WC31", "RodE"),
            ("Q4WE22", "RodD"),
            ("Q4WEK0", "RodF"),
            ("Q4X055", "RodG"),
        )
    ],
    "PF25312": [("O60024", "Asp f 4", "WHO/IUIS Asp f 4")],
}
SCER = {
    "PF10182": [("YIR019C", "FLO11", "SGD name")],
    "PF10528": [
        (o, n, "SGD name, lectin-like flocculin")
        for o, n in (
            ("YAR050W", "FLO1"),
            ("YHR211W", "FLO5"),
            ("YAL063C", "FLO9"),
            ("YKR102W", "FLO10"),
        )
    ],
    "PF07691": [
        (o, n, "SGD name, flocculin with a PA14 lectin domain")
        for o, n in (
            ("YAR050W", "FLO1"),
            ("YHR211W", "FLO5"),
            ("YAL063C", "FLO9"),
            ("YKR102W", "FLO10"),
        )
    ],
    "PF00399": [
        (o, n, "SGD name, Pir protein")
        for o, n in (
            ("YKL164C", "PIR1"),
            ("YKL163W", "PIR3"),
            ("YJL160C", "PIR5"),
            ("YJL159W", "HSP150"),
            ("YJL158C", "CIS3"),
        )
    ],
}
SETS = {
    "Afum_Af293_UniProt": AF293,
    "Afum_A1163": AF293,
    "Afum_W72310": AF293,
    "Scer_S288C": SCER,
    "Scer_S288C_nodubious": SCER,
}
ALIASES = {
    "PF06766": "PF01185",
    "PF28987": "PF01185",
    "PF22354": "PF01185",
}  # other hydrophobin models share the member list


def ids(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    table = SETS.get(a.proteome, {})
    all_ids = ids(a.fasta)
    members, why, missing = [], [], []
    for fam, entries in table.items():
        fams = [fam] + [k for k, v in ALIASES.items() if v == fam]
        for key, label, note in entries:
            hits = [
                i for i in all_ids if i == key or i.split("|")[1:2] == [key] or i.split()[0] == key
            ]
            if not hits:
                missing.append((fam, key, label))
                continue
            for f in fams:
                members.append((f, hits[0]))
                why.append((f, hits[0], label, note, "draft, not reviewed by an expert"))
    with open(f"{a.out_dir}/members.tsv", "w") as f:
        f.write("pfam_acc\tprotein_id\n")
        f.writelines(f"{x}\t{y}\n" for x, y in members)
    with open(f"{a.out_dir}/members.rationale.tsv", "w") as f:
        f.write("pfam_acc\tprotein_id\tlabel\tsource\tstatus\n")
        f.writelines("\t".join(r) + "\n" for r in why)
    print(
        f"{a.proteome}: {len(members)} member rows; not found in FASTA: {missing}", file=sys.stderr
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
