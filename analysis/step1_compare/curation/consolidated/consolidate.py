"""Merge the H99 pass-2 and U. maydis rows, apply the owner rulings and reviewer fixes.

Rulings (owner, 2026-10-03): predictor-selected rows form a report-only sub-stratum;
strong-promoter rows stay with an `overexpressed` flag (option B); Cpl1 is dropped.
Source files stay unchanged. Every edit is listed in changes.tsv and dropped_rows.tsv.
"""

import csv

H99 = "h99_pass2/draft_rows.tsv"
UM = "umay_pilot/draft_rows.tsv"
DATE = "2026-10-03"

# (source, source file line) -> reason. Lines are file line numbers (header = line 1).
DROP = {
    ("H99", 22): "SOD1 N-int: GOA already gives N-sec (PM and raft, PMID 16524904); a curated row cannot override GOA (reviewer)",
    ("H99", 24): "CIG1 plasma membrane: C-terminal mCherry on a GPI protein, ACT1 overexpression, no membrane marker (reviewer)",
    ("H99", 33): "VCX1: sentence reports vacuole shape, not Vcx1 location (reviewer)",
    ("UM", 31): "Cpl1: permeabilisation not stated and a faint supernatant signal; owner ruling 2026-10-03",
}  # fmt: skip

# symbol-based edits: (source, symbol) -> dict(field -> value or ('append', text))
NOTE = "note"
EDITS = []  # (source, symbol, line or None, field, new value, reason)


def edit(source, symbol, line, field, value, reason):
    EDITS.append((source, symbol, line, field, value, reason))


edit("H99", "PLB1", 10, "selected_by_predictor", "no", "reviewer: should be no, not unknown")
edit(
    "UM",
    "Pep1",
    None,
    "selected_by_predictor",
    "yes",
    "reviewer: the paper chose genes for a secretion signal (SignalP, Fig. 1A)",
)
for gene in ("Sir2", "Hst4", "Hst5", "Hst6"):
    edit(
        "UM",
        gene,
        None,
        NOTE,
        "strain stated in the Methods: all strains used derive from the haploid SG200 strain (reviewer, full text)",
        "sirtuin strain",
    )
for gene, line in (("Xyn2", 20), ("Xyn11A", 21)):
    edit(
        "UM",
        gene,
        line,
        NOTE,
        "Table 1 shows the fusion is expressed from the otef promoter (overexpression)",
        "reviewer: promoter",
    )
edit(
    "UM",
    "Stp2",
    None,
    NOTE,
    "panel a (AB33-derived strains, otef), stated in the Methods",
    "reviewer: panel map",
)
edit(
    "UM",
    "Stp3",
    6,
    NOTE,
    "panel a (AB33-derived strains, otef), stated in the Methods",
    "reviewer: panel map",
)
edit(
    "UM",
    "Stp4",
    7,
    NOTE,
    "panel a (AB33-derived strains, otef), stated in the Methods",
    "reviewer: panel map",
)
edit(
    "UM",
    "Stp1",
    4,
    NOTE,
    "panel b (SG200 delta kex2); this panel mapping is the reviewer's reading",
    "reviewer: panel map",
)
edit(
    "UM",
    "Stp5",
    None,
    NOTE,
    "debatable negative: membrane protein in a surface-exposed complex; panel c or d is the reviewer's inference, the legend does not say",
    "reviewer",
)
edit(
    "UM",
    "Stp6",
    None,
    NOTE,
    "debatable negative: membrane protein in a surface-exposed complex; surface speckles map to GO:0009986 (no label change); panel c or d is the reviewer's inference",
    "reviewer",
)
edit(
    "UM",
    "Pit2",
    None,
    NOTE,
    "UMAG ID comes from UniProt REST (UMAG_01375); proteome FASTA headers carry no UMAG ID",
    "reviewer: accession note",
)
edit(
    "UM",
    "Cmu1",
    None,
    NOTE,
    "UMAG ID comes from UniProt REST (UMAG_05731); proteome FASTA headers carry no UMAG ID",
    "reviewer: accession note",
)
edit(
    "UM",
    "Pdi1",
    None,
    NOTE,
    "risk: PMID 34947062 neighbour paper (Cpl1) says Pdi1 is also in the cell wall; Pdi1 S1 Table (cell-wall extract) not read; if listed, this row becomes a conflict",
    "reviewer",
)
edit(
    "UM",
    "Rsp3",
    2,
    NOTE,
    "wall term rests on the figure title plus 'around the outside of fungal hyphae'; strong for P-ext, weak for the wall subset (reviewer)",
    "reviewer",
)
edit(
    "H99",
    "CDA1",
    None,
    NOTE,
    "CDA1 is P-ext only through GOA homology (ISS) terms; the curated row adds GO:0016020 (membrane), which counts as neither surface nor secretory; not a direct P-ext gene (reviewer)",
    "reviewer",
)
for gene in ("BIM1",):
    edit(
        "H99",
        gene,
        None,
        NOTE,
        "parent strain of the tagged mutants is named only in Supplementary Dataset 1 (not read); wild-type H99 is the control in every figure (reviewer)",
        "reviewer",
    )
for gene in ("MNS1", "MNS101"):
    edit(
        "H99",
        gene,
        None,
        NOTE,
        "tagged in 'the WT strain'; H99 not named in the sentence (reviewer)",
        "reviewer",
    )

# overexpressed flag: yes = strong constitutive/heterologous promoter (otef, ACT1, pit2 promoter per reviewer);
# no = native promoter stated; unknown = not judged. Keyed by (source, file line).
OVER = {
    ("UM", 2): "no", ("UM", 3): "yes", ("UM", 4): "yes", ("UM", 5): "yes", ("UM", 6): "yes", ("UM", 7): "yes",
    ("UM", 8): "no", ("UM", 9): "no", ("UM", 10): "unknown", ("UM", 11): "unknown", ("UM", 12): "no",
    ("UM", 14): "no", ("UM", 15): "no", ("UM", 18): "yes", ("UM", 20): "yes", ("UM", 21): "yes",
    ("UM", 22): "yes", ("UM", 28): "unknown", ("UM", 29): "unknown", ("UM", 30): "no",
    ("H99", 32): "yes",
}  # fmt: skip


def read(path, tag):
    with open(path) as fh:
        for line_no, row in enumerate(csv.DictReader(fh, delimiter="\t"), 2):
            row["source_pass"] = tag
            row["source_line"] = str(line_no)
            yield row


rows = list(read(H99, "H99")) + list(read(UM, "UM"))
cols = list(rows[0].keys())
cols = [c for c in cols if c not in ("source_pass", "source_line")] + [
    "overexpressed",
    "source_pass",
    "source_line",
]
out, dropped, changes = [], [], []
for row in rows:
    key = (row["source_pass"], int(row["source_line"]))
    if key in DROP:
        dropped.append([row["source_pass"], row["source_line"], row["symbol"], DROP[key]])
        continue
    for source, symbol, line, field, value, reason in EDITS:
        if source != row["source_pass"] or symbol != row["symbol"]:
            continue
        if line is not None and line != int(row["source_line"]):
            continue
        old = row["evidence_note"] if field == NOTE else row[field]
        if field == NOTE:
            row["evidence_note"] += f" | consolidation {DATE}: {value}"
        else:
            row[field] = value
        changes.append(
            [
                row["source_pass"],
                row["source_line"],
                row["symbol"],
                field,
                old[-60:] if field == NOTE else old,
                value,
                reason,
            ]
        )
    row["overexpressed"] = OVER.get(key, "unknown")
    out.append(row)

with open("consolidated/draft_rows.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(out)
with open("consolidated/dropped_rows.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["source_pass", "source_line", "symbol", "reason"])
    w.writerows(dropped)
with open("consolidated/changes.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["source_pass", "source_line", "symbol", "field", "old", "new", "reason"])
    w.writerows(changes)
print(len(rows), len(out), len(dropped), len(changes))
