# Rubric for the evidence-based curation of the 61 unlabelled hydrophobin calls (tier T4)

*2026-10-08. Written before any decision is made. The owner cannot curate the sheets and asked for a curation by evidence and literature review. These decisions are made by the assistant from the
evidence below, are tier **T4** (curator judgement from evidence, not owner-reviewed, not experimental), are stored apart from the owner's columns in `evidence_sheets.tsv`, are never merged with T1/T2 in any
recall or specificity table, and can be audited and overruled by the owner. Every decision lists the evidence lines that support it.*

## Evidence lines

| Line | Meaning | Independent of the relaxed Pfam call? |
|---|---|---|
| A | **Doublet architecture.** The sequence has 8 cysteines (or 8 to 10) in a window with the 2nd and 3rd cysteine adjacent and the 6th and 7th adjacent (gap 0), a signal peptide (R0 called), and a length of at most 300 aa. The doublet pattern is the qualitative feature in Linder 2005 (p. 878, read in full) and the other reviews. | yes (no Pfam, no model score) |
| B | **Strict Pfam.** A hit to a hydrophobin-class Pfam model at the gathering cutoff (`pfam_hydrophobin` hit). | no (same models, higher cutoff) |
| C | **Curated name or literature.** A paper, or a UniProt, FungiDB or other curated record, names this gene or its orthologous gene product a hydrophobin. The source (PMID, accession or record) is recorded. | yes |
| D | **Ortholog.** Reciprocal best hit of the protein with a T1 or T2 hydrophobin at 40% identity or more over 70% or more of both sequences, or a one-to-one ortholog in a species where hydrophobins are annotated. | partly (uses sequence similarity to the training set) |
| E | **Other function.** The protein has a strong hit (40% identity or more over 70% or more of the query) to a Swiss-Prot entry annotated with a non-hydrophobin function, or a full Pfam-A scan at the gathering cutoff finds a domain of a non-hydrophobin family (for example cutinase, glycoside hydrolase, receptor kinase) that accounts for most of the sequence. | yes |

## Decision rules (applied in this order)

1. **hydrophobin (T4)** if C holds, or if A holds together with at least one of B or D.
2. **not hydrophobin (T4)** if E holds and C does not hold.
3. **unresolved** otherwise. The protein stays unresolved when only A holds, when only the relaxed score supports it, or when evidence conflicts.

A protein that is `unresolved` is not counted as a false positive and not as a hydrophobin. It is reported in its own column.

## Rules for the evidence

- Each line is recorded as `yes`, `no` or `not checked`, with the source. No line is filled from memory.
- Literature (line C): PubMed search by species and "hydrophobin" for genome-wide hydrophobin lists or characterisation papers; any gene ID is mapped to the proteome ID by sequence (BLAST, at least 95% identity over at least 90% of the length). Abstract-only reading is allowed for line C only if the abstract names the gene.
- UniProt names of TrEMBL entries inferred by rule or similarity are not literature; they count for line C only if the record cites a paper or an experimental annotation.
- Precision for the ship rule (spec 6.5) is computed with T4 outcomes, reported in separate columns from T1 and T5. The result is "not owner-reviewed". The owner can overrule any row.
- T4 outcomes are not used to train any model and are not added to the truth set.
