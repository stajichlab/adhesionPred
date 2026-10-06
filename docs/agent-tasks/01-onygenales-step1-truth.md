# Task 01. Step 1 truth set for Onygenales (*Coccidioides* and relatives)

*Read `COMMON-RULES.md` first. Written 2026-10-06.*

## Goal

Build positive and negative controls for rule R0 (`signal_peptide_protein`: "SignalP 6 calls a signal
peptide") for Onygenales, with *Coccidioides immitis* RS first.

## Why this matters

- R0 has a status only in six species (*S. cerevisiae*, *C. albicans*, *A. fumigatus*, *A. nidulans*,
  *C. neoformans*, *U. maydis*). Only *A. nidulans* is `estimated`.
- **No Onygenales truth set exists.** In *Coccidioides* the status of R0 is `unvalidated`.
- Whether SignalP under-calls in *C. immitis* is untested. SignalP calls a signal peptide in 460 of
  9,910 RS proteins (4.6%). In the spherule table, 26 of the 27 extreme *Coccidioides*-specific
  genes have no signal peptide and no GPI call (`docs/reports/2026-10-03-cocci-spherule-surface-table.md`).
- Every *Coccidioides* result of the sorting tool (antigen ranking, repeat call, family scans) depends
  on step 1. The paper needs this set.

## Module and call calibrated

`step1_rule@R0` through `calibrate truth --call signal_peptide_protein --variant R0 --module step1_rule@R0 --taxa <ONE species taxon>`.
One call per species. A genus or higher rank is refused.

## Positive controls

Proteins that are outside the plasma membrane (cell wall, outer face of the membrane, or secreted)
with **experimental** evidence in an Onygenales species:
- GO cell component terms *cell wall*, *extracellular region*, *cell surface* with evidence codes IDA,
  HDA, IMP, IGI, EXP (not IEA, not IBA, ISS, ISO, ISA).
- Literature: secretome or cell wall proteomics, surface shaving, immunolocalisation, secreted
  enzymes with a measured activity in culture filtrate. Each needs a PMID and a quoted sentence.
- Known candidates to check, not to assume: SOWgp (CIMG_04613), Ag2/PRA, CTS1, BAD1, Yps3.
  Mark SOWgp, PRA3, Ag2/PRA and PRA2 `tuned=yes` (they tuned the ranking and repeat work).

`stratum`: `cell_wall`, `secreted`, `outer_pm`.

## Negative controls

Proteins that are **inside** the cell or in a membrane or organelle, with experimental evidence of
that location, and no surface evidence:
- `N-int`: cytosol, nucleus, mitochondrion (experimental GO codes or proteomics of a purified
  fraction).
- `N-sec`: secretory pathway, vacuole or membrane without a wall or extracellular term.
- `PM-TM`: plasma membrane with transmembrane helices and no GPI evidence.
A protein with both a surface and an internal term is `ambiguous` and is **not** a control.
An unlabelled protein is never a negative.

## Size target

- At least 20 clusters of each class per species; about 61 positive clusters for a sensitivity
  half-width of 0.10 at sensitivity 0.8. Expect to fall short for *Coccidioides*. Report the gap.
- Species in order: *C. immitis* RS (taxon check with names.dmp), *C. posadasii*, *Histoplasma
  capsulatum*, *Blastomyces dermatitidis*, *Paracoccidioides brasiliensis*, *Uncinocarpus reesii*.
  Each species is its own entry. Do not pool.

## Where to start

- `analysis/step1_compare/01_extract_go_truth.py` and `labels.py` (the D1 rule that turns GO into
  labels), `analysis/step1_compare/README.md`.
- GOA files per species (EBI) and UniProt. `docs/step1-plain-language-summary.md` for the label
  definitions.
- RS proteome with 9,910 proteins: `_workdir/cocci_spherule/ref/GCF_000149335.2_ASM14933v2_protein.faa.gz`.
  The RS ranking IDs `CIMG_*` equal the RefSeq annotation. Do not use the Fungi_5k RS file (a
  different annotation).
- `analysis/cocci_spherule/` for the spherule table (do not use expression as a label).

## Deliverables

As in `COMMON-RULES.md`, in `data/controls/onygenales-step1/`. Add the column `go_evidence_code`
and `annotation_source` to `controls.tsv`. The sequences must come from the **same annotation** as
the IDs. Report how many controls have a protein that differs between annotations.

## Acceptance checks

1. Every positive and negative has a source line and a quote.
2. No IEA or homology-code evidence.
3. Counts per species and class for rows, sequences and clusters.
4. A check that no control was chosen because R0 calls it (say how you checked).

## Do not

- Do not use SignalP, TMHMM or any predictor to choose a label.
- Do not copy labels from orthologs in yeast.
- Do not use an Onygenales protein as a negative because it is "not secreted" in a prediction.
