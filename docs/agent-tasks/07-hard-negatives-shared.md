# Task 07. Shared hard-negative set (look-alikes that are not the target)

*Read `COMMON-RULES.md` first. Written 2026-10-06. Related: issue #14.*

## Goal

Build a curated set of **hard negatives**: proteins that look like the target in sequence or
architecture but have a characterised function that is **not** the target function. Several modules
need them: the repeat call, the family domain call, and step 1.

## Why this matters

- Random negatives make the problem easy. The first classifier reached ROC-AUC 0.995 without a
  language model, and over-called on look-alikes (precision about 12% on *S. cerevisiae* S288C).
- A false-positive rate on random negatives under-estimates the rate on look-alikes.
- Today the seed list has 31 rows (`data/curated/adhesins/hard_negative_seeds.tsv`). In
  `adhesins.tsv`, 42 `hard_negative` rows form 33 clusters; 28 of 42 are *S. cerevisiae* and 11 are
  *A. fumigatus*. Issue #14 (curated hard negatives) is open.

## Classes of look-alikes to cover

| Stratum | Description | Why it is a look-alike |
|---|---|---|
| `gpi_wall_enzyme` | GPI-anchored wall remodelling enzymes (GEL, PHR, GAS families; yapsins) | Signal peptide, GPI, Ser/Thr-rich; function is wall synthesis, not adhesion |
| `wall_hydrolase` | Chitinases, glucanases | Secreted, wall-associated; not adhesins |
| `mucin_sensor` | Mucin-like sensors (Msb2, Hkr1, Wsc family) | Ser/Thr-rich, glycosylated, surface; signalling function |
| `st_linker_enzyme` | Secreted enzymes with an S/T-rich linker | Look like adhesin repeats |
| `laccase_etc` | Multicopper oxidases (AA1) and other abundant secreted enzymes | Over-called in the Fungi_5k screen (12.5k calls) |
| `repeat_non_adhesin` | Proteins with tandem repeats and a known non-adhesive function | Repeat detectors call them |
| `domain_non_member` | Proteins with a wall-family domain but unrelated function | Family scans call them |

## Task

1. List what exists in `data/curated/adhesins/hard_negative_seeds.tsv`,
   `data/curated/adhesins/adhesins.tsv` (`cls=hard_negative`), and
   `analysis/step1_compare` (N-sec and PM-TM strata; MSB2, HKR1, SAP9 are labelled negative there).
2. For each stratum, find candidates in **non-yeast** species: *Aspergillus*, other Pezizomycotina,
   Onygenales, Basidiomycota. Each needs a source that shows the non-adhesive function (PMID and
   quote). "No adhesion evidence" alone is a weak reason: mark it `evidence_level=N2`. A
   characterised non-adhesive function is `N1`.
3. Cluster with the existing adhesin positives. Remove or flag a hard negative that is in the same
   30% cluster as a positive (2 clusters hold both today).
4. Note which module each row can serve (`serves`: `step1`, `repeat`, `family_domain`).

## Size target

At least 20 clusters per stratum **and** per clade you plan to report. Expect to fall short. Report
the table of clusters per stratum and clade.

## Do not

- Do not call a protein a hard negative because a tool does not call it.
- Do not edit `adhesins.tsv`. Record changes in `manual_overrides.tsv`.
- Do not use a protein as a negative for adhesion if its function includes adhesion in another
  organism (for example moonlighting; check Hsp60 and enolase).
