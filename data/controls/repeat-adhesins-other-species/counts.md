# Counts for task 08: repeat-mechanism controls outside Saccharomycotina

*Generated 2026-10-06 by agent-C. Clustering was with MMseqs2 17-b804f (`easy-cluster`, `--min-seq-id 0.3 -c 0.5`), 204 sequences (the 186 existing rows from `data/controls/repeat-mechanism/curation_table.agent-A.tsv` plus the 18 new rows), 204 sequences clustered, 204 cluster representatives.*

## Rows, sequences and new clusters per group

| Group | New rows | New sequences | New clusters | Class of new rows |
|---|---|---|---|---|
| Cryptococcus | 3 | 3 | 3 | 2 adhesin (Cfl1 `other`, Cda3/MP84 `2d`); 1 negative (Cpl1, repeat_non_adhesin) |
| Botrytis | 4 | 4 | 3 (Bhp2 merged) | 1 adhesin (BcLysM1 `other`); 3 negative (Bhp1/Bhp2/Bhp3, hydrophobin_non_adhesin) |
| Talaromyces | 2 | 2 | 1 (HSP60 merged) | 1 adhesin (GAPDH `2d`); 1 negative (HSP60, moonlighting_non_adhesin) |
| Pneumocystis | 1 | 1 | 1 | 1 adhesin (PCINT1 `other`) |
| Histoplasma | 1 | 1 | 1 | 1 adhesin (YPS3 `other`) |
| Mycosarcoma | 1 | 1 | 1 | 1 adhesin (Lep1 `other`) |
| Fusarium | 1 | 1 | 1 | 1 adhesin (Hyd3 `2c`) |
| Penicillium | 1 | 1 | 0 (HfbA merged) | 1 negative (HfbA, hydrophobin_non_adhesin) |
| Blastomyces | 1 | 1 | 1 | 1 negative (Eng2, wall_hydrolase) |
| Pyricularia | 1 | 1 | 0 (MoMsb2 merged) | 1 negative (MoMsb2, mucin_sensor) |
| Blumeria | 1 | 1 | 1 | 1 negative (Lip1, wall_hydrolase) |
| Zymoseptoria | 1 | 1 | 1 | 1 negative (ZtGT2, gpi_wall_enzyme) |
| **Total** | **18** | **18** | **14 new + 4 merged** | |

## Rows, clusters by class

| Class | Rows | New clusters (independent) |
|---|---|---|
| Adhesin positive (all) | 8 | 8 |
|  `2a` repeat/avidity | **0** | **0** |
|  `2b-i` CFEM | 0 | 0 |
|  `2b-ii` Cys-knot | 0 | 0 |
|  `2b-iii` Bys1 | 0 | 0 |
|  `2c` hydrophobin | 1 | 1 |
|  `2d` moonlighting | 2 | 2 |
|  `other` (novel/integrin/LysM/chitin-binding/effector) | 5 | 5 |
| Hard negative | 10 | 6 new + 4 shared |
|  repeat_non_adhesin | 1 | 1 |
|  hydrophobin_non_adhesin | 4 | 2 new + 2 shared |
|  moonlighting_non_adhesin | 1 | 0 (shared) |
|  mucin_sensor | 1 | 0 (shared) |
|  wall_hydrolase | 2 | 2 |
|  gpi_wall_enzyme | 1 | 1 |

## 2a / stated_in_paper

- New `2a` (repeat-mediated avidity) adhesins: **0**.
- New `2a` clusters resting on `stated_in_paper`: **0** (vacuously).
- This means the repeat (`tandem_repeat_protein`) **2a** call was NOT widened by this work. Every new positive belongs to `2c`, `2d` or `other`. The only ductile-repeat adhesins already in the wider set for these clades are the pre-existing rows (BAD1 `A4D962` Ajellomyces, Msg `O74670` Pneumocystis, Mad1/Mad2 Metarhizium, MPG1 Pyricularia, Rep1 Mycosarcoma). Do NOT read the 8 new positives as evidence the repeat call works outside Saccharomycotina.

## New clusters that are shared with an existing row (flag)

Four new rows land in clusters (30% id) that already hold an existing
`data/controls/repeat-mechanism/curation_table.agent-A.tsv` row. They are not independent
clusters and, for the two with a same-cluster POSITIVE, they conflict:

- Bhp2 (Botrytis hydrophobin, negative) shares a cluster with G4MWK2 (MHP1, Pyricularia class II
  hydrophobin, `E3` adhesin). Conflict: one hydrophobin is negative, one is an E3 adhesin.
- HfbA (Penicillium hydrophobin, negative) shares a cluster with E9QT94 and P41746 (hydrophobin
  `adhesin` rows). Conflict.
- HSP60 (Talaromyces, negative) shares a cluster with P50142 (Histoplasma Hsp60 `2d` adhesin).
  Conflict: same moonlighting chaperonin, adhesive in Histoplasma, tested non-adhesive in T. marneffei.
- MoMsb2 (Pyricularia mucin sensor, negative) shares a cluster with Q5AXD9
  (`surface_other_adhesion_phenotype`, not a positive).

These 4 negatives are flagged for the owner. The `2a` calibration should exclude them (or they need
a person to decide, since the negative/positive boundary sits inside the same homology family).

## Size targets

The `estimated` status needs >=20 positives and >=20 negatives, >=20 independent clusters each,
and a 95% half-width of <=0.10 for sensitivity and specificity. This set is far from that:

- repeat-mediated (`2a`) positives: **0 new** (target 20 per clade). Not met.
- hard negatives: 10 rows (6 independent new clusters). Not met (target 20 per clade).
- No clade reaches 20 clusters.

## Groups that returned nothing (protein-level evidence)

- **Rhodotorula** – none. `Rhodotorula[Title/Abstract] AND adhesin` = 0 hits;
  `Rhodotorula AND (adhesion OR adherence OR biofilm OR flocculation)` = 29 hits (all ~0 abstracts
  with a defined protein adhesin; the only related report is strain-surface attachment / whole-cell
  flocculation harvesting, no protein named).
- **C. gattii** – none. Hits are host-cell-adhesion gene expression in dual-transcriptomes and
  antifungal compounds; no defined adhesin protein.
- **Exophiala** – none. Only whole-cell melanin/chitin-synthase/biofilm work.
- **Wallemia** – none for protein-level adhesion. Closest is the WI1 class IB hydrophobin
  (PMID 32866612) which self-assembles but shows no substrate/host adhesion; listed as a candidate.
- **Hortaea** – none. PMID 8720202 shows strong cell-surface hydrophobicity and whole-cell adhesion
  (98.5% MATH) but identifies no adhesin protein.
- **Knufia** – none. The EPS is a polysaccharide (pullulan, galactofuromannan), not a protein.
- **Colletotrichum** – no row. The 110-kDa spore-coat glycoprotein has E1 evidence (mAb UB20
  blocking, microsphere binding) but no resolvable sequence; moved to `rejected.tsv`.
- **Blumeria** – no positive adhesin; only the Lip1 (wall_hydrolase) negative.

## Leakage

All new rows are from taxa outside the tier-T5 reference proteomes and outside the Saccharomycotina
training set, so `leakage=none` was assigned. The owner should confirm none of these proteins (or a
30% homolog) was used to tune the repeat detector or the reference set; we could not verify this and
did not find any in the detector scripts' reference sources.
