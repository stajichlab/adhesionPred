# `analysis/cocci_repeats/`: repeat surface proteins in *Coccidioides*, and SOWgp

This directory holds the class-2a curation work (repeat/avidity surface proteins, see
`docs/TOOL-ARCHITECTURE.md`) for *Coccidioides*. It then follows one protein family, SOWgp,
across the 496-genome pangenome.

The work has three parts:

| Part | Scripts | Question | Written up in |
|---|---|---|---|
| A. Class-2a candidates | 01-04 | Which secreted proteins in *Coccidioides* have tandem repeat arrays? Which are expressed in spherules? | `docs/reports/2026-09-27-cocci-repeat-surface-proteins.md` |
| B. SOWgp across the pangenome | 05-08 | Is SOWgp present in every strain? Why was it not called in three long-read strains? | This README (section 3) and the report below, section 5 |
| C. SOWgp repeat structure | 09-13 | How do the repeat arrays differ in length and unit structure? Did each length arise once? | `REPORT_2026-09-29_sowgp_repeat_structure.md` |

The model retrain test that uses the class-2a candidates is in
`analysis/model_review/repeat_structure_transfer_2a_retrain.py` (report 2026-09-27, section 6a).

> **Status:** computational results only. No experimental validation.

---

## 1. Environment

- Python scripts: run with `/usr/bin/python3.12`. They need `pandas`, `numpy`, `biopython`,
  `matplotlib`.
- Modules: `signalp/6-gpu` (01), `diamond/2.1.24` (04), `mafft` (06), `ncbi-blast/2.14.1` (07).
- Site paths (`COCCI_PANGENOME`, `COCCI_LONGREAD`, `SPHERULE_RNASEQ`) come from
  `analysis/_common/paths.py` and `paths.sh`, which read `config/site.yaml`.
- Scripts 07 and 08 hardcode the long-read path
  `/bigdata/stajichlab/shared/projects/Onygenales/Coccidioides/UArizona_strains/For_Marc`.
- Script 13 hardcodes the strain tree path (see section 4).

## 2. Part A: class-2a candidates (01-04)

| Script | What it does | Run | Main outputs |
|---|---|---|---|
| `01_signalp.sh` | SignalP 6 on the 5 UArizona long-read proteomes and the RS and Silveira references | `sbatch 01_signalp.sh` (short_gpu) | `signalp/<proteome>/prediction_results.txt` (written in the submit directory; not in this folder) |
| `02_repeat_profile.py` | Periodicity-based tandem repeat detector. For each period *p*, scores positions where s[i] = s[i+p], picks the best *p*, and delimits the repeat region. Also computes composition. | `02_repeat_profile.py <fasta…> --out <tsv>` | `repeat_profile_longread.tsv`, `repeat_profile_reference.tsv` (git-ignored, regenerable) |
| `03_repeat_surface_candidates.py` | Joins repeat profile + SignalP + composition class. Thresholds: coverage ≥ 0.25, ≥ 2.5 copies. | `--profiles … --signalp-dir … --out …` | `class2a_candidates.tsv` (41 candidates) |
| `04_procys_rs_expression.py` | Maps the 20 Pro/Cys-rich candidates to *C. immitis* RS (DIAMOND) and joins spherule/mycelium TPM | on a compute node, see script header | `procys_rs_map.tsv`, `procys_rs_tpm.tsv` |

Main results (details in the 2026-09-27 report):

- 71,044 proteins profiled. 58 have a substantial tandem repeat. 41 of these are also secreted.
- Composition: 20 Pro/Cys-rich (SOWgp/BAD1 type), 10 Ser/Thr-rich (FLO/ALS type), 11 other.
- The detector recovers known periods: SOWgp 47 aa, BAD1 24 aa, ALS1 36 aa.
- Of the Pro/Cys-rich candidates, only SOWgp (`CIMG_04613`) is spherule-induced: about
  15,000 TPM at 48 h, against about 13 TPM in mycelium.

## 3. Part B: SOWgp across the pangenome (05-08)

| Script | What it does | Run | Outputs |
|---|---|---|---|
| `05_sowgp_pangenome.py` | Stage `seed`: 8 seed sequences (5 genome copies + published SOWgp58/66/82). Stage `extract`: every member of the SOWgp orthogroups in the OrthoFinder pangenome (`Cocci_496_OG2_5_5`), profiled with the 02 detector. | `--stage seed`, then `--stage extract` | `sowgp_seed.fa`, `sowgp_pangenome.fa`, `sowgp_pangenome.tsv`, `sowgp_pangenome_summary.tsv`, `sowgp_pangenome_other_families.fa` |
| `06_sowgp_msa.sh` | MAFFT `--localpair` alignment of the seed and pangenome sets | `bash -lc 'module load mafft && bash 06_sowgp_msa.sh'` | `sowgp_seed.msa.fa`, `sowgp_pangenome.msa.fa`, `msa_seed.log` |
| `07_sowgp_locus_tblastn.sh` | tblastn of the seed proteins against the VFC140, Cpos1038 and Cpos3700 assemblies (SEG on) | see script header | `sowgp_tblastn_{VFC140,Cpos1038,Cpos3700}.out` |
| (by hand) | tblastn with SEG off. Command not recorded. | — | `sowgp_tblastn_noseg_*.out` |
| `08_sowgp_loci_extract.py` | Clusters SEG-off hits into loci, extracts a 2 kb window, translates 6 frames, keeps ORFs ≥ 50 aa | `/usr/bin/python3.12 08_sowgp_loci_extract.py` | `sowgp_loci_orfs.fa` (130 ORFs) |
| (by hand) | DIAMOND of the locus ORFs against the seed set. Command not recorded. | — | `sowgp_loci_orfs_diamond.tsv` |

Main results:

- The pangenome has 489 SOWgp gene models (orthogroup OG0001318) in 449 strains, and 464
  models in a second orthogroup (OG0007110, the family of the Silveira fragment QVM10648.1).
- Only 205 of the 489 SOWgp models are ≥ 250 aa. The rest are fragments (72-249 aa).
  Short-read assemblies and gene models can truncate the repeat array.
- A SOWgp locus is present in all three long-read strains with no protein call (VFC140,
  Cpos1038, Cpos3700). The locus ORFs match only short parts of SOWgp (for example aa 1-36 at
  100% identity). The missing call is probably a gene-model failure. See the 2026-09-29 report,
  section 5.

## 4. Part C: SOWgp repeat structure (09-13)

All scripts in this part use `sowgp_units.py`. That module defines a repeat unit as the 47 aa
starting at each `PTDCYGDC` anchor, and gives each unit a class from two diagnostic sites.

| Script | What it does | Outputs |
|---|---|---|
| `sowgp_units.py` | Shared code: anchors, unit classes, distinct alleles, unit differences, colours | — |
| `09_sowgp_repeat_viz.py` | Full-length copies (≥ 250 aa) with anchored units. Architecture, unit-count bars, per-species logos, raster against the modal unit. | `sowgp_repeat_viz.png`, `.copies.tsv`, `.units.tsv` |
| `10_sowgp_unit_map.py` | One row per distinct sequence; units drawn as coloured blocks. Length against unit count. | `sowgp_unit_map.png`, `.tsv` |
| `11_sowgp_unit_identity.py` | Unit-by-unit aa difference matrices, and longer alleles against the 3-unit allele | `sowgp_unit_identity.png`, `.tsv` |
| `12_sowgp_dotplot.py` | Self dot-plots, one sequence per species and unit count | `sowgp_dotplot.png` |
| `13_sowgp_tree_units.py` | Unit arrays on the 488-strain CDS ML tree, and a parsimony test with a within-species permutation null | `sowgp_tree_units.png`, `.pdf`, `.tsv`, `.stats.tsv` |

Run 09 first. 10-13 read `sowgp_repeat_viz.copies.tsv` and can run in any order. Each finishes
in under 20 s.

Tree used by 13:
`/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Phylogeny/results/msa_filter_cds_ascomycota-buildtree/Cocci_cds.488taxa_ascomycota.fa.part.aicc.contree`

Main results (details in `REPORT_2026-09-29_sowgp_repeat_structure.md`):

- Protein length = base length + 47 aa × units. *immitis*: 277 / 324 / 371 aa = 3 / 4 / 5 units.
  *posadasii*: 281 / 328 / 375 / 422 aa = 3 / 4 / 5 / 6 units.
- Anchored counts are whole numbers and match the published copy numbers (SOWgp58: 4,
  SOWgp82: 6). The profiler's fractional `n_copies` (used in the 2026-09-27 report, section 4.2)
  is 0.2-0.6 lower.
- All three published SOWgp alleles (58/66/82) match *posadasii* sequences exactly.
- Some alleles contain identical units. Examples: *posadasii* 5-unit units 1, 3 and 4;
  *immitis* 5-unit units 3 and 4.
- In *posadasii*, unit count is not grouped on the strain tree (P = 0.38), but the flanking
  sequence is (P = 0.001). This fits repeated independent changes in unit count, but
  gene-model errors are not excluded. *immitis* cannot be tested this way.
- Two strains are in the clade of the other species: `Coccidioides_posadasii_B3476` and
  `Coccidioides_immitis_485B-1_L_OLD_CPA0023`.

## 5. Known issues

1. `06_sowgp_msa.sh` and `07_sowgp_locus_tblastn.sh` run `cd "$(dirname "$0")"`. Under
   `sbatch`, `$0` is the spooled copy of the script, so this goes to the wrong directory. Run
   them with `bash` from this directory, or change them to use `$SLURM_SUBMIT_DIR`.
2. `07_sowgp_locus_tblastn.sh` writes BLAST databases to `/tmp`. Inside a SLURM job, `$SCRATCH`
   is the right place.
3. The SEG-off tblastn runs and the DIAMOND run on the locus ORFs are not scripted.
4. Scripts 05 and 09 take species from the strain name. Two strain names disagree with the tree
   (see section 4).
5. The long-read seed rows added in 09 get their species from the sequence ID. That rule labels
   the published SOWgp58/66 as *immitis*. They are not added as rows because their sequences
   already occur in the pangenome, so this does not change any count now.

## 6. File index

| File | From | Tracked in git |
|---|---|---|
| `README.md`, `REPORT_2026-09-29_sowgp_repeat_structure.md` | documentation | no (new) |
| `class2a_candidates.tsv` | 03 | yes |
| `procys_rs_map.tsv`, `procys_rs_tpm.tsv` | 04 | yes |
| `repeat_profile_longread.tsv`, `repeat_profile_reference.tsv` | 02 | no (git-ignored) |
| `sowgp_seed.fa`, `sowgp_pangenome*.fa`, `sowgp_pangenome*.tsv` | 05 | no |
| `sowgp_*.msa.fa`, `msa_seed.log` | 06 | no |
| `sowgp_tblastn_*.out`, `sowgp_tblastn_noseg_*.out` | 07 / by hand | no |
| `sowgp_loci_orfs.fa`, `sowgp_loci_orfs_diamond.tsv` | 08 / by hand | no |
| `rs.dmnd` | DIAMOND database of the RS proteome. Its timestamp matches the 04 run, but 04 writes its own database to a temp directory, so the origin is not recorded. | no |
| `sowgp_repeat_viz.*`, `sowgp_unit_map.*`, `sowgp_unit_identity.*`, `sowgp_dotplot.png`, `sowgp_tree_units.*` | 09-13 | no |
