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
| D. Detector divergence floor | 14-18 | At what unit-to-unit identity does the detector in 02 stop finding a repeat? Does a similarity-scored detector move that floor? | `REPORT_2026-09-29_repeat_detector_divergence.md` |

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
| `01_signalp.sh` | SignalP 6 on the 5 UArizona long-read proteomes and the RS and Silveira references | `sbatch analysis/cocci_repeats/01_signalp.sh` (short_gpu) | `signalp/<proteome>/prediction_results.txt` **in this folder** (gitignored), plus the tracked `signalp_summary.tsv` |
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
- **Reinterpreted 2026-09-30** (`REPORT_2026-09-30_anchored_family_search.md` section 5).
  This previously said only 205 of 489 models are >= 250 aa and "the rest are fragments",
  attributed to short-read truncation. An anchored search that classifies each hit shows
  **81% (383/474) cover the reference completely** and are genuine alleles with a whole number
  of units missing; only 19% are model problems (65 truncated, 26 split). The 178 complete
  alleles below 250 aa are **3- and 4-unit alleles, not fragments**. The 205 figure is right;
  the "mostly fragments" reading is not supported. Copy number is still not safe from
  assembly collapse -- a collapsed array looks like a short allele.
- **Corrected 2026-09-30.** This section previously said VFC140, Cpos1038 and Cpos3700 have a
  SOWgp locus with *no protein call*. **All three do have SOWgp protein models** in their own
  proteomes; the earlier search went through the pangenome orthogroups and never searched the
  long-read proteomes for the anchor directly. Run `20_sowgp_anchor_search.py` to reproduce.
  - VFC140: one model, `VFC140_004201-T1`, 234 aa, **99.6% identity over the whole of SOWgp58**
    with 2 anchored units. A short allele, not a failure. It sits below the 277 aa / 3-unit
    minimum previously recorded for *immitis*.
  - Cpos1038: two consecutive locus tags covering **overlapping halves** (aa 1-164 and
    aa 127-328) — one gene called as two. This is the gene-model failure originally predicted,
    but it yields two short calls, not zero.
  - Cpos3700: a 202 aa model (aa 127-328) plus a separate split pair `010247`/`010248`
    (aa 1-117, aa 221-328). Two loci or one split across contigs is **not established**.
  - The tblastn work in the 2026-09-29 report section 5 stands; only the "no protein call"
    claim was wrong.

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

## 4b. Part D: detector divergence floor and a family-agnostic detector (14-18)

`02_repeat_profile.py` matches residues exactly. That must fail once repeat units diverge,
and its fractional `n_copies` is 0.2-0.6 lower than the anchored counts (section 4). Part D
measures both, and tests a replacement. `02_repeat_profile.py` is unchanged: 03, 05 and 09
still read its output.

| Script | What it does | Run | Main outputs |
|---|---|---|---|
| `14_repeat_detect_general.py` | The new detector. BLOSUM62 similarity instead of exact match; a z-score of the chosen period against unrelated periods as the significance test (applied over the whole protein and again inside the region); integer copy counts from a consensus PSSM scanned back over the protein, with no family-specific motif. Importable (`detect(seq, mode=)`) and a drop-in CLI. | `14_repeat_detect_general.py <fasta…> --out <tsv>` | a TSV whose columns are a superset of 02's |
| `15_repeat_benchmark.py` | Builds the benchmark: 74 real SOWgp alleles (truth from the `PTDCYGDC` anchors), 2880 synthetic arrays over a divergence sweep, 800 non-repeat controls, and the 41 class-2a candidates | `/usr/bin/python3.12 15_repeat_benchmark.py` | `repeat_benchmark.fa`, `repeat_benchmark_truth.tsv` |
| `16_repeat_benchmark_eval.py` | Scores three arms (old, new with exact matching, new) on the benchmark | `/usr/bin/python3.12 16_repeat_benchmark_eval.py` (~4 min) | `repeat_benchmark_calls.tsv`, `repeat_benchmark_metrics.tsv`, `repeat_benchmark.png` |
| `17_repeat_detect_real.sh` | Runs 14 on the same 71,044 proteins 02 was run on, 7 proteomes in parallel in one job | `sbatch 17_repeat_detect_real.sh` (short, ~5 min) | `repeat_general_longread.tsv`, `repeat_general_reference.tsv` |
| `18_repeat_real_compare.py` | Old calls against new calls on the real proteins; prints the gained and lost calls with region complexity for inspection | `/usr/bin/python3.12 18_repeat_real_compare.py` | `repeat_real_compare.tsv` |
| `19_zthreshold_probe.py` | **Review probe (2026-09-30).** Calibrates 14's `Z_MIN` against its true null. The threshold is applied to the best of ~77 periods, so the per-protein null is a maximum, not a single N(0,1) draw: **1.3% of random sequences clear z >= 4**, against 0.02% per period. Prints only. | `/usr/bin/python3.12 19_zthreshold_probe.py` (~90 s) | stdout |
| `20_sowgp_anchor_search.py` | Direct `PTDCYGDC` search of all 7 proteomes, each hit aligned to SOWgp58. Found the SOWgp models the orthogroup search missed (see the 2026-09-30 correction in section 3). | `/usr/bin/python3.12 20_sowgp_anchor_search.py --out sowgp_anchor_hits.tsv` | `sowgp_anchor_hits.tsv` |
| `21_pfam_sets.py` | Builds 4 protein sets for the Pfam test: gained / shared / lost / **length-matched uncalled background** | `/usr/bin/python3.12 21_pfam_sets.py` | `pfam_sets.fa`, `pfam_sets.tsv` |
| `22_pfam_hmmsearch.sh` | Pfam-A at the gathering threshold (`--cut_ga`); also dumps model lengths | `sbatch 22_pfam_hmmsearch.sh` | `pfam_sets.domtbl`, `pfam_models.tsv` |
| `23_pfam_confirm.py` | **Independent confirmation of the 14 calls.** Gained are 13x enriched for Pfam repeat families vs background (31.7% vs 2.4%, P=1e-25); detector period matches Pfam hit spacing (median +0.0 aa). See `REPORT_2026-09-30_pfam_confirmation.md`. | `/usr/bin/python3.12 23_pfam_confirm.py` | `pfam_confirmation.tsv`, `pfam_period_check.tsv`, `pfam_repeat_families.tsv` |
| `30_anchor_family_search.py` | **General motif-anchored family search.** Derives the anchor from the units, measures it against a per-protein shuffle null, and classifies hits as full_length / short_allele / truncated / split. Membership decided by alignment (>=60% id), not by the anchor. | `./30_anchor_family_search.py family --seed sowgp_seed.fa --proteomes longread --out sowgp_anchor_family.tsv` | `sowgp_anchor_family.tsv` |
| `31_anchor_pangenome.sh` | The same search over 493 pangenome proteomes. **SOWgp prevalence 90.5%, union with orthogroups 92.3% — the 92% is NOT an annotation artifact.** See `REPORT_2026-09-30_anchored_family_search.md`. | `sbatch 31_anchor_pangenome.sh` | `sowgp_anchored_pangenome.tsv` |
| `32_class2a_anchor_sweep.py`, `33_anchor_sweep.sh` | Class-2a family sweep. **Written but not run** (the earlier partial output is pre-member-filter and was discarded). | `sbatch 33_anchor_sweep.sh` | — |
| `34_newfam_pangenome.sh` | Anchored search for the PTGIPTEWP family (type protein `CIMG_04070`) over 493 pangenome proteomes. 99.4% prevalence. See `REPORT_2026-09-30_ptgiptewp_family.md`. | `sbatch 34_newfam_pangenome.sh` | `newfam_anchored_pangenome.tsv` |
| `35_unit_genomes.py` | Genome list for the unit search: 6,313 proteomes with taxonomy | `srun -p short -c 2 --mem 4G -t 15 /usr/bin/python3.12 35_unit_genomes.py` | `unit_genomes.tsv` |
| `36_unit_calibrate.py`, `.sh` | Null (shuffled proteomes) and sensitivity calibration of `sowgp_unit.hmm` | `sbatch 36_unit_calibrate.sh` | `unit_calibration_null.tsv`, `unit_calibration_sensitivity.tsv` |
| `37_unit_search.sh` | HMM, anchor and architecture searches over 6,313 proteomes. The unit HMM scores no hit outside *Coccidioides* above 19.7. See `REPORT_2026-09-30_sowgp_architecture_relatives.md`. | `sbatch 37_unit_search.sh` | `unit_hmm_hits.tsv.gz`, `unit_anchor_hits.tsv.gz`, `unit_architecture_onygenales.tsv.gz` |
| `38_unit_bad1_control.py`, `.sh` | BAD1 control. **Written, not yet run.** | `sbatch 38_unit_bad1_control.sh` | — |
| `39_unit_distribution.py`, `.sh` | Copy-number summary. **Written, not yet run.** | `sbatch 39_unit_distribution.sh` | `unit_copynumber.tsv`, `unit_hits_by_protein.tsv` |
| `40_sowgp_array_depth.py` | Read depth over the SOWgp repeat array vs unit count (needs Genotyping mosdepth output) | `/usr/bin/python3.12 40_sowgp_array_depth.py` | `sowgp_array_depth.{tsv,png}` |
| `41_sowgp_status_depth.py` | Pangenome SOWgp status (full-length / fragment / no model) vs read depth at CIMG_04613 | `/usr/bin/python3.12 41_sowgp_status_depth.py` | `sowgp_status_depth.tsv`, `sowgp_missing_model_depth.tsv` |
| `42_sowgp_gene_depth_vs_units.py` | Whole-gene depth vs unit count | `/usr/bin/python3.12 42_sowgp_gene_depth_vs_units.py` | stdout |
| `run_depth_stats.sh` | Runs 41, 42 and 40 | `./run_depth_stats.sh` | `*.log` |
| `43_longread_sowgp_miniprot.sh` | miniprot of SOWgp seed proteins to the 5 long-read assemblies | `./43_longread_sowgp_miniprot.sh` (needs `module load`) | `longread_miniprot/*.sowgp.gff` |
| `44_longread_unit_tblastn.sh` | tblastn of one RS unit to the 5 assemblies (units counted in DNA) | `./44_longread_unit_tblastn.sh` (needs `$SCRATCH`) | `longread_unit_tblastn/*.tsv` |
| `45_longread_locus_summary.py` | DNA unit count vs funannotate and miniprot models at each long-read locus | `/usr/bin/python3.12 45_longread_locus_summary.py` | `longread_locus_summary.tsv` |

Run order: 15, 16, then `sbatch 17`, then 18. Results and limits are in
`REPORT_2026-09-29_repeat_detector_divergence.md`.

## 5. Known issues

1. `06_sowgp_msa.sh` and `07_sowgp_locus_tblastn.sh` run `cd "$(dirname "$0")"`. Under
   `sbatch`, `$0` is the spooled copy of the script, so this goes to the wrong directory. Run
   them with `bash` from this directory, or change them to use `$SLURM_SUBMIT_DIR`.
   **`01_signalp.sh` had the same class of bug and was fixed on 2026-09-30**: it wrote
   `signalp/` relative to `$SLURM_SUBMIT_DIR`, so its output landed wherever the job was
   submitted from and never reached the repository. It now writes to an absolute path under
   `$PROJ_ROOT/analysis/cocci_repeats/signalp/`. The missing output is why the secreted subset
   could not be recomputed for the `14_repeat_detect_general.py` calls (see Part D).
2. `07_sowgp_locus_tblastn.sh` writes BLAST databases to `/tmp`. Inside a SLURM job, `$SCRATCH`
   is the right place.

2b. **SignalP 6 runs only on the GPU build on this cluster (found 2026-09-30).** The CPU
   builds (`signalp/6`, `signalp/6.0i`) ship an **empty** `model_weights/` directory, so
   `--mode fast` fails with `FileNotFoundError: Fast mode requires model to be installed at
   .../distilled_model_signalp6.pt`. Only `signalp/6-gpu` (= `signalp/6.0i-gpu`) has the
   1.6 GB distilled model. `01_signalp.sh` previously fell back to the CPU build when the GPU
   module was unavailable, which could not work; it now checks for the weights and aborts with
   a message instead. **There is no CPU path** — the job must wait for `short_gpu`.
   `--constraint=gpu_latest` was also dropped: it is unnecessary for a small distilled model
   and pushed the scheduled start out by about an hour on a congested queue.
3. The SEG-off tblastn runs and the DIAMOND run on the locus ORFs are not scripted.
4. Scripts 05 and 09 take species from the strain name. Two strain names disagree with the tree
   (see section 4).
5. The long-read seed rows added in 09 get their species from the sequence ID. That rule labels
   the published SOWgp58/66 as *immitis*. They are not added as rows because their sequences
   already occur in the pangenome, so this does not change any count now.

## 6. File index

Large text outputs are tracked as `.gz` copies. The scripts write and read the plain files,
which are not tracked. To use a tracked copy, decompress it here first
(for example `gunzip -k sowgp_repeat_viz.copies.tsv.gz`).

| File | From | In git |
|---|---|---|
| `README.md`, `REPORT_2026-09-29_sowgp_repeat_structure.md`, `REPORT_2026-09-29_repeat_detector_divergence.md` | documentation | yes |
| `class2a_candidates.tsv` | 03 | yes (also `.gz`) |
| `procys_rs_map.tsv`, `procys_rs_tpm.tsv` | 04 | yes (also `.gz`) |
| `repeat_profile_longread.tsv`, `repeat_profile_reference.tsv` | 02 | `.gz` only (plain files are git-ignored) |
| `sowgp_seed.fa`, `sowgp_pangenome.fa`, `sowgp_pangenome_other_families.fa`, `sowgp_pangenome.tsv`, `sowgp_pangenome_summary.tsv` | 05 | `.gz` only |
| `sowgp_seed.msa.fa`, `sowgp_pangenome.msa.fa` | 06 | `.gz` only |
| `msa_seed.log` | 06 | no |
| `sowgp_tblastn_*.out`, `sowgp_tblastn_noseg_*.out` | 07 / by hand | yes (plain) |
| `sowgp_loci_orfs.fa`, `sowgp_loci_orfs_diamond.tsv` | 08 / by hand | `.gz` only |
| `signalp_summary.tsv` | 01. Per-proteome counts: proteins profiled and proteins with a signal peptide. | yes |
| `signalp/<proteome>/prediction_results.txt` | 01. One line per protein. Read by 03. | no (gitignored; regenerable, `.gz` copy written beside each) |
| `rs.dmnd` | DIAMOND database of the RS proteome. Its timestamp matches the 04 run, but 04 writes its own database to a temp directory, so the origin is not recorded. | no |
| `sowgp_repeat_viz.copies.tsv`, `sowgp_repeat_viz.units.tsv`, `sowgp_unit_map.tsv`, `sowgp_unit_identity.tsv`, `sowgp_tree_units.tsv`, `sowgp_tree_units.stats.tsv` | 09-13 | `.gz` only |
| `sowgp_repeat_viz.png`, `sowgp_unit_map.png`, `sowgp_unit_identity.png`, `sowgp_dotplot.png`, `sowgp_tree_units.png`, `sowgp_tree_units.pdf` | 09-13 | yes |
| `repeat_benchmark.fa`, `repeat_benchmark_truth.tsv` | 15 | `.gz` only |
| `repeat_benchmark_calls.tsv`, `repeat_benchmark_metrics.tsv` | 16 | `.gz` only (metrics also plain: it is small) |
| `repeat_benchmark.png` | 16 | yes |
| `repeat_general_longread.tsv`, `repeat_general_reference.tsv` | 17 | `.gz` only |
| `repeat_real_compare.tsv` | 18 | yes (also `.gz`) |
| `logs/17_repeat_detect_real.*.log` | 17 (SLURM) | no |
| `sowgp_unit.hmm` | Profile HMM from 34 SOWgp units, model length 47. Built by `34_sowgp_unit_hmm.py` (which shares the number 34 with `34_newfam_pangenome.sh`). | yes (22 KB) |
| `sowgp_unit_set.fa`, `sowgp_unit.msa.fa`, `sowgp_unit_modal.fa` | 34. Unit set, its alignment, and the modal-length units. `sowgp_unit_modal.fa` is also read by 36, 38 and 39. | yes (plain, small) |
| `class2a_anchor_families.tsv`, `class2a_anchor_hits.tsv`, `class2a_candidates_general.tsv` | 32 / 33 | `.gz` only |
| `newfam_CIMG_04070.faa`, `newfam_anchored_pangenome.tsv` | 34_newfam_pangenome.sh | `.faa` yes (plain, small); `.tsv` `.gz` only |
| `bad1_A4D962.fa` | 38. BAD1 sequence used as the control. | yes (plain, small) |
| `sowgp_anchored_longread.tsv`, `sowgp_anchor_motifs.tsv` | **Origin not recorded.** No script in this folder names them. Probably `30_anchor_family_search.py` run by hand; not verified. | yes (plain, small) |
| `rs.dmnd` | see above | no (gitignored) |
| `unit_architecture_onygenales.tsv.gz` | 37. Repeat and composition columns (read by 38 and 39) for every Onygenales protein searched (644,399 rows, 15 MB). | **no** (gitignored on purpose; regenerate with `37_unit_search.sh`, a SLURM search over all proteomes) |

