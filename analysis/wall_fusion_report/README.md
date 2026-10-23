# Cell wall remodeling and cell fusion family report (report only, option c-i)

This folder holds a script that reports hits to Pfam models of four model sets: `wall_remodeling` (enzymes), `wall_synthesis` (chitin and glucan synthases), `wall_structural_other` (wall models that are not an enzyme family by themselves) and `cell_fusion_candidate_models`. It makes no call. It writes no status file. It has no threshold beyond each model's Pfam GA cutoff. It claims no function. It follows the lectin precedent (`analysis/lectin_class_report/`, option c-i). Plan: `docs/reports/cell-wall-fusion/04-repo-coverage-and-gaps.md` section 4.

**A hit is NOT evidence of adhesion, of cell fusion or of a remodeling function.** It is a hit to a Pfam family model. Wall enzymes are known false-positive sources for surface-glycoprotein detection: in `docs/reports/cell-wall-fusion/07-model-hit-check.md` the ESM classifier calls 6 of 13 S288C wall enzymes above 0.5 (CRH1, UTR2, SCW11, CTS1, SUN4, KRE9), and none of the 13 hits an active row of `data/sorting_hat/family_table.tsv`.

No file under `src/`, `data/sorting_hat/family_table.tsv` or `categories.yaml` is read or changed.

## Files

| File | Purpose |
|---|---|
| `wall_fusion_report.py` | The report script. `--set all` (default) or a comma list of the four set names. |
| `wall_fusion_models.tsv` | The model list: `id`, `name`, `family`, `set`, `supporting_report`, `verification_level` (`05_interpro` or `07_hmmscan_only`), `verification_status` (note). Reviewable. The tests check it against `tests/wall_fusion_report/data/verified_pfam_ids.tsv`, extracted from reports 05 and 07. |
| `*.NOTICE` (in `--out-dir`) | A sidecar file with the no-call notice, one next to each table. | |
| `submit_wall_fusion_report.sbatch` | Optional SLURM wrapper. Not submitted by the author. |
| `results/` | One small run on the 23 proteins of `analysis/cell_wall_check/test_set.faa` (2026-10-17). |

## Models

All IDs are Pfam. Names are the Pfam names. The script reads the NAME of each model from the Pfam file and compares it to the `name` column. It stops with an error on a different name, a missing accession or a missing GA line (exit code 2).

One inclusion rule: a model is listed if it is verified in `05-pfam-id-verification.md` (level `05_interpro`: ID and name resolved in InterPro) or measured in `07-model-hit-check.md` (level `07_hmmscan_only`: seen in the hmmscan of the 23 S288C proteins, name from Pfam-A 38.2, not checked in InterPro). It is applied the same way to wall and fusion models. Set totals never mix enzymes and non-enzymes.

| Set | Models |
|---|---|
| `wall_remodeling` (8, enzymes that remodel) | GH72 PF03198, GH16 PF00722, GH17 PF00332, GH76 PF03663, GH18 PF00704, GH5 PF00150, chitin deacetylase PF01522 (Polysacc_deac_1), BtgE-BtgC PF29664 (`07_hmmscan_only`; Pfam name BtgE-BtgC_glyc-hydro; catches SCW11 266.9, SCW4 283.2, BGL2 224.2) |
| `wall_synthesis` (4) | chitin synthase PF01644 / PF03142 / PF08407, glucan synthase PF02364 |
| `wall_structural_other` (7, not an enzyme family by themselves) | Kre9/Knh1 PF10342, Kre9_KNH1_C PF05390 (`07_hmmscan_only`; KRE9 166.2), SUN PF03856, Ecm33 PF12454, Tos1 PF10287 and PF10290, X8 PF07983 |
| `cell_fusion_candidate_models` (2, `07_hmmscan_only`) | Prm1 PF27758, Fig1 PF12351 |

All others are `05_interpro` (MATCH in 05). SUN is in the non-enzyme set because its nature is unclear: the InterPro name is "Beta-glucosidase (SUN family)", but the reports show no enzyme role. X8 is a companion domain of GH72, also found in plant glucanases. PF01522: 05 verifies the Pfam ID and name (Table A: "MATCH on Pfam name"). What 05 marks UNVERIFIABLE is the CE4 mapping (Table B), not the ID. PF01522 is not measured in 07 (no test protein) and is a large family that also holds non-chitin deacetylases.

Left out on purpose (scope, not verification):

| ID or group | Reason |
|---|---|
| PF00026 (Asp, yapsins) | generic family, no wall-specific tie in the reports |
| PF10333 | the record is Pga1, an ER GPI enzyme (also the 05 Pig-X claim is UNVERIFIABLE) |
| PF00187, PF01607 (CBM18, CBM14), PF01476 (LysM) | in the lectin report |
| PF00018 (SH3_1) | the Fus1 SH3 hit is weak (23.5) and a generic domain (07 row 26) |
| PF06423, PF05024, PF04113, PF04114, PF06432 (GPI anchor synthesis) | ER enzymes of anchor synthesis, not wall remodeling |
| PF01822 (WSC) | a wall-integrity sensor, also in secreted enzymes |
| PF00399, PF22799 and other adhesion rows | already rows of `family_table.tsv` |
| Fus1, Fus2, Kar5, Hap2, Ham proteins, Ags, Kre6, Och1 and others | no Pfam ID verified or measured in reports 01 to 07 |

Caveats from the reports that stay with the rows: GH16, GH17, GH18 and GH5 are large families with intracellular and non-wall members. The Dfg5/Dcw1 (GH76), Bgl2 (GH17), Crh/Utr2 (GH16) and Exg1 (GH5) links come from 04 and 07, not from InterPro or CAZy records. PF01644 also matches animal hyaluronan synthase. S. cerevisiae Ecm33 itself has no hit to PF12454 (07 row 22, none even at E 10), so the label `Ecm33` does not mean an Ecm33 protein. Fusion proteins are membrane proteins outside the signal-peptide population of the sorting hat.

## Outputs

In `--out-dir`:

- `wall_fusion_hits.tsv`: one row per protein and model. Columns: `proteome`, `protein`, `set`, `model_acc`, `model_name`, `family`, `label`, `seq_evalue`, `seq_score`, `domain_score`, `domain_ge_ga`, `verification_level`, `supporting_report`. The `label` names the model only, as in the lectin report (for example `GH72 (PF03198) domain hit`, `Prm1 (PF27758) domain hit`). The `set` column is separate.
- `wall_fusion_counts.tsv`: per proteome and per model: `n_proteins`, `n_hit_proteins`, `best_score`, `n_hit_domain_ge_ga`.
- `wall_fusion_report.md`: the report text.

The two tables have no call or status column. The first two lines of each table start with `#` and hold the no-call notice and the false-positive note. Reading with `pandas.read_csv(..., sep="\t", comment="#")` drops those lines. So the notice is ALSO written to a sidecar file next to each table: `wall_fusion_hits.tsv.NOTICE` and `wall_fusion_counts.tsv.NOTICE`, and it is the first text of `wall_fusion_report.md`.

Repeated IDs: an ID is the first word of the FASTA header. `n_proteins` counts unique IDs. Records that repeat an ID inside a file are merged (one row per protein and model, best sequence score kept). The script prints a warning per proteome and lists it in the report; IDs shared by two proteomes also give a warning (hits stay per proteome). A protein with several domain hits of one model is counted once (best sequence score). The cutoff is `hmmsearch --cut_ga` (sequence GA). No signal peptide, TM or other condition is applied.

## Versions used (2026-10-17, `ssh hpcc` login node)

- Python: `/usr/bin/python3.12`.
- HMMER 3.4: `/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin/hmmsearch` (`module load hmmer/3.4`).
- Pfam 38.2: `/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`.

## Run

Check the queue first. Use a few proteomes on a login node.

```
/usr/bin/python3.12 analysis/wall_fusion_report/wall_fusion_report.py \
  --out-dir OUT_DIR --tmpdir TMP_DIR --cpu 4 \
  --models-tsv analysis/wall_fusion_report/wall_fusion_models.tsv \
  --hmmsearch /opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin/hmmsearch \
  [--set wall_remodeling|cell_fusion|all] [--pfam-hmm PFAM_A_HMM] \
  --fasta-list LIST_OF_FASTA_PATHS
```

Inputs: plain, `.gz` or `.zst` FASTA. The proteome name is the file name without the compression and FASTA suffixes. Two inputs must not give the same name. `TMP_DIR` must exist or be creatable. Without `--tmpdir` the script uses `$SCRATCH` and stops if it is not set. For a SLURM job see `submit_wall_fusion_report.sbatch` (all paths come from the environment).

## Check run (23 S288C proteins)

`results/` holds the output for `analysis/cell_wall_check/test_set.faa`. The scores agree with the measured values of report 07 (for example GAS1 GH72 504.4 with X8 102.4, DFG5 GH76 530.8, DCW1 547.4, BGL2 GH17 353.2, FIG1 192.3). PRM1 gives sequence score 475.9 and domain score 475.5 here. Report 07 lists 475.5 because its column is the hmmscan best-domain score (`full_pfam_ga_hits.tsv`, `best_domain_score`; the hmmscan line itself gives 475.9 and 475.5). Same HMMER 3.4 and Pfam 38.2; only the E-values differ (search space). The two scan-only wall models hit as in 07: PF29664 SCW4 283.2, SCW11 266.9, BGL2 224.4 (sequence) and 224.2 (domain); PF05390 KRE9 166.2. The run of 2026-10-17 was repeated after the review fixes and passes the model-name check against Pfam 38.2. Ecm33 (PF12454) and the chitin synthase and glucan synthase models do not hit this set, as expected from report 07 (none of the 23 proteins is a synthase). This run shows that the script works. It is not a measurement of specificity.

## Tests

```
uv run --no-project --with pytest pytest tests/wall_fusion_report -q
```

The tests use a stub `hmmsearch` and a fake Pfam file built from its own fixture (`tests/wall_fusion_report/data/fake_pfam_models.tsv`, not from the model TSV). They need neither HMMER nor Pfam. They check parsing, grouping, the label, the sets, repeated IDs (warning, unique count), `.zst` input, a model name that differs from Pfam (error), a model without a GA line (error, exit code 2), the no-call guarantee (no call or status column; only the three reports and the two NOTICE files) and the inclusion rule. The verified-ID list `tests/wall_fusion_report/data/verified_pfam_ids.tsv` is extracted from reports 05 (Table A) and 07 (with the Pfam names of `analysis/cell_wall_check/full_pfam_ga_hits.tsv`); a test parses 05 again and compares.

## Custom Ecm33 model (report only, unvalidated; added 2026-10-23)

Owner decision of 2026-10-23: the custom Ecm33 HMM is integrated into this report, REPORT ONLY, until more testing. Nothing here touches `data/sorting_hat`, `family_table.tsv`, `family_measure.tsv`, any call, tier, status or `narrow_ok`.

Why: Pfam PF12454 (Ecm33) does not reach GA on true Ecm33s (S. cerevisiae Ecm33 has no hit; see the notes of the report). The custom model is the `all` model of the ecm33 build (trained on all 481 members; job 30287890). Files: `data/custom_models/ecm33_all_2026-10-23/` (`model.hmm`, `cutoff.json`, `provenance.txt`, `SHA256SUMS`, README). model.hmm sha256 `9a6e2c812540b0b1892f081611bf4044ec674bc524d68f5c6fc3d16001c50937`.

How it is searched: the optional table `wall_fusion_custom_models.tsv` (columns `id, name, model_path, cutoff_bits, sha256, kind, use, validation, supporting_record`; `kind` must be `custom`, `use` `report_only`, `validation` `unvalidated_beyond_circular_holdout`) is passed with `--custom-models-tsv`. Each model is searched with `hmmsearch -T <cutoff_bits>` (sequence bit score, NO `--cut_ga`; the model has no GA line). The Pfam search is unchanged (`--cut_ga`). The script refuses to run if the model file's sha256 differs from the table, or if the NAME in the file differs from the table. Without the option the outputs are exactly as before.

Outputs when the option is given: a separate section "Custom models" at the end of `wall_fusion_report.md` (model sha256 and cutoff in its table), and separate tables `wall_fusion_custom_hits.tsv` and `wall_fusion_custom_counts.tsv` (each with a NOTICE). Custom hits are never merged into the Pfam tables or counts, and never pooled with them. `results/` was NOT regenerated; the section appears on the next run.

Cutoff: 29.6 bits, set from the false-hit side only under the signed freeze (`analysis/wall_fusion_truth/pr6_ecm33_freeze.json`, sha256 `a7553076d4963e4daedf3d7fdf0eba28b8e8349847b2644482fc2e854d7ce915`): 4 extra calls over the four kept proteomes (limit 16.2455), 0 hard-negative calls.

Evidence and its limits (`docs/reports/2026-10-23-pr6-ecm33-build.md`): finds all 469 held-out members in all 11 clades, and both true Ecm33s among the shared negatives (P38248, Q4WNS8); 0 hard-negative calls; 4 to 5 extra calls. BUT the held-out members come from the same phmmer search that supplied the training members (circular), the iteration comparator (jackhmmer) ties it, and no unit was testable under the u rule (u = 0 everywhere). So the model is UNVALIDATED beyond that circular test. A hit is not an Ecm33 call.

Integration: the draft `submit_wall_fusion_report.sbatch` accepts `CUSTOM_MODELS_TSV` and `CUSTOM_MODEL_HMM`, pins the model sha256 and exits (code 3) if it differs. DRAFT: nothing was submitted or run.
