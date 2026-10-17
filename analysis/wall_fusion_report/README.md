# Cell wall remodeling and cell fusion family report (report only, option c-i)

This folder holds a script that reports hits to Pfam models of two gene sets: cell wall remodeling and synthesis, and cell fusion. It makes no call. It writes no status file. It has no threshold. It claims no function. It follows the lectin precedent (`analysis/lectin_class_report/`, option c-i). Plan: `docs/reports/cell-wall-fusion/04-repo-coverage-and-gaps.md` section 4.

**A hit is NOT evidence of adhesion, of cell fusion or of a remodeling function.** It is a hit to a Pfam family model. Wall enzymes are known false-positive sources for surface-glycoprotein detection: in `docs/reports/cell-wall-fusion/07-model-hit-check.md` the ESM classifier calls 6 of 13 S288C wall enzymes above 0.5 (CRH1, UTR2, SCW11, CTS1, SUN4, KRE9), and none of the 13 hits an active row of `data/sorting_hat/family_table.tsv`.

No file under `src/`, `data/sorting_hat/family_table.tsv` or `categories.yaml` is read or changed.

## Files

| File | Purpose |
|---|---|
| `wall_fusion_report.py` | The report script. `--set wall_remodeling`, `cell_fusion` or `all` (default). |
| `wall_fusion_models.tsv` | The model list: `id`, `name`, `family`, `set`, `supporting_report`, `verification_status`. Reviewable. |
| `submit_wall_fusion_report.sbatch` | Optional SLURM wrapper. Not submitted by the author. |
| `results/` | One small run on the 23 proteins of `analysis/cell_wall_check/test_set.faa` (2026-10-17). |

## Models

All IDs are Pfam. Names are the Pfam names. `wall_fusion_models.tsv` says which report supports each row and how the ID was checked.

Wall remodeling and synthesis set (16 models, all marked MATCH in `05-pfam-id-verification.md`): GH72 PF03198, GH16 PF00722, GH17 PF00332, GH76 PF03663, GH18 PF00704, GH5 PF00150, chitin synthase PF01644 / PF03142 / PF08407, glucan synthase PF02364, Kre9/Knh1 PF10342, SUN PF03856, Ecm33 PF12454, Tos1 PF10287 and PF10290, X8 PF07983 (a companion domain of GH72, not an enzyme family by itself).

Cell fusion set (2 models): Prm1 PF27758 and Fig1 PF12351. These two are NOT in the 53 IDs of report 05. Their names (PRM1, Fig1) come from the Pfam-A.hmm 38.2 scan in report 07, and are not checked against InterPro. The tests pin the names.

Left out on purpose:

| ID or group | Reason |
|---|---|
| PF01522 (chitin deacetylase) | 05: CE4 mapping UNVERIFIABLE |
| PF10333 | 05: the Pig-X claim is UNVERIFIABLE (record name is Pga1) |
| PF00187, PF01607 (CBM18, CBM14) | 05: CBM and Avr4 links UNVERIFIABLE; both are in the lectin report |
| PF01476 (LysM) | in the lectin report |
| PF00018 (SH3_1) | the Fus1 SH3 hit is not a fusion family (07) |
| PF00026 (Asp, yapsins) | generic family, no wall-specific tie in the reports |
| PF06423, PF05024, PF04113, PF04114, PF06432 (GPI anchor synthesis) | ER enzymes of anchor synthesis, not wall remodeling |
| PF01822 (WSC) | a wall-integrity sensor, not a synthesis or remodeling family |
| PF00399, PF22799 and other adhesion rows | already rows of `family_table.tsv` |
| PF29664 (BtgE-BtgC), PF05390 (Kre9_KNH1_C) | seen in the 07 scan, not in the 05 verification |
| Fus1, Fus2, Kar5, Hap2, Ham proteins, Ags, Kre6, Och1 and others | no Pfam ID verified in reports 01 to 05 |

Caveats from the reports that stay with the rows: GH16, GH17, GH18 and GH5 are large families with intracellular and non-wall members. The Dfg5/Dcw1 (GH76), Bgl2 (GH17), Crh/Utr2 (GH16) and Exg1 (GH5) links come from 04 and 07, not from InterPro or CAZy records. PF01644 also matches animal hyaluronan synthase. S. cerevisiae Ecm33 has no hit to PF12454 (07 row 22). Fusion proteins are membrane proteins outside the signal-peptide population of the sorting hat.

## Outputs

In `--out-dir`:

- `wall_fusion_hits.tsv`: one row per protein and model. Columns: `proteome`, `protein`, `set`, `model_acc`, `model_name`, `family`, `label`, `seq_evalue`, `seq_score`, `domain_score`, `domain_ge_ga`, `supporting_report`. The `label` is family-only, for example `wall_remodeling_family_hit: GH72` or `cell_fusion_family_hit: Prm1`.
- `wall_fusion_counts.tsv`: per proteome and per model: `n_proteins`, `n_hit_proteins`, `best_score`, `n_hit_domain_ge_ga`.
- `wall_fusion_report.md`: the report text.

The two tables have no call or status column. The first two lines of each table start with `#` and hold the no-call notice and the false-positive note. Skip them when you read the table (for example `pandas.read_csv(..., sep="\t", comment="#")`).

A protein with several domain hits of one model is counted once (best sequence score). The cutoff is `hmmsearch --cut_ga` (sequence GA). No signal peptide, TM or other condition is applied.

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

`results/` holds the output for `analysis/cell_wall_check/test_set.faa`. The scores agree with the measured values of report 07 (for example GAS1 GH72 504.4 with X8 102.4, DFG5 GH76 530.8, DCW1 547.4, BGL2 GH17 353.2, FIG1 192.3). PRM1 gives 475.9 here and 475.5 in report 07; the cause was not checked. Ecm33 (PF12454) and the chitin synthase and glucan synthase models do not hit this set, as expected from report 07 (none of the 23 proteins is a synthase). This run shows that the script works. It is not a measurement of specificity.

## Tests

```
uv run --no-project --with pytest pytest tests/wall_fusion_report -q
```

The tests use a stub `hmmsearch` and a fake Pfam file. They need neither HMMER nor Pfam. They check parsing, grouping, the label, the no-call guarantee (no call or status column, no file other than the three reports) and that the model TSV holds only the verified IDs (committed list, no network).
