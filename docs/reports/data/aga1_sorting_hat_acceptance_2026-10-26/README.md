# Acceptance of the AGA1 modules and calls on ten proteomes, 2026-10-26

The sorting-hat wrappers (`signalp`, `gpi`, `disorder`, `repeat14`, `aga1_model`) and the six AGA1 calls were run on the tool outputs already in hand for the
ten proteomes of the proteome tests, and the engine's calls were compared **protein by protein** with the flags of `analysis/aga1_class/`
(`all_proteins_flags.tsv` of `aga1_proteome_test_2026-10-25`, `aga1_nonyeast_test_2026-10-26`, `aga1_cpar_test_2026-10-26`). Result: **0 differences in 68,534 proteins.**
Modules without a table (Pfam, TM, allergen, repeat02) stay absent and give `not_assessable`; they are not part of this check.

| Species | Taxon | Proteins | Broad flag (engine / analysis) | Narrow flag | Model call | Differences |
|---|---|---|---|---|---|---|
| Calb_SC5314 | 237561 | 6041 | 25 / 25 | 6 / 6 | 30 / 30 | 0 |
| Caur_UniProt | 498019 | 5433 | 14 / 14 | 5 / 5 | 15 / 15 | 0 |
| Cgla_CBS138 | 284593 | 5200 | 25 / 25 | 4 / 4 | 45 / 45 | 0 |
| Kmar_FIM1 | 4911 | 5070 | 10 / 10 | 2 / 2 | 22 / 22 | 0 |
| Scer_S288C | 559292 | 6067 | 11 / 11 | 5 / 5 | 12 / 12 | 0 |
| Afum_Af293 | 330879 | 9647 | 7 / 7 | 3 / 3 | 13 / 13 | 0 |
| Cimm_RS | 246410 | 9712 | 8 / 8 | 2 / 2 | 19 / 19 | 0 |
| Cneo_H99 | 235443 | 7427 | 2 / 2 | 0 / 0 | 1 / 1 | 0 |
| Rtor_NP11 | 1130832 | 8160 | 2 / 2 | 1 / 1 | 1 / 1 | 0 |
| Cpar_CDC317 | 578454 | 5777 | 28 / 28 | 9 / 9 | 23 / 23 | 0 |

(`broad` = `ser_thr_gpi_wall_candidate`, `narrow` = `flocculin_like_candidate`, `model` = `model_surface_candidate`, all at the default gate; the analysis side is `broad_flag`,
`narrow_flag` and "at least 2 of 3 models and a signal peptide". `acceptance_counts.tsv` holds the same table.)

## Commands (local, no SLURM)

```
python3 -I analysis/sorting_hat_run/aga1_acceptance.py --taxdump nodes.dmp --out-root OUT \
  --set yeast5 <proteome-test-in>/chunks <proteome-test-in>/proteins.tsv <proteome-test-merged> <scratch>/rep14_proteomes.tsv <scratch> docs/reports/data/aga1_proteome_test_2026-10-25/all_proteins_flags.tsv \
  --set nonyeast4 <nonyeast-test-in>/chunks ... docs/reports/data/aga1_nonyeast_test_2026-10-26/all_proteins_flags.tsv \
  --set cpar <cpar-test-in>/chunks ... docs/reports/data/aga1_cpar_test_2026-10-26/all_proteins_flags.tsv
python3 -I analysis/sorting_hat_run/aga1_status_truth.py --scoping docs/reports/data/aga1_class_scoping/proteins.tsv --proteome-proteins <proteome-test-in>/proteins.tsv \
  --species Scer_S288C --out-class truth_class_Scer_S288C.tsv --out-gpi truth_uniprot_gpi_Scer_S288C.tsv        # and Calb_SC5314
# status files, per species (taxon 4932 and 5476), in the species work directory of the acceptance run:
cellsurface_sorting_hat_calibrate truth --workdir OUT/<species>/wd --call-status --calls-long OUT/<species>/out/calls.long.tsv.gz --call gpi_long_disorder_protein \
  --truth truth_class_<species>.tsv --calibration-set aga1_850_class_<species> --taxa <4932|5476> --nodes-dmp nodes.dmp --leakage tuned_on_truth --negative-basis assumed --notes "..."
cellsurface_sorting_hat_calibrate truth --workdir OUT/<species>/wd --module gpi --calls-long OUT/<species>/out/calls.long.tsv.gz --call gpi_anchored_protein \
  --truth truth_uniprot_gpi_<species>.tsv --calibration-set aga1_850_uniprot_gpi_<species> --taxa <4932|5476> --nodes-dmp nodes.dmp --leakage unknown --negative-basis assumed --notes "..."
```

## Status files (`docs/reports/data/sorting_hat/call_status/aga1_850/<species>/`)

Both are `smoke` (few positives; negatives assumed; leakage caps the status). They are valid for `categories.yaml` sha256 `a9dc86ce51fec8d5b68e6ea4e98da0718d9177be4c3fd4ff5e100c63c47f8239` only
(a test fails when the packaged config changes, because the files are then stale and must be re-measured).

* `calls/gpi_long_disorder_protein.json` (a call status): *S. cerevisiae* (taxon 4932) 8 tier A against 265 negatives, sensitivity 0.75 [0.33, 1.00], specificity 0.985 [0.956, 0.996];
  *C. albicans* (5476) 10 against 117, sensitivity 0.90 [0.00, 1.00], specificity 0.940 [0.868, 0.981]. **It measures agreement of the property combination (GPI anchor and a disordered run of 200 or more residues)
  with the AGA1 class labels, not that every such protein is in the class.** The 200-residue cut was chosen on these proteins (leakage `tuned_on_truth`); the confidence intervals are wide.
* `gpi.json` (a module status, read by `gpi_anchored_protein`): agreement of the PredGPI call with the UniProt GPI-anchor feature on the labelled 850 proteins that are in the proteome:
  *S. cerevisiae* 58 positives, sensitivity 0.966, specificity 0.965; *C. albicans* 93 positives, sensitivity 0.903, specificity 0.962. Absence of the feature is an assumed negative, the labelled set is
  enriched for wall proteins, and whether these proteins were in PredGPI's training data is unknown (leakage `unknown`).
* Not measured: `long_disordered_protein` (no truth for the property), `model_surface_candidate` (the classifier was trained on the labelled proteins) and the composites (derived from their leaf calls).

## Inputs not in the repository (sha256)

```
171b54b2621d6751eb45c3498036fc3c39ef5b628a9767c0708a415064bf7c4a  nodes.dmp (NCBI taxonomy dump downloaded 2026-10-26, taxdump.tar.gz sha256 d1e5bcb1f5184919f3a3b700ce369bfaa479bc9de79695965901c3de684134fd)
```
The merged raw tool outputs, the detector-14 tables and the score tables are those recorded in the READMEs of `aga1_proteome_test_2026-10-25`, `aga1_nonyeast_test_2026-10-26` and `aga1_cpar_test_2026-10-26`.
Taxon of each proteome: `analysis/sorting_hat_run/aga1_acceptance_taxa.tsv` (the thresholds table selects its row by lineage: 147537 Saccharomycotina for the six yeasts, the fallback row for the four non-yeasts).

## sha256 of the committed files (a test checks every line)

```
2ca4841fd4ff4744e601eb621fb96b3d27685377d0ad2841ae9ac9d6c258a38f  acceptance_counts.tsv
944e73be9f412410597e86d9f0f38cbef4df731c535e2f0943fef46e7576d6cc  truth_class_Calb_SC5314.tsv
69978932b25f221748d68f5f7bd17255c6c62453f034dffb61ff13a667590583  truth_class_Scer_S288C.tsv
61d2ec1de32b6475536dc225c649f2fe9a4a3d22092c77e12947cd2d4a186e8d  truth_uniprot_gpi_Calb_SC5314.tsv
1a05aeb3e43994f1e701932bdcd7332cc3fa89728074f36c009502b8bdb344de  truth_uniprot_gpi_Scer_S288C.tsv
```

**Update 2026-10-26 (internal-TM rule):** the call `gpi_long_disorder_protein` now also requires no transmembrane helix other than the C-terminal one
(`docs/reports/data/aga1_internal_tm_2026-10-26/README.md`). Its two status files were re-measured for the new config (sha256 above): sensitivity is unchanged
(Scer 0.75, Calb 0.90), Scer specificity is 0.985 (one false positive, P41809, is removed), Calb is unchanged. The acceptance equality test above (engine calls equal the analysis flags for 68,534 proteins)
holds for the config of that day; with the new rule the engine's broad flag differs from the analysis flag for the proteins listed in the new record.
