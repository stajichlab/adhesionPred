# AGA1-class features-only classification study, 2026-10-25

Report: `docs/reports/2026-10-25-aga1-feature-classifier.md`. Code: `analysis/aga1_class/feature_classifier.py`
(`features`, `analyse`, `apply`). Tests: `tests/aga1_class/test_feature_classifier.py`. Everything is keyed by UniProt accession.
No residues, FASTA, per-residue scores or repeat-unit text are in this folder.

## Provenance

`features.tsv` (850 proteins, 30 features in 8 groups, plus role, tier, clade and cluster) was built by

```
python3 -I analysis/aga1_class/feature_classifier.py features \
  --proteins docs/reports/data/aga1_class_scoping/proteins.tsv \
  --predictions docs/reports/data/aga1_disorder_2026-10-24/aga1_signalp_predgpi_aiupred.tsv \
  --stalk docs/reports/data/aga1_broad_flag_2026-10-24/stalk_features_w30_f0.5.tsv \
  --rep14 rep14_850.tsv --raw disorder_raw.txt --out features.tsv
```

Inputs that are NOT in the repository (sha256):

| input | sha256 |
|---|---|
| `rep14_850.tsv` (analysis/cocci_repeats/14_repeat_detect_general.py on the 850; the `rep_unit` column is dropped on read) | ee715a4f1e1f26927a09d13c7ede2363c0a5e450d913a6253e2e60e4660449a8 |
| `disorder_raw.txt` (raw AIUPred output with residue letters; used for composition, sequons, entropy and Ser/Thr in disordered residues) | 5810381d06e5fd89e76e9f04b4b18a7e7283428712dbb15ed90ebda9e97b60b3 |
| `aga1_proteins.faa` (staged FASTA; not read by this script, listed for the record) | 708384fbdbe4a2b24714c8547a12ef23671faff7e1bd75fd5c9cf8475ac63a0e |

Missing values are kept missing, never 0: GPI is missing for the 3 proteins of 40 aa or less (`gpi_missing` = 1), repeat features
are missing for the 18 proteins under 80 aa (`rep_missing` = 1), `st_in_disordered` is missing where no residue is disordered.
A cluster-less protein has the cluster id `SINGLETON_<accession>`.

The analysis files come from `analyse --features features.tsv --out-dir . --seeds 10 --n-boot 2000` (seeds 0-9; ablation and
permutation importance use seeds 0-2). Models are scikit-learn defaults named in `feature_classifier.py`.

## Files and sha256 (checked by a test)

| file | sha256 |
|---|---|
| `ablation.tsv` | 008dcad43fdf1146d5e0e17e68e8213ce60150a7504c221629e7f1e59426d64a |
| `clade_holdout.tsv` | 15bf17753edc7b4a6fa09b85eefdfbe86e907833938a7ccdb9655a5350df6238 |
| `extra_roles_scores.tsv` | 4de3e796a050084111bb2b1833a19310a288a457e69cac6c8a66c8e6809c118a |
| `features.tsv` | d60c70d22bdbb80fe186f185663dd8b568b77c1a5cabd4426c467650c24d0f61 |
| `headline_metrics.tsv` | 9ea1c515c1afb804962d1d60758a93299d3dba398afeb581048fbb5b1e3b84f0 |
| `headline_setaside_sed1_dan4.tsv` | c8be371381067d91a2fee67f9f86ff7b779a5e3d804ea6fe2bc99bbf65b51eb5 |
| `length_check.tsv` | 726ae7103944bb42fb8681a34b4bed3a858a39dcb6a5b69b43674f4057f84dd4 |
| `matched_fpr_descriptive.tsv` | 02de3ab298d8455f90c38f9bf4bd7afaec39c2b0a7cdb3e36839514024118c82 |
| `oof_scores.tsv` | 560f257369b6bc33f6c692623c68e94699c04b2b1d9652931dce3ced9b27be63 |
| `paired_vs_rules.tsv` | 91f4c9acb40ad6a7127498edb63e5c20b73349fee0477ebdf4e9bd0e2921bbce |
| `pca_coords.tsv` | ea5b748ae30f6cc4c4d72fce12e33c531464e41ac34b89faf28811ee2f5e7cce |
| `pca_group_loadings.tsv` | 93c360c27a1f33a2a652417dfee5f2f087039ff7b4a65f71351b2ee0dbf392d2 |
| `pca_loadings.tsv` | bbfaa8972d77aadff1dfc351bedb4c3543cb8c35a6eb48abc24d11dba7f045ef |
| `pca_roles.png` | 490b4c17c716f46fe3bfa01b261a1fc2eb141c2b47dca18ef467b2deddccad51 |
| `pca_variance.tsv` | 267c6428d97953b8b7f574e1d542e55911690556fc452dbd7feb36339b2b1a7e |
| `permutation_importance.tsv` | 47d772831a8ca5dee19cea2f7941abf3e6f43e1b15cef055def89c28120c8d06 |
| `single_feature_auc.tsv` | e475864ef7ec57c59e3a66a46e32cdccf9b9faf059bd68e591d6a7d6cc6fcfde |
| `tierB_recall.tsv` | ea6301098a5768b5c2dfdcb5d6771872cc25c78dc6e00d0d8c6beef0416b2663 |
