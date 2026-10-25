# AGA1 class: can a features-only classifier beat the hand rules? (2026-10-25)

Code `analysis/aga1_class/feature_classifier.py`; tests `tests/aga1_class/test_feature_classifier.py`; all tables, the PCA plot and
hashes in `docs/reports/data/aga1_feature_classifier_2026-10-25/` (README there). Everything runs locally; nothing was sent to
HPCC. Proteins are keyed by UniProt accession; no sequence or repeat-unit text is committed.

## Short answer

A small model built only from sequence composition, predictor calls, disorder, repeat and stalk features separates the 35 tier A
proteins from the 561 negatives well (ROC-AUC about 0.97, interval roughly 0.93 to 0.99) when no cluster is shared between
training and test. At the same false-positive rate it finds somewhat more tier A proteins than the broad hand rule, but the
intervals overlap, so **this study does not show an improvement over the rules**. Almost all of the signal sits in amino-acid
composition (low hydrophobic and charged content, high Thr/Ser, low entropy). Predictor calls, disorder, repeats and stalk
features add little once composition is present. The result separates the labelled sets; it does not show that the model finds new
wall proteins.

## What was done

- **Feature table** (`features.tsv`, 850 proteins x 30 features, 8 groups): composition (11), length (1), predictors (SignalP yes
  and probability, PredGPI yes, GPI confidence class; GPI NA for the 3 proteins of 40 aa or less stays missing with an indicator,
  not 0) (5), disorder from AIUPred (fraction, longest run, mean, Ser/Thr fraction inside disordered residues) (4), repeats from
  detector 14 (coverage, copies, z; the 18 proteins under 80 aa are missing with an indicator) (4), stalk and C-tail (3), N-glycan
  sequon density N-x-S/T, x not Pro, per 100 aa (1), amino-acid entropy (1). Length, longest run, stalk run and repeat copies are
  log-transformed inside the pipeline.
- **Labels** (owner decisions): positives = tier A (35 proteins, 19 clusters); negatives = hard_negative (172) + background (389) =
  561. Tier B (32), excluded roles, adhesin_other_mechanism, kiss_row_positive and uncurated_repeat_bearing are never fitted on and
  not in headline metrics; they are scored descriptively, by models whose training part did not contain their cluster. SED1 and
  DAN4 are negatives in the headline; a set-aside run drops them. Cluster-less proteins are singletons.
- **Models**: L2 logistic regression (C 0.1, balanced), random forest (300 trees, balanced), gradient boosting (100 trees, depth
  2). Median imputation and standardisation live in the Pipeline, so they are refit in every fold.
- **Validation**: 5-fold StratifiedGroupKFold by MMseqs cluster, repeated for 10 seeds (a test asserts no cluster crosses). AUC is
  pooled out-of-fold, averaged over seeds, with a 2000-draw cluster bootstrap 95% interval. Recall at 1%, 3.5% and 5% false-positive
  rate uses thresholds read from an inner 3-fold grouped CV on the training part only; the achieved false-positive rate on the test
  part is reported beside each recall. Hand rules: broad = SignalP yes and GPI yes and longest disordered run >= 200; narrow = broad
  and the project repeat call (coverage >= 0.25 and copies >= 2.5; missing counts as not called).

## Headline numbers (cluster-bootstrap 95% intervals)

| method | ROC-AUC | recall at FPR 1% | recall at 3.5% | recall at 5% | achieved FPR (1% / 3.5% / 5% targets) |
|---|---|---|---|---|---|
| logistic regression | 0.973 [0.936, 0.992] | 0.623 [0.403, 0.776] | 0.851 [0.643, 0.969] | 0.869 [0.670, 0.978] | 1.4% / 3.7% / 5.4% |
| random forest | 0.975 [0.945, 0.990] | 0.551 [0.253, 0.741] | 0.774 [0.569, 0.904] | 0.849 [0.656, 0.957] | 1.3% / 3.6% / 5.0% |
| gradient boosting | 0.965 [0.913, 0.988] | 0.586 [0.388, 0.714] | 0.794 [0.645, 0.885] | 0.874 [0.779, 0.934] | 1.3% / 4.6% / 6.8% |
| broad rule (no threshold) | | recall 0.657 [0.350, 0.844] at FPR 2.7% [1.3, 4.2] | | | |
| narrow rule (no threshold) | | recall 0.400 [0.111, 0.576] at FPR 0.4% [0.0, 0.9] | | | |

Tier B (never trained on), recall at the same thresholds: logistic regression 0.371 / 0.713 [0.590, 0.844] / 0.811; random forest
0.336 / 0.662 [0.543, 0.830] / 0.770; gradient boosting 0.444 / 0.669 [0.538, 0.854] / 0.713; broad rule 0.562 [0.423, 0.750];
narrow rule 0.125 [0.025, 0.333]. Across seeds the AUC spread is 0.941-0.982 (logistic), 0.945-0.975 (forest), 0.952-0.971 (boosting).

Set-aside run (SED1 and DAN4 left out): AUC 0.982 [0.959, 0.994], 0.977 [0.952, 0.991], 0.978 [0.960, 0.989]; recall at 3.5%
0.837, 0.746, 0.783; the broad rule's false-positive rate falls to 2.5%. The two proteins are scored as strong positives in the
headline run (logistic scores 0.99 for DAN4 and 0.97 for SED1; forest calls both at 3.5%), as expected for Ser/Thr-rich,
GPI-anchored wall proteins, so they cost the headline some precision.

## Comparison with the hand rules

Paired cluster-bootstrap differences (`paired_vs_rules.tsv`), model minus broad rule:

- At the 1% target the model's false-positive rate is below the rule's, and recall is not different (logistic -0.03 [-0.22, 0.21]).
- At the 3.5% target, logistic regression has +0.19 [0.05, 0.40] more recall with a false-positive rate +1.0 point higher
  [-0.5, 2.6]. Random forest (+0.12 [-0.06, 0.34]) and boosting (+0.14 [-0.03, 0.38]) include zero.
- At 5% the models flag more negatives (+2.4 to +4.2 points, intervals exclude zero), so recall is not comparable.
- Read at the broad rule's own false-positive rate from the pooled out-of-fold scores (`matched_fpr_descriptive.tsv`, optimistic,
  descriptive): 0.86 (logistic), 0.71 (forest), 0.69 (boosting) against 0.66 for the rule.

I do not claim an improvement: the intervals overlap for the two tree models, the logistic gain at 3.5% is one of nine comparisons
with a slightly higher false-positive rate, and the rule thresholds (200 aa, 0.25 / 2.5) were chosen earlier by looking at these
same proteins, so the rules are also not blind. The narrow rule is a high-precision, low-recall rule (0.4% false positives); the
models cannot be compared with it at equal false-positive rate except near the 1% target, where logistic is +0.22 [0.00, 0.48]
with +1.0 point more false positives.

## Held-out clade (descriptive, wide intervals)

Decision: hold out all of one coarse clade (its positives and negatives) and train on the rest, thresholds from inner grouped CV
on the training clades (`clade_holdout.tsv`).

- Saccharomycotina held out (train 14 positives in 9 clusters; test 21 positives in 10 clusters, 417 negatives): AUC 0.967
  [0.925, 0.986] logistic, 0.973 [0.918, 0.992] forest, 0.873 [0.713, 0.931] boosting. The thresholds do not transfer: the
  realised false-positive rate at the 3.5% target is 12% (logistic), 14% (forest) and 36% (boosting) because yeast background
  proteins score higher than the non-yeast negatives used to set the cut. Recall of 0.93 at that point is therefore not a fair
  number. The broad rule on the same proteins: recall 0.81 [0.46, 0.96] at FPR 2.6%.
- Non-Saccharomycotina held out (train 21 positives in 10 clusters; test 14 positives in 9 clusters, 144 negatives): AUC 0.945
  [0.868, 1.000], 0.950 [0.901, 0.986], 0.969 [0.932, 0.994]; recall at the 3.5% target 0.71 / 0.37 / 0.67 with false-positive
  rates 1.4% / 3.1% / 4.7%. The broad rule: recall 0.43 [0.09, 0.80] at FPR 2.8%.

Ranking carries over to a new clade (AUC above 0.87); fixed thresholds do not.

## Principal components

PCA on standardised features of the 596 training-eligible proteins (others projected; `pca_*.tsv`, `pca_roles.png`).

- PC1 explains 29% of the variance (PC2 11%, PC3 9%, first ten 81%). It is the "Ser/Thr-rich, long disordered stalk,
  GPI-anchored" axis: highest loadings are Ser+Thr fraction (0.29), stalk run length (0.27), disorder fraction (0.27), disorder mean
  (0.26), Thr (0.25), GPI call and confidence (0.23), entropy (-0.24, low complexity), hydrophobic fraction (-0.19). By group,
  composition 29%, disorder 25%, stalk 17%, predictors 13%.
- PC2 is mostly size and surface-versus-globular composition: length (0.37), hydrophobic (-0.32), entropy and charged
  residues (0.28, 0.25), signal peptide (0.20). PC3 is repeat content (coverage, copies, z, about 0.36) with a signal-peptide
  term of the opposite sign.
- Tier A separates mainly on PC1 (median 6.3 against -1.7 for background proteins and 0.09 for hard negatives; PC1 alone gives
  AUC 0.94). Tier B, which was projected and never fitted, sits next to tier A (median 5.2). A tail of hard negatives, the
  adhesin_other_mechanism proteins and the projected surface/adhesion-phenotype exclusions overlaps the lower tier A edge. The kiss
  rows and the indirect regulators sit with the background. PC1 vs PC3 shows a repeat-rich subgroup of tier A far from nearly all other
  group.

## Which feature groups carry the signal

Ablation (3 seeds, `ablation.tsv`) and grouped permutation importance (`permutation_importance.tsv`):

| group | alone, AUC (logistic / forest / boosting) | drop it, AUC | permutation AUC drop |
|---|---|---|---|
| composition (11) | 0.973 / 0.969 / 0.974 | 0.911 / 0.938 / 0.897 (from 0.976 / 0.976 / 0.957) | 0.19 / 0.07 / 0.18, intervals exclude 0 |
| disorder (4) | 0.868 / 0.865 / 0.835 | no change | about 0 |
| repeats (4) | 0.854 / 0.820 / 0.801 | no change (boosting 0.976 without) | 0.011 / 0.016 / 0.027, small, intervals exclude 0 |
| stalk (3) | 0.811 / 0.852 / 0.804 | no change | about 0 |
| amino-acid entropy (1) | 0.851 / 0.817 / 0.830 | no change | about 0 |
| length (1) | 0.841 / 0.803 / 0.822 | no change | about 0 |
| predictors (5) | 0.704 / 0.607 / 0.671 | no change | about 0 |
| sequon density (1) | 0.53 / 0.43 / 0.44 (no signal) | no change | about 0 |

Composition alone matches the full model. Single-feature AUCs (`single_feature_auc.tsv`) point to Thr fraction (0.88 over all, 0.97
among proteins of 500 aa or more), hydrophobic fraction (low in tier A; 0.08 means 0.92 inverted), charged fraction (0.91
inverted at 500 aa or more) and entropy (0.86 inverted overall, 0.94 at 500 aa or more). The other groups are largely redundant
with composition (Ser/Thr-rich disordered stalks are composition), which is why dropping any of them does nothing. That the
predictors group is weak alone is expected: signal peptide and GPI are common among the cell-surface hard negatives.

Length proxies: length alone gives AUC 0.85 (tier A median 973 aa against 395 for all negatives) and the longest disordered run
0.89, but entropy (Spearman 0.45 with length among negatives), hydrophobic fraction (-0.38), longest run (0.38) and Asn fraction
(0.34) are partly length-linked. Composition features do not collapse among proteins of 500 aa or more (their single AUCs
rise), so they are not length proxies.

Length check (`length_check.tsv`): the model beats both one-feature baselines. All proteins: logistic AUC minus length-alone 0.122
[0.041, 0.255] and minus longest-run-alone 0.079 [0.023, 0.171]; forest and boosting give the same picture. Among proteins of 500 aa
or more (27 tier A in 15 clusters, 166 negatives; models retrained on that subset): logistic 0.964 [0.909, 0.988] against length
0.871 (difference 0.093 [0.026, 0.215]) and longest run 0.869 (0.095 [0.028, 0.228]); forest similar; boosting 0.910, difference
0.039 [-0.052, 0.184], so boosting does not clearly beat length there.

## Excluded roles (descriptive only; fraction called at the 3.5% target, logistic / forest / boosting)

adhesin_other_mechanism (49) 0.31 / 0.32 / 0.31; uncurated_repeat_bearing (7) 0.31 / 0.24 / 0.18; excluded_surface_other_adhesion_phenotype
(36) 0.006 / 0.03 / 0.04; excluded_indirect_regulator (60) 0.001 / 0.006 / 0.04; kiss_row_positive (70) 0.014 / 0.012 / 0.009. The
model treats the repeat-bearing other-mechanism adhesins like tier A, which the owner may or may not want; the excluded phenotype
groups and the kiss rows are scored as ordinary proteins.

## Limitations

- Only 19 tier A clusters and 35 positives. The intervals are wide (recall intervals about 0.3 wide) and the pooled out-of-fold AUC
  mixes models from different folds.
- Fixed-FPR thresholds come from inner CV on the training part, so they were not read off the test set. Two tables are optimistic
  and labelled so: the matched-FPR recall and the ablation recall (thresholds from the pooled out-of-fold scores). The hand rules are
  not blind either (their cut-offs were set looking at these proteins).
- Hard negatives are cell-surface proteins by construction and background negatives are assumed, not verified; the real false-positive
  rate on whole proteomes is unknown (the earlier S288C finding of about 12% adhesin precision for the old model warns that wall
  glycoproteins are the difficult negatives; SED1, DAN4 and several PGA proteins are the top-scoring negatives here).
- The result shows that the labelled sets separate, not that the model finds new wall proteins. Tier A proteins in one cluster
  (for instance the nine ALS proteins) count as one cluster in the bootstrap but are 9 of 35 positives, so the headline leans on
  a few families.
- Seed-averaged out-of-fold scores are used for AUC and recall; thresholds transfer poorly across clades (shown above).
- ESM-2 embeddings are a separate comparator and were not used here. Composition-based signal may partly describe sequence
  families rather than function.

## Next steps

1. Apply the frozen logistic-regression pipeline to the full proteomes behind the sorting-hat truth sets (needs SignalP, PredGPI,
   AIUPred and detector 14 run on them) and count what it flags beyond the rules; review a sample by hand. The CLI
   (`feature_classifier.py apply`) fits on all tier A plus negatives and scores a new feature table.
2. Choose a threshold per clade or recalibrate on a proteome-level negative set, given the clade-shift result.
3. Decide how the adhesin_other_mechanism proteins should be handled before any model goes into the orchestrator.
4. Compare with the ESM-2 embedding classifier on the same cluster folds and the same 561 negatives.
5. Add more tier A clusters (non-Saccharomycotina) before relying on any recall figure.
