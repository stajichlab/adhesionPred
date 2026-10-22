# PR 6 step 6.0: the 44 missed clusters against the whole of Pfam

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Run record: `docs/reports/data/pr6_pfam_scan_2026-10-22/` (job 30235946, code `ebb29a9`, HMMER 3.4, Pfam 38.2). Spec: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md`, decision 1. Baseline: `2026-10-22-kiss-baseline-recall.md`.*

The 39-model panel missed 44 clusters of the Kiss control set (adhesion 12, wall 28, fusion 4). This step asks whether any other Pfam family hits them. Each of the 66 member sequences was searched with `hmmscan --cut_ga` against all of Pfam-A. The job and this note change no call, no status, no tier and no panel. Every number below is recomputed from the committed tables.

## 1. Result

| class | clusters |
|---|---|
| `hit_panel_model` | 2 |
| `hit_other_pfam` | 38 |
| `no_pfam_hit` | 4 |
| total | 44 |

By slice (`scan_counts.tsv`):

| slice | clusters | panel model | other Pfam | no Pfam hit |
|---|---|---|---|---|
| adhesion | 12 | 0 | 9 | 3 |
| wall | 28 | 2 | 25 | 1 |
| fusion | 4 | 0 | 4 | 0 |

By function group of the spec (the groups are annotation text, not homology):

| group | clusters | panel model | other Pfam | no Pfam hit |
|---|---|---|---|---|
| A1 hydrophobin-like | 2 | 0 | 1 | 1 |
| A2 yeast surface adhesin-like | 4 | 0 | 2 | 2 |
| A3 *Neurospora* lectins | 4 | 0 | 4 | 0 |
| A4 adhesion singletons | 2 | 0 | 2 | 0 |
| W1 glycosylation and glycan supply | 7 | 0 | 7 | 0 |
| W2 chitin synthase regulation | 5 | 0 | 5 | 0 |
| W3 beta-1,6-glucan Kre family | 5 | 0 | 5 | 0 |
| W4 alpha-1,3-glucan and other synthases | 2 | 0 | 2 | 0 |
| W5 GPI wall proteins and signalling | 6 | 2 | 3 | 1 |
| W6 processing and other | 3 | 0 | 3 | 0 |
| F1 pheromone receptors | 3 | 0 | 3 | 0 |
| F2 so (fusion) | 1 | 0 | 1 | 0 |

**The two panel hits reproduce the baseline.** KS0082 (pir1, `PF22799` PIR1-like_C) and KS0116 (cfmA, `PF05730` CFEM) are the two clusters that the baseline found hit by a family model of another slice (they were missed by the wall models of their own slice). The scan finds exactly these two and no other.

**The four clusters with no Pfam hit** (`custom_targets.tsv`): hyd1_like (KS0017, Eurotiomycetes), aga1 (KS0080, Saccharomycotina), ece1 (KS0107, Saccharomycotina) and ecm33 (KS0086, Saccharomycotina). Three are adhesion and one is wall. All four are clade-only clusters. The weak table (E up to 1e-3, descriptive) gives KS0017 one family (`PF22354` Eas), KS0086 two (`PF01030` Recep_L_domain, `PF23952` LRR_EndoS), and none for KS0080 and KS0107. It moves no cluster out of `no_pfam_hit`.

## 2. The 38 clusters that hit another Pfam family

The provisional cover rule (`proposal_fraction` 0.5: the best domain hit covers at least half of the sequence) splits them:

| | adhesion | wall | fusion | total |
|---|---|---|---|---|
| `add_to_panel_proposal` (cover at least 0.5) | 5 | 6 | 3 | 14 |
| `partial_domain_for_owner` (cover under 0.5) | 4 | 19 | 1 | 24 |

The 14 clusters of the first row belong to 12 Pfam families (`panel_proposals.tsv`), listed with the number of clusters and the range of cover:

| family | clusters | slice | cover |
|---|---|---|---|
| `PF01828` Peptidase_A4 | 2 | adhesion | 0.649 to 0.687 |
| `PF02076` STE3 | 2 | fusion | 0.648 to 0.832 |
| `PF02116` STE2 | 1 | fusion | 0.650 |
| `PF03151` TPT | 1 | adhesion | 0.719 |
| `PF03388` Lectin_leg-like | 1 | adhesion | 0.703 |
| `PF29797` Rep1 | 1 | adhesion | 0.598 |
| `PF02055` Glyco_hydro_30 | 1 | wall | 0.700 |
| `PF03200` Glyco_hydro_63 | 1 | wall | 0.607 |
| `PF03935` SKN1_KRE6_Sbg1 | 1 | wall | 0.682 |
| `PF12271` Chs7 | 1 | wall | 0.949 |
| `PF13641` Glyco_tranf_2_3 | 1 | wall | 0.506 |
| `PF17056` KRE1 | 1 | wall | 0.508 |

These are proposals for the owner. No panel model was added.

## 3. What the tables show, stated plainly

- **A cover rule hides strong hits.** The best hits of the 24 partial-domain clusters cover between 0.044 and 0.485 of the sequence, and five of them have a best score (435.1 to 1058.5) above the median best score of the 14 proposals (280.7). iff4 (KS0127) hits `PF11765` Hyphal_reg_CWP with score 435.1 and E-value 2.2e-130 at cover 0.215. That family is not in the panel. Other strong partial-domain hits are so_soft (KS0138, the one fusion cluster of that kind) `PF19050` PhoD_2 (score 1058.5, cover 0.43), KS0115 (ags1 and others) `PF26127` 12TM_Mok13 (684.6, 0.179) and KS0009 (kre5) `PF18404` Glyco_transf_24 (462.9, 0.185). Domain cover is a weak guide to how informative a hit is.
- **glfB is a transporter.** KS0120 (glfB) hits `PF03151` TPT (triose-phosphate transporter) with score 90.6, E-value 1.9e-25, cover 0.719. This agrees with the doubt already listed in the spec (decision 12) about glfB as an adhesion row. This note does not move it.
- **ecm33 (KS0086, Q5AGC4): to check, not explained.** No family hits it at the GA threshold, although `Ecm33` (`PF12454`) is in the panel (`analysis/wall_fusion_report/wall_fusion_models.tsv`, row 17). That file already records, from report 07 row 22, that *S. cerevisiae* Ecm33 itself has no hit to this model, so the miss is not new for its namesake. Whether Q5AGC4 behaves the same way is not established here.
- **A Pfam hit is not specificity.** A hit says that a family covers part of the sequence. It does not say that the family separates adhesion or wall proteins from the rest of a proteome. Broad families such as `Glyco_tranf_2_3` and `Glyco_hydro_30` need the per-family specificity test (D19, issue #69) before anyone considers adding them to the panel.

## 4. What this note does not do

- It changes no tier, no `narrow_ok`, no status file and no panel.
- It makes no claim of recall or specificity. These are class counts of 44 clusters.
- The decision on how to use the 38 hits (which to treat as panel additions, how to rank them, whether to run the specificity test first, and what happens to the four clusters without a hit) is pending with the owner.
