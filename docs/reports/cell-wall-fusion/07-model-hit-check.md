# Do S. cerevisiae agglutinins, wall enzymes and fusion proteins hit the sorting-hat Pfam families?

*2026-10-17. Measurement only. No call, threshold, `src/` or `family_table.tsv` change. Outputs: `analysis/cell_wall_check/`. Branch `cell-wall-fusion-review`.*

## 1. What was run

| Item | Value |
|---|---|
| Test set | 23 S. cerevisiae S288C proteins fetched from `rest.uniprot.org` (reviewed entries, `gene_exact` + taxon 559292), 2026-10-17. `test_set.tsv`, `test_set.faa`, `uniprot_raw/`. Gene name, systematic ID and length come from the UniProt entry. Every UniProt length equals the SGD ORF length in `analysis/model_review/s288c_scores.tsv`. The sequences were not compared residue by residue to SGD. |
| Tool | HMMER 3.4 (Aug 2023), `module load hmmer/3.4` on `ssh hpcc` login node. 6 s run. No slurm job. |
| Pfam | Release 38.2, `/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm` (the path `analysis/lectin_validation/scan_models.py` uses). Versions in `scan.versions.tsv`. |
| Models | All 24 rows of `data/sorting_hat/family_table.tsv`, incl. PF11766, PF05792, plus 12 "context" models that are not in the table (GH72 PF03198, GH16 PF00722, GH17 PF00332, GH76 PF03663, GH5 PF00150, GH18 PF00704, SUN PF03856, Ecm33 PF12454, Kre9_KNH PF10342, X8 PF07983, WSC PF01822, Asp PF00026). `hmmsearch --cut_ga`, as in `lectin_validation`. Script `scan_test_set.py`. |
| Extra scans | (a) same 36 models, no `--cut_ga`, E <= 10: no sub-GA hit exists, so no near misses (`scan_below_ga.tsv`; the two weak GA passes are listed below). (b) `hmmscan --cut_ga` against all of Pfam-A, to name what the proteins hit outside the 36 (`full_pfam_ga_hits.tsv`). |
| Sorting hat | **Not run.** The CLI needs SignalP, TMHMM, repeat, BLAST and Pfam module tables (`analysis/sorting_hat_run/run_converters.sh`, GPU job for SignalP). That is not cheap. Instead I joined the committed S288C run (`docs/reports/data/sorting_hat/classifier_check/curated_label_check.rows.tsv`, made when only PA14 was active) and derived labels from `categories.yaml` rules and my Pfam hits. Derived labels are marked "derived". |
| Shipped classifier | **Not re-run.** Scores are the committed table `analysis/model_review/s288c_scores.tsv`, column `p_shipped_t6_8M` (shipped 8M model, call = p > 0.5 as in `02_cv_and_proteome_eval.py:154`). |
| S288C proteome | SGD `orf_trans_all.fasta.gz`, 6,722 ORFs, sha256 `17e8b47e1ae23178c6000fbc4ab548f102d1b250ef9dff5d811feb3f03dd2c5b` (record: `docs/reports/data/sorting_hat/provenance/Scer_S288C.json`). Path on HPCC only (`_workdir/step1_compare/downloads/`). I did not scan the whole proteome; the committed specificity files agree with my hits (SAG1 244.9, FLO11 240.1, FLO1 PA14 114.0 and GLEYA 67.0). |

## 2. Result table

Scores are hmmsearch sequence scores (best domain in parentheses only when it differs). "Active" is the state of the row in `family_table.tsv` in this worktree. R0 and repeat are from the committed run; "n/m" = protein not in that table.

### Agglutinins and flocculins

| Protein | UniProt (ID) | Family-table hits [status] | Other Pfam hits (full scan) | UniProt GPI site | Shipped ESM p | R0 SP / repeat | Hat label (derived) |
|---|---|---|---|---|---|---|---|
| SAG1 YJR004C | P20840 (SAG1_YEAST, alt. AGAL1) | Candida_ALS_N PF11766 244.9 [active] | none | 627 | 0.959 | called / not_called | wall_family_domain, so cell_wall_adhesion_candidate |
| AGA1 YNR044W | P32323 | none | none | 699 | 0.9999 | called / not_called | other_surface_no_mechanism |
| AGA2 YGL032C | P32781 | none | AGA2 PF17366 95.8 | none | 0.406 (not called) | called / not_called | other_surface_no_mechanism |
| FIG2 YCR089W | P25653 | Flocculin_t3 PF13928 64.8 (5 domains, 24-28 each) [active] | none | 1588 | 0.999 | n/m | wall_family_domain; candidate only if R0 called (not measured) |
| FLO1 YAR050W | P32768 | Flocculin PF00624 388.1 (18 domains); PA14 PF07691 115.7 [active, needs SP]; GLEYA PF10528 67.0 [active] | none | 1514 | 0.998 | called / called | cell_wall_adhesion_candidate (repeat and family) |
| FLO11 YIR019C | P08640 | Flo11 PF10182 240.1 [active] | none | 1346 | 0.9998 | called / not_called | wall_family_domain, so cell_wall_adhesion_candidate |

Als models: PF11766 hits SAG1 only. **PF05792 (Candida_ALS) hits none of the 23 proteins**, and no sub-GA hit either. Neither ALS_M (PF30910, inactive) nor Hyphal_reg_CWP, PIR, PIR1-like_C, CFEM, Hyr1, hydrophobin or HsbA models hit anything in the set. AGA1 and AGA2 do not hit any Als, flocculin, PA14 or other family-table model.

### Wall enzymes and wall proteins

| Protein | UniProt | Family-table hits | Context / other Pfam hits | UniProt GPI | Shipped ESM p | R0 SP / repeat | Hat label (derived) |
|---|---|---|---|---|---|---|---|
| GAS1 YMR307W | P22146 | none | GH72 PF03198 504.4; X8 PF07983 102.4 | 528 | 0.030 | called / not_called | other_surface_no_mechanism |
| CRH1 YGR189C | P53301 | none | GH16 PF00722 170.8 | 482 | 0.713 (called) | called / not_called | other_surface_no_mechanism |
| UTR2 YEL040W | P32623 (CRH2_YEAST) | none | GH16 PF00722 159.5 | 445 | 0.634 (called) | called / not_called | other_surface_no_mechanism |
| SCW4 YGR279C | P53334 | none | GH17 PF00332 45.6 (GA 28.2); BtgE-BtgC PF29664 283.2 | none | 0.033 | n/m | not run |
| SCW11 YGL028C | P53189 | none | BtgE-BtgC PF29664 266.9; no GH17 hit | none | 0.826 (called) | n/m | not run |
| EXG1 YLR300W | P23776 | none | Cellulase PF00150 43.7 (GA 23.5) | none | 0.011 | called / not_called | other_surface_no_mechanism |
| CTS1 YLR286C | P29029 (CHIT_YEAST) | none | GH18 PF00704 81.6; CBM_19 PF03427 32.4 | none | 0.598 (called) | called / not_called | other_surface_no_mechanism |
| SUN4 YNL066W | P53616 | none | SUN PF03856 330.1 | none | 0.935 (called) | called / not_called | other_surface_no_mechanism |
| ECM33 YBR078W | P38248 | none | none in all of Pfam-A at GA; Ecm33 PF12454 no hit even at E 10 | 406 | 0.033 | called / not_called | other_surface_no_mechanism |
| DFG5 YMR238W | Q05031 | none | GH76 PF03663 530.8 | 437 | 0.020 | n/m | not run |
| DCW1 YKL046C | P36091 | none | GH76 PF03663 547.4 | 428 | 0.020 | n/m | not run |
| KRE9 YJL174W | P39005 | none | Kre9_KNH PF10342 46.7; Kre9_KNH1_C PF05390 166.2 | none | 0.605 (called) | called / not_called | other_surface_no_mechanism |
| BGL2 YGR282C | P15703 | none | GH17 PF00332 353.2; BtgE-BtgC PF29664 224.2 | none | 0.017 | n/m | not run |

"UniProt GPI none" means no GPI lipidation feature is annotated in the entry. It is not a test that the protein lacks a GPI anchor.

### Fusion proteins

| Protein | UniProt | Family-table hits | Other Pfam hits | TM helices (UniProt) | Shipped ESM p |
|---|---|---|---|---|---|
| FUS1 YCL027W | P11710 | none | SH3_1 PF00018 23.5 | 1 (72-96) | 0.048 |
| FUS2 YMR232W | Q05670 | none | none | 0 | 0.002 |
| PRM1 YNL279W | P53835 | none | PRM1 PF27758 475.5 | 5 (17-37, 106-126, 297-317, 425-445, 630-650) | 0.011 |
| FIG1 YBR040W | P38224 | none | Fig1 PF12351 192.3 | 4 | 0.007 |

Full table: `analysis/cell_wall_check/summary.tsv`.

## 3. Claim-by-claim check

Source tags: 02 = `02-cell-fusion.md`, 04 = `04-repo-coverage-and-gaps.md`.

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| 1 | 02 s2b: alpha-agglutinin, 650 aa, signal peptide 19 aa, GPI signal; symbol "SAG1" unverified | CONFIRMED | UniProt P20840 gene SAG1 (also AGAL1), 650 aa, signal 1-19, GPI site 627 |
| 2 | 02 s2b: Aga1 is 725 aa with a GPI wall anchor | CONFIRMED | P32323 725 aa, GPI site 699 |
| 3 | 02 s2b: Aga2 is 69 aa | CONFIRMED for the mature chain only | P32781 is 87 aa; chain 19-87 = 69 aa |
| 4 | 02 s3 / 04: Sc alpha-agglutinin is a candidate hit for Candida_ALS_N PF11766 | CONFIRMED | 244.9, domain 50-300 |
| 5 | 02 s3: alpha-agglutinin may hit Candida_ALS PF05792 ("both inactive") | REFUTED | no hit at GA, none even at E 10 |
| 6 | 02 / 04: PF11766 and PF05792 are inactive in `family_table.tsv` | REFUTED for the current table | both `active=yes`, owner 2026-10-09 (commit d94a8ab). This was true on `main` before that commit. Report 04 describes the pre-10-09 state |
| 7 | 04 s1.4/s3.2: Aga1 and Aga2 Pfam "unverified"; no family-table row covers them | CONFIRMED (no table hit). New fact | Aga1: no Pfam-A hit at all. Aga2 hits PF17366 (AGA2), which is not in the table. Neither hits an Als, flocculin or PA14 row |
| 8 | 04 s1.4: SAG1 hits PF11766 (`pfam-proposals.md:41`) | CONFIRMED | same 244.9 as the committed S288C specificity file. That file counts it as a false positive for the Candida-derived member list (`Scer_S288C.PF11766.specificity.tsv`: hits 1, members 0, fp 1), so the specificity test and the biology disagree |
| 9 | 04 s1.4: "FIG2 is a curated adhesin row with no family" | REFUTED | FIG2 hits Flocculin_t3 PF13928 (64.8, 5 weak domains, 24-28 each), an active row |
| 10 | 04 s1.1: no table row is an enzyme family | CONFIRMED for this set | none of the 13 wall proteins hits any of the 24 rows. PA14, noted as also in beta-glucosidases, does not hit BGL2, EXG1, SCW4 or SUN4 |
| 11 | 04 s2: wall housekeeping enzymes (GAS1, CRH1, UTR2, EXG1, CTS1, SUN4, ECM33, KRE9) are not caught by the Pfam families, only by R0 | CONFIRMED | no family-table hit; R0 `called` for all 8 (committed run); repeat `not_called` for all 8. Derived label: `other_surface_no_mechanism` |
| 12 | 04 s2 (spec): "a sorter with 'adhesin' as its only wall label would call these 'other'" | CONFIRMED (derived) | as row 11. Aga1 and Aga2 would also land in `other_surface_no_mechanism`, i.e. sexual agglutinins are missed by the Pfam and repeat calls, while SAG1 is not |
| 13 | 04 s2 (model review): shipped model calls CRH1, UTR2, SUN4, SCW11, KRE9, CTS1 | CONFIRMED | p = 0.713, 0.634, 0.935, 0.826, 0.605, 0.598 (all > 0.5) |
| 14 | 04 s2: shipped model calls known adhesins incl. AGA1, SAG1, FIG2, FLO1, FLO11 and misses 87-aa AGA2 | CONFIRMED | 0.9999, 0.959, 0.999, 0.998, 0.9998; AGA2 0.406 |
| 15 | 04 s2: hard negatives like GAS1, EXG1, ECM33 | NOT all called | GAS1 0.030, EXG1 0.011, ECM33 0.033 are not called. Of my 8 listed hard negatives 5 are called (CRH1, UTR2, SUN4, KRE9, CTS1). The "12/20" figure is not testable on this set |
| 16 | 04 s2: R0 keeps nearly every secreted wall protein (27 of 28 hard negatives) | CONFIRMED on the 8 present | committed R0 `called` 8/8 (the other 5 wall enzymes are not in that table) |
| 17 | 04 s3.1: GH72 PF03198 (Gas1), often with X8 PF07983 | CONFIRMED | Gas1 504.4 + X8 102.4 |
| 18 | 04 s3.1: GH16 PF00722 (Crh1, Utr2) | CONFIRMED | 170.8, 159.5 |
| 19 | 04 s3.1: GH76 PF03663 for Dfg5/Dcw1 (ID and link "my recall, unverified") | CONFIRMED | model name Glyco_hydro_76 (PF03663.21); Dfg5 530.8, Dcw1 547.4. UniProt names both mannan endo-1,6-alpha-mannosidase |
| 20 | 04 s3.1: GH17 PF00332 for Bgl2 (name unverified) | CONFIRMED | model name Glyco_hydro_17; Bgl2 353.2; also SCW4 45.6 |
| 21 | 04 s3.1: GH5 PF00150 (Exg1), GH18 PF00704 (Cts1), Kre9_KNH PF10342 (Kre9), SUN PF03856 (Sun4) | CONFIRMED | Exg1 43.7 (weak, GA 23.5), Cts1 81.6, Kre9 46.7, Sun4 330.1 |
| 22 | 04 s3.1: Ecm33/Sps2 family PF12454, ECM33 a member | REFUTED | S. cerevisiae Ecm33 has no hit to PF12454 (none at E 10) and no Pfam-A hit at all |
| 23 | 04 s3.1: Scw4/Scw10 soluble glucanases have no Pfam in repo | Superseded | SCW4 hits PF00332 (45.6) and PF29664 (283.2); SCW11 hits PF29664 (266.9) only, no GH17 |
| 24 | 02 open Q2: are Scw4 and Scw11 GPI anchored? | NOT TESTABLE here | UniProt has no GPI feature for either. Absence of annotation, not a result |
| 25 | 02 s2c: Fus1 is a 1-pass TM protein; Prm1 and Fig1 are 4-pass (Prm1 stated 4) | PARTLY | Fus1 1 TM (72-96) confirmed; Fig1 4 TM confirmed; Prm1 has 5 UniProt TM features (not 4) |
| 26 | 02 s3 / 04 s1.4: fusion proteins (Fus1, Prm1, Fig1, Fus2) are not covered by any family-table row | CONFIRMED | no hit in any of the 36 models. Pfam does have models for Prm1 (PF27758), Fig1 (PF12351); Fus1 only a weak SH3_1 (23.5); Fus2 none. This closes the "no ID, unverified" gap for two of them |
| 27 | 02 s3: Prm1/Fig1 are not visible to a surface classifier or Pfam adhesion call | CONFIRMED | shipped p 0.011, 0.007; no table hit |
| 28 | 02 s4: FLO1/FLO5/FLO9/FLO10/FLO11 in the table | CONFIRMED for FLO1, FLO11 | hits above; FLO5, FLO9, FLO10 not tested |
| 29 | 02 s3: Aga2 below a length filter | NOT TESTABLE | no length filter was applied here. The shipped model p 0.406 misses it |
| 30 | 02 s3: Fus1 removed by a no-TM rule | NOT TESTABLE | no TM gate on any family it hits (it hits none) |
| 31 | 04 s0: Chitin_synth PF01644/PF03142 swap, GPI-anchor families, WSC, Asp | NOT TESTABLE | no test protein carries them (WSC, Asp context models: no hit) |
| 32 | 04 s2: Scer R0/repeat call fires for PIR1, HSP150, TIR1, TIR4, FIT1 only | NOT TESTABLE | none are in the set; for my 8 proteins repeat is `not_called` except FLO1 |

## 4. What this means (measured facts only)

- With the Als rows active, **SAG1 is the only agglutinin that reaches `wall_family_domain` through an Als model**, and only through PF11766. FLO1, FLO11 and FIG2 get there through flocculin-type rows. **AGA1 and AGA2 hit no table row** and would be `other_surface_no_mechanism`.
- **None of the 13 wall enzymes hits an active family-table row**, so the Pfam families do not mislabel them as adhesion candidates. The shipped ESM classifier does: 6 of 13 are above 0.5 (CRH1, UTR2, SCW11, CTS1, SUN4, KRE9) while GAS1, EXG1, ECM33, BGL2, SCW4, DFG5, DCW1 are not.
- In the sorting hat the wall enzymes would be `other_surface_no_mechanism`; the hat has no enzyme label. Report 04 is right that the table has no enzyme row. The context models show GH72, GH16, GH17, GH76, GH18, SUN, Kre9_KNH all hit at scores far above GA (except weak Exg1, Scw4, Kre9_KNH), so a report-only enzyme list would be easy to build. Ecm33 would not be found by PF12454.
- The Prm1 and Fig1 Pfam models exist and hit strongly (475.5, 192.3). Report 04 s3.2 says no ID is verified. Here they are measured, not yet used.

## 5. Limits

- Single species, one protein per gene, 23 proteins. No negatives beyond the set itself.
- Sorting-hat R0 and repeat values come from an older committed run (family table with only PA14 active). They are not re-measured with the current config, and 9 proteins are absent from it. Derived labels combine them with a fresh Pfam scan and `categories.yaml`; the real CLI was not run.
- UniProt feature annotations are database text (untrusted, not independently verified). Missing GPI features do not show a protein is not GPI anchored.
- Shipped-model scores are read from `s288c_scores.tsv`, produced by the earlier run; I did not repeat the embedding.
- The proteome FASTA exists only on HPCC; the test FASTA was built from UniProt, not from SGD.

## 6. Files

`analysis/cell_wall_check/`: `scan_test_set.py`, `scan_below_ga.sh`, `fetch_features.py`, `make_summary.py`; `test_set.tsv`, `test_set.faa`, `scan.hits.tsv` (per domain), `scan.tbl`, `scan.versions.tsv`, `scan_below_ga.tsv`, `full_pfam_ga_hits.tsv`, `uniprot_features.tsv`, `summary.tsv`, `uniprot_raw/`.
Rebuild: scp `scan_test_set.py`, `test_set.faa` and the family table to HPCC, `module load hmmer/3.4`, run `scan_test_set.py family_table.tsv test_set.faa scan` with python3.12.
