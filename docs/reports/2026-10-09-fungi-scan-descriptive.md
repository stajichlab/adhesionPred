# Descriptive scan of 56 Fungi_5k proteomes

*2026-10-09. Claude Code (claude-sonnet-5-5). Branch `fungi-scan-and-specs`. This is a description, not a measurement. There is no truth set. The counts show what the calls report. They do not show sensitivity, specificity or the true number of cell wall proteins.*

## Method

- 56 proteomes from `/bigdata/stajichlab/shared/projects/Fungi_5k`, chosen by `analysis/fungi_scan/select_proteomes.py` (`selection.tsv`): 24 Onygenales, 16 Eurotiales, plus Serinales (4), Saccharomycetales (2), Tremellales (2), Agaricales (2) and one proteome each of Ustilaginales, Malasseziales, Sordariales, Hypocreales, Magnaporthales and Mucorales.
- Modules: SignalP 6, TMHMM, four Pfam modules (`hmmsearch --cut_ga`), repeat02, repeat14, IUIS allergen BLAST. Core command with the default gate `step1_rule@R0`. All 15 SLURM jobs ended `COMPLETED` (pfam 1 h 24 min, TMHMM 0.9 to 1.5 h, repeats 0.6 to 0.7 h, SignalP 0.4 h, allergen 0.7 h).
- Summary tables: `docs/reports/data/sorting_hat/fungi_scan/` (`proteome_summary.tsv`, `family_summary.tsv`, `attachment_basis.tsv`, `attachment_candidates.tsv`, 1,053 candidate proteins).
- All status values in these runs are those of the shipped config. Calls with status `unvalidated` or `smoke` stay that way here. No new status was measured.

## Counts per proteome (min, median, max over proteomes)

| Order (proteomes) | Signal peptide | Tandem repeat | Wall family domain | Hydrophobin | HsbA | Cell wall adhesion candidate | Surface attachment candidate |
|---|---|---|---|---|---|---|---|
| All (56) | 156, 463, 1324 | 5, 22, 70 | 0, 5, 20 | 0, 2, 12 | 0, 2, 9 | 1, 12, 40 | 1, 18, 49 |
| Onygenales (24) | 192, 438, 738 | 7, 19, 45 | 2, 5, 7 | 0, 1, 3 | 0, 1, 9 | 5, 10, 20 | 5, 12, 21 |
| Eurotiales (16) | 588, 822, 1071 | 16, 33, 70 | 2, 4, 6 | 3, 6, 8 | 0, 7, 9 | 9, 15, 40 | 19, 27, 49 |

Single-proteome orders (Hypocreales, Magnaporthales, Sordariales, Mucorales, Ustilaginales, Malasseziales) are in `proteome_summary.tsv`. With one proteome per order, no range exists and no comparison between those orders is made here.

## Families found (proteins, proteomes), all 56 proteomes

| Family | Module | Proteins | Proteomes |
|---|---|---|---|
| PF05730 (CFEM) | pfam_adhesion | 235 | 51 |
| PF12296 (HsbA) | pfam_hsba | 148 | 33 |
| PF01185 (hydrophobin) | pfam_hydrophobin | 121 | 43 |
| PF07691 (PA14) | pfam_adhesion | 26 | 7 |
| PF22354 | pfam_hydrophobin | 20 | 16 |
| PF28987 | pfam_hydrophobin | 10 | 10 |
| PF06766 | pfam_hydrophobin | 6 | 4 |
| PF29785 | pfam_hydrophobin | 1 | 1 |

Active families with no hit in any of the 56 proteomes do not appear in this table. I did not check whether their absence is expected for these species.

## What the attachment basis shows

`attachment_basis.tsv` counts, per proteome, the distinct proteins held by each evidence type of `surface_attachment_candidate`. In the two *Coccidioides* proteomes the basis is: *C. immitis* RS wall family domain 4, hydrophobin 1, tandem repeat 6; *C. posadasii* Silveira wall family domain 5, hydrophobin 1, tandem repeat 6. Other proteomes are in the table. A protein can have more than one basis.

## SOWgp unit HMM (descriptive)

`sowgp_unit.hmm`, `hmmsearch --nobias`, domain score 26 or more, on the 56 proteomes. Result: hits only in the two *Coccidioides* proteomes (*C. immitis* RS: 1 protein, 4 domains; *C. posadasii* Silveira: 1 protein, 4 domains). The best domain score below 26 in any other hit is 23.0. This matches the earlier 6,313-proteome search (all hits at 26 or more are *Coccidioides*). It adds no new evidence: these proteomes are a subset of that search's space. The score here is a domain score.

## Limits

1. No truth. A count of `called` proteins is not a count of real cell wall or adhesion proteins.
2. The `tandem_repeat_protein` call has smoke status in two species only (S288C and *C. albicans*; precision 0.435 and 0.462 with wide intervals). Its precision in the 56 proteomes is not known.
3. CFEM, PA14 and hydrophobin hits come from Pfam GA cutoffs. Whether a protein in these families works as an adhesin is not tested here.
4. The proteomes are a convenience selection. They do not sample the fungal tree evenly.
5. The sample sizes per order differ (24, 16, then 1 to 4). Medians for orders with 1 to 4 proteomes are not comparable with the Onygenales and Eurotiales rows.

## Reproduce

```
bash analysis/fungi_scan/run/scan_converters.sh        # needs PROJ_ROOT, RUN_LIST, the csh venv
/usr/bin/python3.12 analysis/fungi_scan/summarize_scan.py --dir analysis/fungi_scan --out docs/reports/data/sorting_hat/fungi_scan
```
