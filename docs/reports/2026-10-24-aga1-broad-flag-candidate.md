# AGA1 broad flag candidate: signal peptide + GPI + AIUPred stalk run

Date: 2026-10-24. Owner decisions: AIUPred alone is the disorder tool; minimum disordered run N = 200 residues, kept as a parameter
(`analysis/aga1_class/broad_flag.py --min-disordered-run`), 100 as the pre-declared sensitivity value; Ser/Thr stalk rule = a
30-residue window, wholly inside the disordered run, that is at least 50% Ser/Thr (parameters `--window`, `--min-st-frac` of
`stalk_features.py`). Tables and settings: `docs/reports/data/aga1_broad_flag_2026-10-24/`.

## The stalk rule: what it does to the counts (flag = yes)

| Group | n | no stalk rule, N=200 | **stalk rule, N=200** | no stalk rule, N=100 | stalk rule, N=100 |
|---|---|---|---|---|---|
| tier A | 35 | 23 | **18** | 26 | 22 |
| tier B | 32 | 18 | **18** | 23 | 22 |
| hard negatives | 172 | 8 | **8** | 23 | 20 |
| background assumed negatives | 389 | 7 | **7** | 31 | 22 |
| other-mechanism adhesins | 49 | 6 | **6** | 12 | 10 |
| Kiss rows (other clusters) | 70 | 1 | **1** | 1 | 1 |

**At N = 200 the rule costs five tier A proteins and removes no negative.** It is not doing the job it was added for.

* **The five tier A proteins it drops.** SOWgp58, SOWgp66 and SOWgp82 (Coccidioides) have no Ser/Thr-rich window anywhere (best
  window 27%): their disordered run is not Ser/Thr-rich, so they are a different composition class from the flocculin type, in line
  with keeping the Pfl/SOWgp/Mad groups tagged separately. ALS9 and EAP1 have a long run (457 and 280) without a stalk window; their
  Ser/Thr-rich stretch sits in a shorter run (156 and 126) after the run is broken by one low-scoring residue, so the rule is
  sensitive to how a run is cut.
* **The flagged negatives keep their flag because they are Ser/Thr-rich.** DAN4 (best window 90%), chiA (97%), TIR4, FIT1, HKR1,
  cfmA, chiA1 and xgeA all pass; the Aspergillus chitinase and glucanase proteins have a Ser/Thr-rich linker, which is real
  composition, not noise. A Ser/Thr window cannot tell them from a flocculin stalk. A domain check (a catalytic domain in the same protein) could.
* **Seven of the flagged "background assumed negatives" are probably not negatives for the broad class:** PGA23, PGA25, PGA32,
  PGA39 (CGD names meaning predicted GPI-anchored protein), IHD1, FGR23 and PRY3 all pass with stalk windows of 53 to 80%. The background was assumed negative for the
  narrow adhesion task, so the broad flag may be right here and the labels the weak part.

Sensitivity (`stalk_grid.tsv`, tier A / hard negatives / background yes at N=200): window 20 and 30 with a 40% cut are the same as no
rule (23/8/7 and 20/8/7); 30 at 50% is 18/8/7; 30 at 60% is 13/6/5; window 50 at 60% is 8/5/5. A stricter stalk rule loses tier A
about as fast as it loses negatives. At N=100 the rule does remove some negatives (hard 23 to 20, background 31 to 22), at a cost of four tier A.

## PredGPI misses: is there another route?

Eight tier A proteins have no PredGPI call. They split in two:

* **HWP1** has a UniProt GPI feature and a hydrophobic C-terminal tail (mean Kyte-Doolittle 2.46 over its best 15-residue window,
  ending at the last residue). PredGPI misses it; the tail shows the GPI signal.
* **The seven fission-yeast proteins (map4, pfl3, pfl4, pfl5, pfl6, pfl7, pfl9)** have no UniProt GPI feature and no hydrophobic tail
  (tail values -0.76 to 0.80). There is no GPI signal in them to recover; they are probably not GPI-anchored the way the flocculins
  are. They need another route or a separate class, which fits the planned yeast-only comparison and the Pfl group tag.

An **exploratory tail rescue** is in `broad_flag.py` (`--tail-rescue-min-hydro`, off by default): a protein with a signal peptide, a
PredGPI no and a tail value at or above X counts as GPI-positive. Over the 850 proteins, among those with a signal peptide, 237 of
240 PredGPI-yes proteins have a tail of at least 1.0, against 106 of 379 PredGPI-no. With X = 1.4 the flag gains HWP1 (tier A 18 to 19) and
two other-mechanism adhesins (CSS1 and the Candida glabrata protein Q6FTA2; 6 to 8), and no hard negative or background protein; at X = 2.0, HWP1 and Q6FTA2. That is three proteins. The rule was
written after reading these same proteins, so it is not evidence; it needs proteins it has not seen before it can be trusted.

## What this is not, and what would test it

Still an unevaluated flag on the same 850 proteins that define tier A and the grid. The stalk-rule cut and the rescue were both read
off this set. The next test needs proteins the flag has not seen. The obvious source is the two best-annotated yeast proteomes
(S. cerevisiae S288C 6,722 and C. albicans SC5314 6,212 proteins, already on HPCC): run SignalP, PredGPI and AIUPred on them and
read what the flag calls against the known wall and flocculin proteins. That is three more HPCC jobs and needs your go-ahead.

## Next decisions (owner)

1. Keep the stalk rule as a gate? I recommend not gating on it (it costs five tier A and removes no negative) but keeping it as a
   column and a parameter; the enzyme linkers need a domain check instead.
2. Tail rescue: leave off until it is tested on new proteins (recommended), or turn it on for HWP1-type misses.
3. Run the three tools on the two yeast proteomes as the test on new proteins.
4. Still open: clade unit for narrow classes, glycosylation predictor, ESM-2 tiling, extra UniProt positives.
