# Task 03. Controls for the tandem-repeat call (`tandem_repeat_protein`)

*Read `COMMON-RULES.md` first. Written 2026-10-06.*

## Goal

Build positive and negative controls for the call `tandem_repeat_protein` (repeat detector 02 OR
repeat detector 14), by mechanism label, so that each detector can be calibrated alone.

## Why this matters

- Status today: `unvalidated`. The only checks are a synthetic series and SOWgp (tuned on SOWgp).
- Owner decision C1: positives are curated class-2a repeat (avidity) adhesins. Negatives are curated
  hard negatives and secreted non-repeat proteins. Homology clusters. SOWgp is marked `tuned` and
  excluded from the metrics.
- `calibrate truth` accepts a call that reads one module. `tandem_repeat_protein` reads two
  (`repeat02`, `repeat14`). Until per-call status entries exist, calibrate each detector on a
  call that reads only that module. Ask the owner how to name that call. Do not guess.

## Count made on 2026-10-06 (the starting point)

Table `data/curated/adhesins/adhesins.tsv`; sequences fetched from UniProt; MMseqs2 30% identity,
50% coverage. Script: `analysis/calibration_truth/c1_truth_count.py`.

| Set | Rows | Sequences | Clusters |
|---|---|---|---|
| `adhesin`, evidence E1 | 55 | 54 | 32 |
| `adhesin`, E1 + E2 | 98 | 96 | 50 |
| `adhesin`, E1 + E2, reviewed | 60 | 60 | 33 |
| `adhesin`, E1, `needs_review=no` | 44 | 44 | 27 |
| `hard_negative`, N1 + N2 | 42 | 42 | 33 |
| `hard_negative`, N1 | 32 | 32 | 25 |
| `surface_other_adhesion_phenotype` | 36 | 36 | 28 |
| E1 + E2 with "repeat" or "tandem" in family, summary or name | 8 | 8 | 5 |

- 2 clusters hold both an E1/E2 adhesin and a hard negative. Resolve them.
- 47 of 98 E1/E2 rows have an empty `family` column.
- **The table has no mechanism label.** Only 5 clusters can be marked repeat-mediated from text.
  The 20-cluster floor is not met for repeat adhesins. This is the main gap.
- Taxonomic skew: *C. albicans* 27, *C. glabrata* 14, *C. auris* 14, *S. pombe* 11, *S. cerevisiae* 9,
  *A. fumigatus* 3, *Rhizopus delemar* 3, *C. immitis* 2 (E1+E2 rows).
- 28 of 42 hard negatives are *S. cerevisiae* and 11 are *A. fumigatus*.

## Tasks

1. **Add a mechanism label** to each E1/E2 adhesin, using the class codes of
   `docs/TOOL-ARCHITECTURE.md` section 2: `2a` repeat/avidity, `2b-i` CFEM, `2b-ii` small Cys-knot,
   `2b-iii` Bys1, `2c` hydrophobin, `2d` moonlighting, or `other`, or `unknown`. Add
   `label_confidence` (`high`, `medium`, `low`). Base the label on the paper that gave the adhesion
   evidence, quoted. Start from `data/controls/repeat-mechanism/curation_table.tsv`. Do not label by
   running a repeat detector. The repeat call is calibrated with class `2a` as positives.
2. **Find more repeat-mediated adhesins** outside Saccharomycotina: Pezizomycotina, Mucoromycota,
   Basidiomycota. Examples to check in the literature, not to assume: Msg family (*Pneumocystis*),
   Mad1/Mad2 (*Metarhizium*), BAD1 (*Blastomyces*), Epa and Hyr/Iff families, CotH (Mucorales),
   Flo11 relatives. Each needs a PMID and quote.
3. **Add negatives that look like adhesins but are not**: secreted proteins with Ser/Thr-rich linkers,
   mucin-like sensors (Msb2, Hkr1 are labelled negative in step 1 truth), GPI wall enzymes, and
   secreted proteins with repeats but no adhesion function (check each for evidence). Expand the
   hard negatives beyond *S. cerevisiae* to the species that the adhesins come from.
4. **Mark leakage.** SOWgp alleles, and any protein used to tune or benchmark detector 02 or 14
   (read the scripts' docstrings and `analysis/cocci_repeats/` and `analysis/model_review/` history),
   are `tuned=yes`. Say in `README.md` how you checked.

`stratum` values: positives `2a` (repeat call) and the other class codes; negatives
`hard_negative_N1`, `hard_negative_N2`, `secreted_other`.

## Size target

At least 20 clusters of repeat-mediated positives, and 20 clusters of negatives, per clade you want
to report. Report each clade separately. The paper needs to say that the detectors were validated
only in Saccharomycotina if that stays true.

## Do not

- Do not use detector output to decide a mechanism label.
- Do not treat `surface_other_adhesion_phenotype` rows as positives (the effect may be indirect).
- Do not edit `adhesins.tsv`. Record changes in `manual_overrides.tsv` so the table stays
  reproducible (see `data/curated/adhesins/README.md`).
