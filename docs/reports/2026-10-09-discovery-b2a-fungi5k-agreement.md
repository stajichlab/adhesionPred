# B2a: do the Fungi_5k SignalP 6 and TMHMM results agree with the scan run?

*2026-10-09. Claude Code (claude-sonnet-5-5). Check required by `docs/superpowers/specs/2026-10-09-secretome-discovery-pass-design.md`, section 4 (owner decision O3).*

## Method

Compared the results in `/bigdata/stajichlab/shared/projects/Fungi_5k/results/function/` with the raw tables of the 56 scan proteomes (`_workdir/fungi_scan/*/raw/`).

- SignalP: the set of proteins with Prediction `SP` in the scan (`prediction_results.txt`, SignalP 6.0h, organism Eukarya) against the proteins with a `signal_peptide` feature in the Fungi_5k GFF3.
- TMHMM: every row of the scan `tmhmm.tsv` (length, ExpAA, First60, PredHel, topology) against `*.tmhmm_short.tsv.gz`.

## Result

| Tool | Proteins compared | Identical | Different | Missing in Fungi_5k |
|---|---|---|---|---|
| SignalP 6, `SP` calls | 31,977 SP proteins in 56 proteomes | 31,977 | 0 | 0 |
| TMHMM, all five fields | 489,473 | 489,473 | 0 | 0 |

The scan found only the labels `SP` (31,977) and `OTHER` (457,496). So the R0 rule and the Fungi_5k SignalP agree on every protein of these 56 proteomes.

## Consequence

The discovery pass reuses the Fungi_5k SignalP and TMHMM results (decision O3). `analysis/fungi_scan/run/fungi5k_to_raw.py` writes the scan-format raw tables from them. One difference remains: the GFF3 holds the probability only for signal peptides, so a protein that is not `SP` gets SP probability 0 in the adapter output. The scan has the true small probability. This changes an evidence field, not a call. Checked on one proteome (*Uncinocarpus reesii* 1704): predictions and TMHMM rows are identical to the scan raw files.

## Limits

The check covers the 56 scan proteomes only. Both runs used SignalP 6.0h. Whether the Fungi_5k run used the same mode and settings for every other proteome is not checked here (it matches on these 56).
