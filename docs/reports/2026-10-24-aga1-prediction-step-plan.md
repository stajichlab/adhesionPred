# AGA1-class prediction step: what will be run and why (plan, DRAFT, NOT SUBMITTED)

*2026-10-24. Author: Claude Code (claude-sonnet-5-5). For the project owner. Drafts only: nothing was run on HPCC or SLURM, no sequence is committed, no call, tier, panel or `narrow_ok` was changed. Files: `analysis/aga1_class/` (README there has the commands).*

## Why this step

The owner chose nested labels for the AGA1-class study: a broad flag ("Ser/Thr-rich GPI cell-wall glycoprotein": SED1, DAN4, TIR family and the narrow positives) and a narrow flag ("flocculin/adhesin-like", tier A). UniProt annotation cannot define the broad flag (GPI and signal annotation cover about half of tier A: at S+T of at least 0.35 only 17 of 35 tier A and 8 of 172 hard negatives pass). So the owner decided to predict first and define the flag from the predictions afterwards. This step produces the predictions for the 850 proteins of `docs/reports/data/aga1_class_scoping/proteins.tsv`. It does **not** define the broad flag; that is the owner's decision after seeing counts.

## What will be run

| Job | Tool (module) | What it answers | Draft resources |
|---|---|---|---|
| `signalp.sbatch` | SignalP 6, eukarya, fast mode, CPU build (`signalp/6`) | Does the N-terminus carry a secretion signal peptide, and where is it cut? | `short`, 8 CPUs, 16 GB, 1.5 h limit |
| `predgpi.sbatch` | PredGPI (`predgpi/202001`) | Is there a C-terminal GPI-anchor signal, and where is the omega (attachment) site? | `short`, 4 CPUs, 8 GB, 2 h limit |
| `disorder.sbatch` | AIUPred (`aiupred/2.x`), optionally IUPred3 (`iupred/3`) | Per-residue disorder score; Ser/Thr-rich wall proteins are expected to be long disordered stretches | `short`, 4 CPUs, 8 GB, 1.5 h limit |

Input: one FASTA of the 850 proteins, built locally by `stage_inputs.py` from the out-of-repo sequence cache (header = bare accession), with a sha256 and a lengths table; the owner copies it to HPCC. Every job runs a pin check (the clone must match the reviewed commit for `analysis/aga1_class`), refuses an existing output directory, checks the FASTA against its sha256 and the committed accession manifest, and writes a provenance file with the tool versions.

## What each output means

`convert_predictions.py` turns the raw outputs into one row per protein (join key `accession`):

* **signal_peptide_pred** yes/no (SignalP says SP), **sp_cleavage_pos** (last residue of the signal peptide), **sp_cleavage_prob**, **sp_prob** (SignalP's probability for SP, also given for "no").
* **gpi_pred** yes/no, **gpi_omega_site** (if the tool gives one), **gpi_fp_rate** and **gpi_score** (PredGPI's own confidence numbers; their scale is not yet known).
* **disorder_mean**, **disorder_frac** (share of residues with score of at least 0.5), **disorder_longest_run** (longest stretch of residues at or above 0.5), and the residue count.
* Not applicable values are `NA`; the converter never writes a zero for a missing value and stops with a clear message on a missing or extra accession, a duplicate, an empty file or a row it cannot read.

## What the broad-flag definition will need from this

The owner will want, per label group (tier A, hard negatives, SED1/DAN4, background): how many have a signal peptide, how many have a GPI call, how many have both, how many both and also high Ser/Thr or long disordered runs, and how these counts shift with the PredGPI false-positive cut-off. The table gives all of these; the choice of cut-offs and the combination is open.

## Expected runtime (a guess)

SignalP: about 2 s per protein on CPU in the repo's earlier measurement, so about 30 min on one core and less with 8 threads (not measured for the CPU module). PredGPI and the disorder tools: seconds to a few minutes per hundred proteins in my expectation; the 74 proteins over 1,022 aa are the slow ones. Queue wait is the real unknown (CPU jobs have waited behind `AssocGrpCpuLimit`).

## Risks

* **SignalP on long proteins.** 74 of 850 proteins exceed 1,022 aa (longest 3,971). Whether SignalP 6 truncates long input is unknown; the signal peptide is at the N-terminus so the call should be unaffected, but this is not tested. Shortest proteins (35 aa) may get unreliable calls.
* **PredGPI false positives on non-GPI proteins.** Predictors of GPI signals fire on Ser/Thr-rich C-terminal hydrophobic stretches that are not anchors; the 172 hard negatives (cell-wall enzymes, mucin-type sensors) are exactly the proteins that could trigger this. The false-positive rate PredGPI reports is a model estimate, not measured on these negatives; the table lets the owner measure it on the hard negatives directly.
* **Organism/model settings.** SignalP is run with `eukarya`, correct for fungi but not tuned for them. PredGPI's default model (and whether it has a fungal one) is not known. AIUPred and IUPred3 have different modes (IUPred3 `long`, `short`, `glob`) and different score scales; the two are not interchangeable, so one tool is used per table and the tool name is a column.
* **Tool limits.** None of the three predicts the GPI anchor in the actual wall or the glycosylation; disorder is not adhesin-specific. Yeast-heavy data (see the scoping report) means any cut-off tuned here describes yeast.
* **Format risk.** All three raw formats are unverified (see the README); the first real output may need a small adapter before conversion.

## Open items

1. Confirm read-only on HPCC: PredGPI, AIUPred, IUPred3 command lines and versions; the `signalp/6` CPU module; whether `short` takes the requested CPUs/memory (README section 4). The placeholder lines must be edited and merged before the pin commit is chosen.
2. Which disorder tool is primary (AIUPred drafted first; IUPred3 as an optional second job) and the disorder threshold (0.5 drafted).
3. Whether the Pfl/SOWgp/Mad group tags matter for the prediction step (no; the table is keyed by accession and joins to `proteins.tsv`).
4. A PredGPI adapter from its raw output to the normalised table, after the first output is seen.
5. The owner's go-ahead for each submission, and the copy of the staged FASTA to HPCC (`scp` command in the README, not run).
6. The broad-flag definition itself, which stays with the owner.
