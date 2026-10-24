# AGA1-class prediction step (DRAFT, NOT SUBMITTED)

Plan and reasoning: `docs/reports/2026-10-24-aga1-prediction-step-plan.md`. Nothing here has run on HPCC or SLURM. Every
submission below needs the owner's go-ahead (D21). No sequence is committed.

| File | Role |
|---|---|
| `stage_inputs.py` | local, deterministic: builds `aga1_proteins.faa` for the 850 accessions from the out-of-repo cache, with `.sha256`, `count.txt`, `lengths.tsv` |
| `manifest_accessions.txt` | committed, accessions only (850, sorted); a test compares it with `proteins.tsv` |
| `signalp.sbatch`, `predgpi.sbatch`, `disorder.sbatch` | one job per tool, `short` partition, pin check, `${SCRATCH:?}`, refuse an existing `OUT`, provenance file |
| `convert_predictions.py` | raw outputs to one table per protein (stdlib); fails loudly, no silent zeros |

## 1. Stage the input (local)

```
python3 -I analysis/aga1_class/stage_inputs.py --out "$HOME/aga1_prediction_in"   # a NEW directory outside the repo
cat "$HOME/aga1_prediction_in/count.txt"                                            # 850
```

The owner copies it to HPCC later (NOT RUN; replace the host and path with the real ones):

```
scp -r "$HOME/aga1_prediction_in" <user>@<hpcc-login-host>:/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/aga1_prediction/in
```

Proteins over 1,022 aa (74 of 850) and short ones (shortest 35 aa) are all included; `lengths.tsv` records every length, the
number of non-standard residues (one protein has 2) and the sequence source.

## 2. Submit commands (DRAFT, NOT SUBMITTED; one go-ahead per submission)

Preconditions: the PR holding these files is merged, the clone on HPCC is at that commit and clean for `analysis/aga1_class`
(`CODE_COMMIT` = the merge commit), the log directory exists, `OUT` does not exist. Every job is committed unedited after running.
The two placeholder commands (section 4) must be edited and merged BEFORE the pin commit is chosen, never edited on HPCC.

```
W=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/aga1_prediction
sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=$W/in,OUT=$W/signalp \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log analysis/aga1_class/signalp.sbatch
sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=$W/in,OUT=$W/predgpi,CMD_CONFIRMED=yes \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log analysis/aga1_class/predgpi.sbatch
sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=$W/in,OUT=$W/disorder_aiupred,TOOL=aiupred,CMD_CONFIRMED=yes \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log analysis/aga1_class/disorder.sbatch
# optional second disorder tool, own go-ahead: OUT=$W/disorder_iupred3,TOOL=iupred3
```

Without `CMD_CONFIRMED=yes` the PredGPI and disorder jobs stop after the pin and input checks and run no tool.

## 3. Convert (local, after the raw outputs are copied back)

```
python3 -I analysis/aga1_class/convert_predictions.py --manifest analysis/aga1_class/manifest_accessions.txt \
  --signalp <OUT>/signalp/prediction_results.txt --predgpi <normalised PredGPI table> \
  --disorder <OUT>/disorder_aiupred/<per-residue file> --disorder-tool aiupred --lengths "$HOME/aga1_prediction_in/lengths.tsv" \
  --out <new table path>
```

Output columns: `accession` (join key to `proteins.tsv`), `signal_peptide_pred`, `sp_cleavage_pos`, `sp_cleavage_prob`, `sp_prob`,
`sp_prediction_raw`, `gpi_pred`, `gpi_omega_site`, `gpi_fp_rate`, `gpi_score`, `disorder_tool`, `disorder_threshold` (0.5),
`disorder_n_res`, `disorder_mean`, `disorder_frac`, `disorder_longest_run`. A value that does not apply is `NA`, never 0. The
converter exits 1 with a `FATAL:` message on a missing or extra accession, a duplicate, an empty file, an unparseable row, a
SignalP header that is not a Eukarya run, a disorder score outside [0, 1] or out of order, or (with `--lengths`) a residue count
that differs from the staged length. It refuses an existing output file.

## 4. To confirm read-only on HPCC before submit (placeholders)

1. PredGPI CLI: `module load predgpi/202001; echo "$PREDGPI_HOME"; ls "$PREDGPI_HOME"; <main script> --help`. The draft line
   `predgpi.py -f ... -o ...` is a guess. Also the PredGPI version string and which model/threshold the default uses.
2. AIUPred: `module avail aiupred`, then `module load aiupred; <script> --help`. The draft `aiupred.py -i ... -o ...` is a guess;
   also the exact module name and whether `--gpu` or a GPU-less run is the default.
3. IUPred3: `module load iupred/3; <script> --help`. The draft `iupred3.py FASTA long` is a guess; IUPred3 normally takes one
   sequence per call (loop needed) and has `long`, `short`, `glob` modes (`long` is drafted).
4. SignalP: `module load signalp/6; signalp6 --help` (flag names are those of the repo's GPU job, which ran on HPCC; the CPU
   module `signalp/6` was not run by the author). Whether the CPU build honours `--torch_num_threads`.
5. Whether `short` allows the resources asked (4 to 8 CPUs, 8 to 16 GB) without waiting behind `AssocGrpCpuLimit`.

## 5. Format assumptions, all UNVERIFIED until the first real output is seen

* SignalP 6 `prediction_results.txt`: two `#` header lines (the first holds `Organism: Eukarya`), then TAB rows
  `ID, Prediction, OTHER, SP(Sec/SPI), CS Position`; `CS Position` looks like `CS pos: 22-23. Pr: 0.9876` (the cleavage
  position stored is the first number, the last residue of the signal peptide) and is empty for `OTHER`. The first four columns
  are those the repo's `modules/signalp.py` already parses on real output; the cleavage field text is from the SignalP 6
  documentation as remembered, not from a file read here. Eukarya, fast mode: only `SP` and `OTHER` are expected; anything else stops the converter.
* PredGPI: the converter reads a NORMALISED table (`id`, `prediction` = `GPI-anchored` or `Not GPI-anchored`, `omega_site`,
  `fp_rate`, optional `score`). The raw PredGPI output layout, the meaning and scale of its false-positive rate, and whether it
  reports an omega site are all unknown; an adapter from the raw file to this table is written after the first real output.
* AIUPred and IUPred3: the converter reads `>ID` record lines followed by `pos aa score` rows, scores in [0, 1], `#` lines
  ignored. Raw layouts of both tools are unknown (IUPred3 prints a `#` comment block then `pos aa score`; AIUPred's layout is not known).
* SignalP 6 may only look at the first part of long sequences or truncate them; whether that applies to the 74 proteins over
  1,022 aa is not known (the SP lies at the N-terminus either way, so the call should not change, but this is not tested).
* The FASTA header is the bare accession, so every tool should return it as the first token; a tool that rewrites ids would trip the
  converter's accession checks.
