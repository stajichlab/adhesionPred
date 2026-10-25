# AGA1-class prediction step (DRAFT, NOT SUBMITTED)

Plan and reasoning: `docs/reports/2026-10-24-aga1-prediction-step-plan.md`. Nothing here has run on HPCC or SLURM. Every
submission below needs the owner's go-ahead (D21). No sequence is committed.

| File | Role |
|---|---|
| `stage_inputs.py` | local, deterministic: builds `aga1_proteins.faa` for the 850 accessions from the out-of-repo cache, with `.sha256`, `count.txt`, `lengths.tsv` |
| `manifest_accessions.txt` | committed, accessions only (850, sorted); a test compares it with `proteins.tsv` |
| `signalp.sbatch`, `predgpi.sbatch`, `disorder.sbatch` | one job per tool, `short` partition, pin check, `${SCRATCH:?}`, refuse an existing `OUT`, provenance file |
| `convert_predictions.py` | raw SignalP and PredGPI gff3 (and, later, a disorder table) to one table per protein (stdlib); fails loudly, no silent zeros |
| `build_record_table.py` | joins that table with the readable columns of `proteins.tsv` (role group, gene, UniProt features) |

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
(`CODE_COMMIT` = the merge commit), the log directory exists, `OUT` does not exist. The sbatch files are final: the owner passes
only values through `--export` and the command line, and the files are committed unedited after running.

```
W=/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/aga1_prediction
sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=$W/in,OUT=$W/signalp \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log analysis/aga1_class/signalp.sbatch
sbatch --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=$W/in,OUT=$W/predgpi \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log analysis/aga1_class/predgpi.sbatch
# AIUPred needs a GPU: name the GPU partition on the command line (not written in the file; look it up with `sinfo`)
sbatch -p <gpu-partition> --gres=gpu:1 \
  --export=ALL,PROJ_ROOT=/path/to/clone,CODE_COMMIT=<sha>,IN=$W/in,OUT=$W/disorder_aiupred,TOOL=aiupred \
  -o $W/logs/%x.%j.log -e $W/logs/%x.%j.log analysis/aga1_class/disorder.sbatch
# optional second disorder tool on CPU, own go-ahead: OUT=$W/disorder_iupred3,TOOL=iupred3 (no -p or --gres)
```

The AIUPred job stops on a node without CUDA unless `ALLOW_CPU=yes` is exported (CPU AIUPred is documented as very slow).

## 3. Convert (local, after the raw outputs are copied back)

PredGPI's raw `predgpi_raw.gff3` is read directly by the converter (no adapter file); `--lengths` is required for it. The
disorder block is explicit: give `--disorder FILE --disorder-tool T`, or `--no-disorder` (the disorder columns are then left out,
never zero). The SignalP + PredGPI run of 2026-10-24 used `--no-disorder` (`docs/reports/data/aga1_prediction_step_2026-10-24/`); the disorder tables of the same day are in `docs/reports/data/aga1_disorder_2026-10-24/`.

```
python3 -I analysis/aga1_class/convert_predictions.py --manifest analysis/aga1_class/manifest_accessions.txt \
  --signalp <OUT>/signalp/prediction_results.txt --predgpi <OUT>/predgpi/predgpi_raw.gff3 \
  --lengths "$HOME/aga1_prediction_in/lengths.tsv" \
  --no-disorder   # later: --disorder <OUT>/disorder_aiupred/<per-residue file> --disorder-tool aiupred
  --out <new table path>
python3 -I analysis/aga1_class/build_record_table.py --predictions <that table> --out <new record table>
```

Output columns: `accession` (join key to `proteins.tsv`), `signal_peptide_pred`, `sp_cleavage_pos`, `sp_cleavage_prob`, `sp_prob`,
`sp_prediction_raw`, `gpi_pred` (yes / no / NA), `gpi_omega_site`, `gpi_confidence_class` (high / probable / weak),
`gpi_fpr_upper_bound` (0.0015 / 0.005 / 0.01), `gpi_note`, then, only with a disorder file, `disorder_tool`, `disorder_threshold`
(0.5), `disorder_n_res`, `disorder_mean`, `disorder_frac`, `disorder_longest_run`. A value that does not apply is `NA`, never 0.
The gff3 holds no false positive rate and no raw score, so none is given. PredGPI never scores a sequence of 40 aa or fewer: those
get `gpi_pred` NA with a note, not `no`. The converter exits 1 with a `FATAL:` message on a missing or extra accession, a
duplicate, an empty file, an unparseable row, a SignalP header that is not a Eukarya run, a gff3 row with the wrong column count,
an unknown type or confidence score, a `GPI-anchor` row with start != end, a `Chain` row that does not span the protein, a
disorder score outside [0, 1] or out of order, or a residue count that differs from the staged length. It refuses an existing
output file.

The `.sha256` file written by `stage_inputs.py` holds the hash only (no file name), so `sha256sum -c` does not work on it. Compare
by hand: `sha256sum < aga1_proteins.faa` (or `shasum -a 256 < aga1_proteins.faa`) must equal `cat aga1_proteins.faa.sha256`.
The sbatch pin checks already do exactly that.


## 4. Command lines: confirmed and not (read-only checks on HPCC, 2026-10-24, from the install directories)

Confirmed and used in the sbatch files:

* PredGPI (`predgpi/202001`): `predgpi.py -f FASTA -o OUTF -m json|gff3` (default gff3; `-c` = conservative omega site). The module
  did not set `PREDGPI_HOME` in a non-interactive shell, so the job exports `/opt/linux/rocky/8.x/x86_64/pkgs/predgpi/202001`
  and stops if the directory or `predgpi.py` is missing. Drafted with `-m gff3`; JSON key names are not known, so JSON is not used.
* AIUPred (`aiupred/2.x`): `aiupred.py -i FASTA [-o OUT] [-g GPU] [--force-cpu]`; the module only prepends the PATH, loads no
  python; needs torch (~=2.0.1), numpy, scipy.
* IUPred3 (`iupred/3`): `iupred3.py SEQFILE {long|short|glob} [-a] [-s no|medium|strong]`, needs scipy; one sequence per file, so the
  job loops over records (`long` mode).
* SignalP: versions on disk 6.0g, 6.0h, 6.0i and their `-gpu` builds; the job uses the default `signalp/6` (CPU) and records
  which build it resolved to (`module -t list`, `readlink -f signalp6`, version) in `provenance.txt`.

Still to check or unknown:

1. Whether `python3` after `module load aiupred/2.x` can import torch, numpy and scipy (the job checks and stops with a message),
   and whether a CUDA device is visible; AIUPred and IUPred3 run times (CPU AIUPred is "very slow").
2. The name of the GPU partition (given by the owner with `-p`) and whether `short` takes the requested CPUs/memory without waiting.
3. Whether PredGPI finishes 850 sequences within the 2 h `short` limit (unmeasured).

## 5. Format assumptions (all four tools VERIFIED on the real outputs of 2026-10-24)

* SignalP 6 (verified) `prediction_results.txt`: two `#` header lines (the first holds `Organism: Eukarya`), then TAB rows
  `ID, Prediction, OTHER, SP(Sec/SPI), CS Position`; `CS Position` looks like `CS pos: 22-23. Pr: 0.9876` (the cleavage
  position stored is the first number, the last residue of the signal peptide) and is empty for `OTHER`. Verified on the real output
  (850 rows; the cleavage field is `CS pos: 24-25. Pr: 0.9770`). Eukarya, fast mode: only `SP` and `OTHER` are expected; anything else stops the converter.
* PredGPI (verified): 9-column gff3, one row per protein; `Chain` (no GPI, start 1, end = length) or `GPI-anchor` (start = end =
  omega site, score = confidence class 1.0 / 0.70 / 0.55, attribute `Ontology_term=GO:0046658`). No false positive rate and no raw
  score are in the file. Sequences of 40 aa or fewer are never scored. The converter reads this file directly.
* AIUPred (verified): a banner of `#` lines, then per record a `#>ID` line, `pos<TAB>aa<TAB>score` rows and a blank line. The converter reads `#>ID` as a record line. Scores overshoot [0, 1] slightly (up to 1.0036).
* IUPred3 (verified): the job writes its own `>ID` line and drops the tool's `#` lines; rows are `pos aa score`, scores from -0.0123 to 1.0096. The converter allows 0.02 outside [0, 1] and uses values as written.
* The converter reads `>ID` records followed by `pos aa score` rows, scores in [0, 1], `#` lines ignored. The sbatch jobs only
  warn (not stop) when the record count of the raw disorder or PredGPI output is off, so the raw file is never lost.
* SignalP 6 and the 74 proteins over 1,022 aa: all 74 got a call (53 SP, 21 OTHER); whether the tool truncated them is still
  not known (see the results report).
* The FASTA header is the bare accession, so every tool should return it as the first token; a tool that rewrites ids would trip the
  converter's accession checks.
