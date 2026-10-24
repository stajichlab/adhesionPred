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

## 5. Format assumptions, all UNVERIFIED until the first real output is seen

* SignalP 6 `prediction_results.txt`: two `#` header lines (the first holds `Organism: Eukarya`), then TAB rows
  `ID, Prediction, OTHER, SP(Sec/SPI), CS Position`; `CS Position` looks like `CS pos: 22-23. Pr: 0.9876` (the cleavage
  position stored is the first number, the last residue of the signal peptide) and is empty for `OTHER`. The first four columns
  are those the repo's `modules/signalp.py` already parses on real output; the cleavage field text is from the SignalP 6
  documentation as remembered, not from a file read here. Eukarya, fast mode: only `SP` and `OTHER` are expected; anything else stops the converter.
* PredGPI: the raw gff3 layout and key names are unknown (the source prints `Estimated False Positive Rate =`,
  `HMM log(probability) =`, `GPI-Anchor length=` and one of Highly Probable / Probable / Weakly GPI-Anchored protein /
  Improbable GPI-Anchored protein; whether these appear in the gff3 attributes is not known). The converter reads a NORMALISED
  table (`id`, `prediction` = `GPI-anchored` or `Not GPI-anchored`, `omega_site`, `fp_rate`, optional `score`); an adapter from
  `predgpi_raw.gff3` to it is written after the first real output, together with the choice of which classes count as `yes`.
* AIUPred: per the README, a `>header` line then `pos<TAB>aa<TAB>score` rows; the exact header echo is unverified.
* IUPred3: `#` comment lines then `pos aa score` rows (from memory of the release). The job writes its own `>ID` line and drops the `#` lines.
* The converter reads `>ID` records followed by `pos aa score` rows, scores in [0, 1], `#` lines ignored. The sbatch jobs only
  warn (not stop) when the record count of the raw disorder or PredGPI output is off, so the raw file is never lost.
* SignalP 6 may only look at the first part of long sequences or truncate them; whether that applies to the 74 proteins over
  1,022 aa is not known (the SP lies at the N-terminus either way, so the call should not change, but this is not tested).
* The FASTA header is the bare accession, so every tool should return it as the first token; a tool that rewrites ids would trip the
  converter's accession checks.
