# Lectin-like family report (spec option c-i)

This folder holds a script that reports hits to the lectin-like Pfam models. It makes no call. It writes no status file. It claims no function. Spec: `docs/superpowers/specs/2026-10-10-lectin-class-design.md` (sections 3.1, 3.4, 6.2).

## Files

| File | Purpose |
|---|---|
| `lectin_report.py` | The report script. |
| `compare_committed.py` | Compares the hit table with `analysis/lectin_class_spec/review2/whole_proteome_hits.tsv`. |
| `submit_lectin_report.sbatch` | Optional SLURM wrapper. |
| `f27_six/` | The six PF07938 proteins of spec F27: sequences, SignalP 6 output, BLASTP check, name list, example report. |

## Versions used for the validation runs (2026-10-10)

- Python: `/usr/bin/python3.12`. Do not use the default `python3`.
- HMMER 3.4: `/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin/hmmsearch` (`module load hmmer/3.4`).
- Pfam 38.2: `/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`.
- SignalP 6.0h (module `signalp/6`, package 6.0i): `/opt/linux/rocky/8.x/x86_64/pkgs/signalp/6.0i/env/bin/signalp6`, fast mode, eukarya.
- BLASTP 2.17.0+: `/opt/linux/rocky/8.x/x86_64/pkgs/ncbi-blast/2.17.0+/bin/blastp`.

## Inputs

- Proteome FASTA files: plain, `.gz` or `.zst`. The proteome name is the file name without the compression and FASTA suffixes. Two inputs must not give the same name.
- `--signalp FILE` (optional). Only two formats are accepted. Any other line stops the run with an error that names the line.
  1. One column: the ID of a protein that has a signal peptide.
  2. SignalP 6 `ID<TAB>Prediction`, with Prediction in `SP`, `OTHER`, `LIPO`, `TAT`, `TATLIPO`, `PILIN`. `OTHER` means no signal peptide.
  The ID is the first word of the FASTA header. For a header `sp|ACC|NAME` the accession also matches. Matching uses the protein ID alone, not the proteome name.
- `--name-list FILE` (optional). One protein per line: an ID (full ID or bare accession), then a free note.
- Without `--signalp` the report has no (a2) section. With `--name-list` and no `--signalp` the report has a "Named proteins" section.

## Run: small test on a login node

Check the queue first. Use a few proteomes only.

```
/usr/bin/python3.12 analysis/lectin_class_report/lectin_report.py \
  --out-dir OUT_DIR --tmpdir TMP_DIR --cpu 4 \
  --hmmsearch /opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin/hmmsearch \
  [--signalp SIGNALP_FILE] [--name-list NAME_LIST] \
  --fasta-list LIST_OF_FASTA_PATHS
```

`TMP_DIR` must exist or be creatable. When `--tmpdir` is not given, the script uses `$SCRATCH` and stops if it is not set.

Outputs in `OUT_DIR`: `lectin_report_counts.tsv`, `lectin_report_hits.tsv`, `lectin_report.md`.

## Validation commands (as run on 2026-10-10)

25 proteomes of review 3, then comparison with the committed table:

```
cut -f2 analysis/lectin_class_spec/review3/rescan_sample.txt > list25.txt
cut -f1 analysis/lectin_class_spec/review3/rescan_sample.txt > names25.txt
/usr/bin/python3.12 analysis/lectin_class_report/lectin_report.py --out-dir OUT --tmpdir TMP \
  --hmmsearch HMMSEARCH --cpu 4 --fasta-list list25.txt
/usr/bin/python3.12 analysis/lectin_class_report/compare_committed.py OUT/lectin_report_hits.tsv \
  analysis/lectin_class_spec/review2/whole_proteome_hits.tsv --names names25.txt
```

Result: PF07938 36 proteins in 14 proteomes, the same set as the committed table, maximum score difference 0.0.

The six proteins of spec F27 (`f27_six/`):

```
/usr/bin/python3.12 analysis/lectin_class_report/f27_six/extract_six.py POOLS_DIR \
  SELECTION_TSV analysis/lectin_class_spec/review2/whole_proteome_hits.tsv six.faa six.tsv
signalp6 --fastafile six.faa --organism eukarya --mode fast --output_dir SP6_OUT --format txt
```

`POOLS_DIR` is `_workdir/discovery/pools`. `SELECTION_TSV` is `analysis/fungi_scan/selection_discovery.tsv`. Then run the report with `--signalp SP6_OUT/prediction_results.txt --name-list analysis/lectin_class_report/f27_six/f27_six_names.txt`. See `f27_six/example_report_six_proteomes.md`.

## Full run (653 discovery proteomes)

This recipe was not run. Check the queue before you submit.

1. Make the list of FASTA paths: `tail -n +2 analysis/fungi_scan/selection_discovery.tsv | cut -f3 > fasta_list.txt`.
2. Make signal peptide data. Run SignalP 6 on the proteomes (a separate job; it is not part of this script). Give one SignalP 6 `prediction_results.txt` or a one-column ID list to `--signalp`. Without it the report has no (a2) section.
3. Submit, with all paths absolute:

```
sbatch --export=ALL,REPORT_SCRIPT=/abs/path/lectin_report.py,OUT_DIR=/abs/out,FASTA_LIST=/abs/fasta_list.txt,SIGNALP_FILE=/abs/sp.txt \
  analysis/lectin_class_report/submit_lectin_report.sbatch
```

The wrapper uses `${SCRATCH:?}`, a `trap`, and copies the results back with a temporary name and `mv`. It does not pass `--name-list`; add it to the wrapper if you need it.

4. Memory: the duplicate-ID check keeps all IDs in memory (about 100 bytes per protein; about 7 million proteins in the 653 proteomes). Use `--no-dup-check` if memory is short. The repeated-ID check inside one file is always on.
5. Check the result with `compare_committed.py` against `review2/whole_proteome_hits.tsv` (PF07938, FB_lectin, LDL, Gal-bind_lectin, CCL2-like).

## Behaviour to know

- A hit is a protein with a sequence score at or above the sequence GA (`--tblout`). The tool's own Pfam table uses `--domtblout`. The column `domain_ge_ga` shows whether the best domain also reaches the domain GA.
- A repeated ID inside one proteome is counted once (best score kept). A warning is printed and written in the report.
- An empty FASTA file gives 0 proteins and a warning. The run continues.
- The same ID in two proteomes gives a warning. The run continues.
- `.gz` and `.zst` files are decompressed into the temp directory for the search and removed after.
- H and U values in the report table are copied from the spec. The script does not measure them.
