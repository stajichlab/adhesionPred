# Full run of lectin_report.py on the 653 discovery proteomes (2026-10-10)

- Input list: column `fasta` of `analysis/fungi_scan/selection_discovery.tsv` (653 rows, all files present).
- Pfam: 38.2 (`/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm`). HMMER 3.4. `hmmsearch --cut_ga`.
- Timing test on the login node: 5 proteomes (every 131st of the list), 7 s wall time. This gave about 1.4 s per proteome, so one job of about 15 to 20 min for all 653. One job was used, not an array.
- Job 29846087 (partition short, 4 CPUs, 16 G requested): COMPLETED, wall time 00:16:04, total CPU 13:02, MaxRSS 0.76 GB. Queue wait: about 1 s. Other jobs of this user were pending with AssocGrpCpuLimit at submit time; this job started anyway.
- Wrapper change: `-c 8` to `-c 4`, `--mem=8G` to `--mem=16G`. No failed jobs, no re-submits.
- Outputs (not committed, in `_workdir/lectin_full_run/out/`): hits 1.6 MB, counts 0.7 MB, report 0.14 MB.
- Comparison: `compare.txt` (compare_committed.py, five models, with --names of all 653 proteomes). `model_summary.tsv` has all 15 models. CBM_14 and Chitin_bind_1 were compared with `analysis/lectin_class_spec/r3/cbm_whole_proteome_hits.tsv` by `cbm_cmp.py`. Result: all seven compared models have the same protein sets and equal scores.
- Duplicate protein ID warnings: none in the log or the report.
- Note: the sbatch calls `module load hmmer/3.4` and it worked in the job. On the login node, a `module load` piped into another command does not change PATH; this is a shell effect, not a wrapper fault.
