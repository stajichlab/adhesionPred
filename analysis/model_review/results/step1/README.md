# Step 1 measurement (plan for #9), SLURM job 29301182, 2026-09-30

Run with `analysis/model_review/step1_pooling_scaler.sbatch` at commit 96e791d, 8M model only.
Training embeddings include the 167 exact duplicates. "Legacy pooling" is the emulation (batch 4,
CPU, file order); the real Feb 2026 batching is unknown.

- `pooling_mismatch_summary.tsv`: shipped 8M pickle on S288C under both poolings.
- `pooling_mismatch_scores.tsv.gz`: per-protein scores.
- `cv_scaler1_results.*`: review CV with StandardScaler. `cv_scaler0_results.*`: without.
- Interpretation: issue #25 comment 2026-09-30.
