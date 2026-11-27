# Model review analyses (2026-09-27)

Report: `docs/model-review/2026-09-27-review-and-framework-plan.md`.

Scripts `01_embed_esm2_8M.py`, `02_cv_and_proteome_eval.py`, `03_pooling_mismatch.py`, `run.sh` and `step1_pooling_scaler.sbatch`
used the `surface_glyco` package, which was removed on 2026-11-26. They and the package are at tag `legacy-2026-11-26`. Their
outputs (`cv_results.tsv`, `results.txt`, `s288c_scores.tsv`, `results/`) stay here as the record the report cites.

The remaining scripts (stage 2 ESM C prototype, `esm_vs_hmm.py`, repeat-structure transfer, clade status) do not import that
package.
