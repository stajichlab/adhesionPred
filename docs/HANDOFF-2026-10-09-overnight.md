# Hand-off: overnight work, 2026-10-08 to 2026-10-09

*Branch `fungi-scan-and-specs` (from `hydrophobin-validation`; local until the owner says otherwise). Owner instructions: scan proteomes from Fungi_5k, 40 to 60 jobs on `short` and `short_gpu`; stop at any owner-decision point and at any failed review gate; commit as tasks complete; write this file. No merges, no deletions, no new PRs.*

## Status log (newest last)

- 2026-10-08 late: `surface_attachment_candidate` call added and pushed to `hydrophobin-validation` (`03dcc23`).
- Descriptive scan set up: 56 proteomes from Fungi_5k (`analysis/fungi_scan/selection.tsv`; 24 Onygenales, 12 Aspergillus, 3 Candida, others). 15 grouped SLURM jobs submitted (`_workdir/fungi_scan/jobs.txt`).
