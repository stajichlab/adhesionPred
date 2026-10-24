# Custom Ecm33 HMM, `all` model, 2026-10-23

Custom (non-Pfam) profile HMM of Ecm33, built by the ecm33 build job 30287890 (unit `all`: all 481 phmmer members of Q5AGC4, UniRef50 2025_03). Versioned by folder name and sha256 (`SHA256SUMS`).

- `model.hmm` (233,514 bytes, NAME `ecm33_all`, sha256 `9a6e2c81...c50937`). It has NO GA line.
- `cutoff.json`: bit-score cutoff 29.6, set from the false-hit side only under `analysis/wall_fusion_truth/pr6_ecm33_freeze.json` (sha256 `a7553076d4963e4daedf3d7fdf0eba28b8e8349847b2644482fc2e854d7ce915`): 4 extra calls over the four kept proteomes (limit 16.2455), 0 hard-negative calls.
- `provenance.txt`: job, commit, tool versions.

Status: custom, REPORT ONLY, UNVALIDATED beyond a circular held-out test (owner decision 2026-10-23). It finds all 469 held-out members in 11 clades and both true Ecm33s among the shared negatives (P38248, Q4WNS8), but the held-out members come from the same search that supplied the training, the iteration comparator ties it, and no unit was testable. Record: `docs/reports/2026-10-23-pr6-ecm33-build.md`.

Use: only through `analysis/wall_fusion_report/` (`wall_fusion_custom_models.tsv`, `hmmsearch -T 29.6`). It changes no call, tier, status, `family_table.tsv`, `family_measure.tsv` or `narrow_ok`.

## Owner decision 2026-10-24

The owner decided to KEEP this model, report-only, and to record that it may need additional development. It is unvalidated beyond a circular held-out test (the 469 held-out members come from the same search that supplied the training) and a PANTHER agreement test (verdict `not testable (u < 3): u = 1`; 14 of 15 headline clusters), and it ties seed-only jackhmmer (14 of 15 headline, 185 of 187 in all). The model, `cutoff.json` and `SHA256SUMS` are unchanged. Options for further development: `docs/reports/2026-10-24-pr6-decision-note.md`.
