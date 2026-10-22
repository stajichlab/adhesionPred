# ecm33 pre-freeze check of the shared hard negatives (job 30284040)

Run record of `analysis/wall_fusion_truth/pr6_ecm33_negatives_precheck.sbatch` (code `6fcdd0d2713eb5881e0d99ffcf4abe8a99cbc453`, HMMER 3.4; partition `short`, 2 cores, 2 GB; elapsed 9 min 19 s, exit 0). The job ran `phmmer` of Q5AGC4 (the Ecm33 query) against the 161 shared hard negatives (`data/controls/hard-negatives-shared/`) at E <= 1e-5 (the freeze's labelled_rule definition of a family member) and, separately, at E <= 1e-3 (descriptive).

| File | What it is |
|---|---|
| `negatives_family_table.tsv` | the two family members among the negatives: P38248 (S. cerevisiae ECM33, E 7.6e-56, already in the freeze) and Q4WNS8 (A. fumigatus ecm33, E 8e-48, not in the freeze) |
| `weak_E1e-3.tsv` | descriptive only, never in the list: gel1 P0C7S9 (A. fumigatus GH72), E 4e-4, score 12.7 |
| `proposed_negatives_ignored.json` | the edit the script proposed |
| `provenance.txt`, `summary.txt` | pin, input hashes, settings, the job's summary |

Result: exactly two family members; no hidden member. The raw `phmmer` tables (`phmmer_E1e-5.tbl`, `phmmer_E1e-3.tbl`) stay on HPCC. No sequences are in this folder. The freeze amendment that follows from this result is recorded in `analysis/wall_fusion_truth/pr6_ecm33_freeze.json` (`amendments`).
