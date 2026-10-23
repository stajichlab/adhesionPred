# PR 6 ecm33 build: run record of jobs 30285974 (extract) and 30287890 (build)

Small tables copied back from HPCC after the extract job and the build array of the Ecm33 custom-HMM experiment, plus the local report tables. Freeze: `analysis/wall_fusion_truth/pr6_ecm33_freeze.json`. Spec: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md` (sections 4.5, 4.6, 4.8). Report: `docs/reports/2026-10-23-pr6-ecm33-build.md`. Earlier record (members come from it): `../pr6_gap_uniref50_2026-10-22/`.

## What the run is tied to

| item | value |
|---|---|
| code pin | `ba15e0c83ed52b2636ef6837ad51e0edf3d4c698` (`ba15e0c8`, the merge that carries the signed freeze and its Q4WNS8 amendment) |
| freeze sha256 | `a7553076d4963e4daedf3d7fdf0eba28b8e8349847b2644482fc2e854d7ce915` (file as read by the report) |
| extract job | 30285974 (`pr6_ecm33_extract`): 481 members of the phmmer set of Q5AGC4 against UniRef50 2025_03, 19 fine clades, 11 units |
| build job | 30287890, `sbatch --array=0-22%4`, partition `batch`, 8 cores, 16 GB, 12 h per task. 23 tasks (11 units x 2 independence modes + `all`), all COMPLETED exit 0, elapsed 56 min to 3 h 04 min, no FATAL in any log |
| settings | MAFFT 7.505 L-INS-i, `hmmbuild --amino`, HMMER 3.4, MMseqs2 `b804fbe3` (cluster30: 30 percent identity, 50 percent coverage), Pfam 38.2 PF12454; `min_u` 3; extra-call limit 5 per 10,000 (16.2455 over 32,491 proteins); hard-negative rate 0.05; negatives ignored P38248, Q4WNS8 |
| report | run locally: `python3 -I analysis/wall_fusion_truth/pr6_ecm33_report.py --results <R>/ecm33_build --freeze analysis/wall_fusion_truth/pr6_ecm33_freeze.json --units <R>/members/units.tsv --out <R>/report` |

## Files

| path | what it is |
|---|---|
| `extract_provenance.txt` | provenance of the extract job (pin, UniRef50 file hash, member counts) |
| `members/` | `members.tsv` (481 members with clade), `excluded.tsv`, `units.tsv` (19 clades, 11 units), `member_ids.txt`, `summary.txt` |
| `report/` | `recall_by_unit.tsv` (22 unit/mode rows), `side_by_side.tsv` (11 units, headline against second column), `extra_calls.tsv`, `hard_negative_calls.tsv` (198 rows), `ignored_negatives_line.tsv` (44 rows), `summary.txt` |
| `ecm33_build/<unit>/<mode>/` | one folder per unit and mode (`cluster30` headline, `none` descriptive); the `all` model is `ecm33_build/all/none/` |
| `.../cutoff.json`, `relaxed_cutoff.json` | frozen-rule cutoffs of the custom model and of the relaxed Pfam comparator (each with the freeze sha256 and `positives_used: false`) |
| `.../manifest.json` | held-out ids, training ids with clade and 30 percent cluster, excluded neighbours (so `train_ids.txt` is not copied) |
| `.../hits_custom.tsv`, `hits_pfam_ga.tsv`, `hits_pfam_relaxed.tsv` | held-out member hits (custom and relaxed Pfam with scores; Pfam GA ids only) |
| `.../hits_iter.tsv.gz` | iteration comparator (jackhmmer seeded by the training alignment) hits with E-values, gzipped (8.4 MB before gzip) |
| `.../neg_hits.tsv`, `neg_hits_relaxed.tsv`, `proteome_hits.tsv` | shared-negative and kept-proteome hits |
| `.../heldout_ids.txt`, `hmmbuild.log`, `provenance.txt` | held-out ids, build log, per-task provenance (pin, tool versions, model.hmm sha256) |

All files hold identifiers, scores and labels only. No sequence is committed (a test checks it). Members are UniRef50 representatives (`UniRef50_<accession>`).

## Left on HPCC (not committed)

Per unit and mode `train.aln.fa` (the training alignment, 63 MB in all) and `model.hmm` (5 MB in all), under the build job work directory (`.../_workdir/` of the HPCC clone; the job's `$TMPDIR` copies were scratch and are gone). `ecm33_members.faa` (the 481 member sequences; sha256 `b5e9110d7e315d4035ab16ab1ed7f0f294d410bf603affb429f6d8de3f3eb8ee`) stays there too. `hits_iter.tsv` is committed gzipped (the whole record is about 2.7 MB), so nothing of it stays only on HPCC.

Local copy of the fetched results (outside any repo): `/Users/jstajich/projects/adhesionPred-ecm33-results/`.
