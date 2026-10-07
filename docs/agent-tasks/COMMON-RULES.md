# Common rules for every evidence-set task

*Read this file first. Each task file in this folder points here. Written 2026-10-06.*

## Who you work for and what you produce

You build **control datasets**: lists of proteins that are known positives or known negatives for
one line of evidence. The sorting tool (`cellsurface_sorting_hat`) uses them to measure how well a
module works. You do not train models. You do not change module code. You do not set thresholds.

## What counts as a control

| Term | Meaning |
|---|---|
| Positive control | A protein for which an independent source shows that the module should call it. |
| Negative control | A protein for which an independent source gives a **reason** that the module should not call it. "Not annotated" or "not found in a database" is **not** a reason. |
| Independent source | A paper, an experiment, a curated database entry, or a GO annotation with an experimental evidence code. A prediction from a tool is not independent of that tool. A prediction from the module under test is never a source. |

## Rules that apply to every task

1. **Every row has a source.** Give a PMID, a database record ID, or a file and row. Quote the
   sentence or field that supports the label. If you cannot find a source, do not add the row.
2. **Do not invent.** Do not fill a gap from memory. A claim from memory that you cannot check goes
   in `unverified_notes.md` with the words "from memory, not checked". It does not go in the table.
3. **Do not use the module under test to choose controls.** Do not pick positives because the module
   calls them, or negatives because it does not. This is circular. If a protein was used to tune the
   module, mark `tuned=yes` (see column list).
4. **Group by homology.** Cluster all controls of one dataset with MMseqs2 (`easy-cluster`,
   `--min-seq-id 0.3 -c 0.5`). Report counts of **rows, sequences and clusters** for each class.
   The independent unit is the cluster, not the protein. On HPCC: `module load MMseqs2/17-b804f`.
5. **One species per calibration.** A status is measured per species. Give the species (NCBI
   scientific name and taxon ID) for every row. Check names against `names.dmp`
   (`/srv/projects/db/taxonomy/names.dmp`). A synonym such as `Ustilago maydis` has no taxon ID; use
   the scientific name (`Mycosarcoma maydis`).
6. **Leakage.** For each row say whether the protein or a homolog (30% identity) was used to build,
   tune or train the module, or is in the module's reference set. Use the values `none`, `partial`,
   `tuned_on_truth`, `in_reference`, `unknown`. Do not hide `unknown`.
7. **Size.** The rule for `estimated` status needs at least 20 positives and 20 negatives, at least 20
   independent clusters of each, and a 95% interval half-width of at most 0.10 for **both**
   sensitivity and specificity. A binomial approximation for the clusters needed for a half-width of
   0.10 is about 61 clusters at sensitivity 0.8 and about 96 at 0.5 (it ignores cluster structure, so
   it is a rough guide). For specificity near 0.95 the 20-cluster floor is the limit. Report how far
   your set is from these targets. Do not pad a set to reach them.
8. **Say what you could not do.** End each report with: rows added, rows rejected and why, target
   met or not, and open questions for the owner.

## Output files

Write to `data/controls/<task-name>/` (create it). Use TSV, UTF-8, one header row, no BOM. Large
tables: `.tsv.gz`. Files:

| File | Content |
|---|---|
| `controls.tsv` | One row per protein (columns below). |
| `sequences.faa` | Protein sequences of every row (from UniProt or the genome source named in `source`). Header = `accession`. |
| `clusters.tsv` | `accession<TAB>cluster_id` from MMseqs2. |
| `counts.md` | Table of rows, sequences and clusters per class and per species, and the size targets. |
| `rejected.tsv` | Candidates that you looked at and did not add, with the reason. |
| `unverified_notes.md` | Claims you could not check. |
| `README.md` | How the set was built, the commands, the date, the versions of the databases. |

Required columns of `controls.tsv`:

`accession`, `gene`, `species`, `taxon_id`, `label` (`positive` or `negative`), `stratum`
(task-specific, see the task file), `evidence_level` (task-specific), `source` (PMID or record),
`quote` (the supporting text, short), `tuned` (`yes` or `no`), `leakage` (one of the five values),
`curator` (who or which agent), `date` (YYYY-MM-DD), `notes`.

## Style and tools

- Text: Simplified Technical English. Short sentences. One idea each. No idioms.
- Python: `/usr/bin/python3.12`. On the login node, `python` and `pip` are Python 3.9; do not use them.
- Compress large text output (`zstd -T0` for internal files; `.gz` for files others will read).
- Do not commit publisher texts (copyright). Keep them in `_workdir/curation_texts_archive/` and
  commit only their hashes.
- On SLURM: use `${SCRATCH:?}` for scratch (node-local), a `trap` to remove it, copy results back to
  `/bigdata` with a temporary name and `mv`. Never `BASH_SOURCE`. Never `/tmp`.
- Work on a new branch. Do not merge. Do not push to `main`. Push the branch over HTTPS and open a
  PR for the owner to review.

## Review gate

A control set is not used for calibration until the owner has read `counts.md` and
`unverified_notes.md` and has accepted it in the PR. State in the PR text which rows are weakest.
