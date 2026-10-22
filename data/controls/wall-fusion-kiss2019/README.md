# Control set: wall, fusion and adhesion genes of Kiss et al. 2019

*Written by `analysis/wall_fusion_truth/build_control_set.py` (DRAFT, plan Task 5.2). Do not edit
the files by hand; run the script again.*

Source: Kiss et al. 2019, Nat Commun (PMID:31501435), Supplementary Data 3, 4 and 8, as read in
`docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv` and mapped to sequences by
`analysis/wall_fusion_truth/map_kiss.py` (Task 5.1; report
`docs/reports/2026-10-18-kiss-feasibility.md`). The owner decisions of 2026-10-22 (section 10 of
that report) bind this set.

## Files

| File | Content |
|---|---|
| `controls.tsv` | 112 positive controls with a public UniProt sequence. Checked by the L1 validator. |
| `controls_local_sequence.tsv` | 20 positive controls whose sequence is not public. Same columns. The validator does not read it. **These two files are one set: anything that reads the set must read both.** |
| `sequences.faa` | UniProt sequences of `controls.tsv`. Header = accession. |
| `clusters.tsv` | `accession<TAB>cluster_id`. Kiss cluster per slice: a proxy until Task 5.3. |
| `family_reps.tsv`, `family_rep_sequences.faa` | 35 family representatives (Afu orthologs). Training seeds and baseline report only; not clustered; each keeps its Kiss `cluster_number`; a flagged group in the report. |
| `rejected.tsv` | 13 Kiss rows that are not controls, with the reason. |
| `counts.md`, `unverified_notes.md` | counts and unchecked claims. |

## Keys

The key (`accession`) is the UniProt accession wherever the row has one (owner preference: NCBI or
UniProt accessions, not JGI IDs). A row with no UniProt accession gets the local key
`jgi:<portal>:<proteinId>` (for example `jgi:Neucr2:4224`). A local key is not a public identifier.
Its sequence is not in this repository (JGI data are never committed); `sequence_sha256` is the
check, and `build_kiss_fasta.py` reads the sequence from the downloaded JGI file.

## Columns

The required columns of `docs/agent-tasks/COMMON-RULES.md`, then: `control_role` (`kiss_row`),
`serves` (the calls the row may be used for), `clade_coarse`, `clade_fine`, `strain_taxid`,
`kiss_cluster_number`, `cluster_id`, `cluster_stratum_fine` and `cluster_stratum_coarse`
(`shared`: the Kiss cluster of this slice has strictly mapped members in another clade at that
level; `clade_only`: it has not), `jgi_accept` (`row_verified`, `locus_confirmed`,
`portal_validated`, or empty for a row mapped by a locus tag of its own), `map_status`,
`protein_id_jgi`, `jgi_seq_sha256`, `sequence_sha256`, `sequence_source`, `sequence_matches_jgi`,
`sequence_in_repo`, `seed_class_slice` and `seed_class_all39` (empty until Task 5.3; `seed-like` when the
Kiss cluster has any member at 30% identity or more to a Pfam seed: of the row's own slice models, or of
all 39 models).

`label` is `positive` for the call named in `stratum` (adhesion, wall or fusion). Read the rules
before you use a row:

- adhesion rows: positives of the adhesion-family recall;
- wall rows: positives for the wall models only, and a descriptive false-hit panel for the
  narrow call. A wall row is never a negative and never a truth row of the adhesion call;
- fusion rows: descriptive only (at most 2 clusters per clade: no floor is met);
- no number is pooled across slices (D20).

## Rules of this set

- All strictly mapped rows are accepted (owner, 2026-10-22). The three `jgi_accept` flags stay;
  filter later if a result looks sensitive.
- hyr1, bhp2 and bhp3 are not controls (`name_conflict`, listed in `rejected.tsv`).
- The 35 family representatives are not in `controls.tsv`. They may be training seeds and appear
  in the baseline recall report; they are never test controls and never count toward a floor.
  They are not clustered (the `--include-family-reps` option of `build_kiss_fasta.py` exists and is not used).
- The hard negatives for this set are the shared set of 161 rows, and only that set (owner decision
  2026-10-22; `data/curated/adhesins/hard_negative_seeds.tsv` is excluded and may be added later).
  Eight proteins are both a control here and a shared negative; they keep both roles (the negative
  copy is `<accession>#neg` in the clustering input) and their clusters are flagged `mixed_pos_neg`.
- This set is two files. Anything that reads it must read `controls.tsv` and `controls_local_sequence.tsv`.
- Two Kiss rows that point to the same protein are one row here (`notes` says which).
- `leakage` is `unknown` until `seed_split.py` has compared each row with the Pfam seeds.

## How it was built

```
python3 -I analysis/wall_fusion_truth/build_control_set.py
```

The script is offline and deterministic (two runs give the same bytes). It reads
`analysis/wall_fusion_truth/results/{row_mapping,organism_resolution}.tsv`,
`results/query_cache.json` (the UniProt sequences) and the Kiss table. It never reads the JGI
files. After Task 5.3 run it again with `--clusters-tsv` and `--split-members` (written by
`seed_split.py`) to replace the proxy clusters and fill `leakage`, `seed_class_slice` and `seed_class_all39`.

Check: `python3 -m cellsurface_sorting_hat.controls data/controls/wall-fusion-kiss2019`
(with `PYTHONPATH=src`).
