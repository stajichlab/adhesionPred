# Notes that are not checked

*Written by `build_control_set.py` (Task 5.2, updated with the Task 5.3 cluster split).*

- **Clusters.** `clusters.tsv` holds MMseqs2 clusters (30% identity, 50% coverage; job 30229366,
  record in `docs/reports/data/kiss_cluster_split_2026-10-22/`). The clusters were made over the
  132 controls and the 161 shared negatives together; a cluster that holds both is
  `mixed_pos_neg` in `clusters_split.tsv`. The independent unit is this cluster, not the Kiss
  `cluster_number` (kept in `kiss_cluster_number`).
- **Leakage.** `leakage` is `none` or `in_reference` from the identity to the nearest Pfam seed
  sequence (`seed_split.py`: any member at 30% identity or more makes a cluster `seed-like`).
  The seed sequences are the ones InterPro served on 2026-10-22; the Pfam release they come from
  was not compared with 38.2 row by row. Five models have very few seed sequences
  (`seeds_summary.tsv`: Hydrophobin_D 3, Glyco_hydro_17 6, Chitin_synth_2 7, Hyd1F 8), so a
  `seed-independent` label for them rests on little.
- **portal_validated rows.** These rows have no name evidence of their own. They are accepted by
  the owner (2026-10-22) because their portal agrees on other rows. The flag stays; filter on it
  if a result looks sensitive.
- **locus_confirmed rows.** The UniProt entry found by the JGI header locus tag has that tag as its
  only gene name. The file protein equals the entry. This does not show that the Kiss row means
  this protein.
- **Sequence differs from JGI.** Where `sequence_matches_jgi` is `no`, the committed sequence is
  the UniProt one, not the Kiss (JGI) one.
- **Accession found by sequence.** For a row with no accession by name, the accession is taken from
  the one cached UniProt entry whose sequence has the same sha256 as the JGI protein (the `notes`
  column says so). The Kiss gene name was not used and may differ from the UniProt gene name.
- **Species names.** `species` and `taxon_id` come from NCBI Taxonomy answers cached by `map_kiss.py`
  (`results/organism_resolution.tsv`). They were not checked against `names.dmp`.
- **Wall rows.** The wall rows are a descriptive panel for the wall models and a false-hit panel for
  the narrow call. They are positives for the wall slice only.
- **Family representatives.** The 35 Afu rows stand for a Kiss family. The Afu protein is the
  ortholog named in the Kiss table, not the protein that the authors measured.
- **Bcine1.** Two Botrytis rows (bhp2, bhp3) are not controls: the 2025 annotation renumbered the
  JGI IDs (see `rejected.tsv`).
