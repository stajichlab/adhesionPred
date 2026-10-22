# Notes that are not checked

*Written by `build_control_set.py` (DRAFT, Task 5.2).*

- **Clusters.** `clusters.tsv` uses the Kiss `cluster_number` per slice (`KISS-<slice>-<n>`). It
  is an orthology and MCL family, not a 30% identity cluster. COMMON-RULES rule 4 asks for
  MMseqs2. Task 5.3 makes them. Until then every count of clusters is a proxy.
- **Leakage.** `leakage` is `unknown` for every row: nobody has compared these proteins with the
  Pfam seed sequences yet. Task 5.3 (`seed_split.py`) sets `none` or `in_reference`,
  `seed_class_slice` and `seed_class_all39`.
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
