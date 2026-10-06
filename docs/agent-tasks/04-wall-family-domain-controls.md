# Task 04. Controls for the wall-family domain call (`wall_family_domain`)

*Read `COMMON-RULES.md` first. Written 2026-10-06.*

## Goal

For each Pfam family in `data/sorting_hat/family_table.tsv`, build a reviewed list of **members** and
**non-members**, so the owner can decide whether to activate the family.

## Why this matters

- The call is `pfam_adhesion.hit`. All 15 families start **inactive**. The owner signs off each
  family (decision C4). Until a family is active, the module writes `unavailable`.
- Domain presence is not function. CFEM, Bys1, hydrophobin, Als and Pir families are linked to
  adhesion or wall function in at least one species. The function is not shown for every member.
- Pth11-like receptors carry CFEM. The `no_tm` rule removes them (helices that start inside the first
  35 residues do not count).
- Specificity is never a number from absence in a database. The plan allows **counts of reviewed
  non-member hits** only.

## The families

Read `data/sorting_hat/family_table.tsv` (15 families; columns include accession, model name,
second condition `no_tm` or `signal_peptide`, `active`). Model names were checked against Pfam 38.2
on 2026-10-05 (`hmmfetch -f` by name gives 15 models). Do not change the table. Propose changes in
your report.

## Positive controls (members)

A protein is a member when it has the full-length domain and the **conserved pattern**, or a curated
function:
- 8 Cys for CFEM and for hydrophobins (check the pattern for each family; cite the paper).
- Members known by function: 0 to 3 per proteome today. Find more from papers (deletion or
  localisation phenotypes). Each needs a PMID and a quote.
- `stratum`: `member_function` (curated function), `member_pattern` (pattern only).

## Negative controls (non-members)

- `domain_without_pattern`: domain hit but the conserved pattern is absent (these are errors, not
  members).
- `pattern_positive_uncharacterised`: pattern present, function unknown. Keep as **their own class**.
  They are neither positive nor negative.
- `functionally_unrelated`: proteins with the domain whose characterised function is unrelated to
  wall or adhesion (cite the paper).

## What to run

1. `hmmsearch --cut_ga` with the family models on proteomes of **four clades** (use proteomes that
   exist locally: the Fungi_5k folder `/bigdata/stajichlab/shared/projects/Fungi_5k` (look for its
   `input` folder and `samples.csv`) or the 831-proteome set of `analysis/pf28404_family/`; check the
   path and the file names before you start). Use model **names**, not accessions, with `hmmfetch -f`.
   Write a SLURM script that follows the rules in `COMMON-RULES.md` and run it as one job per
   proteome set. Do not run on the login node.
2. For each family, list every hit, mark members and non-members by the rules above, and count.
3. A person must read the non-member hits. Prepare a review sheet (`review_<family>.tsv`) with
   gene, species, domain coordinates, pattern check, signal peptide, TM count, and a blank
   `reviewer_call` column for the owner.

## Size target

Per family: how many reviewed members, reviewed non-members and uncharacterised pattern-positive
hits, per clade. Most families will not reach 20 clusters of members. Say so. Do not pad.

## Acceptance checks

1. Every member and non-member row has a source or a rule (pattern).
2. No row is labelled by the module's own call.
3. Hits are listed with the HMM name and Pfam release path (`Pfam38.2/`).
4. Separate counts for each of the four clades.

## Do not

- Do not activate a family. Do not edit `family_table.tsv`.
- Do not report specificity as a percentage.
- Do not treat every domain hit as a member.
