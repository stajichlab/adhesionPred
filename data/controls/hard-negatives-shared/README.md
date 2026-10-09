# Task 07: shared hard-negative set

*Built 2026-10-08 by the task 07 agent (claude-opus-5-5). Task file:
`docs/agent-tasks/07-hard-negatives-shared.md`. Rules: `docs/agent-tasks/COMMON-RULES.md`. Not
committed. Not reviewed by the owner. Do not use for calibration before the review gate.*

## Files

| File | Content |
|---|---|
| `controls.tsv` | 161 hard negatives. Required columns, then extra columns (see below). |
| `sequences.faa` | 161 sequences. Header = accession. UniProt, except one NCBI protein (XP_045276491.1). |
| `clusters.tsv` | `accession`, `cluster_id`, `set`. MMseqs2 clusters of the negatives **and** 104 adhesin positives. |
| `counts.md` | Rows, sequences and clusters per stratum, clade, species and module; size targets; inventory of existing data. |
| `rejected.tsv` | Candidates that I looked at and did not add, with the reason and source. |
| `unverified_notes.md` | Claims I could not check, open questions, weakest rows. |
| `manual_overrides.tsv` | **Proposed** changes to `adhesins.tsv` hard-negative rows. Not applied. |
| `leakage_pfam_seed_hits.tsv` | Rows that are Pfam seed members or 30% homologs of a seed domain of a `family_table.tsv` family. |
| `candidates_input.tsv` | The curated input list: accession, stratum, origin, manual PMID and quote, notes. |
| `scripts/` | The scripts that made the files (paths inside point to the agent scratch folder; see Commands). |

## Extra columns in `controls.tsv`

`serves` (`step1`, `repeat`, `family_domain`, `adhesion_level`), `clade`, `strain_taxid` (the
UniProt organism taxon), `protein_name`, `origin` (`new`, `adhesins.tsv`, `task08`), `cluster_id`,
`mixed_cluster_with_positive`, `leakage_step1`, `leakage_repeat`, `leakage_family_domain`,
`family_domain_pfam`, `gpi_feature`, `signal_peptide`, `max_ST30` (highest Ser+Thr fraction in a
30-residue window), `uniprot_repeat_features`, `uniprot_status`.

## How the set was built

1. **Inventory.** I read `hard_negative_seeds.tsv`, the 42 `hard_negative` rows of `adhesins.tsv`,
   the 10 task 08 look-alike rows and the Phase C eval table (`counts.md`, last section).
2. **Candidates.** I queried UniProt REST (`rest.uniprot.org`, 2026-10-08) for reviewed fungal
   entries outside Saccharomycotina, per stratum: GPI-anchor keyword KW-0336; GH72 PF03198; EC
   3.2.1.14, 3.2.1.39, 3.2.1.58, 3.2.1.59, 3.2.1.6, 3.5.1.41, 3.2.1.84, 3.2.1.75; WSC PF01822 and
   msb2/hkr1/mid2 gene names; CBM1 PF00734 and CBM20 PF00686; EC 1.10.3.2, 1.16.3.1, PF07731, PF00394;
   the `family_table.tsv` Pfam families; Repeat features with a signal peptide; aspartic protease
   PF00026; reviewed secreted Onygenales (taxon 33183) and Basidiomycota pathogens. I kept entries
   whose FUNCTION comment has ECO:0000269 (experimental) evidence. I added literature rows from
   PubMed (E-utilities and the PubMed tools): Gel2, Rsp3, Sod3, Cts1, Mp1p.
3. **Selection.** I chose rows by the source, not by any tool call. I did not run or read the output
   of any `cellsurface_sorting_hat` module. I rejected a candidate when the protein or an ortholog has
   a reported adhesion function (task rule), when it is an adhesin positive, or when it shares a
   cluster with a positive and is a paralog of that positive (`rejected.tsv`).
4. **Evidence level.** `N1` = a characterised non-adhesive function (UniProt FUNCTION or CATALYTIC
   ACTIVITY with ECO:0000269, or a PMID abstract that states the function). `N2` = adhesin-like
   architecture and no experimental function (UniProt record only), or a function text that only
   says "may" or "appears to" (10 rows set by hand: Q1E3R8, XP_045276491.1, the 4 *Drepanopeziza*
   CFEM rows, G2WTI6, G2WSY1, G2WYA4, Q09920).
5. **Quote.** The PMID quote, or the UniProt FUNCTION text (cut at 400 characters). The `source`
   column names the experimental PMIDs that UniProt cites.
6. **Species.** UniProt organism taxon, mapped up to rank species with `nodes.dmp`; the name is the
   scientific name in `names.dmp` (`/srv/projects/db/taxonomy`, files of 2026-05-12). Example:
   *A. fumigatus* Af293 (330879) becomes *Aspergillus fumigatus* (746128). *Ustilago maydis* appears
   as *Mycosarcoma maydis* (5270).
7. **Clade.** From the NCBI lineage: Saccharomycotina 147537, Taphrinomycotina 451866, Onygenales
   33183, Eurotiales 5042, other Pezizomycotina 147538, Basidiomycota 5204, other fungi.
8. **Strata.** As in the task file. Rule for GPI enzymes: a GPI-anchored remodelling enzyme or GPI
   chitinase goes to `gpi_wall_enzyme`; a secreted non-GPI enzyme goes to `wall_hydrolase`.
   `laccase_etc` holds AA1 multicopper oxidases and other abundant secreted enzymes and proteins
   (proteases, peroxidases, LPMOs, Mp1p, CBP1, SOD3). Two strata are additions:
   `secretory_non_surface` (step 1 look-alikes) and `wall_structural_other` (carried *S. cerevisiae*
   wall proteins that fit no task stratum). See `unverified_notes.md`, question 3.
9. **Clusters.** MMseqs2 17-b804f, `easy-cluster --min-seq-id 0.3 -c 0.5` (default coverage mode),
   on the negatives plus 104 positives: `adhesins.tsv` rows with `cls=adhesin`, E1 or E2 and an
   accession (96), plus the 8 task 08 adhesin rows. UniProt sequences of 2026-10-08.
10. **Leakage.** `leakage_step1`: `none` for all rows. R0 has no tuned cutoffs and no reference set
    (decision 1 in `docs/paper/03-status-and-validation-rules.md`); SignalP 6 training overlap is
    not measured. `leakage_repeat`: the repeat call was tuned on SOWgp; `easy-search` against SOWgp58,
    SOWgp66, SOWgp82 (Q8NK60, Q8NK61, Q96V71) at 30% identity and 50% coverage found no hit, so
    `none`. I did not audit other sources of the repeat detectors. `leakage_family_domain`:
    `in_reference` when the accession is a Pfam 38.2 seed member of a `pfam_adhesion` family;
    `partial` when the row hits a seed domain sequence of such a family at 30% identity over 50% of
    the seed domain (`--cov-mode 1`); else `none`. Seeds came from the InterPro API
    (`/api/entry/pfam/PFxxxxx/?annotation=alignment:seed`). `leakage` is the most severe value over
    the modules in `serves`. `tuned` is `no` for all rows: no row was used to tune a module, as far as
    the repository documents show.

## Commands

```
# scripts live in the agent scratch folder; copies are in scripts/
python3.12 uq.py NAME 'UNIPROT QUERY' dl            # candidate queries (UniProt REST TSV)
python3.12 -I abs.py PMID ...                       # PubMed abstracts (E-utilities efetch)
python3.12 -I build.py candidates.tsv work          # UniProt JSON, features, species, clade
python3.12 -I fetch_fasta.py positives.list positives.faa
module load MMseqs2/17-b804f
mmseqs easy-cluster comb2.faa clu tmpc --min-seq-id 0.3 -c 0.5 -v 1 --threads 4
mmseqs easy-search sequences.faa seed_domains.faa seed_hits.m8 tmps --min-seq-id 0.3 -c 0.5 --cov-mode 1
mmseqs easy-search sequences.faa sowgp.faa sow_hits.m8 tmps2 --min-seq-id 0.3 -c 0.5
python3.12 -I assemble.py work out "P82476,J9VPD7,J9VND2,C8VIP3" \
  "Q1E3R8=N2,XP_045276491.1=N2,K1XVG1=N2,K1WG73=N2,K1XW16=N2,K1XT82=N2,G2WTI6=N2,G2WSY1=N2,G2WYA4=N2,Q09920=N2"
python3.12 -I counts.py out
```

## Versions

| Resource | Version or date |
|---|---|
| UniProtKB (REST) | queried 2026-10-08 |
| Pfam seeds | Pfam 38.2 via InterPro API, 2026-10-08 (matches `family_table.tsv` release 38.2) |
| NCBI taxonomy dumps | `/srv/projects/db/taxonomy`, names.dmp and nodes.dmp of 2026-05-12 |
| PubMed | E-utilities, 2026-10-08 |
| MMseqs2 | 17-b804f (module `MMseqs2/17-b804f`) |

## Result in short

161 rows, 91 clusters. Strata with 20 or more clusters: `laccase_etc` only. Clades with 20 or more
clusters: Eurotiales, other Pezizomycotina, Saccharomycotina. Details in `counts.md`.
