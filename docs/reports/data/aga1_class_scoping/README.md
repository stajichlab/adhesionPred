# AGA1-class scoping study: tables and scripts (2026-10-24)

Local scoping study for a model of the Ser/Thr/Asn/Gln/Gly/Pro-rich, repeat-bearing, GPI-anchored adhesin-like class
(the class `AGA1` belongs to). Report: `docs/reports/2026-10-24-aga1-class-scoping.md`. Descriptive only: no call, tier,
panel, `narrow_ok` or freeze was changed. Nothing here ran on HPCC or SLURM. The tables hold ids, labels and numbers;
**no sequence is committed** (a test checks this).

## Tables (all tab-separated, one row per protein unless stated)

| File | Content |
|---|---|
| `proteins.tsv` | 850 proteins: accession, gene, organism, clade (coarse, fine, family), `role`, `pos_tier`, curated class, evidence tier, mechanism label, truth_v3 repeat label, hard-negative stratum, candidate narrow class, curated Pfam list, UniProt GPI / signal / repeat-feature annotation, sequence source, source tables |
| `features.tsv` | local sequence features and the numbers of repeat detectors 02 and 14, the existing `tandem_repeat_protein` rule and the exploratory z_seq 3.5 rule |
| `clusters.tsv` | MMseqs2 cluster (30 percent identity, 50 percent coverage) of every protein |
| `source_counts.tsv`, `sequence_availability.tsv`, `class_clade_counts.tsv`, `lectin_overlap.tsv`, `kiss_adhesion_rows.tsv` | inventory counts |
| `baseline_lco.tsv`, `baseline_lcl.tsv`, `single_feature_auc.tsv`, `repeat_call_vs_features.tsv`, `matched_operating_point.tsv`, `oof_scores.tsv` | baseline separation results |
| `repeat_rule_reproduction.tsv` | the existing repeat call recomputed on truth_v3 (reproduces S288C 0.435 / 0.987, C. albicans 0.462 / 0.976) |
| `aga1_vs_named.tsv`, `narrow_slices.tsv` | AGA1 beside named proteins; candidate narrow-class slices |

## Definitions fixed before any model was run

* **role = class_positive**: curated class `adhesin` AND repeat-bearing, meaning the curators' mechanism label is `2a`
  (repeat/avidity; `data/controls/repeat-mechanism/adjudication/curation_table.consensus.tsv`, or `mechanism_label` of the
  other-species table) OR the truth_v3 repeat label is 1. **Tier A** = independent repeat evidence (truth_v3 label 1, or
  mechanism basis `stated_in_paper` / `database_annotation`); **tier B** = mechanism 2a by `family_inference` only.
  GPI and signal peptide are UniProt feature annotations, reported but not part of the rule.
* **role = hard_negative**: curated hard negative in any source table (hard-negatives-shared, other-species, repeat-mechanism,
  truth_v3). Subsets: `hn_st30` (some 30-residue window at least 50 percent S+T), `hn_stfrac25` (whole-protein S+T at least 25
  percent), `hn_repeat` (at least 2 UniProt repeat features). **bg**: S288C and C. albicans proteome proteins with no curated
  class (negatives ASSUMED). Roles `adhesin_other_mechanism`, `excluded_*`, `uncurated_repeat_bearing` and Kiss rows are not
  used as positives or negatives.
* Narrow classes (candidates): by gene name first (AGA1, ALS, FLO11/MUC1, FLO1/5/9/10, HWP/EAP/RBT, HYR/IFF/RBR, EPA, SOWgp,
  PFL), then by curated Pfam (Als domains, Flo11, GLEYA, Flocculin/Flocculin_t3/PA14, Hyphal_reg_CWP/Hyr1). First match wins; the
  order is an assumption. Rule for "could support a class-specific model" (provisional): strict = at least 5 positives in each
  of at least 2 clades; weak = at least 5 in total over at least 2 clades.
* Clade = fine clade of `analysis/wall_fusion_truth/map_kiss.py` (Saccharomycotina; class for Pezizomycotina and Basidiomycota;
  phylum otherwise), read from the UniProt lineage. Family level is reported as an extra, outside the rule.

## Exact commands

Run from the repo root of the worktree. `PY` = the repo environment python
(`/Users/jstajich/projects/adhesionPred/.venv/bin/python`; it has numpy, scikit-learn and biopython). Scripts 00, 01, 03, 05, 07, 08 are stdlib only.

```
D=docs/reports/data/aga1_class_scoping/scripts
python3 $D/00_fetch_uniprot.py           # API calls, see below
python3 $D/01_build_proteins.py
$PY     $D/02_run_detectors.py           # detectors 02 and 14 from analysis/cocci_repeats, local CPU, about 3 s
python3 $D/03_features.py   (or $PY)
python3 $D/04_cluster.py                 # needs mmseqs on PATH (v18-8cc5c used)
python3 $D/05_inventory.py
(cd $D && $PY 06_baseline.py)            # about 1 minute
python3 $D/07_aga1_and_slices.py
python3 $D/08_repeat_rule_reproduction.py
```

Order note: 01 reads the cache from 00; 02 and 04 write only into the cache (FASTA files stay OUTSIDE the repo); 01 must be re-run
before 05 if proteins change. Cache: `~/.cache/adhesionPred/aga1_class_scoping` (override with `AGA1_SCOPING_CACHE`).

## API calls (read-only public UniProt REST, 2026-10-24, UTC about 16:39 to 16:41)

`GET https://rest.uniprot.org/uniprotkb/accessions?accessions=<100 ids>&format=tsv&fields=accession,reviewed,gene_primary,organism_id,organism_name,length,sequence,lineage,ft_signal,ft_lipid,ft_repeat&size=100`,
1 second between requests, 10 requests in 3 runs of the script (the first run, 549 accessions, used an earlier accession list;
later runs fetched what was still missing). Only accessions already named in committed tables were requested
(truth_v3 S288C and C. albicans, repeat-mechanism curation table, hard-negatives-shared, other-species tables): 778 UniProt
accessions. The two ids that are not UniProt accessions (`ACE80261`, `XP_045276491.1`) were not fetched; they have local
sequences, and their clades were set by hand from the species named in the table (`MANUAL_CLADE` in script 01).
Sequences of the 850 proteins come from the committed `data/controls/*/sequences.faa` where present (see
`sequence_availability.tsv`), otherwise from the cache. Kiss JGI-only rows have no committed sequence and are not in the study.

## Features that need HPCC or an external tool later (not computed here)

SignalP 6 (signal peptide), a GPI predictor (PredGPI or similar; here only the UniProt annotation is used), IUPred or AIUPred
(disorder), glycosylation predictors (here only a crude N-X-S/T sequon count), ESM-2 embeddings, a repeat-unit profile, OrthoFinder
(orthology within a clade).
