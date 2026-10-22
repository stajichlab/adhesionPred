# PR 6 gap search, round 1 on Swiss-Prot (job 30260197)

*2026-10-22. Author: Claude Code (claude-sonnet-5-5). Run record: `docs/reports/data/pr6_gap_swissprot_2026-10-22/`. Spec: `docs/superpowers/specs/2026-10-22-pr6-custom-hmms-design.md` (decisions 3, 5, 7 and S3). Scripts: `analysis/wall_fusion_truth/pr6_gap_search.sbatch` and `pr6_gap_summary.py`. All numbers below are recomputed from the committed tables.*

## 1. What was run

For the four clusters that have no Pfam hit (step 6.0), `phmmer` (one pass) and `jackhmmer -N 3` (`--incE 1e-5`) searched Swiss-Prot 2025_03 (573,661 reviewed sequences, sha256 `6a9a6c21…`). Every hit was assigned to the held-out clade, outside the clade, non-fungal or the query itself, using the sorting-hat taxdump. The rule of decision 7 (provisional) says a target is `buildable` when the one-pass set has at least 5 outside fungal sequences in total over at least 2 clades. The job also checked the Ecm33 model `PF12454` (section 3). No model was built.

## 2. The four queries

Hits per method, counted from `gap_homologs.tsv`. "Inside" means a hit in the held-out fine clade (descriptive only). "Self" means the query itself is in Swiss-Prot.

| Cluster | Gene | Query | Held-out clade | Mode | Method | Outside | Inside | Self | Non-fungal | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| KS0017 | hyd1_like | A0ACE4KFX7 | Eurotiomycetes | leave-cluster-out (partial) | phmmer | 0 | 1 | 0 | 0 | |
| | | | | | jackhmmer | 1 | 1 | 0 | 0 | `not_buildable`: only 1 outside sequence (minimum 5) |
| KS0080 | aga1 | P32323 | Saccharomycotina | leave-cluster-out (partial) | phmmer | 0 | 0 | 1 | 0 | |
| | | | | | jackhmmer | 0 | 0 | 1 | 0 | `not_buildable`: no outside fungal homolog |
| KS0107 | ece1 | Q07730 | Saccharomycotina | leave-cluster-out (partial) | phmmer | 0 | 2 | 1 | 0 | |
| | | | | | jackhmmer | 0 | 2 | 1 | 0 | `not_buildable`: no outside fungal homolog |
| KS0086 | ecm33 | Q5AGC4 | Saccharomycotina | leave-clade-out | phmmer | 4 | 8 | 1 | 0 | |
| | | | | | jackhmmer | 4 | 8 | 1 | 37 | `not_buildable`: only 4 outside sequences (minimum 5) |

- **hyd1_like.** The only outside hit is `A0A6V8R0V1` (Sordariomycetes, taxid 101201, E 1.3e-21, score 81.7) and it appears only after iteration. The one-pass set has no outside hit. The one inside hit is `D4APE0` (Eurotiomycetes).
- **aga1.** Swiss-Prot holds the query itself and nothing else that matches.
- **ece1.** The 2 inside hits are in Saccharomycotina (`B9WG30` E 1.2e-113 and `C5M337` E 1.1e-31, one-pass). Both belong to the held-out clade, so they are not outside homologs.
- **ecm33.** The 4 outside hits are the same four sequences with both methods: `D4AZK9` (taxid 663331) and `Q4WNS8` (taxid 330879) in Eurotiomycetes, and `O13960` and `O74411` (both taxid 284812) in Taphrinomycotina. That is 2 clades and 3 species, 4 sequences, one short of the minimum of 5. The 37 non-fungal hits of `jackhmmer` are counted apart and are not part of the rule.

Result: all four targets are `not_buildable` in this round; `buildable` 0 of 4.

## 3. The Ecm33 check

`PF12454` (`Ecm33`, GA 23.2 for the sequence and for the domain) was searched with `hmmsearch --max -E 10` against both Ecm33 proteins:

| Sequence | Length | Sequence score | E-value | Best domain score | Model coverage | Verdict |
|---|---|---|---|---|---|---|
| Q5AGC4 (*C. albicans*) | 413 aa | 5.5 | 0.0022 | 4.0 | 0.610 | `below_GA` |
| P38248 (*S. cerevisiae*) | 429 aa | 0.6 | 0.069 | 3.4 | 0.756 | `below_GA` |

So the earlier "no GA hit, not explained" for Q5AGC4 (note `2026-10-22-pr6-pfam-scan.md`) is a real below-GA result: neither Ecm33 protein reaches the GA threshold of the shipped model,: Q5AGC4 has a sequence score of 5.5 against a GA of 23.2, and the *S. cerevisiae* protein has 0.6. The panel row for `PF12454` therefore cannot detect either protein it is named for. This agrees with what is already recorded: `analysis/wall_fusion_report/wall_fusion_models.tsv` (row 17) and report 07 (row 22, and the table row for ECM33 YBR078W) say that *S. cerevisiae* Ecm33 has no hit to the model, "so a hit does not mean Ecm33", with sensitivity for its namesake 0 of 1. This record adds the second species (Q5AGC4).

Not explained: why the model fails on its namesake. The files do not say how many seed sequences `PF12454` has or from which taxa; I found no committed file that records it, so I state nothing about it. One detail to look at: report 07 lists the length of P38248 as 406 aa, while the sequence in the committed shared negatives (the one searched here) is 429 aa. The files do not explain the difference.

## 4. What this round does and does not show

- Swiss-Prot holds reviewed entries only, a small fraction of the fungal proteins in UniProt. "No outside fungal homolog" for aga1 and ece1 is therefore a weaker finding than it would be against a full database. Small, fast-evolving, species-specific proteins are the expected case, but this round cannot tell a true absence from an absence in the reviewed set.
- ecm33 is one sequence short of the minimum; the number is sensitive to what the database holds. hyd1_like is four short.
- The owner decided on a second round on UniRef50 2025_03 (`/srv/projects/db/Uniprot/2025_03/uniref50.fasta`, 28.8 GB) for queries that find nothing in Swiss-Prot. Which queries go into it is still to be chosen: read strictly, only aga1 and ece1 found nothing (0 outside hits); hyd1_like (1) and ecm33 (4) are below the minimum but did find something. The script takes one database and all four queries, so a subset needs either a small change or a re-run of all four.
- No model was built. No call, status, tier, `narrow_ok` or panel row changed. The training rule is provisional (decision 7), and the counts are sequences, not independent clusters (a proxy).
