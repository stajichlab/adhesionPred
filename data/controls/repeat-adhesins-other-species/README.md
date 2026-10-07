# Task 08: repeat-mechanism controls in fungi outside Saccharomycotina

*Curation by agent-C, 2026-10-06. Task file: `docs/agent-tasks/08-repeat-adhesins-other-species.md`.
Common rules: `docs/agent-tasks/COMMON-RULES.md`.*

## What this directory holds

| File | Content |
|---|---|
| `curation_table.other_species.tsv` | 8 new adhesin rows in the curation-table format (with `evidence_basis`). |
| `hard_negatives_other_species.tsv` | 10 new look-alike negatives (task-07 columns). |
| `search_log.tsv` | Every PubMed query run (query, database, date, hits, PMIDs kept), including zero-hit queries. |
| `candidates.tsv` | Repeat-adhesin-lookalike proteins with no adhesion paper (leads for the owner). |
| `rejected.tsv` | Items looked at and not added, with the reason. |
| `sequences.faa` | Protein sequences of every new row (UniProt, or NCBI where no UniProt entry). Header = accession. |
| `clusters.tsv` | `accession<TAB>cluster_id` (MMseqs2, 30% identity, 50% coverage). |
| `counts.md` | Rows, sequences and clusters per group and per class, and the size targets. |
| `unverified_notes.md` | Claims we could not check. |

## What was done

1. Searched PubMed (NCBI E-utilities `esearch`) for each group with varied words: adhesin, adhesion,
   adherence, agglutinin, flocculation, hydrophobin, cell wall protein, surface protein, biofilm,
   appressorium, attachment. Every query is logged in `search_log.tsv`. Candidate PMIDs were
   fetched (`efetch`) and cached in `_workdir/abstracts/`.
2. Qualified each protein against the task rule: an adhesin row needs a paper showing adhesion or
   binding for that protein (deletion, heterologous expression, antibody blocking, or purified
   protein). `E1` direct, `E2` weaker/family-level; negatives `N1` (characterized non-adhesive) or
   `N2` (adhesin-like, no adhesion evidence).
3. Labelled the mechanism with `docs/TOOL-ARCHITECTURE.md` section 2 class codes (`2a`, `2b-i`,
   `2b-ii`, `2b-iii`, `2c`, `2d`, `other`, `unknown`) plus `label_confidence` and `evidence_basis`,
   each with a `mechanism_quote`. No repeat detector was used; repeats were never inferred from a
   database annotation or from an ortholog in another species.
4. Resolved sequences. UniProt accession where one exists; otherwise the NCBI protein accession
   (`XP_045276491.1` Eng2, `ACE80261` Lip1, `CCD56554` rejected). A protein without a resolvable
   sequence (the Colletotrichum 110-kDa glycoprotein) went to `rejected.tsv`.
5. Clustered the new rows with the existing 186 rows
   (`data/controls/repeat-mechanism/curation_table.agent-A.tsv`) using MMseqs2 17-b804f
   (`easy-cluster`, `--min-seq-id 0.3 -c 0.5`). Report in `counts.md`.

## Commands (reproducibility)

```
# PubMed searches (all logged in search_log.tsv)
python3.12 _workdir/search.py _workdir/queries.json
python3.12 _workdir/fetch_abstracts.py _workdir/queries.json
# sequences
python3.12 _workdir/fetch_build.py       # writes sequences.faa + all_rows.faa
# clustering
MM=/opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs
$MM easy-cluster all_rows.faa _workdir/clum/clu _workdir/clum/tmp --min-seq-id 0.3 -c 0.5 -v 1 --threads 4
python3.12 _workdir/analyze_clusters.py  # cluster per group / class counts
```

Databases accessed 2026-10-06: PubMed (NCBI E-utilities), UniProt (rest.uniprot.org), NCBI protein
(efetch). Taxid values from `data`/NCBI names (checked against the taxonomy; e.g. `Ustilago maydis`
used as the scientific name `Mycosarcoma maydis`, taxid 5270).

## Key findings (point the owner at `counts.md`)

- **No new `2a` repeat adhesins.** All 8 new positives are `2c`, `2d` or `other`. The repeat call
  was not widened outside Saccharomycotina.
- New negatives across 7 clades (hydrophobins that are non-adhesive in their species, wall enzymes,
  mucin sensors, a secreted lipase, a glycosyltransferase).
- 4 negatives share a homology cluster with an existing row (flagged); 2 of these conflict with a
  same-cluster positive.
- Several named groups returned no protein-level evidence (Rhodotorula, C. gattii, Exophiala,
  Wallemia, Hortaea, Knufia, Colletotrichum-as-a-positive, Blumeria-positive).
