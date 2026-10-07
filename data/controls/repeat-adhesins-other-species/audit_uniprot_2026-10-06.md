# Audit of agent-C's result against UniProt (2026-10-06)

*Script: `analysis/calibration_truth/uniprot_nonyeast_repeat_candidates.py`. Output:
`uniprot_nonyeast_repeat_candidates.tsv`. UniProt is queried live, so numbers can change.*

**Question.** Agent-C (task 08) found no new class-2a repeat adhesin outside Saccharomycotina. Did it
miss proteins?

**Method.** All UniProt fungal entries outside Saccharomycotina (taxon 147537) with at least one
`Repeat` feature and a cell-adhesion or cell-wall annotation (keyword KW-0130 or KW-0134, GO:0007155
and children, or "adhesin", "agglutinin", "flocculin", "hydrophobin" in the protein name). Query and
fetch date are in the script and the table.

**Result.**

| Item | Count |
|---|---|
| Entries returned | 61 (20 reviewed) |
| With 2 or more repeat features | 49 |
| With an adhesion keyword or GO term | 36 |
| With an adhesion GO term that has an experimental evidence code | 7 |
| With a signal peptide / a GPI feature | 22 / 5 |
| Already in the curation tables | 11 |

Most of the 49 are annotation noise: automatic (IEA) cell-adhesion terms on phospholipases, kinases,
ankyrin-repeat and FG-GAP-repeat proteins. 12 are *Schizosaccharomyces pombe* (a yeast, not
Saccharomycotina); 8 of them are already in the tables.

**Leads that are not noise** (all candidates, none is a control):

| UniProt | Protein | Organism | Repeat features | Signal peptide / GPI | New cluster at 30% identity? | Literature check |
|---|---|---|---|---|---|---|
| G1XAA9 | Adhesin-like cell surface protein MAD1 | *Arthrobotrys oligospora* | 9 | yes / yes | **No.** Same cluster as *Metarhizium* MAD1 (Q2LC49) | PubMed (PMID 26941065, DOI 10.1038/srep22609) says AoMad1 "has been functionally characterized" in adhesion production; the primary paper was not identified |
| Q4WW98, B0Y8Y8 | Subtelomeric hrmA-associated cluster protein | *Aspergillus fumigatus* | 29, 23 | no / no | **Yes**, one cluster for the two | No PubMed hit for the query. Keyword "cell adhesion" only |
| Q9P403 | Intracellular hyphae protein 1 | *Colletotrichum lindemuthianum* | 10 | yes / no | **Yes** | No PubMed hit for the query. Cell-wall keyword only |

Also present and already known: *Metarhizium* MAD1 (8 repeat features) and MAD2 (3), both with a
signal peptide and a GPI feature.

**Reading.**
- Agent-C's result for controls holds. No non-yeast protein with a repeat feature has protein-level
  adhesion evidence in the tables beyond what is already there, and the one new adhesin-annotated
  protein (*Arthrobotrys* MAD1) joins an existing cluster.
- Agent-C's candidate list (2 items) was too short. Nematode-trapping fungi were not in its genus
  list. The three leads above are new.
- Best case if the literature supports both new clusters: 2 more clusters (19 to 21). Neither has a
  protein-level paper yet.
- The *Metarhizium* MAD1 repeat features (and *Arthrobotrys* MAD1) contradict agent-A's label `2b-i`
  (CFEM) for MAD1.

**Not done.** The leads were not read in full text. The primary AoMad1 paper needs identifying. The
query does not find repeat adhesins that lack a UniProt repeat feature (for example SOWgp and BAD1,
which have none recorded).
