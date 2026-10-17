# Review 2: lectin validation (branch `lectin-validation`, tip d9ecb86)

*2026-10-11. Reviewer: Claude Code (claude-opus-5-5), independent review. Narrow scope: check that the fixes to review 1 are correct and that they added no new error. Read: `docs/reports/2026-10-11-lectin-validation.md` (revision 2), `analysis/lectin_validation/` (code, `results/`, `controls/`, `LOG.md`, `unverified_notes.md`, rejected tables), `tests/lectin_validation/`. Governing documents: lectin spec (owner decisions: option (c), report only, no call, no status, no threshold) and `docs/agent-tasks/COMMON-RULES.md`. I did not edit validation files, the spec or the scripts. All runs were on the login node. No SLURM job was used. My scratch scripts are in `/tmp/lvrev2/` and are not committed.*

## Verdict

**Pass with changes.**

The review 1 fixes are correct in the data and mostly correct in the text. All headline numbers reproduce. The pipeline re-runs to the same outputs. The new tests catch the M1 error. One new text error must be fixed: the "seed-independent" result of 2 of 13 clusters counts G1X5Y2 as a hit, but G1X5Y2 entered the set through the PF07938 cross-reference. Its hit is by construction. The fully independent figure is still **1 of 12 clusters** (CNL), the same as in review 1. The report does not say this. Section 6 gives the plain statement on use.

## 1. What I re-derived

I wrote my own script (`/tmp/lvrev2/rederive.py`). It reads only `controls/controls.tsv` and `out/scan.hits.tsv`. It does not import `lectin_metrics` or `measure.py`.

| Claim in the report | My result | Agree |
|---|---|---|
| 59 primary positives, 32 clusters, 38 species | 59, 32, 38 | yes |
| Union of 15: 42/59 (0.59-0.81), 21/32 (0.48-0.80) | same | yes |
| Independent route: 32/49 (0.51-0.77), 17/28 (0.42-0.76) | same | yes (see N2) |
| Model-gated route: 10/10, 5/5 | same | yes (see N2) |
| T1 rows: 17/18, 9/10 | same | yes |
| PF07938 only: 6/59 (0.05-0.20), 4/32 (0.05-0.28) | same | yes |
| Clusters with an `in_reference` member: 40/40 sequences, 19/19 clusters | same | yes |
| Clusters with no `in_reference` member: 2/19 sequences, 2/13 clusters (0.04-0.42) | same; hits are B2ZRS9 (CNL) and G1X5Y2 | numbers yes; text no (N1) |
| PF14200 family: 6/12 sequences, 5/8 clusters (0.31-0.86); 6/6 and 5/5 without the trefoil rows | same; trefoil rows in 3 clusters (R4GRU5 x4, B0D650, Q7Z8V1) | yes |
| PF07938 scores 525.1, 489.1, 300.7, 158.6, 132.6, 60.7 | Q2UNX8, Q4WW81, A0A084FYP2, P18891, G1XA82, G1X5Y2: same | yes |
| Floors 50 / 100 / 200 on PF07938 family: 6/6, 5/6, 3/6; union 42, 41, 39 of 59; clusters 21, 20, 19 of 32 | same | yes |
| GA values: PF07938 25.3, PF14200 25.6 | `out/scan.models.hmm`: same | yes |
| Stratum clusters 23 + 3 + 5 + 4 + 1 = 36, 4 shared clusters | same (F2ST24, P18891, Q2UNX8, Q9P4R8) | yes |
| Negatives: N1 0/136 (0/74), ER 0/18 (0/8), enzyme 8/12 (1/2), all reasoned 8/166 (1/84), 54 species | same | yes |
| `assumed` 0/497 (0/430), 148 species | same | yes |
| Report-script check: 0 mismatches, 18 of 150 cells non-empty | same (`results/report_script_check.tsv`) | yes |
| Five smallest beta-trefoil margins: 0.7 to 4.4 bits below GA | 25.6 − 24.9 = 0.7; 25.6 − 21.2 = 4.4 | yes |

**Re-run.** I copied the tree to `/tmp/lvrev2/tree`, added the two Swiss-Prot raw files from the author's worktree, and ran `run_all.sh`. Exit 0. The report-script check printed `mismatches 0`. `diff -rq` against the commit: only three files differ (`out/cluster/clu_cluster.tsv`, `out/cluster/seed_hits.m8`, `out/scan_swissprot.tbl`). After sorting and removing comment lines they are identical. So the difference is line order and time stamps only. Every file in `results/` and `controls/` is identical.

**Tests.** `/usr/bin/python3.12 -m pytest -q tests/lectin_validation`: **34 passed**, none skipped (the taxonomy test ran).

**Mutation checks** (on the `/tmp` copy):

| Mutation | Tests that fail |
|---|---|
| `lectin_metrics.cluster_seed_split` made row-level | 3: `test_cluster_level_leakage_split_does_not_double_count`, `test_family_recall_seed_independent_is_cluster_level`, `test_cluster_seed_split_unit` |
| `measure.py` leakage split (line 388) bypasses `cluster_seed_split` and splits rows | 1: `test_cluster_level_leakage_split_does_not_double_count` |

The new tests run `measure.py` as a subprocess on a toy table. They do exercise the grouped numbers. The test of the route split (`test_independent_route_group`) checks sequence counts only, not clusters (see N2).

## 2. Review 1 items: status

| Item | Status | Evidence |
|---|---|---|
| M1 seed-independent recall | **Fixed in data and code. One new text error (N1).** | `cluster_seed_split` is correct. `leakage_split.tsv` labels row-level and cluster-level rows. `family_recall.tsv` splits per model at the cluster level. Cluster-level groups add to 32 (19 + 13). No double count in `leakage_split.tsv` or `family_recall.tsv`. |
| M2 selection route | **Fixed.** Route values correct for all 59 rows. One cluster is shared by the two route groups (N2). | `data/candidate_screen.tsv`: the 10 `model_gated_pfam_xref` rows have no Lectin keyword, are not in UniLectin3D, and carry a PF01476 or PF07938 UniProt cross-reference. B2ZRS9 and Q06100 have the Lectin keyword. G4N906 (Slp1) and Q00363 (Avr4) are named in the spec (L5 and M1 rows) and are in the `EXTRA` list of `make_acc_list.py`. C1G842 is the spec's paracoccin. G3JGU5 came from the review 1 PubMed search. The 43 UniLectin3D rows come from structures. No `independent` row is gated by a model cross-reference. The report says plainly which numbers are partly by construction (summary item 3, sections 5.1 and 8). |
| M3 `assumed` | **Fixed.** | `results/false_hits.tsv` has no `assumed` group (test asserts it). `construction_check.tsv` holds the 0/497. Appendix A says "by construction and has no power". |
| M4 beta-trefoil mapping | **Fixed.** | Both PF14200 figures and the cluster count are given; the mapping choice is stated. |
| M5 one T2 rule | **Fixed and consistent.** | See below. |
| M5 D4AL26 evidence | **Fixed.** | Section 5.2 says protein-level (secretome) evidence and no binding evidence. |
| M6 UniParc sentence | **Fixed.** | Removed; UPI000150AC21 is listed as rejected. |
| m1 to m11 | **Fixed**, except two small new points (N5, N6). | m3: N1-only 16 and 15. m5: "28.8, above its GA of 25.6". m7: "18 of the 150 cells". m8: `run_all.sh` now makes `scan_swissprot.hits.tsv`. m10: exception recorded in section 9. |

**T2 rule applied independently.** `data/candidate_screen.tsv` has 159 candidates. 35 are Swiss-Prot reviewed with an ECO:0000269 sugar- or chitin-binding statement. All 35 are either primary positives, PA14 positives (P32768), or rejected with a valid reason (five duplicates of another strain or accession: A0A1X7S147, A0A1X7S1V5, A0A1X7S690, Q8NJT4, Q8TGE0; two fragments: P86922 at 8 aa, P84526 at 13 aa). I found no inconsistency. Note: UniProt gives G1X5Y2 protein existence "inferred from homology" (PE 3), even with the ECO:0000269 function line. This supports the report's note that the abstract shows no binding assay.

**New rows.**

| Row | Check | Result |
|---|---|---|
| CMRBL G3JGU5 | UniProt ORF name | CCM_03832, *Cordyceps militaris* CM01. UniProt cites only the genome paper (PMID 22112802), not 36632960. The link to the paper is by gene ID only, as `unverified_notes.md` item 20 says. |
| | Abstract of PMID 36632960 | "The Cordyceps militaris gene CCM_03832 encodes a ricin-B like lectin ... CMRBL possesses very high specific binding activity of glycans terminated with type II LacNAc". Quote in the table is correct. T1 is correct. |
| | Pfam | UniProt xref PF00652. Scan: PF00652 40.2 (GA 29.0) and PF14200 43.6. |
| | Taxon | 73501, species rank (strain 983644). Correct. |
| Mgx1LysM P9WEM2 | Abstract of PMID 33797163 | "Mgx1LysM binds chitin, suppresses the chitin-induced ROS burst, and is able to protect fungal hyphae against chitinase hydrolysis". Quote correct. T1 correct. |
| | Taxon | 1047171 (*Zymoseptoria tritici*, species). Correct. Same cluster as Mg1LysM and Mg3LysM (F9XHX3). |
| G1X5Y2 | Abstract of PMID 35142828 | Deletion phenotype only (trap formation, nematocidal activity). No binding assay. The report says this. |
| | Taxon | 2813651 (*Orbilia oligospora*, species). Correct. |

## 3. Findings that must be fixed

### N1. The "seed-independent" hit G1X5Y2 is a model-gated row. The report does not say so.

The 13 clusters with no `in_reference` member are hit 2 of 13: CNL (B2ZRS9) and G1X5Y2. G1X5Y2 has `selection_route = model_gated_pfam_xref`. It entered the candidate list only because UniProt gives it a PF07938 cross-reference. That cross-reference is an InterPro match of PF07938 at GA. So its PF07938 hit here is by construction.

Without G1X5Y2, the seed-independent result is **1 of 12 clusters** (CNL; Wilson 0.01 to 0.35). That is the same as the corrected review 1 figure. The text frames 2 of 13 as the current result in five places and gives no caveat:

- Summary item 4 ("hit 2 of 13 ... CNL (B2ZRS9) and G1X5Y2").
- Section 5.1, "Seed-independent test" ("PF07938 has 1 cluster with no seed-similar row (G1X5Y2, hit)").
- Section 5.3 ("G1X5Y2 has no seed hit, but one positive does not set a cutoff either").
- Section 5.4 ("The two hits are CNL and G1X5Y2").
- Section 8 ("Seed-independent clusters: 2 of 13").

Fix: in each place, say that G1X5Y2 was selected through PF07938, so its hit is by construction. Give "1 of 12 clusters (0.01 to 0.35) when the model-gated row is removed". In section 5.3, say that the only PF07938 positive without a seed hit is also model-gated, so no PF07938 positive is both seed-independent and found without the model. No row change is needed.

## 4. Minor findings

**N2. The two route groups share one cluster.** Cluster F9XHX3 holds independent-route rows (Ecp6 B3VBK9, Mg1LysM F9XHX3, Slp1 G4N906) and model-gated rows (Mg3LysM, Mgx1LysM, ELP1, ELP2). So 28 + 5 = 33 clusters, not 32. This is a row-level split, the same kind of split that M1 fixed for leakage. The independent number (17/28) is still valid as a group statistic, because the cluster is hit by its own independent-route rows. But the model-gated "5 clusters" is not a disjoint count. Only 4 clusters hold model-gated rows alone. Fix: add one sentence in summary item 3 and section 3.1, and label the route rows in `recall_by_group.tsv` as row-level (as `leakage_split.tsv` does). The per-model route rows show the overlap too: PF01476 hits 2 of 28 independent and 4 of 5 gated clusters, but only 5 clusters in all.

**N3. The count "52 non-lectin names" has no script.** I re-derived it from `data/raw/uniprot_fungi_pf07938_xref.tsv.gz`. The 949 entries, 7 reviewed entries (all named "Fucose-specific lectin") and 69 multi-domain entries reproduce. With the rule "name has no 'lectin', 'fucose' or 'uncharacterized'", I get **53** non-lectin names, all TrEMBL (51 "inferred from homology", 2 transcript level). The report says 52 (50 and 2). The difference depends on the name rule. Some of these names are placeholders, not functions (two UUID-like names, "Unplaced genomic scaffold", "hypothetical protein"). No script or LOG line records the fetch query or the rule. Fix: give the rule and the query in `LOG.md`, or say "about 50". The conclusion does not change.

**N4. Section 6 claim is justified, with one nuance.** "No PF07938 specificity claim" is correct and is not overstated. One nuance: "every reviewed protein with a PF07938 cross-reference is a fucose-specific lectin (7 of 7)" is a statement about names. Of the 7, three have no own binding evidence (D4AL26 by similarity; Q8NJT4 and Q8TGE0 are duplicate accessions of AFL and AOL), and G1X5Y2 has no binding assay in its abstract. Suggest "7 of 7 are named fucose-specific lectins; 3 of the 7 have binding data of their own (AAL, AOL, AFL)".

**N5. "0.25% of proteins" mixes units.** Section 5.2 and appendix A say "97 protein-model hits (0.25% of proteins)". 97 is protein-model pairs. The 15 models hit **90 distinct proteins** of 38,451 (0.23%). The argument in appendix A does not change.

**N6. Species total.** Section 7 says "38 of 188 species have one" positive. 188 includes the 148 species of the `assumed` set. Without `assumed`, the set has 83 species. 42 species have a positive if the PA14 rows count. Say which set the 188 refers to.

**N7. G1X5Y2 PE level.** Add to the row note that UniProt gives protein existence "inferred from homology". It is a weak T2.

## 5. Fungal lectins with experimental evidence not in the set

PubMed searches on 2026-10-11 (only abstracts read; I did not add rows):
1. `lectin[Title] AND (fungus OR fungal OR Aspergillus OR ... OR mushroom) AND (recombinant OR glycan array OR crystal structure) AND binding`, 2019 to now (332 hits; 20 read by title and abstract).
2. `"fungal lectin"[Title] OR "mushroom lectin"[Title] OR (lectin[Title] AND <genus list>[Title])`, 2020 to now (223 hits; 20 read).

| Candidate | Species | Evidence in the abstract | PMID | Accession | Note |
|---|---|---|---|---|---|
| PPL | *Pleurocybella porrigens* (Basidiomycota) | "We have cloned PPL and solved its structure by X-ray crystallography and cryo-EM. We report the sugar binding properties of this β-trefoil lectin" | 41342409 | B6ZK51 (UniProt cites this PMID and 19270381; UniProt Pfam PF14200) | High priority. A new beta-trefoil lectin with a structure. It tests PF14200 directly on the open question of section 5.5 / decision 10. Not checked: the PDB entry and whether UniLectin3D lists it. |
| FFSL | *Trichoderma harzianum* (Ascomycota) | "Molecular docking and isothermal titration calorimetry confirmed FFSL's high-affinity ... binding to fucose" | 41534670 | Not verified. A0A2T3ZYI0 is the only PF07938 protein of *T. harzianum* CBS 226.95 in UniProt (291 aa). The paper's strain and gene ID are not checked. | If it is a PF07938 protein, it is an Ascomycota fucose lectin outside the seeds. Medium priority. |
| *Histoplasma* β-galactose lectin | *H. capsulatum* | Haemagglutination inhibited by β-galactose; partly purified | 41805922 | none (no gene) | Cannot be added. Same as the existing note on *Histoplasma*. |
| Rbd proteins | a human fungal pathogen (genus not in the PubMed title field) | Ricin B lectin-like domain proteins; cell wall integrity, virulence. **No binding assay in the abstract.** | 42016327 | not checked | Not a positive on the abstract. |
| Pleurotus galectin-like gene | *Pleurotus ostreatus* | Overexpression gives mite resistance. **No binding assay in the abstract.** | 37180275 | not checked | Not a positive on the abstract. |

The review 1 candidates (BbLec1, CmLec4, CmRlec, RsLysM, Cg2LysM, TAL6, MfAvr4, Cercospora Avr4) are recorded in `controls/rejected.tsv` and `unverified_notes.md` item 22. I did not search those again.

## 6. Owner decisions (report section 10)

The eleven decisions are the right ones and are stated clearly. Each has a default and a reversal. Three additions:

- **Decision 3 (T2 rule)** should add: G1X5Y2 is also model-gated (N1). If it leaves, the seed-independent result is 1 of 12 and PF07938 has no positive without a seed hit.
- **Decision 1 (model-gated rows)** should add that one cluster (F9XHX3) mixes both routes (N2).
- **Decision 6 (new candidates)** could add PPL (B6ZK51) and FFSL (section 5 above).

No decision is missing that would change a call, a status or a threshold. The report creates none of these, as option (c) requires. The validation commits change no file in `src/`, the spec or any status file (checked with `git diff --name-only 0d91255~1 d9ecb86`).

## 7. Can the set be used as a report-grade table?

**Yes, after N1 is fixed in the text.** The rows, sequences, taxa, clusters, leakage labels, selection routes and hit tables are correct and reproduce. N1 changes sentences, not rows. N2 to N7 are small text or provenance fixes.

It is **not** a calibration set. It cannot give a PF07938 specificity. The only test of recall that is both seed-independent and free of model selection is 1 of 12 clusters, one lectin (CNL). The owner review gate of COMMON-RULES still applies.

## Sources

- PubMed abstracts: 36632960 ([doi:10.1016/j.ijbiomac.2023.123207](https://doi.org/10.1016/j.ijbiomac.2023.123207)), 33797163 ([doi:10.1111/mpp.13055](https://doi.org/10.1111/mpp.13055)), 35142828 ([doi:10.1093/femsle/fnac013](https://doi.org/10.1093/femsle/fnac013)), 41342409 ([doi:10.1093/glycob/cwaf082](https://doi.org/10.1093/glycob/cwaf082)), 41534670 ([doi:10.1016/j.biortech.2026.133985](https://doi.org/10.1016/j.biortech.2026.133985)), 41805922 ([doi:10.1371/journal.pntd.0014105](https://doi.org/10.1371/journal.pntd.0014105)), 42016327 ([doi:10.1016/j.tcsw.2026.100173](https://doi.org/10.1016/j.tcsw.2026.100173)), 37180275 ([doi:10.3389/fmicb.2023.1191500](https://doi.org/10.3389/fmicb.2023.1191500)).
- UniProt REST (2026-10-11): G3JGU5, P9WEM2, G1X5Y2, B6ZK51; PF07938 query for *T. harzianum* (taxon 5544).
