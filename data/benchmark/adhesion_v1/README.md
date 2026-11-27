# adhesion_v1 — frozen benchmark for fungal adhesion prediction (built 2026-11-27)

Built by `analysis/curation/build_benchmark.py` (owner decisions D3 and D4 of 2026-11-26; plan
`docs/superpowers/plans/2026-11-26-adhesion-round2-benchmark-cruft-and-layout-study.md`).

**Status: NOT YET FROZEN.** The fold list below is provisional because it counts proteins, and the
real unit is a 30%-identity MMseqs2 cluster. Cluster counts are computed on the HPCC, which was
unreachable when this was built. `VERSION` and `SHA256SUMS` are written so the inputs are pinned, but
the set should not be used to decide anything until the fold list is recomputed on clusters and the
owner signs it off.

## Why this exists

Every model the project has trained reduces to amino-acid composition plus length
(`docs/model-review/`, `docs/reports/data/aga1_feature_classifier_2026-10-25/ablation.tsv`). This set
exists so that a component has to beat that baseline against proteins it cannot separate on
composition, with whole clades held out.

## Files

| File | Content |
|---|---|
| `proteins.tsv` | every protein considered, one row, with its tier and (if excluded) the reason |
| `counts.tsv` | positives and hard negatives per class |
| `folds.tsv` | what each clade is scored as, and which clades cannot be scored |
| `VERSION` | the version name, input checksums and tier counts |
| `SHA256SUMS` | checksums of the files above |
| `sources/*.tsv` | the literature-curation passes, one per clade group |
| `sources/*.unresolved.md` | candidates that failed verification, and why |
| `sources/owner_excluded.md` | rows removed by owner decision, kept in full so it is reversible |

## Tiers

| Tier | n | Meaning |
|---|---|---|
| `positive` | 76 | direct experimental adhesion evidence AND a secretory architecture to predict from |
| `hard_negative` | 90 | a surface or wall protein with a demonstrated non-adhesion function |
| `positive_moonlighting` | 8 | direct evidence, but a cytoplasmic protein found on the surface |
| `excluded` | 143 | failed a rule; `exclude_reason` says which |

Direct evidence means at least one of: `ko_assay` (deletion changes adhesion in a quantitative
assay), `het_expr` (expression in a non-adherent host confers adhesion), `binding` (purified protein
or domain binds a ligand or cell, measured directly). GO annotations, expression induction, homology
and virulence attenuation without an adhesion assay are **not** sufficient.

**`indirect=yes`** marks a hard negative whose deletion *does* reduce adhesion, but downstream of its
enzymatic or structural function — *Cryptococcus* Plb1 (restored by exogenous palmitic acid),
*Botrytis* BcSUN1, and the Ecm33/glucanase class from the curated table. These are the trap a
predictor must not fall into, so they are negatives rather than omissions.

**Moonlighting** proteins are kept on record but stay out of both contrasts and out of the layout
study: no architecture method can predict a protein with no signal peptide and no anchor, so scoring
them would only add a fixed error to every component alike.

## Contrasts

- **C1, positives vs hard negatives** — decides whether a component is kept. Composition-matched, so
  it cannot be won on composition alone.
- **C2, positives vs all other secreted or GPI proteins of the same proteomes** — descriptive;
  assembled at scoring time, not here.

## Folds (provisional — protein counts, not clusters)

A class is its own fold when it has at least 5 positives and 5 hard negatives; otherwise it is pooled
into its subphylum, and a group that is still too thin is reported rather than quietly merged.

| Fold | Level | Positives | Hard negatives | Scored |
|---|---|---|---|---|
| Saccharomycetes | class | 32 | 57 | yes |
| Pezizomycotina | group | 13 | 7 | yes |
| Eurotiomycetes | class | 5 | 11 | yes |
| Taphrinomycotina | group | 11 | 0 | no |
| Microsporidia | group | 9 | 3 | no |
| Mucoromycota | group | 3 | 2 | no |
| Agaricomycotina | group | 2 | 6 | no |
| Ustilaginomycotina | group | 1 | 3 | no |
| Chytridiomycota | group | 0 | 1 | no |

Three folds can be scored. That is the honest state of the field, and it is the main limit on any
claim this project makes about cross-clade transfer.

## What the curation found, and did not

Three literature passes (Basidiomycota; early-diverging fungi; filamentous Ascomycota outside the
yeasts and *A. fumigatus*) added 54 verified rows. Every PMID was checked against PubMed and every
accession against UniProt or NCBI.

Verified absences matter as much as the additions:

- **Basidiomycota is close to empty.** Outside *Schizophyllum* hydrophobins and four already-curated
  proteins, no Basidiomycota adhesin has direct evidence. *Rhodotorula*, *Sporobolomyces* and
  *Malassezia* have none at all — the *M. sympodialis* genome paper states the genus lacks the wall
  proteins that mediate adhesion elsewhere. Any Basidiomycota prediction is therefore essentially
  unvalidatable against current literature.
- **Dothideomycetes adhesion was characterised and never sequenced.** *Parastagonospora nodorum* and
  *Venturia inaequalis* both have well-measured, protein-dependent adhesion with no protein
  identified. *Alternaria* and *Pyrenophora* are genuinely empty (all hits screened by title).
- **Rust fungi** have the strongest functional evidence of any group — adding purified enzymes back
  to autoclaved spores restores adhesion to the bean cuticle — but no sequence was ever deposited.
- **Colletotrichum** has no adhesion gene with a deletion mutant.
- **Microsporidia** supplies 9 of the 76 positives, but it is effectively one heparin-binding-motif
  panel: effect sizes are 10–20% adherence reduction and there is not one knockout among them.

Two pairs are deliberately kept because they test the same architecture on both sides: *U. maydis*
Hum3 (17 Rep1-like repeats, adhesion unchanged on deletion) against the Rep1 positive, and
*S. commune* SC3 (positive) against SC4 (negative). *M. oryzae* Cut2 and *Z. tritici* ZtGT2 are the
cleanest negatives in the set: both had adhesion measured and unchanged.

## Known weaknesses

1. Three scorable folds, and Saccharomycetes holds 42% of the positives and 63% of the negatives.
2. Protein counts, not cluster counts. The Saccharomycetes positives are largely FLO/ALS families and
   will collapse substantially on clustering.
3. Hard negatives come overwhelmingly from *S. cerevisiae* and *C. albicans*, so C1 outside
   Saccharomycotina rests on very few proteins.
4. Several rows are one experimental unit rather than independent observations (both EnP1
   orthologues; SWP9 and SWP7 from one paper; the two *Beauveria* hydrophobin papers disagree about
   which paralogue carries the defect, and both are included with the conflict recorded).
5. `sources/*.unresolved.md` lists what could not be verified. The largest recoverable item is the
   closed-access ricin-B-lectin paper (PMID 34808329); nothing was paid for and no paywall bypassed.

## Rebuilding

```bash
python3 -I analysis/curation/build_benchmark.py --version adhesion_v1 --out data/benchmark/adhesion_v1 --force
```

A change to the rules or the sources is a **new version directory**, not an edit to this one.
