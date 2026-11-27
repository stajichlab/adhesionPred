# adhesion_v1 — frozen benchmark for fungal adhesion prediction (built 2026-11-27)

Built by `analysis/curation/build_benchmark.py` (owner decisions D3 and D4 of 2026-11-26; plan
`docs/superpowers/plans/2026-11-26-adhesion-round2-benchmark-cruft-and-layout-study.md`).

**Status: FROZEN 2026-11-27 by Jason Stajich.** Membership, tiers and folds are fixed. A change to
the rules, the sources or the clustering is a new version directory (`adhesion_v2`), never an edit
to this one.

The fold list below is computed on **clusters**, not proteins: MMseqs2 17-b804f,
`--min-seq-id 0.3 -c 0.5 --cov-mode 0`, 166 sequences to 113 clusters (`clustering.json`,
`clusters.tsv`). Every scored protein has a sequence; none was left unclustered. `VERSION` and
`SHA256SUMS` pin the inputs.

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
| `sequences.faa` | one sequence per scored protein, as fetched from UniProt/NCBI |
| `clusters.tsv` | accession to cluster representative (MMseqs2, 30% identity) |
| `clustering.json` | clustering parameters, counts and MMseqs2 version |

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

## Folds (cluster counts, MMseqs2 30% identity)

A class is its own fold when it has at least 5 positive and 5 hard-negative **clusters**; otherwise
it is pooled into its subphylum, and a group that is still too thin is reported rather than quietly
merged. Counts are clusters, so a family of paralogues counts once.

| Fold | Level | Positives | Hard negatives | Scored | Proteins before clustering |
|---|---|---|---|---|---|
| Saccharomycetes | class | 17 | 38 | yes | 32 / 57 |
| Pezizomycotina | group | 11 | 7 | yes | 13 / 7 |
| Eurotiomycetes | class | 5 | 11 | yes | 5 / 11 |
| Microsporidia | group | 7 | 3 | no | 9 / 3 |
| Taphrinomycotina | group | 7 | 0 | no | 11 / 0 |
| Agaricomycotina | group | 2 | 5 | no | 2 / 6 |
| Mucoromycota | group | 1 | 2 | no | 3 / 2 |
| Ustilaginomycotina | group | 1 | 3 | no | 1 / 3 |
| Chytridiomycota | group | 0 | 1 | no | 0 / 1 |

Three folds can be scored. That is the honest state of the field, and it is the main limit on any
claim this project makes about cross-clade transfer.

**What clustering changed.** Saccharomycetes lost 47% of its positives (32 to 17) and a third of its
negatives: those were FLO/ALS/EPA paralogues, not independent observations. Mucoromycota collapsed
from 3 positives to 1 — the *Rhizopus* CotH proteins are one family. Eurotiomycetes did not move at
all, so its small set is genuinely non-redundant. The three scorable folds survived, which was not
guaranteed.

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

1. Three scorable folds, and Saccharomycetes still holds 47% of the positive clusters after
   clustering.
2. 113 clusters in total. Any claim about cross-clade transfer rests on three folds, one of which
   (Eurotiomycetes) has 5 positive clusters — the minimum the rule allows.
3. Hard negatives come overwhelmingly from *S. cerevisiae* and *C. albicans*, so C1 outside
   Saccharomycotina rests on very few proteins.
4. Several rows are one experimental unit rather than independent observations (both EnP1
   orthologues; SWP9 and SWP7 from one paper; the two *Beauveria* hydrophobin papers disagree about
   which paralogue carries the defect, and both are included with the conflict recorded).
5. `sources/*.unresolved.md` lists what could not be verified. The largest recoverable item is the
   closed-access ricin-B-lectin paper (PMID 34808329); nothing was paid for and no paywall bypassed.

## Rebuilding

```bash
# 1. membership and tiers (laptop or HPCC)
python3 -I analysis/curation/build_benchmark.py --version adhesion_v1 --out data/benchmark/adhesion_v1 --force
# 2. sequences and cluster-based folds (HPCC: needs the network, then MMseqs2)
python3 -I analysis/curation/cluster_benchmark.py fetch   --benchmark data/benchmark/adhesion_v1
module load MMseqs2/17-b804f
python3 -I analysis/curation/cluster_benchmark.py cluster --benchmark data/benchmark/adhesion_v1 --tmpdir "$SCRATCH/mmseqs_benchmark"
```

Step 1 rewrites `folds.tsv` with protein counts, so step 2 must be re-run after it.

A change to the rules or the sources is a **new version directory**, not an edit to this one.
