# Adhesion Protein Molecular Properties & Functional Domain Survey

## Scope

- 749697 adhesion-predicted proteins across the matched species set.
- 749697 background (non-adhesion) proteins sampled 1:1 per species
  for comparison (seed=42, reproducible).

**Regenerating this report:** `protein_universe.csv` and `protein_domain_flags.csv`
(read by this script) are not tracked in git — they are large
(~163MB/~78MB), fully reproducible per-protein intermediate tables.
Regenerate them first with `01_build_protein_universe.py` and
`02_query_functional_domains.py` (and `03_compute_sequence_properties.py`
/ `04_aggregate_and_test.py` for the other inputs this script and its
upstream tables depend on) before re-running this script.

## Domain/topology presence: adhesion vs. background

```
           feature  adhesion_fraction  background_fraction
          has_pfam           0.401193             0.710767
          has_cazy           0.198748             0.038163
        has_merops           0.005957             0.029559
has_signal_peptide           0.706665             0.062824
      has_tm_helix           0.115930             0.179031
```

## Length and Ser/Thr/Pro content: adhesion vs. background, by phylum

Mann-Whitney U test (adhesion vs. background, within each well-powered
phylum), BH-corrected across phyla:

**Length:**
```
            phylum  n_adhesion  n_background  median_adhesion  median_background       u_stat       p_value    p_value_bh
        Ascomycota      555861        555861            345.0              400.0 1.416686e+11  0.000000e+00  0.000000e+00
     Basidiomycota      135340        135340            310.0              362.0 8.326952e+09  0.000000e+00  0.000000e+00
Blastocladiomycota         938           938            410.0              433.0 4.339380e+05  6.100173e-01  6.100173e-01
   Chytridiomycota       14817         14817            432.0              364.0 1.266646e+08 1.769702e-116 5.309106e-116
      Cryptomycota          40            40            303.0              319.5 8.850000e+02  4.161325e-01  4.681490e-01
     Microsporidia          65            65            471.0              346.0 2.423500e+03  1.482347e-01  1.905874e-01
      Mucoromycota       23786         23786            302.0              322.0 3.005379e+08  4.621646e-32  1.039870e-31
  Sanchytriomycota          28            28            445.5              323.0 5.250000e+02  2.991182e-02  4.486773e-02
     Zoopagomycota       18647         18647            338.0              357.0 1.709071e+08  4.567340e-03  8.221211e-03
```

**Ser+Thr+Pro %:**
```
            phylum  n_adhesion  n_background  median_adhesion  median_background       u_stat       p_value    p_value_bh
        Ascomycota      555861        555861        30.331126          18.723404 2.739534e+11  0.000000e+00  0.000000e+00
     Basidiomycota      135340        135340        26.992825          20.041754 1.421713e+10  0.000000e+00  0.000000e+00
Blastocladiomycota         938           938        27.333661          19.499717 6.892610e+05 3.013677e-100 4.520516e-100
   Chytridiomycota       14817         14817        26.528926          17.902542 1.751937e+08  0.000000e+00  0.000000e+00
      Cryptomycota          40            40        27.520142          17.291516 1.311000e+03  9.001945e-07  1.012719e-06
     Microsporidia          65            65        26.229508          16.167665 3.006500e+03  3.176683e-05  3.176683e-05
      Mucoromycota       23786         23786        23.868313          17.802028 3.770963e+08  0.000000e+00  0.000000e+00
  Sanchytriomycota          28            28        28.011540          16.142568 7.380000e+02  1.499221e-08  1.927570e-08
     Zoopagomycota       18647         18647        28.200371          18.371608 2.847865e+08  0.000000e+00  0.000000e+00
```

![Length by phylum](figures/box_length_by_phylum.png)

![Ser/Thr/Pro composition by phylum](figures/box_pct_ser_thr_pro_by_phylum.png)

## Domain presence by phylum

![Pfam presence by phylum](figures/bar_has_pfam_by_phylum.png)

![CAZy presence by phylum](figures/bar_has_cazy_by_phylum.png)

![MEROPS presence by phylum](figures/bar_has_merops_by_phylum.png)

![Signal peptide presence by phylum](figures/bar_has_signal_peptide_by_phylum.png)

![TM helix presence by phylum](figures/bar_has_tm_helix_by_phylum.png)

## Top enriched Pfam domains among adhesion proteins

```
 domain_id  adhesion_count  adhesion_rate  background_count  background_rate  enrichment_ratio
 Flocculin            1493       0.001991                 0              0.0               inf
      Hyr1             586       0.000782                 0              0.0               inf
    Msg2_C             198       0.000264                 0              0.0               inf
   Kazal_1              97       0.000129                 0              0.0               inf
Cadherin_4              94       0.000125                 0              0.0               inf
   DUF5077              81       0.000108                 0              0.0               inf
   DUF3246              57       0.000076                 0              0.0               inf
    FIT1_2              53       0.000071                 0              0.0               inf
  Cornifin              52       0.000069                 0              0.0               inf
     DUF11              51       0.000068                 0              0.0               inf
```

![Top enriched Pfam domains](figures/top_domains_pfam.png)

![Top enriched CAZy families](figures/top_domains_cazy.png)

![Top enriched MEROPS hits](figures/top_domains_merops.png)

## Correlation with model confidence

![Length vs. adhesion probability](figures/scatter_length_vs_probability.png)

![Ser/Thr/Pro vs. adhesion probability](figures/scatter_pct_ser_thr_pro_vs_probability.png)

![Species mean length vs. adhesion fraction](figures/scatter_species_mean_length_vs_fraction.png)

## Headline findings

Adhesion-predicted proteins are consistently enriched for serine+threonine+proline
content relative to their species-matched background across every one of the 9
phyla tested (`tables/clade_phylum_pct_ser_thr_pro_mannwhitney.csv`), and this
holds up after BH correction in all 9: Ascomycota (median 30.3% adhesion vs.
18.7% background, n=555,861 each, p_bh≈0), Basidiomycota (27.0% vs. 20.0%,
n=135,340 each, p_bh≈0), Chytridiomycota (26.5% vs. 17.9%, n=14,817 each,
p_bh≈0), Mucoromycota (23.9% vs. 17.8%, n=23,786 each, p_bh≈0), Zoopagomycota
(28.2% vs. 18.4%, n=18,647 each, p_bh=0), and even the small-n phyla —
Blastocladiomycota (27.3% vs. 19.5%, n=938, p_bh≈4.5e-100), Cryptomycota
(27.5% vs. 17.3%, n=40, p_bh≈1.0e-06), Microsporidia (26.2% vs. 16.2%, n=65,
p_bh≈3.2e-05) and Sanchytriomycota (28.0% vs. 16.1%, n=28, p_bh≈1.9e-08).
The direction and roughly 6-10 percentage-point magnitude are uniform across
the whole fungal tree sampled here — a pattern consistent with Ser/Thr/Pro-rich,
mucin-like or GPI-anchor-proximal sequence stretches that are a recognized
structural hallmark of fungal cell-wall adhesins (e.g. FLO/EPA/Als-family
proteins), since Ser/Thr residues are the attachment sites for the
O-mannosylation that decorates these surface proteins.

Length shows no such uniform pattern (`tables/clade_phylum_length_mannwhitney.csv`):
of the 9 phyla tested, 6 remain significant after BH correction. In the two
best-powered phyla, adhesion proteins are actually *shorter* than
background (Ascomycota: median 345 vs. 400 aa, n=555,861 each, p_bh≈0;
Basidiomycota: 310 vs. 362 aa, n=135,340 each, p_bh≈0), as are Mucoromycota
(302 vs. 322 aa, p_bh≈1.04e-31) and Zoopagomycota (338 vs. 357 aa, p_bh≈0.0082),
while two phyla go the other way — adhesion proteins *longer* than
background — Chytridiomycota (432 vs. 364 aa, n=14,817 each, p_bh≈5.3e-116)
and, at much lower power, Sanchytriomycota (445.5 vs. 323 aa, n=28 each,
p_bh≈0.0449, the weakest of the six significant results). The remaining 3
phyla — Blastocladiomycota,
Cryptomycota, and Microsporidia — show no significant difference after
correction (p_bh=0.61, 0.47, 0.19 respectively). So unlike Ser/Thr/Pro
content, length is not a lineage-independent signature of the
adhesion-predicted set, and its direction is not even consistent in sign
across the phyla where it is significant.

`figures/scatter_length_vs_probability.png` and
`figures/scatter_pct_ser_thr_pro_vs_probability.png` are heavily overplotted
at n=749,697 points each (drawn at alpha=0.3, s=10) and cannot reliably show
a density trend by eye at this scale — the dense region of the point cloud
saturates into a solid block regardless of whether an underlying trend
exists. We do not claim "no relationship" from these plots; assessing
whether model confidence tracks length or Ser/Thr/Pro content would require
a binned-median or hexbin view, which is left for a follow-up (see Follow-up
section). `figures/scatter_species_mean_length_vs_fraction.png` is plotted
at the species level (one point per species, not per protein) and so is not
subject to the same overplotting problem; it is a diffuse cloud centered
around 400-500 aa mean length and a 1-2% species adhesion fraction, with no
visible trend as mean length increases toward 1000-1400 aa.

Domain/topology presence also diverges sharply from length in how consistent
it is across phyla: adhesion proteins carry a signal peptide far more often
than background in every phylum checked (e.g. Ascomycota 74.3% vs. 6.7%,
Basidiomycota 61.4% vs. 4.9%; `tables/clade_phylum_domain_fractions.csv`),
consistent with these being predicted secreted/cell-surface proteins, while
carrying a MEROPS peptidase hit and a TM helix *less* often than background
(Ascomycota MEROPS 0.5% vs. 3.1%, TM helix 10.8% vs. 18.8%) and an overall
Pfam hit less often too (Ascomycota 39.4% vs. 74.7%) — i.e. a large fraction
of the adhesion set carries no recognized Pfam domain at all, consistent with
many fungal adhesins being fast-evolving, repeat-rich, poorly conserved
sequences that escape standard domain databases.

Among proteins that *do* have a domain hit, `tables/top_domains_pfam.csv`
shows the most frequent Pfam domains among adhesion-predicted proteins that
have zero occurrences in the size-matched background sample: 826 distinct
Pfam domains met the ≥5-adhesion-hit threshold for inclusion in this table,
of which 95 have `background_count=0` (an infinite, tie-classed enrichment
ratio); ties are broken by adhesion-hit count (`domains.py`'s
`top_domain_table`, sorted on `["enrichment_ratio", "adhesion_count"]` with a
stable mergesort), so the table is deterministic and reproducible but the
particular 20 domains shown are the most *frequent* background-absent
domains, not necessarily the most *biologically distinctive* ones among the
full set of 95. With that framing, the top hits are Flocculin (1,493
adhesion hits, 0 background) — a domain named for and found in the
fungal flocculin/FLO-family adhesins that this whole analysis is built
around — and Hyr1 (586 hits), a domain associated with the *Candida*
Hyr1/Iff GPI-anchored cell-wall protein family, both consistent with the
adhesin narrative; also present are Cadherin_4 (94 hits) and Ig_3 (32
hits), both cell-adhesion-associated repeat domains, Msg2_C (198 hits, the
C-terminal domain of *Pneumocystis* major surface glycoprotein, a surface
adhesin family in another fungal lineage), and several domains of unknown
function (DUF5077, DUF3246, DUF11, DUF4573, DUF6209, DUF642) that recur in
this kind of screen without a resolvable functional annotation. In
`tables/top_domains_cazy.csv`, the most
enriched families by finite ratio are AA1_3 (12,493 adhesion hits vs. only 2
background, ratio≈6246.5) and AA1_1 (9,279 vs. 5, ratio≈1855.8) — both
"Auxiliary Activity family 1" multicopper-oxidase/laccase domains, an enzyme
family sometimes discussed in the context of fungal cell-wall remodeling and
melanization, though see Caveat 7 below on why an enrichment of this
magnitude should not be treated as settled functional evidence — plus a
cluster of
carbohydrate-binding modules (CBM29, CBM6, CBM43, CBM24, CBM63, CBM20) that
have no catalytic activity of their own but mediate binding to cell-wall or
extracellular polysaccharides, exactly the kind of adhesive function this
analysis is targeting. In `tables/top_domains_merops.csv`, one accession,
MER0003426, stands out with by far the largest adhesion count (1,867 hits vs.
95 background, ratio≈19.65) of any entry in that table, well above the two
next most common finite-ratio accessions MER0038531 (40 vs. 1, ratio≈39.99)
and MER0157380 (20 vs. 1, ratio≈19.99); per Caveat 1 below these are
individual reference-sequence accessions rather than resolvable peptidase
family codes, so no functional family name can be attached to them here.

The three top-domain enrichment charts should not be read as showing the
same phenomenon for the same reason. In `figures/top_domains_pfam.png`, all
20 plotted Pfam domains have `background_count=0` — every bar is legitimately
tied at the plot's display cap because there is truly no background
occurrence to divide by (enrichment is infinite, not merely large), so this
chart cannot distinguish relative magnitude among its 20 domains. In
`figures/top_domains_cazy.png`, only 9 of the 20 domains are background-absent
(`inf`); the other 11 (AA1_3, AA1_1, AA13, AA1_2, AA1, GH132, AA11, CBM24,
CBM63, CBM20, GH43_23) have finite but very large enrichment ratios (from
about 73.5 up to about 6246.5) that still exceed the plot's display cap of 50
and so are also drawn at the cap — the chart looks uniform, but for a mixed
reason: some domains are genuinely background-absent, others are merely
extremely (but finitely) enriched. `figures/top_domains_merops.png` is the
one chart of the three that shows real, visible differentiation: while 17 of
its 20 domains are at `inf`/capped, MER0038531 (≈39.99), MER0157380 (≈19.99),
and MER0003426 (≈19.65) all fall below the display cap and are drawn as
visibly shorter bars than the rest — so, unlike the other two charts, the
MEROPS chart is not "every bar tied at the cap."

## Caveats

1. **MEROPS results are accession-level, not family-level** — this
   database has no peptidase family/clan lookup, so MEROPS hit identities
   are individual sequence accessions (e.g. `MER0080922`), not the
   family codes (e.g. "S08") domain experts would recognize; treat the
   MEROPS enrichment table as "has a hit in this specific reference
   sequence's neighborhood," not a family-level functional claim.
2. **`net_charge_ph7` is a coarse approximation** (simple Lys+Arg minus
   Asp+Glu residue count), not a full pKa-based charge model.
3. **Background sampling is random** (seeded, reproducible) — a different
   seed would produce a slightly different background set; the 1:1
   per-species matching mitigates but doesn't eliminate this.
4. **Non-independence across clades** (species within a genus/family
   share ancestry) applies here exactly as it does in the kingdom-wide
   survey — this pass does not repeat that survey's genus-averaging/
   mixed-model correction; treat clade-level patterns here as
   descriptive/hypothesis-generating, not confirmatory.
5. **Multiple testing across features and clades**: the Mann-Whitney
   tables are BH-corrected *within* each (rank, feature) combination, not
   jointly across every feature/rank combination in this report — a
   stricter joint correction would raise the significance bar further.
6. **Prediction-set circularity**: the "adhesion" set analyzed throughout
   this report is model-*predicted*, not experimentally validated — it is
   the output of a classifier trained on ESM-2 protein embeddings, which
   themselves encode amino acid composition. So the headline finding that
   "adhesion-predicted proteins are Ser/Thr/Pro-rich" partly reflects the
   classifier's own decision function reading back a composition signal
   correlated with its training set (FLO11/ALS-family sequences), not
   necessarily an independent biological discovery about a held-out
   ground truth. Treat this analysis as characterizing what the model
   learned to call "adhesion-like," not as an independent confirmation
   that these proteins are adhesins.
7. **CAZy AA1 (laccase) enrichment is a red flag, not just supporting
   evidence** — summed across the AA1/AA1_1/AA1_2/AA1_3 subfamily entries
   in `tables/top_domains_cazy.csv`, 52,697 adhesion-predicted proteins
   carry an AA1-family domain versus only 33 in the size-matched
   background (a ratio of roughly 1,597-fold; see the table for
   per-subfamily counts). Laccases (multicopper oxidases) are secreted
   **enzymes**, the opposite of the non-enzymatic surface-adhesin category
   this analysis is targeting per the design spec. An enrichment of this
   magnitude is at least as consistent with a systematic classifier
   false-positive mode (e.g. shared secretion-signal or compositional
   features between laccases and true adhesins driving spurious positive
   calls) or an artifact of the training-set composition as it is with
   "laccases are adhesive" — it should not be presented as settled
   functional evidence without this caveat.
8. **P-values at this sample size are not effect sizes**: several
   Mann-Whitney tables compare hundreds of thousands of proteins per group
   (e.g. n=555,861 for Ascomycota) — at this scale, p≈0 is close to
   inevitable even for practically negligible effects. Readers should
   judge practical magnitude from the reported medians (already shown in
   each table) and, where available, the `epsilon_squared` column already
   computed in the `*_omnibus_*.csv` files (e.g.
   `clade_phylum_length_omnibus_adhesion.csv` reports
   `epsilon_squared`≈0.0067 for the phylum effect on length — a negligible
   clade effect despite p≈0 in the pairwise tests) rather than treating a
   tiny p-value alone as evidence of a large or important effect.

## Follow-up

Whether an ESM-2 embedding-space clustering reveals distinct sub-types
within the adhesion-predicted set (e.g. FLO11-like vs. Als-like vs. other
architectures) is deferred to a follow-up phase, informed by which
domain/property patterns found here are worth validating against.
