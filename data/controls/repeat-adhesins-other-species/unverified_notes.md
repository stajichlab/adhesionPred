# Unverified notes / open questions (task 08, agent-C)

1. **Leakage is `unknown`-ish.** I assigned `leakage=none` to every new row because they come from
   taxa outside the tier-T5 reference proteomes and the Saccharomycotina training set. I did NOT read
   the repeat-detector scripts' reference sources to confirm no 30% homolog was used to tune
   detector 02 or 14. The owner should confirm. `analysis/cocci_repeats/` and
   `analysis/model_review/` history were not exhaustively audited this session.

2. **YPS3 accession.** Yps3p was resolved to `A0ACF1AYZ2` (gene `I7I48_01120`), from the G217B
   lineage now classified *Histoplasma ohiense* (taxid 2902605, formerly *H. capsulatum*). I could
   not independently re-derive the locus tag from the 2005 paper's figures; the subagent mapped it.
   A reviewer should double-check this accession maps to *H. capsulatum* YPS3.

3. **BcLysM1 adhesion.** The PubMed abstract of PMID 39655398 confirms chitin binding and plant
   immunity suppression but the mycelial-adhesion claim ("required for mycelial adhesion to
   hydrophobic surfaces, including bean leaves") is from the paper's full text, which I did not read
   in full. The adhesion sentence family is consistent with the title ("A LysM Effector Mediates
   Adhesion...") but should be confirmed against the article text/reviewer.

4. **FgHyd3 vs FgHyd2/FgHyd4.** The PMID 31031728 study reports defects in multiple hydrophobin
   single/triple mutants. I attributed the class I hydrophobin adhesion to FgHyd3 (I1RXJ5); the
   relative roles of FgHyd2/FgHyd4 are not fully separated. Reviewer should confirm the single-gene
   attribution.

5. **Cpl1 strain.** The Cpl1 negative (PMID 32805526) deletion phenotype was measured in *C.
   neoformans* strain JEC21, while the accession J9VQE5 is from the H99 genome. The two are both
   *C. neoformans var. grubii* but not the same strain. Noted; the negative evidence is strain-
   contextual.

6. **Hydrophobin family boundary.** The 2c hydrophobin family straddles adhesive and non-adhesive
   members (e.g. MHP1 E3 vs the B. cinerea/P. expansum hydrophobins that are non-adhesive). This is
   expected but means a family-domain call cannot separate them; the split must come from the
   per-species function, not the domain. Owner decision needed.

7. **Zero-state queries.** For Rhodotorula, Exophiala, Wallemia, Hortaea, Knufia, the "none found"
   conclusions rest on the limited PubMed abstract hits and the keyword scan; a deeper full-text
   search of closed OA repositories might surface more. Recorded as not found rather than
   fabricated.

8. **MMseqs2 coverage mode.** Clustering used MMseqs2 default coverage mode as in
   `analysis/calibration_truth/c1_truth_count.py` (`--min-seq-id 0.3 -c 0.5`). This matches the
   existing set so counts are comparable.
