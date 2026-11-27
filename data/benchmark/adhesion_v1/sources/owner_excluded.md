# Rows removed from the sources by owner decision (2026-11-26)

Both rows were curated and verified, then removed because they fail a rule the owner set. They are
recorded here in full so the decision is reversible and auditable; neither is in `proteins.tsv`.

## J9VXZ9 — Cryptococcus neoformans H99 Mpr1 (was: hard_negative)

**Why removed.** The curation rule admits a protein as an adhesin on `het_expr`: expression in a
non-adherent host confers adhesion. Mpr1 meets exactly that — expressing it in *S. cerevisiae* makes
yeast attach to and cross brain endothelium (PMID 24895304). Labelling it a hard negative on the
grounds that its mechanism is proteolysis would let curator judgement override the evidence rule,
which is the failure mode the rule exists to prevent. The accession is also an inference: the paper
names no locus tag, and CNAG_04735 was chosen as the only peptidase M36 in the H99 genome.

Excluding it costs nothing measurable: Tremellomycetes has no positives and so cannot be a scored
fold either way.

Row as curated:

```
J9VXZ9	uniprot	MPR1	Extracellular metalloproteinase Mpr1 (M36 fungalysin, CNAG_04735)	Cryptococcus neoformans (strain H99)	235443	Tremellomycetes	hard_negative	m36_protease	function_shown	Mpr1 is a secreted M36 fungalysin with a predicted prodomain and a conserved Zn-binding catalytic motif whose proteolytic activity was measured directly for the recombinant protein (fluorescent protease assay, activity comparable to a trypsin control), and although mpr1 deletion in H99 strongly reduces CFU of cryptococci adhering to brain endothelium at 20 min, 2 h and 8 h (restored in mpr1 complemented strain) and expression of CnMPR1 alone lets S. cerevisiae attach to and cross the endothelium, the proposed mechanism is proteolytic remodelling of host endothelial surface proteins rather than an adhesive domain.	24895304	831	pmid+acc	adhesion phenotype on deletion but indirect mechanism; owner decision 2026-11-26 to label hard_negative. CONTESTED - this is the strongest of the three adhesion phenotypes because the heterologous S. cerevisiae result would otherwise read as het_expr; it is labelled negative on mechanism, not on weak evidence, and should be excluded from any sensitivity estimate that treats negatives as untested. Accession assignment is an inference: the primary paper names no locus, and CNAG_04735 is the only peptidase M36 (Pfam PF02128) protein in the H99 genome, which is the background used throughout that study.
```

## Q8SWH2 — Encephalitozoon cuniculi Hsp70 ECU02_0100 (was: hard_negative)

**Why removed.** Two rules. Its negative status rests on a null adherence result beside EnP1, and
"a null adhesion result alone does not make a hard negative" (the same rule that keeps *Beauveria*
Adh1/Adh3 out). It is also an intracellular chaperone, so under the moonlighting rule it does not
belong in a surface-protein negative set; its Hsp70 function is inferred from homology rather than
assayed in PMID 20811483.

Row as curated:

```
Q8SWH2	uniprot	ECU02_0100	Heat shock protein 70 family member ECU02_0100 ("C1" Hsp70)	Encephalitozoon cuniculi (strain GB-M1)	284813	Microsporidia	hard_negative	intracellular_hsp70_chaperone_no_adherence	function_shown	Immunoelectron microscopy placed the protein on internal structures and not the spore surface, and in the same spore adherence inhibition assay in which anti-EnP1 antibodies blocked adherence, anti-C1 Hsp70 antibodies did not inhibit spore adherence to host cells.	20811483	683	pmid+acc	Explicit experimental negative control run side by side with the EnP1 positive (Q8SWL3) in one assay, which is why it is worth keeping despite being intracellular rather than a wall protein.
```

## Also held out, never added

- ***Beauveria* Adh1/Mad1 and Adh3** — named adhesins whose papers report a null adhesion result,
  with no alternative function shown. Neither label fits.
- ***Vairimorpha ceranae* NcSWP8** (EEQ83062.1, 189 aa, verified) — antibody blocking lowered the
  *infection rate*, not adherence. Calling it an adhesin overreaches; calling it a hard negative
  could bury a real adhesin. Ready to promote if a later version adds an indirect tier.
