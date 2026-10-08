# Adjudication of the repeat-mechanism labels (2026-10-06)

Script: `analysis/calibration_truth/adjudicate_mechanism_labels.py`. Inputs (sha256, first 16 characters):
agent-A table `49bdc7ac7dacb14a`, agent-B table (version of 2026-10-06, after its last edit) `3b83b7b895c0942b`, UniProt repeat features `6a358dcca6df8627`.
The result is the same with agent-B's earlier committed version (checked). Re-run the script when the agent
tables change.

**Rule.** The evidence tier decides: paper statement for that protein (3), UniProt repeat features (2), family
inference (1), background knowledge (0). `unknown` is no vote. Two labels at the same top tier go to an expert.

**Two questions are kept apart.** `repeat_evidence_tier` answers "has the protein a tandem repeat region?" and
only `2a` votes count. `final_label` answers "what mechanism does the evidence give?". A vote of `other` is
not evidence against repeats: agent-B used it when the paper it read gives another mechanism or none.

**Outcome.** 96 rows: 74 agreement, 15 settled by tier, 5 no claim, 2 conflicts. Clusters with a repeat region:
8 on paper statements, 19 with UniProt, 35 with family inference. Needed for the cluster floor: 20. For
`estimated` the half-width also matters (about 61 clusters at sensitivity 0.8).

**MAD1 (Q2LC49).** Agent-A `2b-i` (family inference) is overruled: agent-B `2a` (paper) and UniProt (8 repeat
features) agree. Final label `2a`, tier paper-stated. MAD2 (Q2LC47) already agreed.

**The 2 conflicts** (`expert_review_packet.tsv`, expert columns blank for a person to fill):
- Q6FTA2, *C. glabrata*, "Agglutinin-like protein N-terminal domain-containing protein".
- Q59TP1, RBT1, *C. albicans*.
Both are evidence level E2 (domain membership, no protein-level paper), UniProt records no repeats, and neither
agent cites a paper about the protein itself (agent-A cites the ALS1 and HWP1 papers; agent-B cites a Pfam
match). **Recommendation from the assistant, not from an expert:** exclude both from the controls, and set the
repeat label to `unknown`. A person should confirm.

## Owner decision of 2026-10-07 on the two conflicts

`owner_decisions.tsv` records the decision. It is a separate file because `adjudicate_mechanism_labels.py` rewrites
`expert_review_packet.tsv` with blank expert columns on every run. The script does not read `owner_decisions.tsv`.

- Q6FTA2 and Q59TP1 (RBT1) are excluded from the repeat controls. Their repeat label is `unknown`. This is an owner
  decision, not an expert review.
- RBT1 keeps its adhesin E2 label. The note cites PMIDs 10978273 (Braun 2000), 19837954 (Ene and Bennett 2009) and
  21414038 (Bonhomme 2011).
- The cluster counts do not change (8 / 19 / 35). Each of the two rows shares its cluster (P20840, P46593) with a
  row at the `database_annotation` tier.
- Braun 2000 was read in full (`to_import/genetics0031.pdf`). The two RBT1 alleles differ by a segment of
  aa 612 to 640 (the genome project sequence lacks it). The paper does not describe a tandem repeat in RBT1.
  The RBT1 repeat label stays at the family-inference tier.
