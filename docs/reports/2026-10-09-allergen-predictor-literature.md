# Allergen prediction literature for an evidence-only model (2026-10-09)

Purpose: inputs for the spec of an allergen predictor in `cellsurface_sorting_hat`. Evidence only,
`unvalidated`. Owner decision of 2026-10-09 lifts the "no allergenicity model" rule for this use.

Checking levels: **read** = abstract or full-text fields checked by fetch; **search** = search
snippet only, not opened; **unverified** = the source page did not load.

## Methods

| Method | Year | Inputs | Training data | Reported result | Availability | Checked |
|---|---|---|---|---|---|---|
| Allerdictor (Dang and Lawrence) | 2014 | overlapping 6-mers, sparse counts; linear SVM | ~3,907 allergens (IUIS, Allergome, SDAP, AllergenOnline, AllerMatch); Swiss-Prot non-allergens filtered for allergens; ~10 negatives per allergen | AUPRC 0.97 / 0.91 / 0.85 on benchmarks A / B / C | Web server does not resolve; no license found | read (paper) |
| Applm (Wong et al., arXiv 2508.10541) | 2025 | xTrimoPGLM (100 B-parameter protein language model) | not stated on the page | "consistently outperforms seven" methods; no numbers on the page | "open-source", no link given; CC BY-NC-ND | read (abstract page) |
| AllergenAI (BMC Bioinformatics) | 2025 | CNN on sequence | not checked | not checked | not checked | unverified (publisher redirect) |
| AllerTrans | 2025 | deep learning on sequence | not checked | sensitivity 97.91%, specificity 97.69% (2-fold evaluation, self-reported) | not checked | search |
| pLM4Alg | 2024 | ESM-2 embeddings, CNN head | not checked | accuracy up to 95.1%, AUC up to 0.99 (self-reported) | not checked | search |
| AlgPred 3.0 | 2024-25 | ML plus IgE epitope mapping | 3,694 allergens and 3,694 non-allergens | no accuracy in snippet | web server | search |
| AllerCatPro 2.0 | 2022-25 | sequence and predicted 3D structure similarity | 4,979 allergens, 162 low-allergenic, 165 autoimmune | "better accuracy" claimed, no numbers in snippet | web server | search |
| Fungal SVM classifiers (Bentham, 2017) | 2017 | annotations of ~1,927 fungal allergens, 194 genera; SVM | fungal allergens | not checked | said to be downloadable from the authors' site; not verified | search |

## What the recent work says about evaluation

- Applm's abstract names three hard cases: novel allergens with no close training example; allergens
  versus non-allergens among close homologs; and single mutations. This is the same concern as our
  leakage rules. The split design is the part to copy. The page gives no numbers.
- A 2026 preprint (arXiv 2606.22181, title only checked) reports that residue-level attributions in
  protein language models do not recover allergen epitopes. Do not claim epitope-level
  interpretation from attributions.
- No head-to-head benchmark of the tools above was found. Self-reported numbers use different
  splits. They are not comparable with each other or with our data.

## Fungal-specific work

- No current Aspergillus or fungal-only deep learning classifier was found in search.
- The 2017 fungal SVM is the closest prior work. It is old, and its availability is not verified.
- A 2026 bioRxiv preprint (AlphaFold-based fungal allergen structures) uses 1,289 WHO/IUIS
  non-fungal allergens for cross-kingdom comparison. It is a structure source, not a classifier.

## Implications for our design

1. Include a composition-only baseline in every run. `AGENTS.md` records that composition alone
   reached ROC-AUC 0.995 under leave-family-out CV on the existing model.
2. Include the 6-mer linear model (the paper's method) as a baseline. It is cheap.
3. Include an ESM-2 embedding with a linear head as the main candidate. It is the method most
   recent tools use. Run it only after the spec passes review.
4. Split by cluster (MMseqs2, 30% identity) and also by family. Report results on high-similarity
   pairs separately, following Applm's third scenario.
5. Weak negatives (Swiss-Prot, not allergen-annotated) train the model. Evaluation negatives need an
   independent reason, as task 05 requires. Report the two sets separately.
6. Report AUPRC and sensitivity at a fixed specificity with cluster-bootstrap intervals. Set gates
   after measurement.
7. Do not use the self-reported numbers above as targets.

## Full-text update (2026-10-09)

Full text read for: DeepPlantAllergy (Briefings Bioinf. 2025, bbaf605, from `to_import/allergen/`),
AllerTrans (Biol. Methods Protoc. 2025, bpaf040, from `to_import/allergen/`), Applm (arXiv 2508.10541v1,
from `to_import/allergen/`), AllergenAI (BMC Bioinf. 2025, via PMC12625376). Fukutomi and Taniguchi
(2015 fungal allergen review, `to_import/allergen/64_5.pdf`) read for background.

Key facts, each checked in the full text:
- **DeepPlantAllergy:** positives 955 plant allergens clustered to 614 at 85% identity. Negatives from
  UniProt reviewed entries with allergen, antigen and cancer words excluded, 1,817 retrieved, 769 kept
  after removing any with >20% identity to a positive. Cluster-level 70/20/10 split; test 61 positives.
  ESM-1B and ProtBert test MCC 0.879 and 0.873; one-hot encoding MCC 0.617. Code: github.com/Lilly-dh/DeepPlantAllergy.
- **AllerTrans:** ESM-2 650M and ProtT5-XL embeddings concatenated; three-layer network; 5-fold CV on the
  AlgPred 2.0 sets. The methods text describes no identity-based split. Independent AlgPred 2.0 test:
  sensitivity 93.5%, specificity 95.38%, MCC 0.900. Recombinant check on 42 + 42 UniProt entries.
  Code: github.com/faezesarlakifar/AllerTrans.
- **Applm:** similarity-aware pipeline with an inter-split threshold (T_S) and an inter-class threshold
  (T_C). Negatives are paired with positives at T_C. Average AUROC 0.844-0.864 with pLM embeddings, 0.742
  with one-hot and 0.684 with BLOSUM62. External sets By Date, Homologs and Mutations; average AUPRC
  0.700 against 0.202 for random. Performance falls as T_S falls and T_C rises (AUROC 0.942 easiest, 0.722
  hardest). Length Control helps only when the test set has uniform length. Code: github.com/brianwongsh/Applm.
  Licence of the code not stated in the text.
- **AllergenAI:** weak negatives are cupin family (PF00190) proteins, 2,775 kept after removing hits to
  SDAP at E <= 0.1. Test accuracy 0.94 and AUROC 0.97. The paper's counts for the training pool do not
  reconcile between Methods and Results.

Open items:
- Whether the 2017 fungal SVM files still exist.
- Licences of the DeepPlantAllergy, AllerTrans and Applm code before any reuse.
- Owner decisions Q1-Q4 in `docs/superpowers/specs/2026-10-09-allergen-predictor-spec.md`.
