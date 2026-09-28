# Class 2b: PRA/Ag2's fold, and who else has it

**2026-09-28.** Queued as class 2b in `docs/HANDOFF-HPCC.md` §4 ("structure + homology search").
Reframed below from a Foldseek structural search over the whole protein to a domain-level
question, because most of PRA is not a stable fold at all.

## 1. PRA/Ag2 is not one fold — most of it is disordered

`analysis/cocci_antigens/literature_antigens.fa` carries PRA's sequence (UniProt Q6QJA6, from
Zhu et al. 1996, *Gene* 181:121-125, [doi:10.1016/S0378-1119(96)00486-6](https://doi.org/10.1016/S0378-1119(96)00486-6)).
AlphaFold DB already has a model (`AF-Q6QJA6-F1`, v6, mean pLDDT 62 — moderate/low overall).
Per-residue pLDDT lines up exactly with the domain boundaries the 1996 paper reported from
Chou-Fasman/hydropathy prediction alone:

| region | aa | mean pLDDT | interpretation |
|---|---|---|---|
| N-terminal signal peptide | 1–18 | 51.4 | cleaved, low confidence expected |
| **CFEM domain (Pfam PF05730, e=1.5e-13)** | 20–84 | **83.7** | **confidently folded** |
| Pro/Thr tetrapeptide-repeat region | 89–141 | 51.8 | disordered/extended, not a fold |
| C-terminal GPI-anchor signal | 180–194 | 51.1 | cleaved, low confidence expected |

So the "fold" question is specifically about the **CFEM domain** (aa 20–84), not the antigenic
Pro/Thr-repeat tail that the classical serology literature focuses on (the tail is what an
antibody actually binds; the CFEM domain is what gives the protein its structure). No structure
prediction had to be run — AlphaFold DB already covers this TrEMBL accession.

## 2. CFEM is a real gene family in Coccidioides, not unique to PRA

Queried the Fungi_5k `functionalDB` duckdb's precomputed Pfam scan (no new HMMER run needed) for
PF05730 across all 71 Onygenales genomes in Fungi_5k — script:
`analysis/model_review/pra_cfem_survey.py`, output: `analysis/model_review/cfem_onygenales.tsv`
(untracked/regenerable, same convention as other `*_profile_*.tsv` outputs).

- **415 CFEM-domain proteins across 71 Onygenales genomes** (~6/genome on average).
- **7 in *C. immitis* RS, 7 in *C. posadasii* Silveira** — a modest paralog family, not a
  single gene.
- PRA is the strongest CFEM hit in both genomes by e-value (`FA2214EC_002110-T1`,
  `FA93A29E_002418-T1`, domain coordinates 20–84/20–84, matching InterPro's call on Q6QJA6
  exactly) — confirming these two records ARE the PRA orthologs, not a different Fungi_5k
  numbering artifact.
- **The other 6 CFEM proteins per genome are new leads**, not previously flagged anywhere in
  this repo's antigen or adhesin candidate lists.
- **9 of the 14 Coccidioides CFEM proteins have a high-confidence SignalP call (>0.999)**; the
  other 5 have no SignalP call. Per this project's standing rule 5 (absence ≠ evidence of
  absence), that is reported as "no call," not "not secreted" — Fungi_5k's SignalP coverage is
  known to be thin (~4% overall, per `docs/HANDOFF-HPCC.md`).

## 3. Not done, flagged rather than guessed

- **No Foldseek structural confirmation.** `foldseek` is not installed on this HPCC node (no
  module, no binary on PATH). The finding above rests on Pfam profile-HMM assignment (CFEM),
  which is a well-established, independently-curated domain family — a reasonable stand-in for
  "same fold," but not the same evidence as an actual structural superposition. If a real
  structural check is wanted later: `foldseek` ships as a single static binary (no root needed
  to install into `~/.local/bin` or a repo-local `bin/`), and could compare AlphaFold models of
  all 14 Coccidioides CFEM proteins pairwise, or against the PDB structures of characterized
  CFEM domains (e.g. *Candida albicans* Rbt5/Csa1-family heme receptors).
- **No functional annotation attempted for the other 6 CFEM paralogs per genome.** CFEM domains
  are broadly associated in the literature with fungal heme/iron-acquisition receptors (e.g.
  *C. albicans* Rbt5, Csa1) and cell-surface signaling (e.g. *Magnaporthe* Pth11) — general
  background, not verified here against these specific Coccidioides paralogs. Whether any are
  known iron-acquisition genes, secreted antigens, or uncharacterized is an open question.
- **Fungi_5k protein IDs (`FA2214EC_002110-T1` etc.) are not yet mapped back to FungiDB locus
  tags** (`CIMG_xxxxx`) or cross-referenced against `cocci_antigen_ranking.tsv` / the shipped
  adhesion classifier's scores. The existing `cocci_antigens/01_build_inputs.sh` idmap
  (pangenome reference ID → Fungi_5k ID) already does exactly this mapping for a different
  purpose and could be reused.

## 4. Suggested next step

Map the 6 non-PRA Coccidioides CFEM proteins to FungiDB locus tags, check whether they score as
antigens/adhesins in the existing pipelines, and check pangenome prevalence/copy-number the same
way `analysis/cocci_antigens/NOTES.md` did for SOWgp and Ag2/PRA — this is the natural
continuation of both class 2a (pangenome variability) and this class 2b survey.
