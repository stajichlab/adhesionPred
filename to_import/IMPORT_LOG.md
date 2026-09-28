# Literature import log

This project has no Mycelium living-repo (`.living/`) structure, so imports are tracked here
instead: one entry per paper dropped into `to_import/`, recording the citation, what was
extracted, and where it landed. Source files stay in `to_import/` for provenance; they are not
deleted after processing.

---

## molecules-23-03145.xml

**Citation:** Duarte Escalante E, Frías De León MG, Martínez García LG, Herrera J, Acosta
Altamirano G, Cabello C, Palma G, Reyes Montes MdR. "Selection of Specific Peptides for
*Coccidioides* spp. Obtained from Antigenic Fractions through SDS-PAGE and Western Blot Methods
by the Recognition of Sera from Patients with Coccidioidomycosis." *Molecules* 2018, 23(12),
3145. DOI: [10.3390/molecules23123145](https://doi.org/10.3390/molecules23123145). PMID:
30513599, PMC6321320.

**Read:** 2026-09-28. Full JATS XML (`to_import/molecules-23-03145.xml`), Table 3
("Identification of proteins that recognize *Coccidioides* spp.") specifically.

**What Table 3 contains:** 3 proteins identified by mass spec from SDS-PAGE/western-blot bands
(100 kDa ×2, 50 kDa ×1) that were most frequently recognized by antibodies in sera from 27
patients with confirmed coccidioidomycosis, with BLAST cross-reactivity against related fungi
(*C. immitis*/*C. posadasii*, *H. capsulatum*, *A. fumigatus*, *P. brasiliensis*, *P. lutzi*).

**Extracted into:**
- `analysis/cocci_antigens/literature_antigens.fa` — full-length sequences for the 3
  accessions (A0A0J6YAG6, A0A0J8R1I4, A0A0J6F383). All three are now **deleted from UniProtKB**
  ("not part of a reference proteome") and were recovered from **UniParc** instead
  (UPI0000D87081, UPI00065EBA54, UPI0001A7E586) — UniProt's archive keeps sequences after
  UniProtKB deletion.
- `analysis/cocci_antigens/literature_antigens_meta.tsv` — per-protein band, source species,
  and the paper's own cross-reactivity BLAST results (species, % identity), transcribed
  directly from Table 3.

**Not yet done:** these 3 proteins are not yet wired into `02_score_antigens.py` as a scored
input (distinct from the IEDB confounder panel in `00_download_iedb_antigens.py` — these are
serum-validated *Coccidioides* antigens, not cross-reactivity confounders from other fungi).

**Why this matters beyond antigen scoring:** A0A0J6YAG6 ("GPI anchored serine/threonine-rich
protein") is architecturally similar to the tandem-repeat surface proteins (class 2a) in the
queued retrain-test task (`docs/HANDOFF-HPCC.md` §4) — worth checking against
`analysis/model_review/repeat_structure_transfer.py` as a possible *Coccidioides*-native
positive example, not just an antigen candidate.
