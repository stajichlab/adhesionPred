# Fungal allergen scoping note (issue #19)

Note, 2026-10-05: the call names `allergen_homolog_hit` and `allergen_candidate` were renamed to `iuis_allergen_similarity` and `iuis_allergen_homolog`. The text below keeps the old names.

**Working note, 2026-10-04.** Input to `docs/superpowers/specs/2026-10-04-orchestrator-design.md`
(category `allergen_candidate`, decision D6). Script and data: `analysis/allergen_scoping/`.

> **Status.** Computational summary of public data. No allergen module exists. No model was
> trained or tested. Numbers marked "measured" come from the downloads on 2026-10-04
> (WHO/IUIS tables as served that day; UniProt release 2026_03, dated 02-September-2026).
> Numbers marked "reported" come from search-result summaries or web pages and were not
> re-derived.

## 1. Questions from issue #19

1. How many distinct fungal allergen families exist beyond *Alternaria*, *Aspergillus*,
   *Cladosporium* and *Malassezia*?
2. How many fungal allergens are surface or secreted (overlap with step 1), and how many are
   intracellular?
3. Is allergenicity predictable beyond homology to known allergens? (Baseline to beat: the FAO/WHO
   rule of 35% identity over 80 aa.)
4. Which databases can we use for training and classification?

## 2. Databases

| Database | Fungal content | Access | Terms | Status for us |
|---|---|---|---|---|
| **WHO/IUIS Allergen Nomenclature** (allergen.org) | **Measured:** 120 fungal allergen molecules from 31 species. Of 1,169 allergens in total, 95 are Ascomycota, 23 Basidiomycota, 2 Mucoromycota. 121 fungal isoallergen rows, 116 with a protein sequence, 104 with a UniProt accession. | CSV tables, no login: `csv.php?table=allergen`, `isoallergen`, `joint`, `idmapping` | "Free to download". Users must cite the site URL and a recent IUIS publication. | **Primary source.** Downloaded and summarised. This is the reference for allergen names. |
| **UniProt keyword KW-0020 (Allergen)**, Fungi (taxon 4751) | **Measured:** 111 entries, all reviewed (Swiss-Prot). 15 have "homolog" in the protein name (13 are *Arthroderma*). | REST API | Not checked by me. | **Second source.** Includes similarity-based "homolog" entries, so not all rows are demonstrated allergens. Use for sequence features (signal peptide, location, Pfam). |
| **AllergenOnline** (FARRP, Univ. of Nebraska) | **Reported:** version 24, released 2026-01-26, 2,373 sequence entries in 986 taxonomic-protein groups from 455 species (all kingdoms). The site lists "aero fungi" and "food fungi" categories. I did not get a fungal count. | Search page and FASTA search tool. A page `databasefasta.shtml` exists. I did not confirm how to download the full set. | Not stated on the page I read. | **Check with FARRP.** Peer-reviewed list with stated inclusion rules: "allergen" (IgE binding plus biological activity) and "putative allergen" (IgE binding only). |
| **COMPARE** (HESI) | **Reported:** the 2026 release (tenth) is announced; an entry needs literature evidence of IgE binding and an expert panel reviews it. No fungal count found. | FASTA plus metadata (xls). The download is a form: files are **emailed to the submitter**. | Not stated. | **Owner action:** submit the form with an email address. I did not submit it. |
| **Fungal Allergen Database** (Iwate University, Bioinformation 2019, 15:820-823) | **Reported:** 2,486 entries from 105 fungal allergen genes, with genomic data from four *Aspergillus* species. Sources are WHO/IUIS, DDBJ and AspGD. Entries have an allergenicity rank (1 to 4) and a sequence-homology rank. | Web portal (`fungusallergen.agr.iwate-u.ac.jp`), full download through TogoDB. | Not stated. | **Maintenance unclear.** I did not open the portal. Treat as a derivative of WHO/IUIS. |
| **AllFam** (Med. Univ. Vienna) | **Reported:** 151 allergen-containing Pfam families (all kingdoms) on Pfam 30.0 (2016). | Browse and search. Bulk export not confirmed. | Not stated. | **Stale.** Last update on the page: 2017-03-07. Local Pfam is release 38. |
| **AlgPred 2.0** (Raghava group) | **Reported:** training set of 10,075 allergens and 10,075 non-allergens, built from COMPARE, AllergenOnline, Swiss-Prot, AllerTOP and AlgPred; no pair of proteins above 40% similarity across the train/validation split. | Web server | Not stated. | **Method reference and baseline**, not a source of fungal truth. |
| **IEDB** | Already used for the antigen table (`data/curated/antigens/`). It has IgE epitope data. | Download | Not checked. | Possible epitope source. Not examined. |
| AllerBase, Allergome, SDAP, AllerTrans | Names appear in search results. **Not opened.** | | | Unchecked. |
| A 2026 bibliographic review of airborne fungal allergens (*Medicina* 62:1186, PMC13302811) | Not opened. | | | A lead for under-sampled genera. |

## 3. Answers from the data (measured)

### 3.1 Size and taxonomic spread (WHO/IUIS)

- 120 fungal molecules in 31 species. Exposure: 103 airway, 17 contact.
- By genus: *Aspergillus* 38, *Penicillium* 17, *Malassezia* 13, *Alternaria* 12, *Cladosporium* 10,
  *Coprinus* 5, *Curvularia* 4, *Fusarium* 4, *Trichophyton* 4, *Candida* 3, and 7 more genera with
  five or fewer.
- By order: Eurotiales 55, Pleosporales 18, Malasseziales 13, Capnodiales 10, Agaricales 8,
  Hypocreales 5, **Onygenales 4** (all *Trichophyton*), Saccharomycetales 3, Mucorales 2,
  Sporidiobolales 2.
- **No *Coccidioides*, *Histoplasma* or *Blastomyces* allergen is listed.** The only Onygenales
  entries are four *Trichophyton* proteases (Tri r 2, Tri r 4, Tri t 1, Tri t 4).

### 3.2 Question 1: families

The 97 IUIS fungal accessions found in UniProt carry 56 distinct Pfam domains; 11 entries have no
Pfam domain. The most common: Peptidase_S8 PF00082 (13), Inhibitor_I9 PF05922 (10), Ribosomal_60s
PF00428 (8), Redoxin PF08534 (6), Enolase_C PF00113 and Enolase_N PF03952 (6 each), Thioredoxin
PF00085 (6), AltA1 PF16541 (3), Pro_isomerase PF00160 (3), HSP70 PF00012 (3). So the answer to
question 1 is **at least 56 Pfam domains, dominated by a few enzyme and housekeeping families**,
and not a short list of fungus-specific allergen families. Only AltA1 (and Allergen_Asp_f_4,
PF25312, in Pfam) look allergen-specific. This counts domains in annotated entries. It is not a
count of independent protein families.

### 3.3 Question 2: surface or intracellular

Using UniProt annotation only (no prediction was run), for the 97 IUIS accessions found:

| Annotation | Count |
|---|---|
| Signal peptide feature | 30 |
| GPI-anchor feature | 1 (Asp f 9, a Crh-like chitinase/transglycosylase) |
| "Secreted" in subcellular location | 21 |
| Intracellular location (cytoplasm, mitochondrion, peroxisome, nucleus, ER, vacuole) | 24 |
| No location annotation | 51 |

For the 111 UniProt KW-0020 entries: 49 have a signal peptide (37 of the 96 that are not marked
"homolog").

Reading: **about one third of fungal allergens have a signal peptide.** Most of the rest are
cytosolic or organellar housekeeping proteins. The protein names in IUIS agree: of 120 molecules,
about 26 are proteases (7 vacuolar serine protease, 6 alkaline serine protease, 5 serine protease
and others), about 20 are enolase, aldolase or dehydrogenase types, about 20 are redox or chaperone
types (thioredoxin, cyclophilin, HSP, superoxide dismutase), and 9 are ribosomal proteins. These
keyword counts are rough and overlap.

Cell-surface members exist in the sets: Asp f 9 (Crh-like, GPI), Asp f 34 / PhiA (cell wall protein,
three UniProt entries), a class I hydrophobin (*Cladosporium*, entry Q8NIN9), cerato-platanin
(Pfam PF07249, which the architecture documents place in the surface class), and cell wall
mannoprotein 1 (O60025). Asp f 2 and Asp f 1 are already in our curated tables (`antigens.tsv`; Asp f 2
as an adhesin hard negative in `eurotiomycetes_seeds.tsv`).

**Consequence for design.** An allergen category cannot rely on the step 1 surface call. Most
allergens are not surface proteins. The orchestrator must let `allergen_candidate` be called
without `surface_glycoprotein` (the draft spec rule "allergen homology or allergen Pfam hit" already
does this). Allergen and surface calls overlap in a minority.

### 3.4 Question 3: is it predictable beyond homology

Not answered by data here. I ran no prediction. The facts that bear on the question:

- Many allergens are members of ordinary enzyme families (proteases, enolases, thioredoxins). A
  homology rule (the FAO/WHO 35% over 80 aa rule) will also hit non-allergen members of those
  families. AllergenOnline states that no simple score boundary makes cross-reactivity certain.
- A published ML method exists (AlgPred 2.0, section 2). Its reported training set is large and
  cross-kingdom. Its performance on fungal proteins is not known to me.
- Our own tools cannot test this yet, because no fungal non-allergen set with matching families
  exists in the repository.

### 3.5 Overlap with other tables

- UniProt KW-0020 fungi and IUIS (UniProt accessions) overlap on 63 accessions. 48 KW-0020 entries
  are not in IUIS (15 of the KW-0020 entries are "homolog" entries) and 41 IUIS accessions are not in
  the KW set (7 of those are not found in UniProt).
- A union of the two sources would give a larger positive set but mixes demonstrated allergens
  with similarity-based ones. The spec should keep the evidence class as a column.

## 4. Recommendations for the allergen module (for decision D6)

1. **Truth source.** Use WHO/IUIS fungal rows as the demonstrated-allergen set (120 molecules,
   116 sequences). Add AllergenOnline and COMPARE sequences once downloaded, keeping the source and
   AllergenOnline's "allergen / putative allergen" class as columns. Keep UniProt "homolog"
   entries out of positives.
2. **Version 1 method: homology plus domain.** Hits to the IUIS sequences (identity and coverage
   reported, not just a yes/no) and allergen-specific Pfam models (AltA1, Allergen_Asp_f_4). Mark
   the module `unvalidated` until a held-out test exists. Do not treat a hit to a large enzyme
   family as an allergen call by domain alone.
3. **Evidence column, not only a call.** Report the best hit (allergen name, identity, coverage,
   source). Reports for the 35%/80 aa rule can be made from the same table.
4. **Validation design (needs its own spec).** Positives: IUIS fungal allergens, grouped by
   homology so that paralogs and isoallergens do not cross the train/test split. Negatives: proteins
   from the same families and the same genomes that are not listed as allergens. Their negative
   status is **not known**: absence from a database means "never tested". Report this limit.
5. **Clade scope.** The set covers Eurotiales, Pleosporales, Malasseziales, Capnodiales and a few
   others. For Onygenales it has only four *Trichophyton* proteases. For *Coccidioides* the module
   can report homology to known allergens only, and the report should say that no *Coccidioides*
   allergen is known to the sources checked.
6. **Licence and citation.** WHO/IUIS asks for citation of allergen.org and a recent IUIS
   publication. Check the terms for UniProt, AllergenOnline and COMPARE before any redistribution
   of derived tables. The files in `analysis/allergen_scoping/` are WHO/IUIS and UniProt extracts
   and are meant for internal analysis.

## 5. Open points

1. Download AllergenOnline and COMPARE (COMPARE needs an email form), then count fungal entries
   and check the overlap with WHO/IUIS.
2. Open AllerBase, Allergome, SDAP and the Iwate portal, and check licences and maintenance.
3. Read the 2026 *Medicina* review for genera with few listed allergens.
4. Run SignalP and PredGPI on the IUIS fungal sequences, to replace the UniProt-annotation counts in
   section 3.3 (51 of 97 entries have no location annotation).
5. Decide the negative set for validation (recommendation 4).
6. Check whether any *Coccidioides* proteins have IgE-binding evidence in the IEDB data already
   used for the antigen table.
