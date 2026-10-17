# Family-matched negative controls

Negatives for the allergen predictor. Each negative is a reviewed (Swiss-Prot) fungal protein. It shares a Pfam domain with at least one allergen positive. It is not an allergen by UniProt annotation.

Date built: 2026-10-10. UniProt release 2026_03 (release date 02-September-2026), from the `X-UniProt-Release` header of the REST responses.

## Files

- `family_negatives.faa`: sequences, IDs `F_<accession>`.
- `family_negatives.tsv`: accession, family, family_name, organism, length, sha256 of the sequence, reviewed, overlap_with_candidates.
- `counts.md`: per-family counts, totals, cluster counts, skipped positives.
- `positive_pfam.tsv`: Pfam domains of each UniProt positive (step 1).
- `positive_ncbi.txt`, `positive_ncbi_map.tsv`: NCBI-ID positives and the failed mapping attempt.
- `cluster/`: MMseqs2 input (`input.faa`: 200 positives plus new negatives), clusters, and `cluster_members.tsv`.
- `scripts/`: stdlib Python 3.12 scripts, run in order: `step1_positives.py`, `step1b_idmap.py`, `step2_families.py`, `step3_build.py`, `step4_summary.py`, `step5_write_docs.py`. `uni.py` holds the REST helper. `step2_state.json` and `step3_state.json` hold the query log and selection state.

## Step 1: positive Pfam domains

Query per UniProt-format positive, in batches of 40:

    (accession:A OR accession:B OR ...)
    fields: accession,xref_pfam,reviewed,organism_name,length,protein_name

The 26 positives with NCBI protein IDs (contain a dot) were not given Pfam domains. See counts.md.

## Step 2: negative query per family

For each Pfam family used by a positive:

    (taxonomy_id:4751) AND (xref:pfam-PFxxxxx) AND (reviewed:true)
    AND (NOT keyword:KW-0020) AND (NOT cc_allergen:allergic) AND (length:[10 TO 850])

Fields: accession,reviewed,organism_name,length,protein_name,xref_pfam.

Filters applied in code after the query:
- Remove accessions in `positives/positives.tsv`.
- Remove accessions in `candidates/candidates.tsv`.
- Remove accessions already picked for an earlier family.

After selection, the final set was checked with a second query. Zero of the 516 final accessions have keyword KW-0020 or cc_allergen:allergic.

## Step 3: cap

At most 10 negatives per family. Sampling uses `random.Random(20261010).sample` on the accessions sorted by ID. There is no separate total cap. The total follows from the per-family cap and the number of families.

## Step 4: clustering

    /opt/linux/rocky/8.x/x86_64/pkgs/mmseqs2/17-b804f/bin/mmseqs easy-cluster \
      input.faa clu tmp --min-seq-id 0.3 -c 0.5 --cov-mode 0

This is the non-AVX2 build. It ran directly on host c07, not under SLURM.

## Sequences and names

Sequences come from `https://rest.uniprot.org/uniprotkb/accessions?format=fasta`. Pfam family names come from the InterPro API (`https://www.ebi.ac.uk/interpro/api/entry/pfam/<PF>`).
