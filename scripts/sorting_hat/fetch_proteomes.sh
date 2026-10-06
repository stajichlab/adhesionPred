#!/bin/bash
# Download the two extra A. fumigatus proteomes (decision D13) (provenance records are written later, Task 12).
# Run on the login node: WORKDIR=/path/work bash fetch_proteomes.sh
# A1163 (CEA10, FGSC A1163, CBS 144.89): UniProt proteome UP000001699 (9,942 proteins on 2026-10-04).
# W72310: NCBI GCA_040167795.1 (UCR_Afum_W72310_1.0).
set -euo pipefail
: "${WORKDIR:?export WORKDIR}"
DEST="$WORKDIR/proteomes"
mkdir -p "$DEST"
trap 'rm -f "$DEST"/.tmp.*' EXIT
curl -fsSL "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000001699&format=fasta" -o "$DEST/.tmp.A1163.faa"
mv "$DEST/.tmp.A1163.faa" "$DEST/Afum_A1163_UP000001699.faa"
BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/167/795/GCA_040167795.1_UCR_Afum_W72310_1.0"
curl -fsSL "$BASE/GCA_040167795.1_UCR_Afum_W72310_1.0_protein.faa.gz" -o "$DEST/.tmp.W72310.faa.gz"
mv "$DEST/.tmp.W72310.faa.gz" "$DEST/Afum_W72310_GCA_040167795.1.faa.gz"
echo "downloaded; now write provenance (Task 12 of the plan)"
