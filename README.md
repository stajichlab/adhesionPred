# Adhesion Protein Predictor

A bioinformatics tool for predicting adhesion proteins from FASTA formatted sequence files.

## Setup

1. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Setup Pre-commit Hooks:**
   ```bash
   pre-commit install
   ```

## Usage

Place your FASTA files (with extensions `.faa`, `.pep`, `.pep.fa`, or `.fasta`) inside the \`input/\` directory.

To run the prediction:

\`\`\`bash
python adhesion_pred.py
\`\`\`