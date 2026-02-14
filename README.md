# Adhesion Protein Predictor

A bioinformatics tool for predicting adhesion proteins from FASTA formatted sequence files.

# Requirements

This will run with pytorch on CPUs but will be much faster if on a GPU system with torchvision installed having CUDA bindings.

# Usage

### Training

I have built a set of FLO, ALS1 related proteins from Saccharomyces and Candida for starters. This seems to have some reasonable power.

```
python scripts/training.py --positive training/positive --negative training/negative
```

### Application 

You can run on a single file at a time and produce a report for each query file
```
mkdir -p query
pushd query
curl -O https://fungidb.org/a/service/raw-files/release-68/CalbicansSC5314/fasta/data/FungiDB-68_CalbicansSC5314_AnnotatedProteins.fasta
curl -O https://fungidb.org/a/service/raw-files/release-68/CneoformansJEC21/fasta/data/FungiDB-68_CneoformansJEC21_AnnotatedProteins.fasta
curl -O https://fungidb.org/a/service/raw-files/release-68/Spombe972h/fasta/data/FungiDB-68_Spombe972h_AnnotatedProteins.fasta
popd
for qorg in $(ls query/*.fasta)
do
   python scripts/predict.py --input $qorg --output $(basename $qorg .fasta).adhesion_predict.csv
done
```

You can run on a single folder and all results will be combined in a single file. It will look for all .fasta, .fa, .pep, .aa with or without .gz extensions.

```
python scripts/predict.py --input query --output Combinedquery_adhesion_predict.csv
```

## Development Setup

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
