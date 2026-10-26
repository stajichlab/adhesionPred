#!/usr/bin/env python3
"""Add the allergen_homology and antigen_lookup modules to existing sorting-hat work directories and re-run the engine.

Per species folder `<runs-root>/<species>/` (with `p.faa` and `wd/`): BLASTP of the proteome against the IUIS fungal allergen sequences
(options of scripts/sorting_hat/blast_allergen.sbatch: -evalue 1 -max_target_seqs 200 -seg no), module `allergen`, module `antigen`
(applicable taxon 246410 only; every other proteome gets `not_applicable`), then the engine into `<species>/out_allergen_antigen`.

The antigen ranking is keyed on C. immitis RS gene IDs (CIMG_nnnnn). The UniProt proteome of C. immitis RS is keyed on UniProt accessions, so a
protein map (accession -> gene) is made from the `GN=CIMG_...` field of the UniProt FASTA header (`--cimm-uniprot-fasta`); accessions without
such a field are `not_in_reference`. Run it with `python3 -I`; blastp and makeblastdb must be on PATH. Nothing is overwritten.
"""

import argparse
import hashlib
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

APPLICABLE = "246410"
RANKING = ROOT / "analysis/cocci_antigens/cocci_antigen_ranking.tsv"
ISOALLERGEN = ROOT / "analysis/allergen_scoping/iuis_fungal_isoallergens.tsv"
ALLERGEN = ROOT / "analysis/allergen_scoping/iuis_fungal_allergens.tsv"
CAP, EVALUE, SEG = 200, "1", "no"
BLAST_COLUMNS = "6 qseqid sseqid pident length qlen slen bitscore evalue"
GN = re.compile(r"\bGN=(CIMG_\d+)\b")


class StepError(Exception):
    pass


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def uniprot_protein_map(fasta, out):
    """accession -> CIMG gene from the UniProt header (`>tr|ACC|NAME ... GN=CIMG_nnnnn ...`); returns (n_mapped, n_proteins)."""
    n = m = 0
    with open(fasta) as fh, open(out, "w") as o:
        o.write("protein_id\tgene_id\n")
        for line in fh:
            if not line.startswith(">"):
                continue
            n += 1
            parts = line[1:].split()[0].split("|")
            acc = parts[1] if len(parts) == 3 else parts[0]
            g = GN.search(line)
            if g:
                o.write(f"{acc}\t{g.group(1)}\n")
                m += 1
    return m, n


def blast_version():
    return subprocess.run(
        ["blastp", "-version"], check=True, capture_output=True, text=True
    ).stdout.splitlines()[0]


def run_species(species, taxon, runs_root, nodes, allergen_fa, pmap, tmp, threads):
    from cellsurface_sorting_hat.cli import main as engine
    from cellsurface_sorting_hat.modules.cli import main as module

    base = Path(runs_root) / species
    if (base / "out_allergen_antigen").exists():
        raise StepError(f"{base / 'out_allergen_antigen'} exists; nothing is overwritten")
    fa, wd = base / "p.faa", base / "wd"
    for name in ("allergen_homology", "antigen_lookup"):
        if (wd / "modules" / f"{name}.json").exists():
            raise StepError(f"{wd / 'modules' / name}.json exists; nothing is overwritten")
    blast_out = tmp / f"{species}.blast.tsv"
    subprocess.run(
        ["blastp", "-query", str(fa), "-db", str(tmp / "allergens"), "-evalue", EVALUE, "-max_target_seqs", str(CAP), "-seg", SEG,
         "-outfmt", BLAST_COLUMNS, "-num_threads", str(threads), "-out", str(blast_out)],
        check=True,
    )  # fmt: skip
    common = ["--fasta", str(fa), "--workdir", str(wd)]
    steps = [
        ["allergen", *common, "--blast", str(blast_out), "--allergen-fasta", str(allergen_fa), "--blast-version", blast_version(),
         "--evalue", EVALUE, "--seg", SEG, "--max-target-seqs", str(CAP)],
        ["antigen", *common, "--taxon", taxon, "--applicable-taxa", APPLICABLE, "--ranking", str(RANKING), "--protein-map", str(pmap)],
    ]  # fmt: skip
    for argv in steps:
        if module(argv) != 0:
            raise StepError(f"{species}: module {argv[0]} failed")
    rc = engine(
        [
            "--fasta",
            str(fa),
            "--taxon",
            taxon,
            "--taxdump",
            str(nodes),
            "--workdir",
            str(wd),
            "--out",
            str(base / "out_allergen_antigen"),
        ]
    )
    if rc != 0:
        raise StepError(f"{species}: the engine failed")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs-root", required=True)
    ap.add_argument("--taxdump", required=True, help="nodes.dmp")
    ap.add_argument("--taxa", required=True, help="TSV with columns species, taxon_id")
    ap.add_argument(
        "--cimm-uniprot-fasta",
        required=True,
        help="UniProt FASTA of C. immitis RS (UP000001261) for the accession-to-gene map",
    )
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args(argv)
    from cellsurface_sorting_hat.modules.allergen import build_allergen_fasta

    taxa = {}
    for line in Path(a.taxa).read_text().splitlines()[1:]:
        f = line.split("\t")
        taxa[f[0]] = f[1]  # species, taxon_id
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        allergen_fa = tmp / "iuis_fungal_allergens.faa"
        n, skipped = build_allergen_fasta(ISOALLERGEN, allergen_fa, ALLERGEN)
        if n >= CAP:
            raise StepError(
                f"{n} allergen sequences; the -max_target_seqs cap {CAP} must exceed the database size"
            )
        print(
            f"allergen database: {n} sequences ({len(skipped)} skipped); sha256 {sha256(allergen_fa)}"
        )
        subprocess.run(
            [
                "makeblastdb",
                "-in",
                str(allergen_fa),
                "-dbtype",
                "prot",
                "-out",
                str(tmp / "allergens"),
                "-logfile",
                "/dev/null",
            ],
            check=True,
        )
        pmap = tmp / "protein_map.tsv"
        m, total = uniprot_protein_map(a.cimm_uniprot_fasta, pmap)
        print(
            f"protein map: {m} of {total} UniProt accessions have GN=CIMG_; sha256 {sha256(pmap)}; ranking sha256 {sha256(RANKING)}"
        )
        for species, taxon in taxa.items():
            run_species(species, taxon, a.runs_root, a.taxdump, allergen_fa, pmap, tmp, a.threads)
            print(f"{species}: done")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except StepError as err:
        print(f"FATAL: {err}", file=sys.stderr)
        sys.exit(1)
