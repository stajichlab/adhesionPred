#!/usr/bin/env python3
"""Add the Pfam family-table, TMHMM and repeat02 modules to existing sorting-hat work directories and re-run the engine.

For each species folder `<runs-root>/<species>/` (with `p.faa` and `wd/`, made by aga1_acceptance.py) the script merges the chunk outputs of the
HPCC proteome array (`pfam_family/<species>_NNN/domtbl.txt`, `tmhmm/<species>_NNN/tmhmm.tsv`), runs repeat02 (02_repeat_profile.py) on the FASTA,
writes the modules `tm`, `pfam_*` (all of them, from data/sorting_hat/family_table.tsv) and `repeat02` into `wd/modules`, and writes the engine
output to `<species>/out_domains` (nothing is overwritten; a species folder with an existing `out_domains` is refused).

Run it with the Python of the project environment, e.g.
    python3 -I analysis/sorting_hat_run/aga1_domain_modules.py --runs-root R --taxdump nodes.dmp --results RESULTS_DIR [--results ...]
Each --results folder holds `pfam_family/`, `tmhmm/` chunk folders. Downloaded tool output is untrusted data; it is only parsed here.
"""

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

CHUNK = re.compile(r"^(?P<species>.+)_(?P<n>\d{3})$")


class DomainError(Exception):
    pass


def chunks(results, tool, species):
    found = []
    for base in results:
        d = Path(base) / tool
        if d.is_dir():
            found += [
                p
                for p in sorted(d.iterdir())
                if (m := CHUNK.match(p.name)) and m["species"] == species
            ]
    if not found:
        raise DomainError(f"{species}: no {tool} chunks in {results}")
    return found


def provenance(path):
    d = {}
    for line in Path(path).read_text().splitlines():
        if ": " in line:
            k, v = line.split(": ", 1)
            d.setdefault(k, v)
    return d


def merge_tmhmm(files, out):
    header, n = None, 0
    with open(out, "w") as fh:
        for f in files:
            lines = (f / "tmhmm.tsv").read_text().splitlines()
            if not lines:
                raise DomainError(f"{f}: empty tmhmm.tsv")
            if header is None:
                header = lines[0]
                fh.write(header + "\n")
            elif lines[0] != header:
                raise DomainError(f"{f}: header differs from the first chunk")
            for line in lines[1:]:
                fh.write(line + "\n")
                n += 1
    return n


def merge_domtbl(files, out):
    """Concatenate the hits of finished chunks (each must end with the `# [ok]` trailer) under one set of comments and one trailer."""
    head, hits = [], []
    for i, f in enumerate(files):
        lines = (f / "domtbl.txt").read_text().splitlines()
        if not any(line.strip() == "# [ok]" for line in lines[-3:]):
            raise DomainError(f"{f}: no '# [ok]' trailer; the hmmsearch run is not finished")
        if i == 0:
            head = [
                x
                for x in lines
                if x.startswith("#") and x.strip() != "# [ok]" and not x.startswith("# Target file")
            ]
        hits += [x for x in lines if x.strip() and not x.startswith("#")]
    Path(out).write_text("\n".join(head + hits + ["# [ok]"]) + "\n")


def run_species(species, results, runs_root, nodes, taxon, tmp):
    from cellsurface_sorting_hat.cli import main as engine
    from cellsurface_sorting_hat.modules.cli import main as module

    base = Path(runs_root) / species
    if (base / "out_domains").exists():
        raise DomainError(f"{species}: {base / 'out_domains'} exists; nothing is overwritten")
    fa, wd = base / "p.faa", base / "wd"
    common = ["--fasta", str(fa), "--workdir", str(wd)]
    tm_chunks, pf_chunks = (
        chunks(results, "tmhmm", species),
        chunks(results, "pfam_family", species),
    )
    prov = provenance(pf_chunks[0] / "provenance.txt")
    for c in pf_chunks[1:]:
        if provenance(c / "provenance.txt").get("pfam_sha256") != prov.get("pfam_sha256"):
            raise DomainError(f"{c}: Pfam release differs between chunks")
    hmmer = re.search(
        r"HMMER (\S+)", provenance(pf_chunks[0] / "provenance.txt").get("hmmer", "")
    ) or re.search(r"hmmer/(\S+)", prov.get("module", ""))
    tm_file, dom_file, rep_file = (
        tmp / f"{species}.tmhmm.tsv",
        tmp / f"{species}.domtbl",
        tmp / f"{species}.rep02.tsv",
    )
    merge_tmhmm(tm_chunks, tm_file)
    merge_domtbl(pf_chunks, dom_file)
    subprocess.run(
        [
            sys.executable,
            "-I",
            str(ROOT / "analysis/cocci_repeats/02_repeat_profile.py"),
            str(fa),
            "--out",
            str(rep_file),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    steps = [
        ["tm", *common, "--table", str(tm_file), "--tmhmm-version", "2.0c"],
        ["repeat02", *common, "--table", str(rep_file), "--script", str(ROOT / "analysis/cocci_repeats/02_repeat_profile.py")],
        ["pfam", *common, "--domtbl", str(dom_file), "--family-table", str(ROOT / "data/sorting_hat/family_table.tsv"),
         "--pfam-release", prov["pfam_release"], "--pfam-sha256", prov["pfam_sha256"], "--hmmer-version", hmmer.group(1) if hmmer else "3.4",
         "--sp-module", "step1_rule@R0", "--tm-module", "tm"],
    ]  # fmt: skip
    for argv in steps:
        if module(argv) != 0:
            raise DomainError(f"{species}: module {argv[0]} failed")
    rc = engine(
        [
            "--fasta",
            str(fa),
            "--taxon",
            str(taxon),
            "--taxdump",
            str(nodes),
            "--workdir",
            str(wd),
            "--out",
            str(base / "out_domains"),
        ]
    )
    if rc != 0:
        raise DomainError(f"{species}: the engine failed")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs-root", required=True)
    ap.add_argument("--taxdump", required=True, help="nodes.dmp")
    ap.add_argument(
        "--results",
        action="append",
        required=True,
        help="folder with pfam_family/ and tmhmm/ chunk folders",
    )
    ap.add_argument(
        "--taxa", required=True, help="TSV with columns species, taxon (the acceptance taxa table)"
    )
    a = ap.parse_args(argv)
    taxa = {}
    for line in Path(a.taxa).read_text().splitlines()[1:]:
        f = line.split("\t")
        taxa[f[0]] = f[1]  # species, taxon_id
    with tempfile.TemporaryDirectory() as t:
        for species, taxon in taxa.items():
            run_species(species, a.results, a.runs_root, a.taxdump, taxon, Path(t))
            print(f"{species}: done")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except DomainError as err:
        print(f"FATAL: {err}", file=sys.stderr)
        sys.exit(1)
