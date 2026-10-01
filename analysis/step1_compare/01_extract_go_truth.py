#!/usr/bin/env python3
"""D1: extract the GO truth set (spec 3.1) into truth_set.tsv.gz and counts.tsv.

Inputs are read from --input-dir (default $STEP1_WORKDIR/downloads). Every input is checked
against manifest.tsv before it is read; a strict-mode SHA-256 mismatch stops the run.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import gaf
import go_obo
import manifest
import paths
import truth_table


def read_species(path: str | Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(
    species_rows: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
    input_dir: Path,
    out_dir: Path,
    obo_name: str = "go-basic.obo",
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    by_file = {r["file"]: r for r in manifest_rows}
    obo_path = input_dir / obo_name
    obo_sha = manifest.verify_against_manifest(obo_path, by_file[obo_name])
    ontology = go_obo.parse_obo(obo_path)
    truth_rows, count_rows, log = [], [], {"obo": obo_name, "obo_sha256": obo_sha, "sources": []}
    for sp in species_rows:
        gaf_path = input_dir / sp["gaf_file"]
        sha = manifest.verify_against_manifest(gaf_path, by_file[sp["gaf_file"]])
        filtered = gaf.filter_gaf(gaf_path, ontology, sp["taxon_filter"] or None)
        info = truth_table.SourceInfo(
            source_id=sp["source_id"],
            species=sp["species"],
            taxon_id=sp["taxon_id"],
            in_clade=sp["in_clade"],
            role=sp["role"],
            source_file=sp["gaf_file"],
            source_sha256=sha,
            source_date=gaf.header_value(gaf_path, "date-generated"),
            obo_sha256=obo_sha,
        )
        rows = truth_table.build_truth_rows(filtered, ontology, info)
        truth_rows.extend(rows)
        count_rows.append(truth_table.count_rows(rows, filtered, sp["source_id"]))
        log["sources"].append(
            {"source_id": sp["source_id"], "sha256": sha, "dropped": dict(filtered.dropped)}
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    truth_table.write_tsv(out_dir / "truth_set.tsv.gz", truth_table.TRUTH_COLUMNS, truth_rows)
    truth_table.write_tsv(out_dir / "counts.tsv", truth_table.COUNT_COLUMNS, count_rows)
    (out_dir / "extract_log.json").write_text(json.dumps(log, indent=2, sort_keys=True) + "\n")
    return truth_rows, count_rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--manifest", default=str(paths.STEP1_DIR / "manifest.tsv"))
    parser.add_argument("--input-dir", default=None, help="default: $STEP1_WORKDIR/downloads")
    parser.add_argument("--out-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--sources", nargs="*", default=None, help="source_id values to run")
    args = parser.parse_args(argv)
    species_rows = read_species(args.species)
    if args.sources:
        species_rows = [r for r in species_rows if r["source_id"] in args.sources]
    input_dir = Path(args.input_dir) if args.input_dir else paths.downloads_dir()
    out_dir = Path(args.out_dir) if args.out_dir else paths.workdir()
    try:
        _, counts = run(species_rows, manifest.read_manifest(args.manifest), input_dir, out_dir)
    except manifest.DownloadError as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for c in counts:
        print(
            "\t".join(f"{k}={c[k]}" for k in ("source_id", "p_ext", "n_int", "n_sec", "ambiguous"))
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
