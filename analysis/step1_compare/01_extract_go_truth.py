#!/usr/bin/env python3
"""D1: extract the GO truth set (spec 3.1) into truth_set.tsv.gz and counts.tsv.

Inputs are read from --input-dir (default $STEP1_WORKDIR/downloads). Every input is checked
against manifest.tsv before it is read; a strict-mode SHA-256 mismatch stops the run.
"""

import argparse
import csv
import sys
from pathlib import Path

import gaf
import go_obo
import manifest
import paths
import runinfo
import truth_table

OUTPUT_NAMES = ("truth_set.tsv.gz", "counts.tsv", "extract_log.json")

EPILOG = """\
Columns are listed in COLUMNS.md. The nohom_* columns of counts.tsv reproduce d1_count.py and
are NOT used as truth. The direct_* columns (label unchanged without homology codes) are the
direct-evidence truth.
"""


SPECIES_COLUMNS = (
    "source_id",
    "species",
    "taxon_id",
    "taxon_filter",
    "in_clade",
    "role",
    "gaf_file",
)


def read_species(path: str | Path) -> list[dict[str, str]]:
    """Read species.tsv. Raise ValueError for a missing column or a row with a wrong field count."""
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        header = reader.fieldnames or []
    missing = [c for c in SPECIES_COLUMNS if c not in header]
    if missing:
        raise ValueError(f"{path}: header lacks {missing}")
    for lineno, row in enumerate(rows, start=2):
        if None in row or None in row.values():
            raise ValueError(f"{path}: row {lineno} does not have {len(header)} fields")
    return rows


def _manifest_row(by_file: dict[str, dict[str, str]], name: str) -> dict[str, str]:
    if name not in by_file:
        raise manifest.DownloadError(f"{name}: no manifest row")
    return by_file[name]


def run(
    species_rows: list[dict[str, str]],
    manifest_rows: list[dict[str, str]],
    input_dir: Path,
    out_dir: Path,
    obo_name: str = "go-basic.obo",
    provenance: dict | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build the truth table and counts. Nothing is written until every source has succeeded."""
    if not species_rows:
        raise manifest.DownloadError("no sources selected")
    by_file = {r["file"]: r for r in manifest_rows}
    obo_path = input_dir / obo_name
    obo_sha = manifest.verify_against_manifest(obo_path, _manifest_row(by_file, obo_name))
    ontology = go_obo.parse_obo(obo_path)
    truth_rows, count_rows = [], []
    log = {
        "obo": obo_name,
        "obo_sha256": obo_sha,
        "input_dir": str(input_dir),
        "git_commit": runinfo.git_commit(),
        "python": runinfo.python_version(),
        "sources": [],
    }
    log.update(provenance or {})
    for sp in species_rows:
        gaf_path = input_dir / sp["gaf_file"]
        sha = manifest.verify_against_manifest(gaf_path, _manifest_row(by_file, sp["gaf_file"]))
        filtered = gaf.filter_gaf(gaf_path, ontology, sp["taxon_filter"] or None)
        date = gaf.header_value(gaf_path, "date-generated")
        info = truth_table.SourceInfo(
            source_id=sp["source_id"],
            species=sp["species"],
            taxon_id=sp["taxon_id"],
            in_clade=sp["in_clade"],
            role=sp["role"],
            source_file=sp["gaf_file"],
            source_sha256=sha,
            source_date=date,
            obo_sha256=obo_sha,
        )
        rows = truth_table.build_truth_rows(filtered, ontology, info)
        truth_rows.extend(rows)
        count_rows.append(truth_table.count_rows(rows, filtered, sp["source_id"]))
        log["sources"].append(
            {
                "source_id": sp["source_id"],
                "sha256": sha,
                "date_generated": date,
                "unknown_term_rows": filtered.unknown_term_rows,
                "dropped": dict(filtered.dropped),
            }
        )
    runinfo.atomic_write_all(
        out_dir,
        {
            "truth_set.tsv.gz": lambda p: truth_table.write_tsv(
                p, truth_table.TRUTH_COLUMNS, truth_rows
            ),
            "counts.tsv": lambda p: truth_table.write_tsv(p, truth_table.COUNT_COLUMNS, count_rows),
            "extract_log.json": lambda p: runinfo.write_json(p, log),
        },
    )
    return truth_rows, count_rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--manifest", default=str(paths.STEP1_DIR / "manifest.tsv"))
    parser.add_argument("--input-dir", default=None, help="default: $STEP1_WORKDIR/downloads")
    parser.add_argument("--out-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--sources", nargs="*", default=None, help="source_id values to run")
    args = parser.parse_args(argv)
    try:
        input_dir = Path(args.input_dir) if args.input_dir else paths.downloads_dir()
        out_dir = Path(args.out_dir) if args.out_dir else paths.workdir()
        species_rows = read_species(args.species)
        all_source_ids = [r["source_id"] for r in species_rows]
        if args.sources is not None:
            for source_id in args.sources:
                if source_id not in all_source_ids:
                    raise manifest.DownloadError(
                        f"unknown source {source_id}; valid ids: {', '.join(all_source_ids)}"
                    )
            species_rows = [r for r in species_rows if r["source_id"] in args.sources]
        provenance = {
            "species_sha256": manifest.sha256_file(args.species),
            "manifest_sha256": manifest.sha256_file(args.manifest),
            "arguments": list(argv) if argv is not None else sys.argv[1:],
            "all_sources": [r["source_id"] for r in species_rows] == all_source_ids,
        }
        _, counts = run(
            species_rows,
            manifest.read_manifest(args.manifest),
            input_dir,
            out_dir,
            provenance=provenance,
        )
    except (manifest.DownloadError, gaf.GafFormatError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for c in counts:
        print(
            "\t".join(f"{k}={c[k]}" for k in ("source_id", "p_ext", "n_int", "n_sec", "ambiguous"))
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
