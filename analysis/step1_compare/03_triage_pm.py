#!/usr/bin/env python3
"""D8: triage P-ext genes with a non-IEA plasma-membrane term (P-gpi, PM-TM, pm-unresolved).

Reads $STEP1_WORKDIR/truth_set.tsv.gz, species.tsv and curated_gpi.tsv (checked by
d8_triage.read_curated_gpi). Queries UniProtKB REST.
Writes d8_triage.tsv, d8_gpi_outside_pext.tsv, d8_curated_conflicts.tsv,
curated_gpi_unmatched.tsv, d8_counts.tsv, truth_set_triaged.tsv.gz and the
raw UniProt JSON pages (d8_uniprot/; pages of older runs are not deleted and can remain) to the
work directory. Stops (exit 2, no output tables) on an invalid curated_gpi.tsv (missing column,
bad value, or a row that matches no truth gene), on an HTTP error, a response without the
expected fields, or when PM candidates exist and none of them matches a UniProt entry. Outputs are written to temp names and moved with os.replace
after all are complete, so a failure leaves earlier outputs untouched. d8_run.json records the
selected sources, `all_sources`, the truth set SHA-256, the UniProt release, the git commit, the
Python version and the arguments. The script refuses a truth_set.tsv.gz whose extract_log.json
does not say `all_sources: true`, unless --allow-partial-truth-set is given. A run with
--sources replaces the full outputs: check `all_sources` in d8_run.json before you use them.
"""

import argparse
import http.client
import json
import sys
import urllib.error
from pathlib import Path

import d8_triage
import manifest
import paths
import runinfo
import truth_table

TRIAGE_COLUMNS = (
    "source_id",
    "gene_id",
    "symbol",
    "d8_class",
    "d8_reason",
    "uniprot_accessions",
    "reviewed",
    "gpi_eco",
    "tm_count",
    "tm_eco",
)
OUTSIDE_COLUMNS = ("source_id", "gene_id", "symbol", "label", "uniprot_accession", "gpi_eco")
OUTPUT_NAMES = (
    "d8_triage.tsv",
    "d8_gpi_outside_pext.tsv",
    "d8_curated_conflicts.tsv",
    "curated_gpi_unmatched.tsv",
    "d8_counts.tsv",
    "truth_set_triaged.tsv.gz",
    "d8_run.json",
)
D8_COUNT_COLUMNS = (
    "source_id",
    "pm_candidates",
    "p_gpi",
    "pm_tm",
    "pm_unresolved",
    "no_uniprot_entry",
    "gpi_feature_no_evidence",
    "organism_curated_gpi_entries",
    "curated_gpi_outside_pext",
    "uniprot_release",
)


def uniprot_fetch(url: str, raw_dir: Path, tag: str) -> tuple[list[dict], str]:
    """Return (entries, release). Every page is saved as JSON for provenance."""
    entries, release = [], ""
    raw_dir.mkdir(parents=True, exist_ok=True)
    for n, (body, headers) in enumerate(manifest.uniprot_pages(url)):
        (raw_dir / f"{tag}_{n:03d}.json").write_bytes(body)
        page = json.loads(body)
        if "results" not in page:
            raise d8_triage.UniprotError(f"{url}: response has no 'results' field")
        entries.extend(page["results"])
        release = release or headers.get("x-uniprot-release", "")
    return entries, release


def triage_source(sp, rows, literature, fetch):
    mapping = sp["id_mapping"]
    candidates = [r for r in rows if r["pm_candidate"] == "yes"]
    entries, release = [], ""
    terms = d8_triage.query_terms(candidates, mapping)
    for i, url in enumerate(d8_triage.search_urls(sorted(terms))):
        got, release = fetch(url, f"{sp['source_id']}_candidates_{i}")
        entries.extend(d8_triage.parse_entry(e) for e in got)
    matched_any = any(
        d8_triage.entries_for_gene(r["gene_id"], mapping, entries) for r in candidates
    )
    if candidates and not matched_any:
        raise d8_triage.UniprotError(
            f"{sp['source_id']}: none of {len(candidates)} candidates matched a UniProt entry"
        )
    triage, no_entry = [], 0
    for r in candidates:
        mine = d8_triage.entries_for_gene(r["gene_id"], mapping, entries)
        no_entry += not mine
        key = (sp["source_id"], r["gene_id"])
        d8_class, reason = d8_triage.classify_pm(
            mine, key in literature, literature.get(key, False)
        )
        triage.append(
            {
                "source_id": r["source_id"],
                "gene_id": r["gene_id"],
                "symbol": r["symbol"],
                "d8_class": d8_class,
                "d8_reason": reason,
                "uniprot_accessions": ",".join(e.accession for e in mine),
                "reviewed": ",".join("yes" if e.reviewed else "no" for e in mine),
                "gpi_eco": ";".join(",".join(e.gpi_eco) for e in mine),
                "tm_count": ";".join(str(e.tm_count) for e in mine),
                "tm_eco": ";".join(",".join(e.tm_eco) for e in mine),
            }
        )
    organism, release2 = fetch(d8_triage.organism_gpi_url(sp["taxon_id"]), f"{sp['source_id']}_gpi")
    curated = [e for e in (d8_triage.parse_entry(x) for x in organism) if e.curated_gpi]
    outside = []
    for r in rows:
        if r["label"] == "P-ext":
            continue
        for e in d8_triage.entries_for_gene(r["gene_id"], mapping, curated):
            outside.append(
                {
                    "source_id": r["source_id"],
                    "gene_id": r["gene_id"],
                    "symbol": r["symbol"],
                    "label": r["label"],
                    "uniprot_accession": e.accession,
                    "gpi_eco": ",".join(e.gpi_eco),
                }
            )
    counts = {
        "source_id": sp["source_id"],
        "pm_candidates": str(len(candidates)),
        "p_gpi": str(sum(t["d8_class"] == d8_triage.P_GPI for t in triage)),
        "pm_tm": str(sum(t["d8_class"] == d8_triage.PM_TM for t in triage)),
        "pm_unresolved": str(sum(t["d8_class"] == d8_triage.PM_UNRESOLVED for t in triage)),
        "no_uniprot_entry": str(no_entry),
        "gpi_feature_no_evidence": str(
            len({e.accession for e in entries if e.gpi_features_without_eco})
        ),
        "organism_curated_gpi_entries": str(len(curated)),
        "curated_gpi_outside_pext": str(len(outside)),
        "uniprot_release": release or release2,
    }
    return triage, outside, counts


def apply_triage(truth_rows, triage_rows):
    """Return truth rows with d8_class and d8_reason; stratum becomes the D8 class."""
    by_gene = {(t["source_id"], t["gene_id"]): t for t in triage_rows}
    out = []
    for r in truth_rows:
        t = by_gene.get((r["source_id"], r["gene_id"]))
        new = {
            **r,
            "d8_class": t["d8_class"] if t else "",
            "d8_reason": t["d8_reason"] if t else "",
        }
        if t:
            new["stratum"] = t["d8_class"]
        out.append(new)
    return out


def main(argv=None, fetch=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--curated-gpi", default=str(paths.STEP1_DIR / "curated_gpi.tsv"))
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--sources", nargs="*", default=None)
    parser.add_argument(
        "--allow-partial-truth-set",
        action="store_true",
        help="use a truth_set.tsv.gz whose extract_log.json does not say all_sources: true",
    )
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    if fetch is None:

        def fetch(url, tag):
            return uniprot_fetch(url, work / "d8_uniprot", tag)

    try:
        species_rows = truth_table.read_tsv(args.species)
        all_source_ids = [r["source_id"] for r in species_rows]
        if args.sources is not None:
            for source_id in args.sources:
                if source_id not in all_source_ids:
                    raise d8_triage.UniprotError(
                        f"unknown source {source_id}; valid ids: {', '.join(all_source_ids)}"
                    )
            species_rows = [r for r in species_rows if r["source_id"] in args.sources]
        if not species_rows:
            raise d8_triage.UniprotError("no sources selected")
        runinfo.require_full(
            work / "extract_log.json", "truth_set.tsv.gz", args.allow_partial_truth_set
        )
        truth_path = work / "truth_set.tsv.gz"
        truth_rows = truth_table.read_tsv(truth_path)
        curated_rows = d8_triage.read_curated_gpi(args.curated_gpi)
        literature = {
            (r["source_id"], r["gene_id"]): r["override_tm"] == "yes" for r in curated_rows
        }
        unmatched = d8_triage.check_curated_gpi(curated_rows, truth_rows)
        no_gene = [u for u in unmatched if u["reason"] == "no_truth_gene"]
        if no_gene:
            names = ", ".join(f"{u['source_id']} {u['gene_id']}" for u in no_gene[:10])
            raise d8_triage.CuratedGpiError(f"curated_gpi.tsv rows match no truth gene: {names}")
        triage, outside, counts = [], [], []
        for sp in species_rows:
            rows = [r for r in truth_rows if r["source_id"] == sp["source_id"]]
            if not any(r["pm_candidate"] == "yes" for r in rows):
                continue
            t, o, c = triage_source(sp, rows, literature, fetch)
            triage += t
            outside += o
            counts.append(c)
        triaged = apply_triage(truth_rows, triage)
        conflicts = d8_triage.tm_conflicts(triage, literature)
        selected = [r["source_id"] for r in species_rows]
        run_log = {
            "sources": selected,
            "all_sources": selected == all_source_ids,
            "truth_set_sha256": manifest.sha256_file(truth_path),
            "uniprot_release": {c["source_id"]: c["uniprot_release"] for c in counts},
            "git_commit": runinfo.git_commit(),
            "python": runinfo.python_version(),
            "arguments": list(argv) if argv is not None else sys.argv[1:],
        }
        columns = (*truth_table.TRUTH_COLUMNS, "d8_class", "d8_reason")
        tables = {
            "d8_triage.tsv": (TRIAGE_COLUMNS, triage),
            "d8_gpi_outside_pext.tsv": (OUTSIDE_COLUMNS, outside),
            "d8_curated_conflicts.tsv": (d8_triage.CONFLICT_COLUMNS, conflicts),
            "curated_gpi_unmatched.tsv": (d8_triage.UNMATCHED_COLUMNS, unmatched),
            "d8_counts.tsv": (D8_COUNT_COLUMNS, counts),
            "truth_set_triaged.tsv.gz": (columns, triaged),
        }
        writers = {
            name: (lambda p, c=c, r=r: truth_table.write_tsv(p, c, r))
            for name, (c, r) in tables.items()
        }
        writers["d8_run.json"] = lambda p: runinfo.write_json(p, run_log)
        runinfo.atomic_write_all(work, writers)
    except (
        d8_triage.UniprotError,
        urllib.error.URLError,
        http.client.HTTPException,
        OSError,
        ValueError,
        KeyError,
    ) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for name, rows_for_review in (
        ("d8_curated_conflicts.tsv", conflicts),
        ("curated_gpi_unmatched.tsv", unmatched),
    ):
        if rows_for_review:
            print(
                f"NOTE: {len(rows_for_review)} row(s) need review; see {name}",
                file=sys.stderr,
            )
    for c in counts:
        print("\t".join(f"{k}={c[k]}" for k in D8_COUNT_COLUMNS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
