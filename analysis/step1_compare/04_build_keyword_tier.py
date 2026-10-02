#!/usr/bin/env python3
"""D10: build keyword_tier.tsv.gz (T-c rows for V-kw) and keyword_tier_removed.tsv.

Reads data/curated/surface/surface.tsv, data/curated/adhesins/eurotiomycetes_seeds.tsv,
species.tsv, $STEP1_WORKDIR/truth_sequences.tsv.gz and (for its SHA-256 only)
$STEP1_WORKDIR/truth_set.tsv.gz. Sequences for surface.tsv and the seed
accessions come from $STEP1_WORKDIR/keyword_sequences.fasta.gz; --fetch downloads that file
from UniProtKB REST if it is missing and records its SHA-256 and UniProt release in
keyword_sequences.json. A cached sequence file without a matching sidecar is refused.

Errors print `STOP: <message>` and exit 2 with no output files written. Outputs (keyword_tier.tsv.gz,
keyword_tier_removed.tsv, keyword_tier_run.json) are written to temp names and moved into place
only when all are complete. The script refuses a truth set or sequence table whose run log
(extract_log.json, sequence_run.json) is missing or does not say `all_sources: true`, unless
--allow-partial-truth-set is given. The script also stops (re-run 02) when the
truth_set_sha256 in sequence_run.json differs from the SHA-256 of the current truth_set.tsv.gz;
--allow-partial-truth-set does not disable this check. keyword_tier_run.json records
all_sources, truth_set_sha256, the SHA-256 of the inputs, the UniProt release, the git commit,
the Python version and the arguments (no timestamps).

Sequence checks (a missing sequence could hide a test protein):
- The sidecar keeps the sorted list of requested accessions and its SHA-256. A cached file is
  refused if the accession list for the current surface.tsv and seeds file differs; rerun with
  --fetch after deleting the cache.
- The run STOPs if any requested accession has no sequence (names up to 20), unless
  --allow-missing-sequences is given. Missing T-c rows are then removed as `no_sequence` and
  appear in keyword_tier_removed.tsv.
- The run STOPs if a literature-seed accession has no sequence, unless
  --allow-missing-seed-sequences is given (this disables the hash rule for those seeds). The run
  JSON lists every seed with `in_keyword_tier` and `has_sequence`.

Limitation: removal is by accession and by exact cleaned-sequence hash only. Near-identical
orthologs or paralogs under other accessions stay in the training table. Homology clustering in
the later dataset plan must remove them.
"""

import argparse
import datetime
import gzip
import hashlib
import http.client
import json
import os
import sys
import urllib.error
import zlib
from pathlib import Path

import keyword_tier
import manifest
import paths
import runinfo
import sequences
import truth_table

NET_ERRORS = (urllib.error.URLError, OSError, http.client.HTTPException)
OUTPUT_NAMES = ("keyword_tier.tsv.gz", "keyword_tier_removed.tsv", "keyword_tier_run.json")


def accessions_digest(accessions) -> str:
    """SHA-256 of the sorted, newline-joined accession list."""
    return hashlib.sha256("\n".join(sorted(accessions)).encode("ascii")).hexdigest()


def fetch_sequences(accessions: list[str], dest: Path, opener=None) -> dict:
    """Download the sequences and write a JSON sidecar with SHA-256 and UniProt release."""
    kwargs = {} if opener is None else {"opener": opener}
    part = dest.with_name(dest.name + ".part")
    meta_part = sidecar(dest).with_name(sidecar(dest).name + ".part")
    release = ""
    try:
        with (
            open(part, "wb") as raw,
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz,
        ):
            pages = (
                page
                for url in keyword_tier.accession_fasta_urls(accessions)
                for page in manifest.uniprot_pages(url, **kwargs)
            )
            try:
                for body, headers in pages:
                    gz.write(body)
                    release = release or headers.get("x-uniprot-release", "")
            except NET_ERRORS as exc:
                raise manifest.DownloadError(f"UniProt download failed: {exc!r}") from exc
        manifest.check_payload(part)
        part.replace(dest)
        returned = {
            key
            for header, _ in sequences.read_fasta(dest)
            if (key := sequences.fasta_key(header, "uniprot"))
        }
        meta = {
            "sha256": manifest.sha256_file(dest),
            "uniprot_release": release,
            "accessions": sorted(accessions),
            "accessions_sha256": accessions_digest(accessions),
            "requested": len(accessions),
            "returned": len(returned & set(accessions)),
            "fetched_on": datetime.datetime.now(datetime.UTC).date().isoformat(),
        }
        meta_part.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
        os.replace(meta_part, sidecar(dest))
    finally:
        part.unlink(missing_ok=True)
        meta_part.unlink(missing_ok=True)
    return meta


def sidecar(seq_path: Path) -> Path:
    return seq_path.with_name("keyword_sequences.json")


def check_cached(seq_path: Path, wanted=None) -> dict:
    """Return the sidecar of a cached sequence file. Refuse a file without a matching sidecar.

    If `wanted` (the accessions needed now) is given, also refuse a cache requested for a
    different accession list."""
    meta_path = sidecar(seq_path)
    if not meta_path.exists():
        raise manifest.DownloadError(f"{meta_path} missing; delete {seq_path} and rerun --fetch")
    meta = json.loads(meta_path.read_text())
    if not meta.get("uniprot_release") or meta.get("sha256") != manifest.sha256_file(seq_path):
        raise manifest.DownloadError(
            f"{seq_path}: SHA-256 or UniProt release does not match {meta_path.name}"
        )
    if wanted is not None:
        if isinstance(meta.get("accessions"), str) or "accessions_sha256" not in meta:
            raise manifest.DownloadError(
                f"{meta_path.name} is in the old format (no accession list); delete {seq_path} "
                "and rerun with --fetch"
            )
        if meta["accessions_sha256"] != accessions_digest(wanted):
            raise manifest.DownloadError(
                f"{seq_path.name} was fetched for a different accession list than the current "
                "surface.tsv and seeds file; delete it and rerun with --fetch"
            )
    return meta


def heldout_sets(species_rows, seq_rows):
    heldout_sources = {r["source_id"]: r for r in species_rows if r["role"] != "train"}
    accessions, hashes = {}, {}
    for r in seq_rows:
        sp = heldout_sources.get(r["source_id"])
        if sp is None or r["label"] == "unlabelled":
            continue
        label = f"{r['source_id']}:{r['gene_id']}"
        if sp["id_mapping"] == "uniprot":
            accessions[r["gene_id"]] = label
        hashes.setdefault(r["seq_sha256"], label)
    return accessions, hashes


def _require_current_truth_set(work: Path) -> str:
    """Stop unless sequence_run.json records the SHA-256 of the current truth_set.tsv.gz.

    Return that SHA-256. A difference means that 01 ran again after 02, so
    truth_sequences.tsv.gz is stale."""
    run_path = work / "sequence_run.json"
    try:
        recorded = json.loads(run_path.read_text()).get("truth_set_sha256")
    except (OSError, ValueError) as exc:
        raise manifest.DownloadError(
            f"cannot read {run_path} ({exc}); re-run 02_attach_sequences.py"
        ) from exc
    current = manifest.sha256_file(work / "truth_set.tsv.gz")
    if recorded != current:
        raise manifest.DownloadError(
            f"truth_sequences.tsv.gz is stale: {run_path.name} records truth_set_sha256 "
            f"{recorded!r}, but the current truth_set.tsv.gz has {current}; "
            "re-run 02_attach_sequences.py"
        )
    return current


def main(argv=None) -> int:
    root = paths.repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface", default=str(root / "data/curated/surface/surface.tsv"))
    parser.add_argument(
        "--seeds", default=str(root / "data/curated/adhesins/eurotiomycetes_seeds.tsv")
    )
    parser.add_argument("--species", default=str(paths.STEP1_DIR / "species.tsv"))
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument(
        "--allow-partial-truth-set",
        action="store_true",
        help="use truth tables whose run logs do not say all_sources: true",
    )
    parser.add_argument(
        "--allow-missing-sequences",
        action="store_true",
        help="continue when requested accessions have no sequence (T-c rows become no_sequence)",
    )
    parser.add_argument(
        "--allow-missing-seed-sequences",
        action="store_true",
        help="continue when a literature seed has no sequence; disables the hash rule for it",
    )
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()

    try:
        surface_rows = truth_table.read_tsv(args.surface)
        seeds = keyword_tier.read_seed_accessions(args.seeds)
        seq_path = work / "keyword_sequences.fasta.gz"
        wanted = sorted({r["accession"] for r in surface_rows} | set(seeds))
        if seq_path.exists():
            meta = check_cached(seq_path, wanted)
        elif args.fetch:
            meta = fetch_sequences(wanted, seq_path)
        else:
            raise manifest.DownloadError(f"{seq_path} missing; rerun with --fetch")
        allow = args.allow_partial_truth_set
        runinfo.require_full(work / "extract_log.json", "truth_set.tsv.gz", allow)
        runinfo.require_full(work / "sequence_run.json", "truth_sequences.tsv.gz", allow)
        truth_set_sha256 = _require_current_truth_set(work)
        all_sources = all(
            runinfo.says_all_sources(work / name)
            for name in ("extract_log.json", "sequence_run.json")
        )
        seq_by_acc = {
            sequences.fasta_key(h, "uniprot"): s
            for h, s in sequences.read_fasta(seq_path)
            if sequences.fasta_key(h, "uniprot")
        }
        missing = [a for a in wanted if a not in seq_by_acc]
        if missing and not args.allow_missing_sequences:
            raise manifest.DownloadError(
                f"{len(missing)} of {len(wanted)} requested accessions have no sequence "
                f"(first {min(20, len(missing))}: {', '.join(missing[:20])}); "
                "use --allow-missing-sequences to continue"
            )
        missing_seeds = [a for a in seeds if a not in seq_by_acc]
        if missing_seeds and not args.allow_missing_seed_sequences:
            raise manifest.DownloadError(
                f"literature seeds without a sequence: {', '.join(missing_seeds[:20])}; the "
                "hash rule needs them; use --allow-missing-seed-sequences to disable it for them"
            )
        tc_accessions = {
            r["accession"] for r in surface_rows if r["source"] == keyword_tier.TC_SOURCE
        }
        truth_seq_path = work / "truth_sequences.tsv.gz"
        species_rows = truth_table.read_tsv(args.species)
        accessions, hashes = heldout_sets(species_rows, truth_table.read_tsv(truth_seq_path))
        kept, removed = keyword_tier.build_keyword_tier(
            surface_rows, seq_by_acc, seeds, accessions, hashes
        )
        run_log = {
            "all_sources": all_sources,
            "truth_set_sha256": truth_set_sha256,
            "inputs_sha256": {
                "surface": manifest.sha256_file(args.surface),
                "seeds": manifest.sha256_file(args.seeds),
                "species": manifest.sha256_file(args.species),
                "truth_sequences": manifest.sha256_file(truth_seq_path),
                "keyword_sequences": meta["sha256"],
            },
            "uniprot_release": meta["uniprot_release"],
            "missing_sequences": missing,
            "seeds": {
                a: {
                    "gene": g,
                    "in_keyword_tier": a in tc_accessions,
                    "has_sequence": a in seq_by_acc,
                }
                for a, g in sorted(seeds.items())
            },
            "kept": len(kept),
            "removed": len(removed),
            "git_commit": runinfo.git_commit(),
            "python": runinfo.python_version(),
            "arguments": list(argv) if argv is not None else sys.argv[1:],
        }
        runinfo.atomic_write_all(
            work,
            {
                "keyword_tier.tsv.gz": lambda p: truth_table.write_tsv(
                    p, keyword_tier.KEYWORD_COLUMNS, kept
                ),
                "keyword_tier_removed.tsv": lambda p: truth_table.write_tsv(
                    p, keyword_tier.REMOVED_COLUMNS, removed
                ),
                "keyword_tier_run.json": lambda p: runinfo.write_json(p, run_log),
            },
        )
    except (
        manifest.DownloadError,
        sequences.MappingError,
        ValueError,
        OSError,
        KeyError,
        EOFError,
        zlib.error,
        http.client.HTTPException,
    ) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    reasons: dict[str, int] = {}
    for r in removed:
        reasons[r["reason"]] = reasons.get(r["reason"], 0) + 1
    print(
        f"kept={len(kept)} removed={len(removed)} {reasons} "
        f"sequences_sha256={meta['sha256'][:12]} uniprot_release={meta['uniprot_release']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
