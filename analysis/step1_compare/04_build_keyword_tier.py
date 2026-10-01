#!/usr/bin/env python3
"""D10: build keyword_tier.tsv.gz (T-c rows for V-kw) and keyword_tier_removed.tsv.

Reads data/curated/surface/surface.tsv, data/curated/adhesins/eurotiomycetes_seeds.tsv,
species.tsv and $STEP1_WORKDIR/truth_sequences.tsv.gz. Sequences for surface.tsv and the seed
accessions come from $STEP1_WORKDIR/keyword_sequences.fasta.gz; --fetch downloads that file
from UniProtKB REST if it is missing and records its SHA-256 and UniProt release in
keyword_sequences.json. A cached sequence file without a matching sidecar is refused.

Errors print `STOP: <message>` and exit 2 with no output files written. Outputs (keyword_tier.tsv.gz,
keyword_tier_removed.tsv, keyword_tier_run.json) are written to temp names and moved into place
only when all are complete. The script refuses a truth set or sequence table whose run log
(extract_log.json, sequence_run.json) is missing or does not say `all_sources: true`, unless
--allow-partial-truth-set is given. keyword_tier_run.json records the SHA-256 of the inputs, the
UniProt release, the git commit, the Python version and the arguments (no timestamps).
"""

import argparse
import datetime
import gzip
import http.client
import json
import os
import platform
import subprocess
import sys
import urllib.error
import zlib
from pathlib import Path

import keyword_tier
import manifest
import paths
import sequences
import truth_table

NET_ERRORS = (urllib.error.URLError, OSError, http.client.HTTPException)
OUTPUT_NAMES = ("keyword_tier.tsv.gz", "keyword_tier_removed.tsv", "keyword_tier_run.json")


def fetch_sequences(accessions: list[str], dest: Path, opener=None) -> dict[str, str]:
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
        meta = {
            "sha256": manifest.sha256_file(dest),
            "uniprot_release": release,
            "accessions": str(len(accessions)),
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


def check_cached(seq_path: Path) -> dict[str, str]:
    """Return the sidecar of a cached sequence file. Refuse a file without a matching sidecar."""
    meta_path = sidecar(seq_path)
    if not meta_path.exists():
        raise manifest.DownloadError(f"{meta_path} missing; delete {seq_path} and rerun --fetch")
    meta = json.loads(meta_path.read_text())
    if not meta.get("uniprot_release") or meta.get("sha256") != manifest.sha256_file(seq_path):
        raise manifest.DownloadError(
            f"{seq_path}: SHA-256 or UniProt release does not match {meta_path.name}"
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


def git_commit() -> str:
    """HEAD commit of the repository, or 'unknown'."""
    try:
        done = subprocess.run(
            ["git", "-C", str(paths.repo_root()), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return done.stdout.strip() or "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _require_full(log_path: Path, what: str) -> None:
    """Stop unless the run log says that the table covers all sources."""
    try:
        all_sources = json.loads(log_path.read_text()).get("all_sources")
    except (OSError, ValueError) as exc:
        raise manifest.DownloadError(
            f"cannot read {log_path} ({exc}); use --allow-partial-truth-set to skip this check"
        ) from exc
    if all_sources is not True:
        raise manifest.DownloadError(
            f"{log_path.name} does not say all_sources: true (value: {all_sources!r}), so {what} "
            "may hold only some sources; rerun without --sources or use --allow-partial-truth-set"
        )


def _write_outputs(work: Path, kept, removed, run_log) -> None:
    """Write to temp names in the output directory, then os.replace all of them."""
    work.mkdir(parents=True, exist_ok=True)
    temps = {name: work / f".tmp.{name}" for name in OUTPUT_NAMES}
    try:
        truth_table.write_tsv(temps["keyword_tier.tsv.gz"], keyword_tier.KEYWORD_COLUMNS, kept)
        truth_table.write_tsv(
            temps["keyword_tier_removed.tsv"], keyword_tier.REMOVED_COLUMNS, removed
        )
        temps["keyword_tier_run.json"].write_text(
            json.dumps(run_log, indent=2, sort_keys=True) + "\n"
        )
        for name in OUTPUT_NAMES:
            os.replace(temps[name], work / name)
    finally:
        for temp in temps.values():
            temp.unlink(missing_ok=True)


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
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()

    try:
        surface_rows = truth_table.read_tsv(args.surface)
        seeds = keyword_tier.read_seed_accessions(args.seeds)
        seq_path = work / "keyword_sequences.fasta.gz"
        if seq_path.exists():
            meta = check_cached(seq_path)
        elif args.fetch:
            wanted = sorted({r["accession"] for r in surface_rows} | set(seeds))
            meta = fetch_sequences(wanted, seq_path)
        else:
            raise manifest.DownloadError(f"{seq_path} missing; rerun with --fetch")
        if not args.allow_partial_truth_set:
            _require_full(work / "extract_log.json", "truth_set.tsv.gz")
            _require_full(work / "sequence_run.json", "truth_sequences.tsv.gz")
        seq_by_acc = {
            sequences.fasta_key(h, "uniprot"): s
            for h, s in sequences.read_fasta(seq_path)
            if sequences.fasta_key(h, "uniprot")
        }
        truth_seq_path = work / "truth_sequences.tsv.gz"
        species_rows = truth_table.read_tsv(args.species)
        accessions, hashes = heldout_sets(species_rows, truth_table.read_tsv(truth_seq_path))
        kept, removed = keyword_tier.build_keyword_tier(
            surface_rows, seq_by_acc, seeds, accessions, hashes
        )
        run_log = {
            "inputs_sha256": {
                "surface": manifest.sha256_file(args.surface),
                "seeds": manifest.sha256_file(args.seeds),
                "species": manifest.sha256_file(args.species),
                "truth_sequences": manifest.sha256_file(truth_seq_path),
                "keyword_sequences": meta["sha256"],
            },
            "uniprot_release": meta["uniprot_release"],
            "kept": len(kept),
            "removed": len(removed),
            "git_commit": git_commit(),
            "python": platform.python_version(),
            "arguments": list(argv) if argv is not None else sys.argv[1:],
        }
        _write_outputs(work, kept, removed, run_log)
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
