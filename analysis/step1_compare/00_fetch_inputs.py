#!/usr/bin/env python3
"""Fetch the pinned step 1 inputs into the work directory and check them against manifest.tsv.

Environment: PROJ_ROOT (repository root) and STEP1_WORKDIR (output root). Files go to
$STEP1_WORKDIR/downloads. A strict-mode SHA-256 mismatch stops the run (exit 2) unless
--update-manifest is given. Small downloads only: this needs no SLURM job.
"""

import argparse
import datetime
import gzip
import http.client
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

import manifest
import paths

NET_ERRORS = (urllib.error.URLError, OSError, http.client.HTTPException)


class NetworkError(RuntimeError):
    """An HTTP error, a URL error or a timeout during a download."""


def fetch_one(row: dict[str, str], dest_dir: Path, opener) -> dict[str, str]:
    dest = dest_dir / row["file"]
    part = dest.with_name(dest.name + ".part")
    headers: dict[str, str] = {}
    try:
        if row["kind"] == "http":
            try:
                body, headers = manifest.http_get(row["url"], opener)
            except NET_ERRORS as exc:
                raise NetworkError(f"{row['url']}: {exc!r}") from exc
            part.write_bytes(body)
        elif row["kind"] == "uniprot_fasta":
            with (
                open(part, "wb") as raw,
                gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz,
            ):
                pages = manifest.uniprot_pages(row["url"], opener)
                while True:
                    try:
                        body, page_headers = next(pages)
                    except StopIteration:
                        break
                    except NET_ERRORS as exc:
                        raise NetworkError(f"{row['url']}: {exc!r}") from exc
                    gz.write(body)
                    headers = headers or page_headers
        else:
            raise ValueError(f"{row['file']}: unknown kind {row['kind']!r}")
    except NetworkError:
        part.unlink(missing_ok=True)
        raise
    try:
        manifest.check_payload(part)
    except manifest.DownloadError:
        part.unlink()
        raise
    sha = manifest.sha256_file(part)
    return {
        "part": str(part),
        "dest": str(dest),
        "sha256": sha,
        "size": str(part.stat().st_size),
        "last_modified": headers.get("last-modified", ""),
        "http_date": headers.get("date", ""),
        "uniprot_release": headers.get("x-uniprot-release", ""),
    }


def main(argv=None, opener=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=str(paths.STEP1_DIR / "manifest.tsv"))
    parser.add_argument("--dest", default=None, help="default: $STEP1_WORKDIR/downloads")
    parser.add_argument("--only", nargs="*", default=None, help="file names to fetch")
    parser.add_argument("--update-manifest", action="store_true")
    args = parser.parse_args(argv)

    opener = opener or urllib.request.urlopen
    dest_dir = Path(args.dest) if args.dest else paths.downloads_dir()
    dest_dir.mkdir(parents=True, exist_ok=True)
    rows = manifest.read_manifest(args.manifest)
    known = {r["file"] for r in rows}
    for name in args.only or []:
        if name not in known:
            print(f"STOP: unknown file {name}", file=sys.stderr)
            return 2
    today = datetime.datetime.now(datetime.UTC).date().isoformat()
    log_rows, failed, changed = [], [], False
    for row in rows:
        if args.only and row["file"] not in args.only:
            continue
        dest = dest_dir / row["file"]
        if (
            dest.exists()
            and row["mode"] == "strict"
            and manifest.sha256_file(dest) == row["sha256"].lower()
        ):
            print(f"ok (cached) {row['file']}")
            continue
        try:
            got = fetch_one(row, dest_dir, opener)
        except NetworkError as exc:
            print(f"STOP: {exc}", file=sys.stderr)
            return 2
        except manifest.DownloadError as exc:
            print(f"FAILED {row['url']} -> {row['file']}: {exc}", file=sys.stderr)
            failed.append(row["file"])
            continue
        # --update-manifest accepts a new strict hash: verify the payload in record mode.
        check_row = {**row, "mode": "record"} if args.update_manifest else row
        try:
            manifest.verify_against_manifest(got["part"], check_row)
        except manifest.DownloadError as exc:
            Path(got["part"]).unlink(missing_ok=True)
            print(
                f"FAILED {row['url']} -> {row['file']}: {exc}; "
                "rerun with --update-manifest to accept it",
                file=sys.stderr,
            )
            failed.append(row["file"])
            continue
        mismatch = got["sha256"] != row["sha256"].lower()
        os.replace(got["part"], got["dest"])
        if mismatch:
            print(f"changed {row['file']}: {row['sha256'][:12]} -> {got['sha256'][:12]}")
        if mismatch and args.update_manifest:
            row.update(
                sha256=got["sha256"],
                size=got["size"],
                last_modified=got["last_modified"] or "not verified",
                verified=today,
            )
            changed = True
        log_rows.append({"file": row["file"], **{k: got[k] for k in got if k not in ("part",)}})
    if changed:
        manifest.write_manifest(args.manifest, rows)
    if log_rows:
        log = dest_dir / "fetch_log.tsv"
        write_header = not log.exists()
        with open(log, "a", encoding="utf-8") as handle:
            cols = [
                "file",
                "dest",
                "sha256",
                "size",
                "last_modified",
                "http_date",
                "uniprot_release",
            ]
            if write_header:
                handle.write("\t".join(["fetched_on", *cols]) + "\n")
            for entry in log_rows:
                handle.write("\t".join([today, *(entry.get(c, "") for c in cols)]) + "\n")
    return 2 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
