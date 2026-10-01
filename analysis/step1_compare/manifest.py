"""Pinned-input manifest: read it, hash files, validate downloads. Standard library only.

Mode `strict`: a SHA-256 that differs from the manifest stops the run.
Mode `record`: the upstream file changes often (PomBase daily, UniProt per release). A different
SHA-256 is logged, not fatal. The payload checks below still apply to both modes.
"""

import csv
import gzip
import hashlib
import re
import urllib.request
import zlib
from pathlib import Path

MANIFEST_COLUMNS = (
    "file",
    "kind",
    "url",
    "sha256",
    "size",
    "last_modified",
    "date_generated",
    "mode",
    "verified",
)
HTML_PREFIXES = (b"<!doctype html", b"<html", b"<?xml", b"<head")


class DownloadError(RuntimeError):
    """A download is empty, HTML, truncated, or has a SHA-256 that differs from the manifest."""


def read_manifest(path: str | Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        missing = [c for c in MANIFEST_COLUMNS if c not in row]
        if missing:
            raise ValueError(f"{path}: manifest row {row.get('file')} lacks {missing}")
    return rows


def write_manifest(path: str | Path, rows: list[dict[str, str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(MANIFEST_COLUMNS), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_payload(path: str | Path) -> None:
    """Refuse a missing file, an empty file, an HTML error page and a truncated gzip file.

    Limits: an empty but valid gzip stream and a truncated plain-text file (for example a cut
    go-basic.obo) pass this check. For `strict` files the SHA-256 comparison catches both. For
    `record` files only the gzip end-of-stream check applies.
    """
    path = Path(path)
    if not path.exists():
        raise DownloadError(f"{path}: file not found")
    if path.stat().st_size == 0:
        raise DownloadError(f"{path.name}: empty file")
    with open(path, "rb") as handle:
        head = handle.read(512)
    if head[:2] == b"\x1f\x8b":
        try:
            with gzip.open(path, "rb") as gz:
                first = gz.read(512)
                while gz.read(1 << 20):
                    pass
        except (EOFError, OSError, zlib.error) as exc:
            raise DownloadError(f"{path.name}: truncated or corrupt gzip ({exc})") from exc
        head = first
    if head.lstrip().lower().startswith(HTML_PREFIXES):
        raise DownloadError(f"{path.name}: HTML page, not data")


def verify_against_manifest(path: str | Path, row: dict[str, str]) -> str:
    """Return the file's SHA-256. Raise DownloadError on a strict-mode mismatch."""
    check_payload(path)
    actual = sha256_file(path)
    if row["mode"] == "strict" and actual != row["sha256"]:
        raise DownloadError(
            f"{Path(path).name}: SHA-256 {actual[:12]} differs from manifest {row['sha256'][:12]}"
        )
    return actual


def http_get(url: str, opener=urllib.request.urlopen):
    """Return (body bytes, headers dict with lower-case keys)."""
    request = urllib.request.Request(url, headers={"User-Agent": "adhesionPred-step1/1"})
    with opener(request, timeout=300) as response:
        body = response.read()
        headers = {k.lower(): v for k, v in response.headers.items()}
    return body, headers


_NEXT = re.compile(r'<([^>]+)>;\s*rel="next"')


def uniprot_pages(url: str, opener=urllib.request.urlopen):
    """Yield (body, headers) for each page of a UniProt REST search, following Link rel=next."""
    while url:
        body, headers = http_get(url, opener)
        yield body, headers
        match = _NEXT.search(headers.get("link", ""))
        url = match.group(1) if match else ""
