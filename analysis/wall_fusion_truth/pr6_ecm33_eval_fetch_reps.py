#!/usr/bin/env python3
"""Fetch the UniRef50 cluster representatives of the Ecm33 independent-set test clusters from the public UniProt REST service (LOCAL, read-only).

Run on a laptop, never on HPCC. Stdlib only; `python3 -I`. NOT run by its author (the job stays DRAFT); one UniRef50 URL was read once to confirm the form.

Why: the 187 cluster ids of candidates.tsv were taken from the LIVE UniProt on 2026-10-23; the HPCC file is the 2025_03 release and holds only
51 of them (job 30324005). The signed freeze says the test sequence is "the cluster representative (UniRef50 FASTA record)"; the live record of
the cluster id satisfies that. The test set is therefore the live PANTHER / UniRef50 view of the fetch date, recorded with date and sha256.

What it does
  - reads the cluster ids with the SAME selection as the extract job (pr6_ecm33_eval_members.select_proteins + cluster_table);
  - GET https://rest.uniprot.org/uniref/<UniRef50_ID>.fasta (one request per id, a user-agent, a pause between requests, retries with
    backoff on 429/5xx/network errors); raw answers are cached under <out-dir>/raw (a cached answer is not fetched again);
  - FAILS LOUDLY, listing every problem and writing no staged file, when an id returns nothing, is HTTP 404, has no or several records, or a
    first header token that is a DIFFERENT id (a merged or obsolete cluster); nothing is skipped silently;
  - writes <out-dir>/reps.faa (header = the cluster id only), reps.faa.sha256 (sha256sum format) and provenance.txt (URL pattern, date,
    ids requested and received, service release from the response headers, sha256 of the candidates table).
  --out-dir is required and must be OUTSIDE this repository (the same guard as build_kiss_fasta.py); no sequence is ever committed.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import importlib.util
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CANDIDATES = ROOT / "docs/reports/data/ecm33_independent_set_candidates/candidates.tsv"
URL = "https://rest.uniprot.org/uniref/{id}.fasta"
USER_AGENT = "adhesionPred-ecm33-eval-fetch/1.0 (jason.stajich@ucr.edu; read-only, one request per cluster id)"
ID_RE = re.compile(r"^UniRef50_[A-Za-z0-9_-]+$")
SEQ_RE = re.compile(r"^[A-Za-z*]+$")
RELEASE_HEADERS = ("x-uniprot-release", "x-uniprot-release-date", "x-api-deployment-date")


def _members():
    spec = importlib.util.spec_from_file_location(
        "pr6_ecm33_eval_members", HERE / "pr6_ecm33_eval_members.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pr6_ecm33_eval_members"] = mod
    spec.loader.exec_module(mod)
    return mod


def refuse_repo_path(path: Path, root: Path = ROOT) -> Path:
    """Return the resolved path; refuse a path inside the repository."""
    p = path.expanduser().resolve()
    if p == root.resolve() or p.is_relative_to(root.resolve()):
        raise SystemExit(f"refusing --out-dir {p}: it is inside the repository ({root})")
    return p


def cluster_ids(candidates: Path) -> list[str]:
    """The test cluster ids: the selection of the extract job (reused, not re-implemented)."""
    m = _members()
    return [c["cluster_id"] for c in m.cluster_table(m.select_proteins(m.read_tsv(candidates)))]


def parse_record(cid: str, text: str) -> tuple[str | None, str | None]:
    """(sequence, problem). Exactly one FASTA record whose first header token is `cid`."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return None, "empty answer"
    heads = [ln for ln in lines if ln.startswith(">")]
    if len(heads) != 1 or not lines[0].startswith(">"):
        return None, f"{len(heads)} FASTA records (expected 1)"
    got = lines[0][1:].split()[0] if lines[0][1:].split() else ""
    if got != cid:
        return None, f"answered with a different cluster id {got!r} (merged or obsolete?)"
    seq = "".join(lines[1:])
    if not seq or not SEQ_RE.match(seq):
        return None, "empty or non-sequence body"
    return seq, None


def fetch_one(
    cid: str,
    opener=urllib.request.urlopen,
    retries: int = 5,
    backoff: float = 2.0,
    sleep=time.sleep,
) -> tuple[str | None, dict, str | None]:
    """(text, response headers (lower-case), problem). HTTP 404/410 is a problem at once; 429, 5xx and network errors are retried."""
    req = urllib.request.Request(
        URL.format(id=cid), headers={"User-Agent": USER_AGENT, "Accept": "text/plain"}
    )
    last = "no attempt"
    for attempt in range(retries):
        try:
            with opener(req, timeout=60) as r:
                hdr = {k.lower(): v for k, v in r.headers.items()}
                return r.read().decode("utf-8"), hdr, None
        except urllib.error.HTTPError as e:
            if e.code in (404, 410):
                return None, {}, f"HTTP {e.code} (not found or obsolete)"
            last = f"HTTP {e.code}"
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last = f"network error: {e}"
        sleep(backoff * (2**attempt))
    return None, {}, f"gave up after {retries} attempts ({last})"


def wrap(seq: str, n: int = 60) -> str:
    return "\n".join(seq[i : i + n] for i in range(0, len(seq), n))


def run(
    ids: list[str],
    out: Path,
    opener=urllib.request.urlopen,
    delay: float = 0.3,
    retries: int = 5,
    sleep=time.sleep,
    candidates_sha: str = "",
) -> int:
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    seqs: dict[str, str] = {}
    problems: list[tuple[str, str]] = []
    releases: dict[str, set] = {h: set() for h in RELEASE_HEADERS}
    bad = [i for i in ids if not ID_RE.match(i)]
    if bad:
        raise SystemExit(f"FATAL: ids that are not UniRef50 ids: {bad[:5]}")
    if len(set(ids)) != len(ids):
        raise SystemExit("FATAL: duplicate cluster ids in the selection")
    for n, cid in enumerate(ids):
        cache = raw / f"{cid}.fasta"
        if cache.exists():
            text, hdr, prob = cache.read_text(), {}, None
        else:
            text, hdr, prob = fetch_one(cid, opener=opener, retries=retries, sleep=sleep)
            if text is not None:
                cache.write_text(text)
            if n + 1 < len(ids):
                sleep(delay)
        if prob is None:
            seq, prob = parse_record(cid, text)
            if seq is not None:
                seqs[cid] = seq
        for h in RELEASE_HEADERS:
            if hdr.get(h):
                releases[h].add(hdr[h])
        if prob:
            problems.append((cid, prob))
    if problems:
        print(
            f"FATAL: {len(problems)} of {len(ids)} cluster ids did not give their own record; no staged file is written:",
            file=sys.stderr,
        )
        for cid, prob in problems:
            print(f"  {cid}\t{prob}", file=sys.stderr)
        return 1
    faa = out / "reps.faa"
    faa.write_text("".join(f">{cid}\n{wrap(seqs[cid])}\n" for cid in ids))
    digest = hashlib.sha256(faa.read_bytes()).hexdigest()
    (out / "reps.faa.sha256").write_text(f"{digest}  reps.faa\n")
    rel = "; ".join(
        f"{h}: {', '.join(sorted(v)) or 'not in the response headers (cached answers keep no headers)'}"
        for h, v in releases.items()
    )
    (out / "provenance.txt").write_text(
        "job: pr6_ecm33_eval_fetch_reps (UniRef50 representatives of the independent-set test clusters)\n"
        f"source: {URL} (UniProt REST, one request per id)\n"
        f"date: {datetime.date.today().isoformat()} (the test set is the live view of this date, not the 2025_03 release)\n"
        f"ids requested: {len(ids)}\nids received: {len(seqs)}\n"
        f"service: {rel}\n"
        f"candidates.tsv sha256: {candidates_sha}\n"
        f"reps.faa sha256: {digest}\n"
    )
    print(f"wrote {faa} ({len(seqs)} of {len(ids)} ids), sha256 {digest}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="OUTSIDE the repository; must be new or empty of reps.faa",
    )
    ap.add_argument("--candidates", type=Path, default=CANDIDATES)
    ap.add_argument("--delay", type=float, default=0.3, help="seconds between requests")
    ap.add_argument("--retries", type=int, default=5)
    a = ap.parse_args(argv)
    out = refuse_repo_path(a.out_dir)  # before anything is read or fetched
    if (out / "reps.faa").exists():
        raise SystemExit(f"refusing: {out / 'reps.faa'} exists; nothing is overwritten")
    ids = cluster_ids(a.candidates)
    sha = hashlib.sha256(a.candidates.read_bytes()).hexdigest()
    return run(ids, out, delay=a.delay, retries=a.retries, candidates_sha=sha)


if __name__ == "__main__":
    sys.exit(main())
