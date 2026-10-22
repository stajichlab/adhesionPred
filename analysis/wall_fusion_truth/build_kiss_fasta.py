#!/usr/bin/env python3
"""Build the FASTA of the Kiss control sequences for Task 5.3. LOCAL ONLY, never run on HPCC.

DRAFT, not run. The output holds JGI protein sequences (for the rows of
`controls_local_sequence.tsv`). JGI data are never committed, so:
  * --out-dir is required and must be OUTSIDE this repository (the script refuses a path inside it);
  * every sequence is checked against the sha256 that the committed control set holds
    (`sequence_sha256`); a mismatch stops the script and nothing is written;
  * every JGI file is checked against the sha256 of the committed file table
    (`analysis/wall_fusion_truth/results/jgi_files.tsv`) before it is read.
The JGI files are untrusted data: they are only read (gzip, text), never executed.

Written to --out-dir:
  kiss_sequences.faa   one FASTA: the control sequences (header = accession; `jgi:<portal>:<id>`
                       for a local-key row), plus the negatives and the family representatives only
                       when asked for (a protein in two classes is an error unless --on-overlap keep-both)
  members.tsv          accession, class, slice, source, cluster_id (Kiss proxy), jgi_accept
  *.sha256             the sha256 of each file, to be checked after the copy to HPCC

Negative sources are never defaulted: the owner has not yet confirmed which two are meant (spec
section 9 item 12). Give each one with --negative-set DIR (a control-set directory whose
controls.tsv has label=negative rows and whose sequences.faa holds their sequences).

Usage:
  python3 -I analysis/wall_fusion_truth/build_kiss_fasta.py --jgi-dir to_import/kiss2019_jgi \\
      --out-dir ~/kiss_5_3_inputs [--negative-set data/controls/hard-negatives-shared] [--include-family-reps]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONTROL_DIR = ROOT / "data/controls/wall-fusion-kiss2019"
JGI_TABLE = HERE / "results/jgi_files.tsv"


def load_map_kiss():
    spec = importlib.util.spec_from_file_location("map_kiss_for_fasta", HERE / "map_kiss.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["map_kiss_for_fasta"] = mod
    spec.loader.exec_module(mod)
    return mod


def refuse_repo_path(path: Path, root: Path = ROOT) -> Path:
    """Return the resolved path; refuse a path inside the repository."""
    p = path.expanduser().resolve()
    if p == root.resolve() or p.is_relative_to(root.resolve()):
        raise SystemExit(f"refusing --out-dir {p}: it is inside the repository ({root})")
    return p


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    cur = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif cur is not None:
            out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


class OverlapError(ValueError):
    """A protein is in two classes (for example a Kiss wall control that is also a hard negative)."""


def gather(
    control_dir: Path,
    jgi_dir: Path,
    jgi_table: Path,
    negative_dirs: list[Path],
    include_reps: bool,
    on_overlap: str = "error",
    mk=None,
) -> tuple[dict[str, str], list[dict]]:
    """Return ({accession: sequence}, members rows); raise on any hash mismatch.

    A negative that is also a positive (or a family representative) is an overlap. With
    on_overlap=error (the default) the call raises OverlapError and lists every overlap. With
    keep-both the negative copy is written as `<accession>#neg` (class negative, `overlap_with`
    names the other entry): the two copies fall in one cluster, which is then a mixed cluster.
    """
    mk = mk or load_map_kiss()
    seqs: dict[str, str] = {}
    members: list[dict] = []
    overlaps: list[str] = []

    def add(acc, seq, cls, sl, source, cluster, accept):
        overlap = ""
        if acc in seqs:
            overlaps.append(acc)
            if on_overlap != "keep-both" or cls != "negative":
                return
            overlap, acc = acc, f"{acc}#neg"
        seqs[acc] = seq
        members.append(
            {
                "accession": acc,
                "class": cls,
                "slice": sl,
                "source": source,
                "cluster_id": cluster,
                "jgi_accept": accept,
                "overlap_with": overlap,
            }
        )

    public = read_fasta(control_dir / "sequences.faa")
    for r in read_tsv(control_dir / "controls.tsv"):
        seq = public[r["accession"]]
        if sha(seq) != r["sequence_sha256"]:
            raise ValueError(f"{r['accession']}: sequences.faa differs from sequence_sha256")
        add(
            r["accession"],
            seq,
            "positive",
            r["stratum"],
            "controls.tsv",
            r["cluster_id"],
            r["jgi_accept"],
        )

    local = read_tsv(control_dir / "controls_local_sequence.tsv")
    if local:
        table = mk.read_jgi_table(jgi_table)
        wanted: dict[str, set[str]] = {}
        target = {}
        for r in local:
            first = r["protein_id_jgi"].split("+")[0]
            pid, code = mk.jgi_protein_id(first), first.split("_", 1)[1]
            portal = mk.jgi_portal(code)
            if r["accession"].startswith("jgi:") and r["accession"] != f"jgi:{portal}:{pid}":
                raise ValueError(f"{r['accession']}: does not match {first}")
            wanted.setdefault(portal, set()).add(pid)
            target[r["accession"]] = (portal, pid)
        loaded = mk.load_jgi_proteins(jgi_dir, table, wanted)  # checks each file's sha256
        for r in local:
            portal, pid = target[r["accession"]]
            if pid not in loaded.get(portal, {}):
                raise ValueError(f"{r['accession']}: protein {pid} not found in the {portal} file")
            seq = loaded[portal][pid][1]
            if sha(seq) != r["sequence_sha256"]:
                raise ValueError(f"{r['accession']}: JGI sequence differs from sequence_sha256")
            add(
                r["accession"],
                seq,
                "positive",
                r["stratum"],
                "controls_local_sequence.tsv",
                r["cluster_id"],
                r["jgi_accept"],
            )

    if include_reps:
        fam = read_fasta(control_dir / "family_rep_sequences.faa")
        for r in read_tsv(control_dir / "family_reps.tsv"):
            if sha(fam[r["accession"]]) != r["sequence_sha256"]:
                raise ValueError(f"{r['accession']}: family sequence differs from sequence_sha256")
            if r["accession"] in seqs:
                continue  # the same protein is already a control; one sequence, one class
            add(
                r["accession"],
                fam[r["accession"]],
                "family_rep",
                r["slice"],
                "family_reps.tsv",
                f"KISS-{r['slice']}-{r['kiss_cluster_number']}",
                "",
            )

    for d in negative_dirs:
        neg = read_fasta(d / "sequences.faa")
        for r in read_tsv(d / "controls.tsv"):
            if r.get("label") != "negative":
                continue
            add(
                r["accession"],
                neg[r["accession"]],
                "negative",
                "",
                d.name,
                r.get("cluster_id", ""),
                "",
            )
    if overlaps and on_overlap == "error":
        raise OverlapError(
            f"{len(overlaps)} sequences are in two classes: {', '.join(sorted(set(overlaps)))}. "
            "Choose --on-overlap keep-both (the negative copy is named <accession>#neg) "
            "or remove them from the negative set."
        )
    return seqs, members


def write_outputs(out: Path, seqs: dict[str, str], members: list[dict]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fa = out / "kiss_sequences.faa"
    fa.write_text("".join(f">{a}\n{seqs[a]}\n" for a in sorted(seqs)))
    mp = out / "members.tsv"
    cols = ["accession", "class", "slice", "source", "cluster_id", "jgi_accept", "overlap_with"]
    with mp.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for m in sorted(members, key=lambda x: x["accession"]):
            w.writerow([m[c] for c in cols])
    for p in (fa, mp):
        (out / f"{p.name}.sha256").write_text(hashlib.sha256(p.read_bytes()).hexdigest() + "\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jgi-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--control-dir", type=Path, default=CONTROL_DIR)
    ap.add_argument("--jgi-table", type=Path, default=JGI_TABLE)
    ap.add_argument("--negative-set", type=Path, action="append", default=[])
    ap.add_argument("--include-family-reps", action="store_true")
    ap.add_argument("--on-overlap", choices=("error", "keep-both"), default="error")
    a = ap.parse_args(argv)
    out = refuse_repo_path(a.out_dir)  # before anything is read
    seqs, members = gather(
        a.control_dir, a.jgi_dir, a.jgi_table, a.negative_set, a.include_family_reps, a.on_overlap
    )
    write_outputs(out, seqs, members)
    by = {}
    for m in members:
        by[m["class"]] = by.get(m["class"], 0) + 1
    print(f"{len(seqs)} sequences written to {out}: {by}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
