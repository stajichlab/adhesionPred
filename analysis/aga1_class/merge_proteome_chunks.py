"""Merge the per-chunk outputs of proteome_array.sbatch into one raw file per tool (stdlib only).

Reads $RESULTS/{signalp,predgpi,aiupred}/<chunk>/ for every chunk of chunks.tsv (staged by stage_proteomes.py), checks that every
chunk directory and its provenance.txt exist, that all provenance files name the same pinned commit and the expected tool, and that
the chunk sha256 in provenance matches chunks.tsv, then concatenates in chunks.tsv order into a NEW directory:
  signalp/prediction_results.txt (the two '#' header lines once), predgpi/predgpi_raw.gff3, aiupred/disorder_raw.txt (the banner
  once; records start with '#>ID'), merge_provenance.txt (commit, per-tool chunk counts, sha256 of every merged file).
The merged AIUPred file carries residue letters and must stay outside the repository.

  python3 -I analysis/aga1_class/merge_proteome_chunks.py --chunks CHUNKS_TSV --results DIR --out NEW_DIR

Refuses an existing --out. FATAL: and exit 1 on a missing chunk or provenance, a commit or tool mismatch, or a chunk hash mismatch.
"""

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

TOOLS = {
    "signalp": ("prediction_results.txt", "signalp"),
    "predgpi": ("predgpi_raw.gff3", "predgpi"),
    "aiupred": ("disorder_raw.txt", "aiupred"),
}


class MergeError(Exception):
    pass


def provenance(path):
    out = {}
    for line in Path(path).read_text().splitlines():
        m = re.match(r"(commit|tool|chunk sha256|chunk): (.+)|(chunk sha256) (\S+)", line)
        if m and m.group(1):
            out[m.group(1)] = m.group(2)
        elif m and m.group(3):
            out["chunk sha256"] = m.group(4)
    return out


def merge(chunks_tsv, results, out):
    out = Path(out)
    if out.exists():
        raise MergeError(f"{out} exists; nothing is overwritten")
    with open(chunks_tsv, newline="") as fh:
        chunks = list(csv.DictReader(fh, delimiter="\t"))
    if not chunks:
        raise MergeError("chunks.tsv has no chunks")
    commits, parts = set(), {}
    for tool, (raw_name, tool_name) in TOOLS.items():
        texts = []
        for c in chunks:
            d = Path(results) / tool / c["chunk"]
            if not (d / raw_name).is_file() or not (d / "provenance.txt").is_file():
                raise MergeError(f"{d}: missing {raw_name} or provenance.txt")
            pv = provenance(d / "provenance.txt")
            if pv.get("tool") != tool_name:
                raise MergeError(
                    f"{d}: provenance tool is {pv.get('tool')!r}, expected {tool_name!r}"
                )
            if pv.get("chunk sha256") != c["sha256"]:
                raise MergeError(f"{d}: provenance chunk sha256 differs from chunks.tsv")
            commits.add(re.sub(r" .*", "", pv.get("commit", "")))
            texts.append((d / raw_name).read_text())
        parts[tool] = texts
    if len(commits) != 1 or "" in commits:
        raise MergeError(f"the provenance files name several commits: {sorted(commits)}")
    merged = {}
    for tool, texts in parts.items():
        if tool == "signalp":
            head = texts[0].splitlines()[:2]
            if not head[0].startswith("# SignalP") or "Eukarya" not in head[0]:
                raise MergeError("signalp: first chunk lacks the Eukarya header line")
            body = [ln for t in texts for ln in t.splitlines() if not ln.startswith("#")]
            merged[tool] = "\n".join(head + body) + "\n"
        elif tool == "predgpi":
            merged[tool] = "".join(t if t.endswith("\n") else t + "\n" for t in texts)
        else:
            merged[tool] = "".join(t if t.endswith("\n") else t + "\n" for t in texts)
    out.mkdir(parents=True)
    lines = [f"commit: {commits.pop()}", f"chunks: {len(chunks)}"]
    for tool, text in merged.items():
        name = TOOLS[tool][0]
        (out / tool).mkdir()
        (out / tool / name).write_text(text)
        lines.append(
            f"{tool}: {len(chunks)} chunks, sha256 {hashlib.sha256(text.encode()).hexdigest()}  {tool}/{name}"
        )
    (out / "merge_provenance.txt").write_text("\n".join(lines) + "\n")
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("chunks", "results", "out"):
        ap.add_argument(f"--{k}", required=True)
    a = ap.parse_args(argv)
    try:
        lines = merge(a.chunks, a.results, a.out)
    except (MergeError, OSError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
