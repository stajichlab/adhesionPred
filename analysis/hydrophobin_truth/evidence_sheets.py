#!/usr/bin/python3.12
"""L8: evidence sheets for the owner (tier T5). One row per protein that has a hydrophobin call (strict Pfam or relaxed level) and no
T1/T2/T3/LP label, in the 12 proteomes. Writes evidence_sheets.tsv. The owner's choices go in owner_decisions.tsv; none is written here.
Needs blastp on PATH (module load ncbi-blast/2.14.0+). Usage: evidence_sheets.py --dir analysis/hydrophobin_truth
"""

import argparse
import csv
import gzip
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SWISSPROT = "/bigdata/operations/pkgadmin/srv/projects/db/Swissprot/2023_03/uniprot_sprot.fasta"


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def cys_gaps(seq):
    pos = [i for i, a in enumerate(seq.upper()) if a == "C"]
    return [b - a - 1 for a, b in zip(pos, pos[1:], strict=False)]


def kinds(strict, relaxed, labelled):
    if labelled or not (strict or relaxed):
        return None
    if strict and relaxed:
        return "strict_and_relaxed_unlabelled"
    return "strict_only_unlabelled" if strict else "relaxed_only_unlabelled"


def parse_blast(lines, top=3):
    out = {}
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        q, s, pid, cov, ev, bits, title = f[:7]
        out.setdefault(q, [])
        if len(out[q]) < top:
            out[q].append(f"{s} {pid}% cov{cov} {ev} {title}")
    return out


def blast(query, db, top, title=True):
    fmt = "6 qseqid sseqid pident qcovs evalue bitscore " + ("stitle" if title else "sseqid")
    res = subprocess.run(
        [
            "blastp",
            "-query",
            str(query),
            "-db",
            str(db),
            "-evalue",
            "1e-3",
            "-max_target_seqs",
            "10",
            "-outfmt",
            fmt,
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    lines = sorted(res.splitlines(), key=lambda x: (x.split("\t")[0], -float(x.split("\t")[5])))
    return parse_blast(lines, top)


def read_module(path):
    with gzip.open(path, "rt") as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    pd = load("prefreeze_data")
    runs = list(csv.DictReader(open(d / "run_list.tsv"), delimiter="\t"))
    pmap = list(csv.DictReader(open(d / "proteome_map.tsv"), delimiter="\t"))
    lp = {
        k.split("|")[1]: v for k, v in pd.read_fasta(d / "literature/jensen_sequences.faa").items()
    }
    rows, qseqs = [], {}
    for r in runs:
        name, w = r["name"], Path(r["workdir"])
        seqs = pd.read_fasta(r["fasta"])
        headers = {}
        for line in open(r["fasta"]):
            if line.startswith(">"):
                headers[line[1:].split()[0]] = line[1:].strip()[:110]
        rel = read_module(w / "modules/hydrophobin_relaxed.tsv.gz")
        strict = read_module(w / "modules/pfam_hydrophobin.tsv.gz")
        r0 = read_module(w / "modules/step1_rule@R0.tsv.gz")
        tm = read_module(w / "modules/tm.tsv.gz")
        cys = read_module(w / "modules/cys8_pattern.tsv.gz")
        labelled = {x["proteome_id"] for x in pmap if x["proteome"] == name} | pd.exact_matches(
            lp, seqs
        )
        for i in seqs:
            is_s = strict.get(i, {}).get("hit") == "1"
            is_r = rel.get(i, {}).get("hit") == "1"
            k = kinds(is_s, is_r, i in labelled)
            if not k:
                continue
            key = f"{name}|{i}"
            qseqs[key] = seqs[i]
            rows.append(
                {
                    "key": key,
                    "proteome": name,
                    "id": i,
                    "kind": k,
                    "header": headers.get(i, ""),
                    "length": len(seqs[i]),
                    "n_cys": seqs[i].upper().count("C"),
                    "cys_gaps": ",".join(map(str, cys_gaps(seqs[i]))),
                    "r0": r0.get(i, {}).get("call", ""),
                    "n_tm_mature": tm.get(i, {}).get("n_tm_mature", ""),
                    "strict_families": strict.get(i, {}).get("families", ""),
                    "relaxed_score": rel.get(i, {}).get("score", ""),
                    "relaxed_model": rel.get(i, {}).get("query", ""),
                    "cys8_spacing_class": cys.get(i, {}).get("spacing_class", ""),
                }
            )
    with tempfile.TemporaryDirectory() as tmp:
        q = Path(tmp) / "q.faa"
        with open(q, "w") as f:
            for k, s in qseqs.items():
                f.write(f">{k}\n{s}\n")
        known = Path(tmp) / "known"
        subprocess.run(
            [
                "makeblastdb",
                "-in",
                str(d / "positives_t2t3.faa"),
                "-dbtype",
                "prot",
                "-out",
                str(known),
            ],
            check=True,
            capture_output=True,
        )
        hk = blast(q, known, 2, title=False)
        sp = blast(q, SWISSPROT, 3)
    for r in rows:
        r["blast_known_hydrophobins"] = " | ".join(hk.get(r["key"], [])) or "no hit at 1e-3"
        r["blast_swissprot"] = " | ".join(sp.get(r["key"], [])) or "no hit at 1e-3"
        r["owner_decision"] = ""
        r["owner_reason"] = ""
    keys = [k for k in rows[0] if k != "key"]
    with open(d / "evidence_sheets.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter

    print(len(rows), "sheet rows;", dict(Counter(r["kind"] for r in rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
