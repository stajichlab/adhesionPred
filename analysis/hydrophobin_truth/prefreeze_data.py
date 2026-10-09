"""Data assembly for the cutoff rule (pure functions plus small readers)."""

import csv
import gzip

NOT_CALLED = (
    -1.0
)  # a hard-negative protein that the conditions do not call; below any floor of 0 bits


def n_cys(seq):
    return seq.upper().count("C")


def read_fasta(path):
    out, k, buf = {}, None, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if k:
                    out[k] = "".join(buf).rstrip("*")
                k, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if k:
        out[k] = "".join(buf).rstrip("*")
    return out


def exact_matches(literature, proteome):
    """Proteome IDs whose sequence equals the sequence of a literature protein."""
    seqs = set(literature.values())
    return {p for p, s in proteome.items() if s in seqs}


def candidates(proteome, scores, meta):
    """Extra-call candidates: a score, R0 called, at least 8 cysteines, not strict, not labelled."""
    out = []
    for pid, s in scores.items():
        x = meta.get(pid)
        if x and x["r0"] == "called" and x["n_cys"] >= 8 and not x["strict"] and not x["labelled"]:
            out.append({"id": pid, "proteome": proteome, "score": s})
    return sorted(out, key=lambda c: c["id"])


def hard_negative_scores(scores, meta, groups):
    """Per group, one value per protein: its score if the conditions call it, else NOT_CALLED."""
    out = {}
    for pid, g in groups.items():
        x = meta[pid]
        called = x["r0"] == "called" and x["n_cys"] >= 8 and pid in scores
        out.setdefault(g, []).append(scores[pid] if called else NOT_CALLED)
    return out


def read_r0(path):
    with gzip.open(path, "rt") as fh:
        return {r["id"]: r["call"] for r in csv.DictReader(fh, delimiter="\t")}


def read_hits(path):
    with gzip.open(path, "rt") as fh:
        return {r["id"] for r in csv.DictReader(fh, delimiter="\t") if r.get("hit") == "1"}


def load_tuning(d):
    """Tuning-data context: tuning proteomes (sizes and per-protein metadata) and the tuning hard-negative parts."""
    from pathlib import Path

    d = Path(d)
    sc = d / "scores_pre_freeze"
    runs = {r["name"]: r for r in csv.DictReader(open(d / "run_list.tsv"), delimiter="\t")}
    split = {
        r["proteome"]: r for r in csv.DictReader(open(d / "proteome_split.tsv"), delimiter="\t")
    }
    tuning = [p for p, r in split.items() if r["part"] == "tuning"]
    pmap = list(csv.DictReader(open(d / "proteome_map.tsv"), delimiter="\t"))
    lp = read_fasta(d / "literature/jensen_sequences.faa")
    lp = {k.split("|")[1]: v for k, v in lp.items()}
    meta, sizes = {}, {}
    for p in tuning:
        r = runs[p]
        seqs = read_fasta(r["fasta"])
        w = Path(r["workdir"])
        r0 = read_r0(w / "modules/step1_rule@R0.tsv.gz")
        strict = read_hits(w / "modules/pfam_hydrophobin.tsv.gz") | read_hits(
            w / "modules/pfam_hsba.tsv.gz"
        )
        labelled = {x["proteome_id"] for x in pmap if x["proteome"] == p} | exact_matches(lp, seqs)
        sizes[p] = len(seqs)
        meta[p] = {
            i: {
                "r0": r0.get(i, ""),
                "n_cys": n_cys(s),
                "strict": i in strict,
                "labelled": i in labelled,
            }
            for i, s in seqs.items()
        }
    hn_seqs = read_fasta(sc / "hn_tuning.faa")
    groups = {
        r["id"]: r["group"]
        for r in csv.DictReader(open(sc / "hn_tuning_groups.tsv"), delimiter="\t")
    }
    r0ref = read_r0(d / "r0_reference/modules/step1_rule@R0.tsv.gz")
    hn_meta = {
        i: {"r0": r0ref.get("HN_" + i.rsplit("-", 1)[0], ""), "n_cys": n_cys(s)}
        for i, s in hn_seqs.items()
    }
    return {
        "tuning": tuning,
        "sizes": sizes,
        "meta": meta,
        "groups": groups,
        "hn_meta": hn_meta,
        "r0ref": r0ref,
    }
