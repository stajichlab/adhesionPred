#!/usr/bin/env python3
"""Phase C step 3: fit every candidate on every outer training set and score (spec 3.3).

Reads $STEP1_WORKDIR/phasec/split_members.tsv.gz and clusters.tsv.gz (09; SHA-256 checked
against splits_run.json), and the Phase B universe: phaseb/features_unique.tsv.gz and
phaseb/unique_sequences.tsv.gz (SHA-256 checked against build_run.json) and the four
embedding matrices (SHA-256 checked against embedding_run.json).

One work unit = one (split, fold, variant). Training rows: part `train` (V-go) or `train` and
`train_tc` (V-kw); y = class pos. Scored rows: parts test, test_tc, test_lit; for FULL every
Phase B unique sequence. The nested protocol is in models.py. Units run in --workers processes.

Writes to $STEP1_WORKDIR/phasec/:

  scores.tsv.gz    one row per unit, candidate and scored sequence: score (decision value;
                   empty for rules), prob (Platt), call (1/0)
  scores_run.json  fitted settings per unit and candidate (C, g, t, threshold, Platt a and b,
                   H variant, inner PR-AUC), the training-row J of every rule grid cell
                   (`rule_grid` of R1 and R2), seed, input hashes, library versions

A hash that is in two scored parts of one unit (for example a GO `test` row and a `test_lit`
row with the same sequence) is scored once; its `part` lists the parts, sorted and joined with
a comma (`test,test_lit`). Readers that index by hash ignore `part`.

STOP (exit 2, no output): stale 08 or 09 outputs; a 09 run made from other 08 outputs than
build_run.json records; changed Phase B files; an outer training set that cannot be fitted (for
example fewer than 3 positive clusters); a missing key in an input or run file.
"""

import argparse
import gzip
import io
import multiprocessing
import sys
from concurrent.futures import FIRST_EXCEPTION, ProcessPoolExecutor, wait
from pathlib import Path

import evalio
import manifest
import models
import paths
import runinfo
import truth_table
import universe

OUTPUT_NAMES = ("scores.tsv.gz", "scores_run.json")
SCORE_COLUMNS = (
    "split_id",
    "fold",
    "variant",
    "candidate",
    "seq_sha256",
    "part",
    "score",
    "prob",
    "call",
)
SCORED_PARTS = ("test", "test_tc", "test_lit")
CLUSTERS = "clusters.tsv.gz"


def unit_order(members, splits_order) -> list[tuple[str, str]]:
    keys = {(m["split_id"], m["fold"]) for m in members}
    rank = {s: i for i, s in enumerate(splits_order)}
    return sorted(keys, key=lambda k: (rank[k[0]], int(k[1])))


def tasks(members, splits_order, all_hashes, candidates, clusters=None) -> list[dict]:
    """One task per (split, fold, variant). Only training rows carry labels into a task.

    `clusters` (hash -> cluster_id from clusters.tsv.gz), when given, must agree with the
    cluster ids of the members. A training hash that is also scored stops the run, except in
    FULL, which scores every sequence."""
    by_unit: dict[tuple, list] = {}
    for m in members:
        by_unit.setdefault((m["split_id"], m["fold"]), []).append(m)
    out = []
    for key in unit_order(members, splits_order):
        rows = by_unit[key]
        if key[0] == "FULL":
            score, parts = list(all_hashes), {}
        else:
            scored = [m for m in rows if m["part"] in SCORED_PARTS]
            found: dict[str, set] = {}
            for m in scored:
                found.setdefault(m["seq_sha256"], set()).add(m["part"])
            score = list(found)  # one entry per hash, in member order
            parts = {h: ",".join(sorted(p)) for h, p in found.items()}
        for variant in evalio.VARIANTS:
            use = ("train", "train_tc") if variant == "V-kw" else ("train",)
            train = [m for m in rows if m["part"] in use]
            bad = [m for m in train if m["class"] not in ("pos", "neg")]
            if bad:
                raise evalio.StopError(
                    f"{key[0]}|{key[1]}: training row {bad[0]['seq_sha256']} has class {bad[0]['class']}"
                )
            if clusters is not None:
                wrong = [m for m in train if clusters.get(m["seq_sha256"]) != m["cluster_id"]]
                if wrong:
                    raise evalio.StopError(
                        f"{key[0]}|{key[1]}: cluster_id of {wrong[0]['seq_sha256']} differs from clusters.tsv.gz"
                    )
            if key[0] != "FULL":
                shared = {m["seq_sha256"] for m in train} & set(score)
                if shared:
                    raise evalio.StopError(
                        f"{key[0]}|{key[1]}: {len(shared)} sequences are in training and scored parts"
                    )
            out.append(
                {
                    "key": f"{key[0]}|{key[1]}|{variant}",
                    "split_id": key[0],
                    "fold": key[1],
                    "variant": variant,
                    "train": [m["seq_sha256"] for m in train],
                    "y": [m["class"] == "pos" for m in train],
                    "groups": [m["cluster_id"] for m in train],
                    "score": score,
                    "parts": parts,
                    "seed": evalio.SEED,
                    "candidates": candidates,
                }
            )
    return out


def run_tasks(task_list, u, workers: int, phaseb: Path):
    if workers <= 1:
        return [models.run_task(t, u) for t in task_list]
    ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=ctx, initializer=models.init_worker, initargs=(str(phaseb),)
    ) as pool:
        futures = [pool.submit(models.run_task, t) for t in task_list]
        wait(futures, return_when=FIRST_EXCEPTION)
        failed = [f for f in futures if f.done() and not f.cancelled() and f.exception()]
        if failed:
            pool.shutdown(wait=True, cancel_futures=True)  # do not start the other units
            raise failed[0].exception()
        return [f.result() for f in futures]


def _fmt(x) -> str:
    return format(float(x), ".10g")


def write_scores(path: Path, task_list, results) -> None:
    raw = open(path, "wb")
    with raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        text = io.TextIOWrapper(gz, encoding="utf-8", newline="")
        text.write("\t".join(SCORE_COLUMNS) + "\n")
        for task, (key, _, scores) in zip(task_list, results, strict=True):
            if key != task["key"]:
                raise evalio.StopError(f"result for unit {key} where {task['key']} was expected")
            head = f"{task['split_id']}\t{task['fold']}\t{task['variant']}\t"
            for cand in task["candidates"]:
                s = scores[cand]
                for i, h in enumerate(task["score"]):
                    part = task["parts"].get(h, "all")
                    score = "" if s["score"] is None else _fmt(s["score"][i])
                    prob = "" if s["prob"] is None else _fmt(s["prob"][i])
                    call = "1" if s["call"][i] else "0"
                    text.write(f"{head}{cand}\t{h}\t{part}\t{score}\t{prob}\t{call}\n")
        text.flush()
        text.detach()


def run(work: Path, workers: int, candidates=models.CANDIDATES, arguments=()):
    work = Path(work)
    out = evalio.out_dir(work)
    build = evalio.require_current(out, "build_run.json", ("eval_table.tsv.gz",), "08")
    split_log = evalio.require_current(
        out, "splits_run.json", ("clusters.tsv.gz", "split_members.tsv.gz"), "09_make_splits.py"
    )
    for name in ("eval_table.tsv.gz", "eval_literature.tsv", "eval_sequences.fasta.gz"):
        if split_log["input_sha256"].get(name) != build["outputs_sha256"].get(name):
            raise evalio.StopError(
                f"splits_run.json was made from another {name} than build_run.json records; "
                "re-run 09_make_splits.py"
            )
    phaseb = work / "phaseb"
    for name in ("features_unique.tsv.gz", "unique_sequences.tsv.gz"):
        if manifest.sha256_file(phaseb / name) != build["input_sha256"][name]:
            raise evalio.StopError(f"phaseb/{name} changed after 08 read it; re-run 08 to 10")
    emb_run = evalio.read_json(phaseb / "emb" / "embedding_run.json")
    if emb_run.get("unique_sequences_sha256") != build["input_sha256"]["unique_sequences.tsv.gz"]:
        raise evalio.StopError(
            "embedding_run.json names another unique_sequences.tsv.gz than 08 read"
        )
    u = universe.load(phaseb, verify=True)
    members = truth_table.read_tsv(out / "split_members.tsv.gz")
    clusters = {r["seq_sha256"]: r["cluster_id"] for r in truth_table.read_tsv(out / CLUSTERS)}
    task_list = tasks(members, split_log["splits"], u.hashes, tuple(candidates), clusters)
    results = run_tasks(task_list, u, workers, phaseb)
    log = {
        "all_sources": build["all_sources"],
        "truth_set_sha256": build["truth_set_sha256"],
        "input_sha256": {
            "split_members.tsv.gz": split_log["outputs_sha256"]["split_members.tsv.gz"],
            "clusters.tsv.gz": split_log["outputs_sha256"]["clusters.tsv.gz"],
            "features_unique.tsv.gz": build["input_sha256"]["features_unique.tsv.gz"],
            "unique_sequences.tsv.gz": build["input_sha256"]["unique_sequences.tsv.gz"],
        },
        "embedding_array_sha256": {
            f"{m}.{w}": emb_run["models"][m][w]["array_sha256"]
            for m in universe.MODELS
            for w in ("nterm", "cterm")
        },
        "seed": evalio.SEED,
        "candidates": list(candidates),
        "c_grid": list(models.C_GRID),
        "inner_folds": models.INNER_FOLDS,
        "units": {key: params for key, params, _ in results},
        "git_commit": runinfo.git_commit(),
        "library_versions": evalio.library_versions(),
        "arguments": list(arguments),
    }
    evalio.write_outputs(
        out,
        {"scores.tsv.gz": lambda p: write_scores(p, task_list, results)},
        "scores_run.json",
        log,
    )
    return log


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--candidates", default=",".join(models.CANDIDATES), help="comma list (tests, golden file)"
    )
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        cands = tuple(args.candidates.split(","))
        unknown = sorted(set(cands) - set(models.CANDIDATES))
        if unknown:
            raise evalio.StopError(f"unknown candidates {unknown}")
        cands = tuple(c for c in models.CANDIDATES if c in cands)
        log = run(work, args.workers, cands, list(argv) if argv is not None else sys.argv[1:])
    except (evalio.StopError, models.ModelError, ValueError, OSError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    except KeyError as exc:
        print(f"STOP: missing key {exc} in an input or run file", file=sys.stderr)
        return 2
    print(f"units={len(log['units'])} candidates={len(log['candidates'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
