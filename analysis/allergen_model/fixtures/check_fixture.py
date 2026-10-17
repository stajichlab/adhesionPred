#!/usr/bin/env python3.12
"""T1.6 checks: assert the dry-run outputs on the fixture, then assert that broken inputs fail cleanly.

Run after run_dryrun.sh:  /usr/bin/python3.12 check_fixture.py
Exit 0 if every check passes; otherwise the failing checks are listed and the exit code is 1.
Identities are read from the MMseqs2 tables (measured), not from the planted values, because the
planted identity is approximate.
"""

import json
import os
import subprocess
import sys
import tempfile

FIX = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(FIX, "dryrun_out")
JOBS = os.path.join(FIX, "..", "jobs")
EXPECTED = json.load(open(os.path.join(FIX, "data", "expected.json")))
COV_MIN = 0.5  # owner decision O7
T_C_VALUES = (0.6, 0.7)
failures = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (f"  [{detail}]" if detail and not cond else ""))
    if not cond:
        failures.append(name)


def zst_rows(path):
    out = subprocess.run(["zstd", "-dc", path], capture_output=True, text=True, check=True).stdout
    lines = [ln.split("\t") for ln in out.strip().split("\n")]
    header, rows = lines[0], lines[1:]
    return [dict(zip(header, r, strict=False)) for r in rows]


# 1. Clustering
clu = {}
for ln in open(os.path.join(WORK, "out", "cluster", "clusters.tsv")):
    acc, cid = ln.rstrip("\n").split("\t")
    clu[acc] = cid
check("cluster: one row per sequence (9)", len(clu) == 9, str(len(clu)))
for a, b in EXPECTED["cluster_same"]:
    check(
        f"cluster: {a} and {b} share a cluster",
        clu.get(a) == clu.get(b),
        f"{clu.get(a)} vs {clu.get(b)}",
    )
for a, b in EXPECTED["cluster_different"]:
    check(f"cluster: {a} and {b} differ", clu.get(a) != clu.get(b), f"{clu.get(a)} vs {clu.get(b)}")

# 2. Pairs, and the pairing rule (T_C plus the O7 coverage condition)
pairs = zst_rows(os.path.join(WORK, "out", "pairs", "pairs_candidates_vs_positives.tsv.zst"))
check("pairs: N_c4 has no hit", not any(r["query"] == "N_c4" for r in pairs))
for t in T_C_VALUES:
    got = {}
    for r in pairs:
        if float(r["fident"]) >= t and float(r["qcov"]) >= COV_MIN and float(r["tcov"]) >= COV_MIN:
            got.setdefault(r["query"], set()).add(r["target"])
    if t == 0.6:
        check(
            "pairing T_C 0.6: N_c1 pairs with P_B1",
            got.get("N_c1") == {"P_B1"},
            str(got.get("N_c1")),
        )
    check(
        f"pairing T_C {t}: N_c2 does not pair (identity 0.33)",
        "N_c2" not in got,
        str(got.get("N_c2")),
    )
    check(
        f"pairing T_C {t}: N_c3 does not pair (coverage 0.23 fails O7)",
        "N_c3" not in got,
        str(got.get("N_c3")),
    )
# Without the coverage rule N_c3 would pair; the check shows the rule is doing work.
n_c3_hi = [r for r in pairs if r["query"] == "N_c3" and float(r["fident"]) >= 0.6]
check(
    "pairs: N_c3 has high identity but low coverage (rule matters)",
    bool(n_c3_hi) and float(n_c3_hi[0]["qcov"]) < COV_MIN,
)

# 3. Cross-split
cs = os.path.join(WORK, "out", "xsplit")
mx = float(open(os.path.join(cs, "cross_split_max.txt")).read().strip())
check(
    "cross-split: planted violation found (max fident >= 0.75)",
    mx >= EXPECTED["cross_split_max_fident_min"],
    str(mx),
)
xrows = zst_rows(os.path.join(cs, "cross_split_pairs.tsv.zst"))
check(
    "cross-split: pair P_B1m vs P_B1 is listed",
    any({r["query"], r["target"]} == {"test_P_B1m", "train_P_B1"} for r in xrows),
)
check("cross-split: each unordered pair written once", all(r["query"] < r["target"] for r in xrows))
done = open(os.path.join(cs, "cross_split.done")).read()
check(
    "cross-split: done marker lists train, val, test counts",
    all(f"  {s}" in done for s in ("2 test", "3 train", "1 val")),
)

# 4. Failure cases: each must exit nonzero and leave no .done marker
JOB_ENV = {"SCRATCH": os.path.join(WORK, "scratch_neg")}


def run_job(job, env_extra, outdir):
    os.makedirs(os.path.join(WORK, "scratch_neg"), exist_ok=True)
    env = dict(os.environ, **JOB_ENV, **env_extra, OUTDIR=outdir)
    p = subprocess.run(
        ["bash", "-l", os.path.join(JOBS, job)], env=env, capture_output=True, text=True
    )
    return p


tmpd = tempfile.mkdtemp(prefix="negcases.", dir=WORK)

# 4a. Candidates without N_ prefix: the cluster job must refuse.
bad = os.path.join(tmpd, "bad_prefix.faa")
with open(bad, "w") as fh:
    fh.write(open(os.path.join(FIX, "data", "positives.faa")).read())
    fh.write(">X_novel\nMKVLAAGIVGLLLAAQPAMA\n")
out = os.path.join(tmpd, "c1")
p = run_job("cluster_mmseqs.sbatch", {"SEQS": bad, "ARM": "shared"}, out)
check(
    "neg: cluster job refuses an input with no N_ candidates",
    p.returncode != 0 and not os.path.exists(os.path.join(out, "cluster.done")),
    p.stderr[-200:],
)

# 4b. Split FASTA with one prefix only: the cross-split job must refuse (it cannot pass vacuously).
one = os.path.join(tmpd, "one_split.faa")
with open(one, "w") as fh:
    fh.write(open(os.path.join(FIX, "data", "positives.faa")).read().replace(">P_", ">train_P_"))
out = os.path.join(tmpd, "x1")
p = run_job("cross_split_check.sbatch", {"SPLIT_FASTA": one, "ARM": "shared"}, out)
check(
    "neg: cross-split job refuses a single-split input",
    p.returncode != 0 and not os.path.exists(os.path.join(out, "cross_split.done")),
    p.stderr[-200:],
)

# 4c. Empty candidate set: the pairs job must refuse before searching.
empty = os.path.join(tmpd, "empty.faa")
open(empty, "w").close()
out = os.path.join(tmpd, "p1")
p = run_job(
    "pairs_search.sbatch",
    {"QUERY": empty, "TARGET": os.path.join(FIX, "data", "positives.faa"), "ARM": "shared"},
    out,
)
check(
    "neg: pairs job refuses an empty candidate file",
    p.returncode != 0 and not os.path.exists(os.path.join(out, "pairs.done")),
    p.stderr[-200:],
)

if failures:
    print(f"\n{len(failures)} check(s) failed: {failures}")
    sys.exit(1)
print("\nall fixture checks passed")
