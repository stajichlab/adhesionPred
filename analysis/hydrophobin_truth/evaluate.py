#!/usr/bin/python3.12
"""L5 evaluation (after the freeze). Reads freeze.json and the post-freeze score files.

Refuses to run unless the freeze commit is an ancestor of HEAD and freeze.json is unchanged. Negatives are assumed (absent
annotation). Leave-cluster-out numbers are labelled `partial`: the author saw the Pfam-missed entries before the freeze.
Usage: evaluate.py --dir analysis/hydrophobin_truth --out docs/reports/data/sorting_hat/hydrophobin_ext
"""

import argparse
import csv
import gzip
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
Z = 1.959964
PFAM_MISSED = set("A0A0A2KI06 A0A6V8R0V1 G9MKB2 I1RDP9 I1RDQ0 O94217 Q00367 Q0KKA0 Q7Z9L4".split())


def wilson(k, n):
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def called_set(scores, cutoff, r0, n_cys):
    """Ids with score at or above the cutoff, R0 called and at least 8 cysteines."""
    return {
        i
        for i, s in scores.items()
        if s >= cutoff and r0.get(i) == "called" and n_cys.get(i, 0) >= 8
    }


def recovered_clusters(called, clusters, pfam_missed):
    """A Pfam-missed cluster is recovered when at least one Pfam-missed member of it is called."""
    return {c for c, members in clusters.items() if members & pfam_missed & called}


def unrecovered(recovered, all_missed_clusters):
    return [c for c in all_missed_clusters if c not in recovered]


def relaxed_recall_pass(recovered, n_clusters, need=3):
    return len(recovered) >= need


def hmm_decision(u, recovered_by_hmm, unrecovered_by_relaxed):
    if u <= 1:
        return "not testable"
    need = math.ceil(u / 2)
    got = len(set(recovered_by_hmm) & set(unrecovered_by_relaxed))
    return "recall rule met" if got >= need else "recall rule not met"


def parse_tblout(path):
    best = {}
    if not Path(path).exists():
        return best
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        c = line.split()
        s = float(c[5])
        if c[0] not in best or s > best[c[0]]:
            best[c[0]] = s
    return best


def read_tsv(p):
    opener = gzip.open if str(p).endswith(".gz") else open
    with opener(p, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check_freeze(d):
    freeze = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", "analysis/hydrophobin_truth/freeze.json"],
        capture_output=True,
        text=True,
        cwd=ROOT,
    ).stdout.strip()
    if not freeze:
        raise SystemExit("freeze.json is not committed")
    if subprocess.run(["git", "merge-base", "--is-ancestor", freeze, "HEAD"], cwd=ROOT).returncode:
        raise SystemExit("the freeze commit is not an ancestor of HEAD")
    if subprocess.run(
        ["git", "diff", "--quiet", freeze, "--", "analysis/hydrophobin_truth/freeze.json"], cwd=ROOT
    ).returncode:
        raise SystemExit("freeze.json changed after the freeze")
    return freeze


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d, out = Path(a.dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    freeze_commit = check_freeze(d)
    fz = json.load(open(d / "freeze.json"))
    sc = d / "scores_post_freeze"
    pre = d / "scores_pre_freeze"
    cut_rel = fz["relaxed_pfam"]["full_sequence_cutoff_bits"]
    fold_cut = fz["custom_hmm"]["per_fold_cutoffs"]
    spec = importlib.util.spec_from_file_location(
        "cys8", ROOT / "src/cellsurface_sorting_hat/modules/cys8.py"
    )
    cys8 = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ROOT / "src"))
    spec.loader.exec_module(cys8)
    spacing = cys8.load_sets(ROOT / "data/sorting_hat/cys8_spacing.yaml")

    truth = {r["accession"]: r for r in read_tsv(d / "truth_all.tsv")}
    seqs = {i: r["sequence"] for i, r in truth.items()}
    t2 = {i for i, r in truth.items() if r["tier"] == "T2"}
    r0ref = {r["id"]: r["call"] for r in read_tsv(d / "r0_reference/modules/step1_rule@R0.tsv.gz")}
    ncys = {i: s.upper().count("C") for i, s in seqs.items()}
    clusters = {}
    for r in read_tsv(d / "clusters_positives.tsv"):
        if r["tier"] == "T2":
            clusters.setdefault(r["cluster"], set()).add(r["id"])
    folds = {int(r["fold"]): r for r in read_tsv(d / "folds.tsv")}
    missed_clusters = [c for c, m in clusters.items() if m & PFAM_MISSED]

    res = {
        "freeze_commit": freeze_commit,
        "label": "leave-cluster-out numbers are partial (the author saw the Pfam-missed entries before the freeze); negatives are assumed",
    }

    # strict Pfam GA on T2
    ga_t2 = parse_tblout(sc / "gagate.t2.tblout")
    # The strict level of the tool reads domain rows at --cut_ga (module pfam_hydrophobin): 122 T2 entries. The sequence-level
    # --tblout at --cut_ga also lists proteins that pass the sequence cutoff but have no domain at the domain cutoff.
    strict_t2 = t2 - PFAM_MISSED
    res["pfam_ga_missed_t2"] = sorted(PFAM_MISSED)
    res["strict_sequence_level_tblout_hits_t2"] = len(ga_t2)
    res["tblout_hits_that_are_not_strict_domain_hits"] = sorted(set(ga_t2) & PFAM_MISSED)
    # relaxed on T2 (scored before the freeze for the option choice; same file)
    rel_t2 = parse_tblout(pre / "t2_positives.nobias.tblout")
    called_rel = called_set(rel_t2, cut_rel, r0ref, ncys)
    rec_rel = recovered_clusters(called_rel, clusters, PFAM_MISSED)
    u_list = unrecovered(rec_rel, missed_clusters)
    res["relaxed"] = {
        "cutoff": cut_rel,
        "option": "nobias",
        "pfam_missed_proteins_called": sorted(called_rel & PFAM_MISSED),
        "pfam_missed_clusters": len(missed_clusters),
        "clusters_recovered": sorted(rec_rel),
        "clusters_recovered_n": len(rec_rel),
        "wilson_clusters": wilson(len(rec_rel), len(missed_clusters)),
        "wilson_proteins": wilson(len(called_rel & PFAM_MISSED), len(PFAM_MISSED)),
        "u": len(u_list),
        "unrecovered": u_list,
        "recall_rule_pass": relaxed_recall_pass(rec_rel, len(missed_clusters)),
    }
    ext = strict_t2 | called_rel
    res["extended_recall_t2"] = {
        "strict_ga": len(strict_t2),
        "extended": len(ext),
        "n": len(t2),
        "wilson_strict": wilson(len(strict_t2), len(t2)),
        "wilson_extended": wilson(len(ext), len(t2)),
    }

    # frozen-spacing rule on T2 (with R0)
    cys_called = set()
    for i in t2:
        c = cys8.classify(seqs[i], spacing)
        if c["pattern_match"] and r0ref.get(i) == "called":
            cys_called.add(i)
    res["cys8_pattern"] = {
        "pfam_missed_called": sorted(cys_called & PFAM_MISSED),
        "t2_called": len(cys_called),
    }

    # HMM leave-cluster-out (frozen option: default filters)
    hmm_called, hmm_called_nobias = set(), set()
    for i, f in folds.items():
        name = f"fold_{i:02d}"
        for opt, bag in (("default", hmm_called), ("nobias", hmm_called_nobias)):
            s = parse_tblout(sc / f"hmm.{name}.test.{opt}.tblout")
            cut = fold_cut[name][opt]
            bag |= called_set(s, cut, r0ref, ncys) & set(f["test_members"].split(","))
    rec_hmm = recovered_clusters(hmm_called, clusters, PFAM_MISSED)
    rec_hmm_nb = recovered_clusters(hmm_called_nobias, clusters, PFAM_MISSED)
    res["hmm_loco"] = {
        "option": "default (frozen)",
        "t2_called": len(hmm_called & t2),
        "t2_n": len(t2),
        "wilson_t2": wilson(len(hmm_called & t2), len(t2)),
        "pfam_missed_proteins_called": sorted(hmm_called & PFAM_MISSED),
        "clusters_recovered": sorted(rec_hmm),
        "clusters_recovered_n": len(rec_hmm),
        "nobias_clusters_recovered_n": len(rec_hmm_nb),
        "nobias_t2_called": len(hmm_called_nobias & t2),
        "decision": hmm_decision(len(u_list), rec_hmm, u_list),
    }
    # T2 lost to the conditions
    lost = sorted(i for i in t2 if r0ref.get(i) != "called" or ncys[i] < 8)
    res["t2_lost_to_conditions"] = {
        "n": len(lost),
        "ids": lost,
        "no_r0_call": sorted(i for i in t2 if r0ref.get(i) != "called"),
        "fewer_than_8_cys": sorted(i for i in t2 if ncys[i] < 8),
    }

    # hard negatives, test part
    thr = json.load(open(d / "hard_negative_thresholds.json"))
    groups = {r["id"]: r["group"] for r in read_tsv(sc / "fasta/hn_test_groups.tsv")}
    hn_seq = {}
    k = None
    for line in open(sc / "fasta/hn_test.faa"):
        line = line.rstrip()
        if line.startswith(">"):
            k = line[1:]
            hn_seq[k] = ""
        else:
            hn_seq[k] += line
    hn_n = {i: s.upper().count("C") for i, s in hn_seq.items()}
    hn_r0 = {i: r0ref.get("HN_" + i.rsplit("-", 1)[0], "") for i in hn_seq}
    rel_hn = parse_tblout(sc / "relaxed.hn_test.tblout")
    ga_hn = parse_tblout(sc / "gagate.hn_test.tblout")
    hmm_hn = parse_tblout(sc / "hmm.fold_all.hn_test.default.tblout")
    cut_all = fold_cut["fold_all"]["default"]
    hn = {}
    for g in sorted(set(groups.values())):
        ids = [i for i, gg in groups.items() if gg == g]
        rel_c = called_set(rel_hn, cut_rel, hn_r0, hn_n) & set(ids)
        hmm_c = called_set(hmm_hn, cut_all, hn_r0, hn_n) & set(ids)
        ga_c = set(ga_hn) & set(ids)
        hn[g] = {
            "n_test": len(ids),
            "pfam_ga": len(ga_c),
            "relaxed": len(rel_c),
            "hmm_all": len(hmm_c),
            "rate_relaxed": round(len(rel_c) / len(ids), 4),
            "rate_hmm_all": round(len(hmm_c) / len(ids), 4),
            "threshold": thr.get(g),
            "relaxed_passes": None if thr.get(g) is None else len(rel_c) / len(ids) <= thr[g],
        }
    res["hard_negatives_test"] = hn

    # LP secondary
    lp_seq = {}
    for line in open(sc / "fasta/lp_scored.faa"):
        line = line.rstrip()
        if line.startswith(">"):
            k = line[1:]
            lp_seq[k] = ""
        else:
            lp_seq[k] += line
    lp_n = {i: s.upper().count("C") for i, s in lp_seq.items()}
    lp_r0 = {i: r0ref.get(i, "") for i in lp_seq}
    ga_lp = set(parse_tblout(sc / "gagate.lp_scored.tblout"))
    rel_lp = called_set(parse_tblout(sc / "relaxed.lp_scored.tblout"), cut_rel, lp_r0, lp_n)
    hmm_lp = called_set(
        parse_tblout(sc / "hmm.fold_all.lp_scored.default.tblout"), cut_all, lp_r0, lp_n
    )
    lp_missed = set(lp_seq) - ga_lp
    cys_lp = {
        i
        for i, s in lp_seq.items()
        if cys8.classify(s, spacing)["pattern_match"] and lp_r0.get(i) == "called"
    }
    res["lp_secondary"] = {
        "scored": len(lp_seq),
        "pfam_ga_hits": len(ga_lp),
        "pfam_ga_missed": sorted(lp_missed),
        "relaxed_on_missed": sorted(rel_lp & lp_missed),
        "hmm_all_on_missed": sorted(hmm_lp & lp_missed),
        "cys8_on_missed": sorted(cys_lp & lp_missed),
        "note": "predictions (Jensen 2010 genome screen); consistency check, not a measurement; the reserved clusters and the 8 unverified proteins are not scored",
    }
    json.dump(res, open(out / "l5_results.json", "w"), indent=1)
    # per-cluster recovery table
    with open(out / "pfam_missed_clusters.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "cluster",
                "members",
                "pfam_missed_members",
                "relaxed_recovered",
                "hmm_loco_recovered",
                "cys8_recovered",
            ]
        )
        for c in sorted(missed_clusters):
            mm = clusters[c] & PFAM_MISSED
            w.writerow(
                [
                    c,
                    ",".join(sorted(clusters[c])),
                    ",".join(sorted(mm)),
                    int(c in rec_rel),
                    int(c in rec_hmm),
                    int(bool(mm & cys_called)),
                ]
            )
    print(
        json.dumps(
            {
                k: v
                for k, v in res.items()
                if k
                in (
                    "relaxed",
                    "hmm_loco",
                    "extended_recall_t2",
                    "t2_lost_to_conditions",
                    "pfam_missed_definition_matches_the_9",
                )
            },
            indent=1,
        )
    )
    print("hard negatives:", json.dumps(hn, indent=1))
    print("LP:", json.dumps(res["lp_secondary"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
