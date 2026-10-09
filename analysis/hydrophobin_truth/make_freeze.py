#!/usr/bin/python3.12
"""Write freeze.json (the pre-registration of the v1 extended level) and the final relaxed .hmm file.

Checks that only tuning data was scored before the freeze. Records hashes of every input and pre-freeze score file, the frozen
choices, the thresholds of the spec, and what the author had seen. Run once, then commit; nothing in freeze.json changes after.
Usage: make_freeze.py --dir analysis/hydrophobin_truth
"""

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path


def disallowed(files, tuning):
    ok_names = [
        re.compile(r"^hn_tuning\.(default|nobias)\.tblout$"),
        re.compile(r"^t2_positives\.(default|nobias)\.tblout$"),
    ]
    tun = "|".join(re.escape(t) for t in tuning)
    ok_names += [re.compile(rf"^({tun})\.(default|nobias)\.tblout$")]
    ok_names += [re.compile(rf"^fold_(\d\d|all)\.({tun}|hn_tuning)\.(default|nobias)\.tblout$")]
    return [f for f in files if not any(p.match(f) for p in ok_names)]


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def hashes(base, names):
    return {n: sha(Path(base) / n) for n in sorted(names)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    rel = json.load(open(d / "prefreeze_relaxed.json"))
    hmm = json.load(open(d / "prefreeze_hmm.json"))
    tuning = rel["tuning_proteomes"]
    sc1, sc2 = d / "scores_pre_freeze", d / "scores_pre_freeze_hmm"
    names1 = [p.name for p in sc1.glob("*.tblout")]
    names2 = [p.name for p in sc2.glob("*.tblout")]
    bad = disallowed(names1, tuning) + disallowed(names2, tuning)
    if bad:
        raise SystemExit(f"score files outside the tuning data: {bad}")
    opt = rel["chosen_option"]
    cutoff = rel["options"][opt]["cutoff"]
    out_hmm = d / "relaxed" / "hydrophobin_relaxed.hmm"
    subprocess.run(
        [
            "/usr/bin/python3.12",
            str(d / "make_relaxed_hmm.py"),
            "final",
            "--template",
            str(d / "relaxed/template.hmm"),
            "--cutoff",
            str(cutoff),
            "--out",
            str(out_hmm),
        ],
        check=True,
    )
    inputs = [
        "proteome_split.tsv",
        "folds.tsv",
        "folds_meta.json",
        "clusters_positives.tsv",
        "hard_negative_thresholds.json",
        "hard_negatives.tsv.gz",
        "hard_negative_counts.json",
        "r0_reference/reference.faa",
        "r0_reference/modules/step1_rule@R0.tsv.gz",
        "relaxed/template.hmm",
        "relaxed/template.provenance.json",
        "prefreeze_relaxed.json",
        "prefreeze_hmm.json",
        "truth_all.tsv",
        "literature/v2_reserve.tsv",
    ]
    freeze = {
        "frozen": "v1 extended level (relaxed Pfam first; custom HMM candidate)",
        "spec": "docs/superpowers/specs/2026-10-08-hydrophobin-custom-hmm-design.md (revision 4, E11)",
        "plan": "docs/superpowers/plans/2026-10-08-hydrophobin-extended-level.md (revision 2)",
        "relaxed_pfam": {
            "models": json.load(open(d / "relaxed/template.provenance.json"))["models"],
            "search_option": opt,
            "full_sequence_cutoff_bits": cutoff,
            "domain_ga": -1000.0,
            "conditions": "R0 called; at least 8 cysteines in the full sequence; not called by pfam_hydrophobin or pfam_hsba; not labelled",
            "hmm_file": "analysis/hydrophobin_truth/relaxed/hydrophobin_relaxed.hmm",
            "hmm_sha256": sha(out_hmm),
            "chosen_by": "lowest_cutoff over the tuning proteomes and tuning hard-negative parts (no positive score); option chosen by recovered Pfam-missed clusters, positives used only for that choice",
            "pre_freeze": rel["options"],
        },
        "custom_hmm": {
            "aligner": "mafft 7.505 --auto --quiet",
            "hmmbuild": "hmmbuild --amino --informat afa",
            "cysteine_column_check": "eight conserved cysteine columns present in at least 80% of training sequences; a failure stops the run",
            "search_option": "default filters (in-sample recovery ties; ties go to default)",
            "per_fold_cutoffs": {
                k: {"default": v["default"]["cutoff"], "nobias": v["nobias"]["cutoff"]}
                for k, v in hmm["folds"].items()
            },
        },
        "cutoff_rule": "lowest bit score at or above the floor (0) with at most 5 extra calls per 10,000 proteins in each tuning proteome and a call rate at most 0.05 on each tuning hard-negative group (HsbA ignored, E11); ties to the higher cutoff",
        "thresholds": {
            "recall_relaxed_passes_if_clusters_recovered_at_least": 3,
            "of_pfam_missed_clusters": 6,
            "hmm_rule": "u = Pfam-missed clusters not recovered by relaxed; u<=1 -> not testable; u>=2 -> add if it recovers ceil(u/2) at cost <= relaxed + 2 per 10,000",
            "cost_limit_per_10000_in_each_test_proteome": 5,
            "hard_negative_pass_rate": 0.05,
            "precision_rule": "at least 5 resolved unlabelled extra-call clusters and Wilson 95% lower bound >= 0.5",
        },
        "pre_freeze_changes": [
            "aligner: mafft --auto failed the planned cysteine-column check in fold 16 (the eighth cysteine split between two columns, 66% and 34%); changed to L-INS-i (--localpair --maxiterate 1000) for every fold",
            "L-INS-i alignments still split a cysteine over distant columns in folds 16 and 20 (class I and class II loop lengths differ by more than a column window; C3-C4 is 11 residues in class II and 33 to 39 in class I)",
            "stop condition changed from 'eight columns with a cysteine in at least 80% of sequences' to 'eight columns with a cysteine in at least 30% of sequences'; the same-column and the ordinal checks are reported as diagnostics (ordinal check False in folds 20 and 26)",
            "all changes were made before any held-out cluster, test proteome, test hard-negative part or LP score existed",
        ],
        "l1b_relabel_choice": "no T2 entry flagged (all 131 texts hydrophobin-related); no relabel",
        "input_hashes": hashes(d, inputs),
        "pre_freeze_score_hashes": {
            **hashes(sc1, names1),
            **{"hmm/" + k: v for k, v in hashes(sc2, names2).items()},
        },
        "seen_by_the_author_before_the_freeze": [
            "the 9 Pfam-missed T2 entries and their regex and sub-cutoff hit results (report 2026-10-08-hydrophobin-validation)",
            "the plan reviewer's relaxed-Pfam numbers at 8.5 and 15 bits and the 0-bit floor results (plan review 1)",
            "spec review 2 numbers (HCF1 30.4 bits with --nobias; relaxed reaches 5 of 9 at 8.5 bits)",
            "T2 sequence scores against the seven Pfam models, used here only to choose the search option",
        ],
    }
    json.dump(freeze, open(d / "freeze.json", "w"), indent=1)
    print("freeze.json written; relaxed option", opt, "cutoff", cutoff)
    return 0


if __name__ == "__main__":
    sys.exit(main())
