#!/usr/bin/python3.12
"""Anchored search applied to every class-2a candidate family, not only SOWgp.

Question: does a motif-anchored search find family members that the periodicity detectors
(02_repeat_profile.py, 14_repeat_detect_general.py) do not call?

Method, one pass over the proteomes per phase so 58 families cost about the same as one:

  1. For each class-2a candidate protein, take its repeat units from 14's detector and
     enumerate candidate anchor k-mers (30_anchor_family_search.derive path, without the
     background counting).
  2. ONE pass over the 7 proteomes counting every candidate motif of every family, in the
     real proteins and again in a per-protein residue shuffle. The shuffled count is the
     chance rate. This is the step that decides whether a family HAS a usable anchor.
  3. Pick one anchor per family: highest conservation among motifs whose chance count is at
     or below --max-null. Families with no such motif are reported as failures, with the
     composition class, because that is the result for Ser/Thr-rich families.
  4. ONE pass collecting hits for all selected anchors, then align each family's hits to its
     own reference protein and classify them with 30's classifier.
  5. Compare the hit set to the proteins 02 and 14 call as repeats.

Outputs
  --out-families  one row per family: anchor, specificity, hits, and how many hits 02/14 miss
  --out-hits      one row per hit, the 30_anchor_family_search columns plus `family`

Run: sbatch 33_anchor_sweep.sh   (about 25 min on one core; see the script header)
"""

import argparse
import importlib.util
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "_common"))


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, HERE / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


A = load("anchor30", "30_anchor_family_search.py")
DET = load("detect14", "14_repeat_detect_general.py")


def read_candidates(paths):
    """(strain, protein, comp_class) rows from the 03-style candidate TSVs."""
    rows = []
    for p in paths:
        if not Path(p).exists():
            print(f"skip missing {p}", file=sys.stderr)
            continue
        with A._open(p) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {k: i for i, k in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rows.append(
                    (
                        f[ix["strain"]],
                        f[ix["protein"]],
                        f[ix["comp_class"]] if "comp_class" in ix else "",
                        Path(p).name,
                    )
                )
    return rows


def detector_calls(paths):
    """{protein_id} for every protein a detector TSV calls a repeat (rep_period > 0)."""
    out = set()
    for p in paths:
        if not Path(p).exists():
            continue
        with A._open(p) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {k: i for i, k in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                try:
                    if float(f[ix["rep_period"]] or 0) > 0:
                        out.add(f[ix["protein"]])
                except (ValueError, IndexError):
                    continue
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "--candidates",
        nargs="+",
        default=["class2a_candidates_general.tsv", "class2a_candidates.tsv"],
    )
    ap.add_argument("--proteomes", nargs="+", default=["longread"])
    ap.add_argument(
        "--profiles-02",
        nargs="+",
        default=["repeat_profile_longread.tsv", "repeat_profile_reference.tsv"],
    )
    ap.add_argument(
        "--profiles-14",
        nargs="+",
        default=["repeat_general_longread.tsv", "repeat_general_reference.tsv"],
    )
    ap.add_argument(
        "--max-null",
        type=int,
        default=5,
        help="max proteins in the shuffled null a usable anchor may hit",
    )
    ap.add_argument("--min-conservation", type=float, default=0.6)
    ap.add_argument("--kmin", type=int, default=A.KMIN)
    ap.add_argument("--kmax", type=int, default=A.KMAX)
    ap.add_argument(
        "--per-family",
        type=int,
        default=25,
        help="candidate motifs carried to the background count, per family",
    )
    ap.add_argument("--out-families", default="class2a_anchor_families.tsv")
    ap.add_argument("--out-hits", default="class2a_anchor_hits.tsv")
    args = ap.parse_args()

    proteomes = A.resolve_proteomes(args.proteomes)
    print(f"reading {len(proteomes)} proteomes", flush=True)
    seqs = {}  # protein id -> sequence
    where = {}  # protein id -> proteome label
    per_proteome = {}  # label -> {id: seq}
    for p in proteomes:
        d = A.read_fasta(p)
        lab = A.label_of(p)
        per_proteome[lab] = d
        for k, v in d.items():
            seqs[k] = v
            where[k] = lab
    print(f"{len(seqs)} proteins", flush=True)

    cands = read_candidates(args.candidates)
    print(f"{len(cands)} candidate rows", flush=True)

    # --- phase 1: per-family candidate motifs -------------------------------------------
    fams = {}
    for strain, prot, comp, src in cands:
        if prot in fams or prot not in seqs:
            continue
        s = seqs[prot]
        d = DET.detect(s)
        p = d.get("period") or 0
        if not p:
            fams[prot] = {
                "protein": prot,
                "strain": strain,
                "comp_class": comp,
                "src": src,
                "period": 0,
                "units": [],
                "cands": [],
                "reason": "no_period",
            }
            continue
        units = [s[i : i + p] for i in range(d["start"], d["end"] - p + 1, p)]
        cons = A.kmer_conservation(units, args.kmin, args.kmax)
        pool = []
        for dd in cons.values():
            pool.extend((v, km) for km, v in dd.items() if v >= args.min_conservation)
        pool.sort(reverse=True)
        fams[prot] = {
            "protein": prot,
            "strain": strain,
            "comp_class": comp,
            "src": src,
            "period": p,
            "units": units,
            "cons": cons,
            "cands": [km for _, km in pool[: args.per_family]],
            "reason": "",
        }
    print(
        f"{len(fams)} families, "
        f"{sum(1 for f in fams.values() if not f['cands'])} with no candidate motif",
        flush=True,
    )

    # --- phase 2: one background pass ----------------------------------------------------
    allmot = sorted({m for f in fams.values() for m in f["cands"]})
    print(f"counting {len(allmot)} motifs in {len(seqs)} proteins + shuffled null", flush=True)
    bg = list(seqs.values())
    obs = A.count_motifs(bg, allmot)
    print("  real done", flush=True)
    null = A.count_motifs(A.shuffled_copy(bg), allmot)
    print("  null done", flush=True)

    # --- phase 3: pick one anchor per family ---------------------------------------------
    for f in fams.values():
        best = None
        for m in f["cands"]:
            if null[m][0] > args.max_null:
                continue
            c = f["cons"][len(m)][m]
            key = (c, -null[m][0], -len(m))
            if best is None or key > best[0]:
                best = (key, m)
        if best is None:
            f["anchor"] = ""
            f["reason"] = f["reason"] or "no_specific_motif"
        else:
            f["anchor"] = best[1]
            f["conservation"] = round(f["cons"][len(best[1])][best[1]], 3)
            f["obs_prot"] = obs[best[1]][0]
            f["null_prot"] = null[best[1]][0]

    usable = {k: f for k, f in fams.items() if f.get("anchor")}
    print(
        f"{len(usable)}/{len(fams)} families have a specific anchor "
        f"(chance hits <= {args.max_null})",
        flush=True,
    )

    # collapse families that selected the same anchor: they are one family
    by_anchor = defaultdict(list)
    for f in usable.values():
        by_anchor[f["anchor"]].append(f)
    reps = {}
    for anchor, fl in by_anchor.items():
        fl.sort(key=lambda f: -len(seqs[f["protein"]]))
        reps[anchor] = fl[0]
        reps[anchor]["members"] = [f["protein"] for f in fl]
    print(f"{len(reps)} distinct anchors after collapsing duplicates", flush=True)

    # --- phase 4: one hit-collection pass, then per-family alignment ----------------------
    anchors = sorted(reps)
    hits_by_anchor = defaultdict(list)
    for lab, d in per_proteome.items():
        for name, s in d.items():
            for m in anchors:
                if m in s:
                    hits_by_anchor[m].append((lab, name, s))

    aligner = A.make_aligner()
    called02 = detector_calls(args.profiles_02)
    called14 = detector_calls(args.profiles_14)
    print(f"02 calls {len(called02)} proteins, 14 calls {len(called14)}", flush=True)

    fam_rows, hit_rows = [], []
    for anchor in anchors:
        f = reps[anchor]
        ref = seqs[f["protein"]]
        rows_by_prot = defaultdict(list)
        for lab, name, s in hits_by_anchor[anchor]:
            pos = [i for i in range(len(s)) if s.startswith(anchor, i)]
            r = {
                "family": f["protein"],
                "proteome": lab,
                "protein": name,
                "length": len(s),
                "n_anchor": len(pos),
                "anchor_pos": ",".join(str(x + 1) for x in pos),
            }
            r.update(A.align_to_ref(aligner, s, ref))
            rows_by_prot[lab].append(r)
        fam_hits = []
        for lab, rs in rows_by_prot.items():
            path = next(p for p in proteomes if A.label_of(p) == lab)
            genes = A.read_gff_genes(A.gff_for(path), {A.gene_id(r["protein"]) for r in rs})
            A.classify(rs, len(ref), f["period"], genes)
            fam_hits.extend(rs)
        hit_rows.extend(fam_hits)
        ids = {r["protein"] for r in fam_hits}
        new02 = ids - called02
        new14 = ids - called14
        calls = Counter(r["call"] for r in fam_hits)
        fam_rows.append(
            {
                "family": f["protein"],
                "strain": f["strain"],
                "comp_class": f["comp_class"],
                "source": f["src"],
                "ref_len": len(ref),
                "period": f["period"],
                "anchor": anchor,
                "anchor_k": len(anchor),
                "conservation": f["conservation"],
                "obs_prot": f["obs_prot"],
                "null_prot": f["null_prot"],
                "n_members_collapsed": len(f["members"]),
                "n_hits": len(fam_hits),
                "n_proteomes": len({r["proteome"] for r in fam_hits}),
                "n_missed_by_02": len(new02),
                "n_missed_by_14": len(new14),
                "missed_by_both": len(new02 & new14),
                "split_model": calls.get("split_model", 0),
                "split_candidate": calls.get("split_candidate", 0),
                "truncated": calls.get("truncated", 0),
                "short_allele": calls.get("short_allele", 0),
                "full_length": calls.get("full_length", 0),
                "fragment": calls.get("fragment", 0),
            }
        )

    for f in fams.values():
        if f.get("anchor"):
            continue
        fam_rows.append(
            {
                "family": f["protein"],
                "strain": f["strain"],
                "comp_class": f["comp_class"],
                "source": f["src"],
                "ref_len": len(seqs.get(f["protein"], "")),
                "period": f["period"],
                "anchor": "",
                "anchor_k": 0,
                "conservation": "",
                "obs_prot": "",
                "null_prot": "",
                "n_members_collapsed": 1,
                "n_hits": 0,
                "n_proteomes": 0,
                "n_missed_by_02": 0,
                "n_missed_by_14": 0,
                "missed_by_both": 0,
                "split_model": 0,
                "split_candidate": 0,
                "truncated": 0,
                "short_allele": 0,
                "full_length": 0,
                "fragment": 0,
                "fail_reason": f["reason"],
            }
        )

    ffields = [
        "family",
        "strain",
        "comp_class",
        "source",
        "ref_len",
        "period",
        "anchor",
        "anchor_k",
        "conservation",
        "obs_prot",
        "null_prot",
        "n_members_collapsed",
        "n_hits",
        "n_proteomes",
        "n_missed_by_02",
        "n_missed_by_14",
        "missed_by_both",
        "full_length",
        "short_allele",
        "truncated",
        "fragment",
        "split_model",
        "split_candidate",
        "fail_reason",
    ]
    A.write_tsv(fam_rows, args.out_families, ffields)
    A.write_tsv(hit_rows, args.out_hits, ["family"] + A.FIELDS)
    print(
        f"\nwrote {args.out_families} ({len(fam_rows)} families) and "
        f"{args.out_hits} ({len(hit_rows)} hits)"
    )

    # Count MEMBERS only. A hit that fails the alignment filter is labelled "unrelated" and
    # is not a family member, so counting it as "found by anchoring, missed by the detector"
    # overstates the gain. Before 2026-09-30 these totals were over all hits: they read 52 and
    # 42 where the member-only numbers are 24 and 14. The per-family n_missed_by_* columns are
    # still all-hit counts -- read them with the "unrelated" column beside them.
    members = {r["protein"] for r in hit_rows if r.get("call") != "unrelated"}
    allhits = {r["protein"] for r in hit_rows}
    tot02 = len(members - called02)
    tot14 = len(members - called14)
    print(
        f"\ndistinct anchored hits: {len(allhits)}  "
        f"(members {len(members)}, rejected as unrelated {len(allhits - members)})"
    )
    print(f"  MEMBERS not called by 02: {tot02}")
    print(f"  MEMBERS not called by 14: {tot14}")
    print(
        f"  (all-hit counts, including unrelated: "
        f"{len(allhits - called02)} and {len(allhits - called14)})"
    )
    cc = Counter(r["call"] for r in hit_rows)
    print("  calls: " + ", ".join(f"{k}={v}" for k, v in sorted(cc.items())))
    bad = Counter(f["comp_class"] for f in fams.values() if not f.get("anchor"))
    print("\nfamilies with no specific anchor, by composition class:")
    for k, v in bad.most_common():
        print(f"  {k or '(none)'}: {v}")


if __name__ == "__main__":
    main()
