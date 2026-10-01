#!/usr/bin/python3.12
"""Motif-anchored search for members of a tandem-repeat protein family.

Generalises 20_sowgp_anchor_search.py. That script hardcoded the SOWgp anchor PTDCYGDC,
which was picked by eye. This one derives the anchor from the family's repeat units,
measures how specific the anchor is against a real null, searches any set of proteomes,
and classifies each hit against the failure modes that hid three SOWgp models
(see the 2026-09-30 correction in REPORT_2026-09-29_sowgp_repeat_structure.md).

WHY A MOTIF ANCHOR AND NOT A REPEAT DETECTOR
--------------------------------------------
02_repeat_profile.py and 14_repeat_detect_general.py both need a detectable periodicity.
A family member loses its periodicity when the gene model is truncated, when one gene is
called as two, or when the allele carries too few units to clear a copy-number threshold.
An anchored search needs ONE occurrence of a short conserved motif, so it survives all three.
It cannot replace a detector: it only finds members of a family you already have units for.

WHAT IT REPORTS
---------------
Per hit: protein id, length, anchor count and positions, the span of the reference the hit
aligns to, identity, and a call:
  full_length   ref span is (nearly) complete and the alignment has no large internal gap
  short_allele  ref span is complete but the query is missing an internal block, sized as a
                whole number of repeat units. This is NOT an error (VFC140 is one).
  split_model   two or more hits in one proteome whose ref spans overlap or abut, whose
                locus tags are consecutive, and which the GFF places adjacent on one contig
                and strand. Reports the intergenic gap in bp.
  truncated     partial ref span with no partner
  fragment      partial ref span, short, no partner, below --min-frac of the reference

USAGE
-----
  # 1. derive an anchor from a family's units and measure its specificity
  ./30_anchor_family_search.py derive --seed sowgp_seed.fa --proteomes longread --out-tsv m.tsv

  # 2. search with the chosen anchor
  ./30_anchor_family_search.py search --motif PTDCYGDC --reference sowgp_seed.fa \
      --reference-id SOWgp58_Cocci_immitis_published --period 47 --proteomes longread \
      --out hits.tsv

  # both at once
  ./30_anchor_family_search.py family --seed sowgp_seed.fa --proteomes longread --out hits.tsv

--proteomes takes `longread` (the 7 proteomes 01/02 use), `pangenome` (the 493-strain
Pangenome/input set), or an explicit list of FASTA paths.

Needs biopython and numpy. The 7-proteome run takes about a minute; the 493-proteome run is
a SLURM job (31_anchor_pangenome.sh).
"""

import argparse
import gzip
import json
import random
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "_common"))
sys.path.insert(0, str(HERE))

# --- defaults ---------------------------------------------------------------------------

KMIN, KMAX = 5, 12
# Fraction of units that must contain the motif. Not 1.0 and not 0.9: the terminal unit of a
# repeat array usually runs into the C-terminus and shares little with the internal units, and
# the first and last tiled units are often partial. For SOWgp the ceiling over all k-mers is
# 0.79, so a stricter cut returns nothing. derive_anchors() relaxes to the observed ceiling
# rather than returning an empty list.
MIN_CONSERVATION = 0.75
SHUFFLE_SEED = 20260930
FULL_FRAC = 0.85  # ref coverage at or above this counts as a complete span
MIN_FRAC = 0.25  # below this a partial hit is called a fragment, not a truncation
DEL_UNIT_FRAC = 0.6  # an internal query deletion of >= this many periods marks a short allele
SPLIT_REF_SLOP = 40  # aa: ref spans this close count as abutting
SPLIT_TAG_SLOP = 3  # locus-tag numbers this far apart still count as consecutive
SPLIT_BP_SLOP = 20000  # bp: genomic distance beyond which a GFF-confirmed pair is rejected
# A motif that is rare in a per-protein-residue shuffle can still occur in UNRELATED real
# proteins: the shuffle destroys homology and convergent sequence, so it under-states the
# false-positive rate for an ordinary-composition motif. Measured here: the anchor LAAKIS has
# zero chance hits and still lands in eleven 1589 aa proteins that align to the reference at
# 44.8% over 11% of it. The alignment, not the anchor, is what rejects those. A hit is a
# family member only when it also passes these.
MIN_IDENTITY = 60.0  # percent over the aligned ungapped columns
MIN_ALN_LEN = 30  # aa


# --- io ---------------------------------------------------------------------------------


def _open(path):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return {k: v.rstrip("*") for k, v in seqs.items()}


def resolve_proteomes(spec):
    """`longread`, `pangenome`, or explicit paths."""
    import paths  # noqa: E402  -- from ../_common

    if len(spec) == 1 and spec[0] == "longread":
        lr = Path(paths.COCCI_LONGREAD)
        pan = Path(paths.COCCI_PANGENOME) / "input_run2"
        out = sorted(lr.glob("*/*.proteins.fa"))
        for f in ("CimmitisRS_FungiDB.fasta", "CposadasiiSilveira2022_FungiDB.fasta"):
            if (pan / f).exists():
                out.append(pan / f)
        return out
    if len(spec) == 1 and spec[0] == "pangenome":
        return sorted((Path(paths.COCCI_PANGENOME) / "input").glob("*.proteins.fa"))
    return [Path(p) for p in spec]


def label_of(path):
    return Path(path).name.replace(".proteins.fa", "").replace(".fasta", "").replace(".fa", "")


def gff_for(proteome_path):
    """The .gff3 beside the proteome FASTA (both datasets keep them in one directory)."""
    p = Path(proteome_path).resolve()
    for cand in (
        p.parent / (p.name.replace(".proteins.fa", ".gff3")),
        p.parent / (p.stem + ".gff3"),
    ):
        if cand.exists():
            return cand
    return None


def read_gff_genes(path, wanted):
    """{gene_id: (contig, start, end, strand)} for the wanted gene ids only."""
    out = {}
    if path is None:
        return out
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"ID=([^;]+)", f[8])
            if not m or m.group(1) not in wanted:
                continue
            out[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
    return out


# --- units ------------------------------------------------------------------------------


def units_from_seed(seed_seqs, mode="similarity"):
    """Repeat units of every seed sequence, from 14_repeat_detect_general.py.

    Units are tiled from the detector's called region start at its called period. The
    detector searches the tiling phase, so the region start is a unit boundary when it
    finds one at all.
    """
    import importlib

    det = importlib.import_module("14_repeat_detect_general")
    units, periods = [], []
    for s in seed_seqs.values():
        d = det.detect(s, mode=mode)
        p = d.get("period") or 0
        if not p:
            continue
        a, b = d["start"], d["end"]
        periods.append(p)
        for i in range(a, b - p + 1, p):
            units.append(s[i : i + p])
    period = int(Counter(periods).most_common(1)[0][0]) if periods else 0
    return units, period


# --- anchor derivation ------------------------------------------------------------------


def kmer_conservation(units, kmin=KMIN, kmax=KMAX):
    """{k: {kmer: fraction of units containing it}}"""
    n = len(units)
    out = {}
    for k in range(kmin, kmax + 1):
        c = Counter()
        for u in units:
            for km in {u[i : i + k] for i in range(len(u) - k + 1)}:
                if "X" in km or "*" in km:
                    continue
                c[km] += 1
        out[k] = {km: v / n for km, v in c.items() if v / n >= 0.5}
    return out


def shuffled_copy(seqs, seed=SHUFFLE_SEED):
    """Per-protein residue shuffle: keeps each protein's own composition and length.

    This is the null the specificity number is measured against. It is stricter than a
    whole-proteome shuffle: a Pro/Cys-rich protein stays Pro/Cys-rich, so a motif that is
    only a composition artifact will hit the shuffled copy about as often as the real one.
    """
    rng = random.Random(seed)
    out = []
    for s in seqs:
        ls = list(s)
        rng.shuffle(ls)
        out.append("".join(ls))
    return out


def count_motifs(seqs, motifs):
    """{motif: (n_proteins containing it, n_total occurrences)}"""
    npro = Counter()
    nocc = Counter()
    for s in seqs:
        for m in motifs:
            c = s.count(m)
            if c:
                npro[m] += 1
                nocc[m] += c
    return {m: (npro[m], nocc[m]) for m in motifs}


def derive_anchors(
    units,
    background_seqs,
    kmin=KMIN,
    kmax=KMAX,
    min_conservation=MIN_CONSERVATION,
    top=25,
    per_k=40,
):
    """Rank candidate anchors by measured specificity, not by eye.

    For each candidate k-mer we report:
      conservation  fraction of the family's repeat units that contain it
      obs_prot      proteins in the background set that contain it (family members included)
      null_prot     proteins in a per-protein-shuffled copy of the same set that contain it
    null_prot is the chance rate. A motif with null_prot in the hundreds is useless however
    conserved it is.
    """
    cons = kmer_conservation(units, kmin, kmax)
    ceiling = max((v for d in cons.values() for v in d.values()), default=0.0)
    cut = min_conservation
    if ceiling < min_conservation:
        cut = 0.95 * ceiling
        print(
            f"note: no k-mer reaches conservation {min_conservation}; ceiling is "
            f"{ceiling:.2f}, using {cut:.2f}",
            file=sys.stderr,
        )
    cands = []
    for d in cons.values():
        keep = sorted((v, km) for km, v in d.items() if v >= cut)
        keep = [km for _, km in sorted(keep, reverse=True)[:per_k]]
        cands.extend(keep)
    cands = sorted(set(cands))
    if not cands:
        return []
    obs = count_motifs(background_seqs, cands)
    null = count_motifs(shuffled_copy(background_seqs), cands)
    rows = []
    for m in cands:
        k = len(m)
        rows.append(
            {
                "motif": m,
                "k": k,
                "conservation": round(cons[k][m], 3),
                "obs_prot": obs[m][0],
                "obs_occ": obs[m][1],
                "null_prot": null[m][0],
                "null_occ": null[m][1],
            }
        )
    # rank: fewest chance hits first, then most conserved, then shortest
    rows.sort(key=lambda r: (r["null_prot"], -r["conservation"], r["k"]))
    return rows[:top]


# --- alignment --------------------------------------------------------------------------


def make_aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices

    al = Align.PairwiseAligner()
    al.mode = "local"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -11
    al.extend_gap_score = -1
    return al


def align_to_ref(aligner, query, ref):
    """Local alignment stats in reference coordinates.

    Returns ref_start (1-based), ref_end, identity%, ref coverage, query coverage, and
    max_query_del: the largest block of the reference that the query skips over with no
    matching insertion. A short allele shows a large max_query_del sized as whole repeat
    units; a truncated model shows a short ref span instead.
    """
    a = aligner.align(query, ref)[0]
    qb, rb = a.aligned[0], a.aligned[1]
    rs, re_ = int(rb[0][0]), int(rb[-1][1])
    qs, qe = int(qb[0][0]), int(qb[-1][1])
    ident = sum(1 for x, y in zip(a[0], a[1], strict=False) if x == y and x != "-")
    alen = sum(1 for x, y in zip(a[0], a[1], strict=False) if x != "-" and y != "-")
    max_del = 0
    for i in range(len(rb) - 1):
        rgap = int(rb[i + 1][0] - rb[i][1])
        qgap = int(qb[i + 1][0] - qb[i][1])
        if rgap - qgap > max_del:
            max_del = rgap - qgap
    return {
        "ref_start": rs + 1,
        "ref_end": re_,
        "ref_cov": round((re_ - rs) / len(ref), 3),
        "query_start": qs + 1,
        "query_end": qe,
        "query_cov": round((qe - qs) / max(len(query), 1), 3),
        "identity_pct": round(100 * ident / max(alen, 1), 1),
        "aln_len": alen,
        "max_query_del": max_del,
    }


# --- classification ---------------------------------------------------------------------


TAG_RE = re.compile(r"^(.*?)(\d+)(?:-T\d+)?$")


def tag_parts(protein_id):
    """('CPOS1038_', 3233) from 'CPOS1038_003233-T1'. None when the id has no number."""
    base = protein_id.split("-T")[0]
    m = TAG_RE.match(base)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


def gene_id(protein_id):
    return protein_id.split("-T")[0]


def spans_overlap_or_abut(a, b, slop=SPLIT_REF_SLOP):
    lo, hi = (a, b) if a["ref_start"] <= b["ref_start"] else (b, a)
    return lo["ref_end"] + slop >= hi["ref_start"]


def classify(hits, ref_len, period, genes, min_identity=MIN_IDENTITY):
    """Set `call` and split-pair fields on every hit of one proteome, in place.

    genes: {gene_id: (contig, start, end, strand)} from the GFF, may be empty.
    """
    for h in hits:
        h["call"] = ""
        h["partner"] = ""
        h["pair_ref_cov"] = ""
        h["gff_contig"] = ""
        h["gff_gap_bp"] = ""
        h["collinear"] = ""
        h["units_missing"] = ""
        h["member"] = h["identity_pct"] >= min_identity and h["aln_len"] >= MIN_ALN_LEN

    # pass 1: split pairs.
    #
    # Two tiers, because the evidence differs in kind:
    #   split_model      the two ref spans OVERLAP or abut. One gene called as two, and the
    #                    duplicated span is direct evidence of it.
    #   split_candidate  the two ref spans leave a hole. Consistent with one gene called as
    #                    two with the middle lost, but NOT separable from two independently
    #                    truncated models at one degraded locus. Reported, not concluded.
    partial = [h for h in hits if h["ref_cov"] < FULL_FRAC and h["member"]]
    used = set()
    for i, A in enumerate(partial):
        if A["protein"] in used:
            continue
        pa, na = tag_parts(A["protein"])
        if pa is None:
            continue
        for B in partial[i + 1 :]:
            if B["protein"] in used:
                continue
            pb, nb = tag_parts(B["protein"])
            if pb != pa or nb is None or abs(nb - na) > SPLIT_TAG_SLOP:
                continue
            union = (
                max(A["ref_end"], B["ref_end"]) - min(A["ref_start"], B["ref_start"]) + 1
            ) / ref_len
            if union < FULL_FRAC:
                continue
            ga, gb = genes.get(gene_id(A["protein"])), genes.get(gene_id(B["protein"]))
            gap, contig, collinear = "", "", ""
            if ga and gb:
                if ga[0] != gb[0] or ga[3] != gb[3]:
                    continue  # different contig or strand: not one gene called as two
                # intergenic bases strictly between the two gene features. The 2026-09-30
                # addendum quotes 81 bp for the Cpos1038 pair; that is the same gap counted
                # as next_start - prev_end, i.e. inclusive of one endpoint.
                gap = max(ga[1], gb[1]) - min(ga[2], gb[2]) - 1
                if gap > SPLIT_BP_SLOP:
                    continue
                contig = ga[0]
                # does genomic order match reference order? On the minus strand the gene
                # with the higher coordinate comes first in the transcript.
                first = A if (ga[1] < gb[1]) == (ga[3] == "+") else B
                other = B if first is A else A
                collinear = "yes" if first["ref_start"] <= other["ref_start"] else "no"
            elif not spans_overlap_or_abut(A, B):
                continue  # no GFF and no span overlap: not enough to call anything
            call = "split_model" if spans_overlap_or_abut(A, B) else "split_candidate"
            for X, Y in ((A, B), (B, A)):
                X["call"] = call
                X["partner"] = Y["protein"]
                X["pair_ref_cov"] = round(union, 3)
                X["gff_contig"] = contig
                X["gff_gap_bp"] = gap
                X["collinear"] = collinear
            used.update({A["protein"], B["protein"]})
            break

    # pass 2: everything else
    for h in hits:
        if h["call"]:
            continue
        if not h["member"]:
            h["call"] = "unrelated"
            continue
        if h["ref_cov"] >= FULL_FRAC:
            if period and h["max_query_del"] >= DEL_UNIT_FRAC * period:
                h["call"] = "short_allele"
                h["units_missing"] = round(h["max_query_del"] / period, 1)
            else:
                h["call"] = "full_length"
        elif h["ref_cov"] >= MIN_FRAC:
            h["call"] = "truncated"
        else:
            h["call"] = "fragment"
    return hits


# --- search -----------------------------------------------------------------------------


FIELDS = [
    "proteome",
    "protein",
    "length",
    "n_anchor",
    "anchor_pos",
    "ref_start",
    "ref_end",
    "ref_cov",
    "query_cov",
    "identity_pct",
    "aln_len",
    "max_query_del",
    "member",
    "call",
    "units_missing",
    "partner",
    "pair_ref_cov",
    "gff_contig",
    "gff_gap_bp",
    "collinear",
]


def search(proteomes, motifs, ref, period, use_gff=True, progress=False, min_identity=MIN_IDENTITY):
    aligner = make_aligner()
    ref_len = len(ref)
    all_rows = []
    for pi, path in enumerate(proteomes):
        label = label_of(path)
        seqs = read_fasta(path)
        hits = []
        for name, s in seqs.items():
            pos = []
            for m in motifs:
                pos.extend(i for i in range(len(s)) if s.startswith(m, i))
            if not pos:
                continue
            pos = sorted(set(pos))
            row = {
                "proteome": label,
                "protein": name,
                "length": len(s),
                "n_anchor": len(pos),
                "anchor_pos": ",".join(str(p + 1) for p in pos),
            }
            row.update(align_to_ref(aligner, s, ref))
            hits.append(row)
        if hits:
            genes = {}
            if use_gff:
                genes = read_gff_genes(gff_for(path), {gene_id(h["protein"]) for h in hits})
            classify(hits, ref_len, period, genes, min_identity)
            all_rows.extend(sorted(hits, key=lambda r: (r["ref_start"], r["protein"])))
        if progress and (pi + 1) % 25 == 0:
            print(
                f"  {pi + 1}/{len(proteomes)} proteomes, {len(all_rows)} hits",
                file=sys.stderr,
                flush=True,
            )
    return all_rows


def write_tsv(rows, path, fields=FIELDS):
    with open(path, "w") as fh:
        fh.write("\t".join(fields) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(f, "")) for f in fields) + "\n")


def print_hits(rows):
    hdr = (
        f"{'proteome':<34} {'protein':<26} {'len':>5} {'anch':>4} "
        f"{'ref span':>12} {'cov':>5} {'id%':>6} {'call':<13} {'note'}"
    )
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        note = ""
        if r["call"].startswith("split"):
            note = (
                f"partner {r['partner']}, gap {r['gff_gap_bp']} bp, "
                f"pair cov {r['pair_ref_cov']}, collinear {r['collinear'] or '?'}"
            )
        elif r["call"] == "short_allele":
            note = f"internal deletion {r['max_query_del']} aa = {r['units_missing']} units"
        span = f"{r['ref_start']}-{r['ref_end']}"
        print(
            f"{r['proteome']:<34} {r['protein']:<26} {r['length']:>5} {r['n_anchor']:>4} "
            f"{span:>12} {r['ref_cov']:>5} {r['identity_pct']:>6} {r['call']:<13} {note}"
        )


# --- commands ---------------------------------------------------------------------------


def cmd_derive(args):
    seed = read_fasta(args.seed)
    if args.units:
        units = list(read_fasta(args.units).values())
        period = int(np.median([len(u) for u in units])) if units else 0
    else:
        units, period = units_from_seed(seed)
    if not units:
        sys.exit("no repeat units found in the seed set -- supply --units instead")
    print(f"{len(units)} units from {len(seed)} seed sequences, modal period {period} aa")

    proteomes = resolve_proteomes(args.proteomes)
    bg = []
    for p in proteomes:
        bg.extend(read_fasta(p).values())
    print(f"background: {len(bg)} proteins from {len(proteomes)} proteomes\n")

    rows = derive_anchors(units, bg, args.kmin, args.kmax, args.min_conservation, args.top)
    hdr = f"{'motif':<14} {'k':>2} {'cons':>5} {'obs_prot':>9} {'obs_occ':>8} {'null_prot':>10} {'null_occ':>9}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(
            f"{r['motif']:<14} {r['k']:>2} {r['conservation']:>5} {r['obs_prot']:>9} "
            f"{r['obs_occ']:>8} {r['null_prot']:>10} {r['null_occ']:>9}"
        )
    print("\nobs_prot  proteins in the real background containing the motif (family included)")
    print("null_prot same count on a per-protein residue shuffle = the chance rate")
    print("A motif is usable when null_prot is ~0 and obs_prot is close to the family size.")
    if args.out_tsv:
        write_tsv(
            rows,
            args.out_tsv,
            ["motif", "k", "conservation", "obs_prot", "obs_occ", "null_prot", "null_occ"],
        )
        print(f"\nwrote {args.out_tsv}")
    return rows, period


def _reference(args, seed):
    if args.reference_id:
        if args.reference_id not in seed:
            sys.exit(f"{args.reference_id} not in the reference FASTA")
        return args.reference_id, seed[args.reference_id]
    name = max(seed, key=lambda k: len(seed[k]))
    return name, seed[name]


def cmd_search(args, motifs=None, period=None):
    seed = read_fasta(args.reference or args.seed)
    rname, ref = _reference(args, seed)
    motifs = motifs or args.motif
    period = period or args.period or 0
    proteomes = resolve_proteomes(args.proteomes)
    print(f"anchor(s): {', '.join(motifs)}")
    print(f"reference: {rname}, {len(ref)} aa, unit period {period} aa")
    print(f"proteomes: {len(proteomes)}\n")
    rows = search(
        proteomes,
        motifs,
        ref,
        period,
        use_gff=not args.no_gff,
        progress=len(proteomes) > 25,
        min_identity=getattr(args, "min_identity", MIN_IDENTITY),
    )
    print_hits(rows)
    n = Counter(r["call"] for r in rows)
    print(
        f"\nfamily members (identity >= {getattr(args, 'min_identity', MIN_IDENTITY)}%): "
        f"{sum(1 for r in rows if r['member'])}/{len(rows)} hits"
    )
    print("\ncalls: " + ", ".join(f"{k}={v}" for k, v in sorted(n.items())))
    print(f"proteomes with >=1 hit: {len({r['proteome'] for r in rows})}/{len(proteomes)}")
    if args.out:
        write_tsv(rows, args.out)
        print(f"wrote {args.out}")
    if args.out_json:
        Path(args.out_json).write_text(
            json.dumps(
                {
                    "motifs": motifs,
                    "reference": rname,
                    "ref_len": len(ref),
                    "period": period,
                    "calls": dict(n),
                },
                indent=2,
            )
        )
    return rows


def cmd_family(args):
    rows, period = cmd_derive(args)
    if not rows:
        sys.exit("no candidate anchor passed the conservation cut")
    best = rows[0]["motif"]
    print(
        f"\nchosen anchor: {best} "
        f"(conservation {rows[0]['conservation']}, chance hits {rows[0]['null_prot']})\n"
    )
    args.reference = args.reference or args.seed
    return cmd_search(args, motifs=[best], period=period)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument(
            "--proteomes",
            nargs="+",
            default=["longread"],
            help="`longread`, `pangenome`, or FASTA paths",
        )
        p.add_argument("--out", default=None)
        p.add_argument("--out-json", default=None)
        p.add_argument(
            "--min-identity",
            type=float,
            default=MIN_IDENTITY,
            help="a hit below this identity is called `unrelated`, not a member",
        )
        p.add_argument(
            "--no-gff", action="store_true", help="skip GFF adjacency confirmation for split calls"
        )

    d = sub.add_parser("derive", help="derive and score candidate anchors")
    d.add_argument("--seed", required=True, help="FASTA of family members")
    d.add_argument("--units", default=None, help="FASTA of repeat units (skips detection)")
    d.add_argument("--kmin", type=int, default=KMIN)
    d.add_argument("--kmax", type=int, default=KMAX)
    d.add_argument("--min-conservation", type=float, default=MIN_CONSERVATION)
    d.add_argument("--top", type=int, default=25)
    d.add_argument("--out-tsv", default=None)
    d.add_argument("--proteomes", nargs="+", default=["longread"])

    s = sub.add_parser("search", help="search proteomes for an anchor")
    s.add_argument("--motif", nargs="+", required=True)
    s.add_argument("--reference", required=True, help="FASTA holding the reference sequence")
    s.add_argument("--reference-id", default=None, help="default: the longest sequence")
    s.add_argument("--period", type=int, default=0, help="repeat unit length in aa")
    s.add_argument("--seed", default=None)
    common(s)

    f = sub.add_parser("family", help="derive an anchor then search with it")
    f.add_argument("--seed", required=True)
    f.add_argument("--units", default=None)
    f.add_argument("--reference", default=None, help="default: --seed")
    f.add_argument("--reference-id", default=None)
    f.add_argument("--period", type=int, default=0)
    f.add_argument("--kmin", type=int, default=KMIN)
    f.add_argument("--kmax", type=int, default=KMAX)
    f.add_argument("--min-conservation", type=float, default=MIN_CONSERVATION)
    f.add_argument("--top", type=int, default=25)
    f.add_argument("--out-tsv", default=None)
    common(f)

    args = ap.parse_args()
    {"derive": cmd_derive, "search": cmd_search, "family": cmd_family}[args.cmd](args)


if __name__ == "__main__":
    main()
