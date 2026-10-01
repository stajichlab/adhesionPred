#!/usr/bin/env python3
"""SOWgp (class-2a) across the Coccidioides pangenome: seed FASTA + all annotated copies.

Motivation (docs/reports/2026-09-27-cocci-repeat-surface-proteins.md section 4.2): the
periodicity detector found SOWgp at 4 loci (RS CIMG_04613, CiB10637_003943, CiB10992_003451,
Silveira QVM09276.1) but no call at VFC140, Cpos1038, Cpos3700. Absence is NOT established:
gene models can fail on repeat arrays and short-read annotations can truncate them. To
separate real absence from annotation/assembly artifact we need the SOWgp ortholog set across
the 496-genome pangenome, aligned, with copy number and extent measured on every member.

Stage 1 (--stage seed): build the seed multiple-sequence file from the taxa the report lists
plus the published SOWgp58/66/82 alleles (UniProt, Hung et al. 2002).

Stage 2 (--stage extract): find the SOWgp orthogroups in the OrthoFinder pangenome
(results/Cocci_496_OG2_5_5/), pull every annotated member (all strains), and profile each
sequence with the same periodicity detector used for the report (02_repeat_profile.py), so
copy-number variation is measured the same way on all copies.

Outputs (analysis/cocci_repeats/):
  sowgp_seed.fa              stage 1: 8 sequences (5 genome copies + 3 published alleles)
  sowgp_pangenome.fa         stage 2: every member of the SOWgp orthogroups, headers = id|strain
  sowgp_pangenome.tsv        stage 2: per-sequence profile (length, period, copies, coverage)
  sowgp_pangenome_summary.tsv stage 2: per-strain count, length range, copy-number range, OG

Caveats: the pangenome is short-read assemblies (repeats collapse/truncate), so copy numbers
there are a LOWER bound relative to long-read calls. Sequences that fell outside the SOWgp
orthogroup (e.g. the 117 aa Silveira fragment QVM10648.1 in OG0007110) are reported but flagged
separately rather than silently merged.
"""

import argparse
import csv
import importlib.util
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from analysis._common.paths import COCCI_LONGREAD, COCCI_PANGENOME  # noqa: E402

OUTDIR = Path(__file__).resolve().parent

# tracked genome copies from the report's section 4.2 table (+ the 5th family member found by
# analysis/model_review/repeat_structure_transfer_2a_retrain.py)
GENOME_COPIES = [
    # (label, fasta path, record prefix)
    (
        "Cocci_immitis_RS_CIMG_04613",
        COCCI_PANGENOME / "input_run2" / "CimmitisRS_FungiDB.fasta",
        "CIMG_04613-t26_1-p1",
    ),
    (
        "COCCIDIOIDES_IMMITIS_CIB10637_003943",
        COCCI_LONGREAD / "CiB10637" / "Coccidioides_immitis_CiB10637.proteins.fa",
        "CIB10637_003943",
    ),
    (
        "COCCIDIOIDES_IMMITIS_CIB10992_003451",
        COCCI_LONGREAD / "CiB10992" / "Coccidioides_immitis_CiB10992.proteins.fa",
        "CIB10992_003451",
    ),
    (
        "Coccidioides_posadasii_Silveira_QVM09276",
        COCCI_PANGENOME / "input_run2" / "CposadasiiSilveira2022_FungiDB.fasta",
        "QVM09276.1",
    ),
    (
        "Coccidioides_posadasii_Silveira_QVM10648_trunc",
        COCCI_PANGENOME / "input_run2" / "CposadasiiSilveira2022_FungiDB.fasta",
        "QVM10648.1",
    ),
]

PUBLISHED = {
    "Q8NK60": "SOWgp58_Cocci_immitis_published",
    "Q8NK61": "SOWgp66_Cocci_immitis_published",
    "Q96V71": "SOWgp82_Cocci_posadasii_published",
}

ORTHO_DIR = COCCI_PANGENOME / "results" / "Cocci_496_OG2_5_5" / "Orthogroups" / "Orthogroups.tsv"
OG_SEQ_DIR = COCCI_PANGENOME / "results" / "Cocci_496_OG2_5_5" / "Orthogroup_Sequences"


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def fetch_uniprot(accessions):
    seqs = {}
    for i in range(0, len(accessions), 90):
        q = urllib.parse.urlencode(
            {"accessions": ",".join(accessions[i : i + 90]), "format": "fasta"}
        )
        with urllib.request.urlopen(
            f"https://rest.uniprot.org/uniprotkb/accessions?{q}", timeout=60
        ) as r:
            for blk in r.read().decode().split(">"):
                if not blk.strip():
                    continue
                h, *rest = blk.split("\n")
                seqs[h.split("|")[1]] = "".join(rest).replace(" ", "")
    return seqs


def make_seed(out, fetch):
    out = Path(out)
    recs = []
    for label, fa, prefix in GENOME_COPIES:
        for pid, seq in read_fasta(fa):
            if pid.startswith(prefix):
                recs.append((label, pid, seq))
                break
        else:
            print(f"  WARN: {prefix} not found in {fa}", file=sys.stderr)
    # published alleles
    if out.exists() and out.stat().st_size:
        # use cached published sequences if the seed file already exists
        for pid, seq in read_fasta(out):
            if "published" in pid:
                recs.append((pid, pid, seq))
    if fetch:
        seqs = fetch_uniprot(list(PUBLISHED))
        for acc, label in PUBLISHED.items():
            if acc in seqs:
                recs.append((label, acc, seqs[acc]))
    # dedupe (published may be cached and fetched again)
    seen, uniq = set(), []
    for label, pid, seq in recs:
        if label in seen:
            continue
        seen.add(label)
        uniq.append((label, pid, seq))
    with open(out, "w") as fh:
        for label, _pid, seq in uniq:
            fh.write(f">{label} [{pid}] len={len(seq)}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    print(f"wrote {len(uniq)} sequences to {out}")
    for label, _pid, seq in uniq:
        print(f"   {label:<55} {len(seq):>4} aa")


def parse_seed_ids(seed_fa):
    """Return {protein_id: label} from the seed FASTA (header "label [pid] len=N")."""
    return dict(parse_seed_fasta(seed_fa))


def parse_seed_fasta(seed_fa):
    """Return {protein_id: sequence} from the seed FASTA (header "label [pid] len=N")."""
    import re

    pat = re.compile(r"^>([^\[\n]+)\s*\[([^\]]+)\]\s*len=\d+")
    seqs, cur = {}, None
    with open(seed_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                mo = pat.match(line.strip())
                cur = mo.group(2) if mo else line.strip()[1:]
                seqs.setdefault(cur, [])
            else:
                seqs[cur].append(line.strip())
    return {pid: "".join(b) for pid, b in seqs.items() if b}


def load_seed_ogs(orotho_tsv, seed_fa):
    """Return {seed_id: [OGs]} using protein IDs extracted from the seed FASTA.

    Matches protein IDs against the OT column of Orthogroups.tsv rows; moving a seed to a
    different OG (e.g. after re-running OrthoFinder) only requires re-running this stage.
    """
    seeds = parse_seed_ids(seed_fa)
    ogs = defaultdict(list)
    with open(orotho_tsv) as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        for row in r:
            for s in seeds:
                if s in row:
                    ogs[s].append(row[0])
    return header, seeds, ogs


def sowgp_identity(member_seq, seed_seqs):
    """Best local identity (%) of a member against any confirmed SOWgp seed.

    Requires the local alignment span >= MEANINGFUL_SPAN (default 50 aa) so that short
    low-complexity near-matches (identical 6-10 aa stretches, e.g. the period-16 NWIKRED
    family vs SOWgp58) are not counted as orthology.
    """
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    best = 0.0
    for seed in seed_seqs:
        aln = aligner.align(member_seq, seed)
        a = aln[0].aligned
        ids = 0
        amt = 0
        for (st, en), (qt, _qe) in zip(a[0], a[1], strict=False):
            ids += sum(1 for i in range(en - st) if member_seq[st + i] == seed[qt + i])
            amt += en - st
        if amt >= 50:
            best = max(best, 100 * ids / amt)
    return best


def extract_pangenome(orotho_tsv, og_seq_dir, og_ids, seed_fa, outfa, outtsv, outsum):
    assert og_seq_dir.exists(), og_seq_dir
    # confirmed SOWgp references = the 4 full-length genome copies + published 58/66/82;
    # the 117 aa period-16 Silveira seed (QVM10648.1) is a DIFFERENT family (see docstring)
    seed_seqs = parse_seed_fasta(seed_fa)
    ref_seqs = [s for pid, s in seed_seqs.items() if len(s) >= 300]

    og_members = {}  # og -> {strain: [gene, ...]} (cells may hold comma-separated paralogs)
    with open(orotho_tsv) as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        ncols = len(header)
        for row in r:
            if row[0] not in og_ids:
                continue
            og_members[row[0]] = {
                header[i]: [g.strip() for g in row[i].split(",") if g.strip()]
                for i in range(1, ncols)
                if row[i].strip()
            }
    # load profiles by reusing the report's detector
    spec = importlib.util.spec_from_file_location(
        "repeat_profile", Path(__file__).with_name("02_repeat_profile.py")
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load 02_repeat_profile.py")
    rp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rp)

    rows = []
    out_other = Path(str(outfa).replace(".fa", "_other_families.fa"))
    with open(outfa, "w") as fh_main, open(out_other, "w") as fh_other:
        for og in sorted(og_ids):
            fa = og_seq_dir / f"{og}.fa"
            if not fa.exists():
                print(f"  WARN missing {fa}", file=sys.stderr)
                continue
            strain2genes = og_members.get(og, {})
            gene2strain = {g: s for s, genes in strain2genes.items() for g in genes}
            n_confirmed = 0
            n_members = 0
            for _gene, seq in read_fasta(fa):
                n_members += 1
                seq = seq.rstrip("*")
                ident = sowgp_identity(seq, ref_seqs)
                if ident >= 35.0:
                    n_confirmed += 1
            family = "SOWgp" if n_confirmed / max(n_members, 1) >= 0.5 else "other"
            print(
                f"  {og}: {n_members} members, {n_confirmed} >=35% identical to SOWgp -> {family}"
            )
            for gene, seq in read_fasta(fa):
                seq = seq.rstrip("*")
                strain = gene2strain.get(gene, f"unknown_in_{og}")
                ident = sowgp_identity(seq, ref_seqs)
                out = fh_main if family == "SOWgp" else fh_other
                out.write(f">{gene}|{strain}|{og}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i : i + 60] + "\n")
                rep = rp.best_repeat(seq)
                rows.append(
                    {
                        "orthogroup": og,
                        "family": family,
                        "strain": strain,
                        "protein": gene,
                        "length": len(seq),
                        "best_ident_sowgp": round(ident, 1),
                        "rep_period": rep["period"],
                        "rep_n_copies": rep["n_copies"],
                        "rep_coverage": rep["coverage"],
                        "rep_unit": rep["unit"],
                        "species": "immitis"
                        if "immitis" in strain
                        else ("posadasii" if "posadasii" in strain else "other"),
                    }
                )
    n_seq = len(rows)
    with open(outtsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    # per-strain summary
    per = defaultdict(list)
    for row in rows:
        per[row["strain"]].append(row)
    with open(outsum, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "strain",
                "species",
                "family",
                "n_copies",
                "length_range",
                "copies_range",
                "max_cov",
                "orthogroups",
            ]
        )
        for strain, rs in sorted(per.items()):
            lens = [r["length"] for r in rs]
            cps = [r["rep_n_copies"] for r in rs]
            w.writerow(
                [
                    strain,
                    rs[0]["species"],
                    rs[0]["family"],
                    len(rs),
                    f"{min(lens)}-{max(lens)}",
                    f"{min(cps)}-{max(cps)}",
                    max(r["rep_coverage"] for r in rs),
                    ",".join(sorted({r["orthogroup"] for r in rs})),
                ]
            )
    strains = len(per)
    print(f"wrote {n_seq} sequences ({strains} strains) from {sorted(og_ids)}")
    if outtsv:
        print(outtsv)


def main():
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--stage", choices=["seed", "extract"], required=True)
    ap.add_argument(
        "--fetch-uniprot",
        action="store_true",
        help="seed: fetch SOWgp58/66/82 from UniProt REST (cached in sowgp_seed.fa)",
    )
    ap.add_argument(
        "--og-input",
        help="extract: path to a file with orthogroup IDs (default: discover from seeds)",
    )
    args = ap.parse_args()

    seed = OUTDIR / "sowgp_seed.fa"
    if args.stage == "seed":
        make_seed(seed, args.fetch_uniprot)
        return

    if args.og_input:
        og_ids = [line.strip() for line in open(args.og_input) if line.strip()]
    else:
        if not seed.exists():
            make_seed(seed, True)
        header, seed_labels, ogs = load_seed_ogs(ORTHO_DIR, seed)
        print("found OGs for seeds:")
        for s, o in sorted(ogs.items()):
            print(f"   {s:<45} -> {o}")
        og_ids = sorted({og for ogs_ in ogs.values() for og in ogs_})
    print("using orthogroups:", og_ids)
    extract_pangenome(
        ORTHO_DIR,
        OG_SEQ_DIR,
        og_ids,
        seed,
        OUTDIR / "sowgp_pangenome.fa",
        OUTDIR / "sowgp_pangenome.tsv",
        OUTDIR / "sowgp_pangenome_summary.tsv",
    )


if __name__ == "__main__":
    main()
