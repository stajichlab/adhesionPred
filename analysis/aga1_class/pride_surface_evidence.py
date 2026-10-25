"""Conidial-surface proteomics evidence (PRIDE PXD031199, PXD044190) joined to the AGA1-class pipeline calls (pandas, stdlib XML).

PXD031199 (Pinzan et al. 2024, PMID 39191887): trypsin shaving of resting and swollen conidia of A. fumigatus, A. fischeri,
A. oerlinghausenensis and A. lentulus, searched with Proteome Discoverer against FungiDB (Af293 ids `Afu6g07910-T-p1`). Runs 1-3 are
Af293, 4-6 A. oerlinghausenensis, 7-9 A. fischeri, 10-12 A. lentulus (from the reference database named in each file); the `s` files are the
second condition of the same replicate. The paper describes resting and swollen conidia; mapping the `s` suffix to "swollen" is an
INFERENCE from the file names and is reported as `run_suffix`, not as a stage. This script reads only the protein lists
(`*.prot.xml`), not spectra, and the Af293 runs are the ones that map to our Af293 proteome (Afu6g07910 -> AFUA_6G07910).
PXD044190 is a MaxQuant search of A. fumigatus A1163 (AFUB_ ids) plus mouse in an infection experiment (samples MO_*); its ids do not map
to Af293 directly; the `maxquant` command maps its A. fumigatus protein groups (about 5% of the groups, the rest is mouse) to the Af293 proteome by the
best shared 10-mers of the two near-identical proteomes (the best Af293 match must hold at least 80% of the A1163 protein's 10-mers and be unique
at that score) and writes only accession-level numbers; it needs the two FASTA files, which are never written anywhere.

  python3 -I analysis/aga1_class/pride_surface_evidence.py detections --prot-dir DIR --out runs.tsv
  python3 -I analysis/aga1_class/pride_surface_evidence.py maxquant --groups proteinGroups.txt --a1163-names N --a1163-fasta F1 --af293-fasta F2 --out mq.tsv
  python3 -I analysis/aga1_class/pride_surface_evidence.py join --runs runs.tsv --uniprot-names NAMES.tsv --flags FLAGS.tsv [--maxquant mq.tsv] --out-dir NEW_DIR

`detections` writes one row per (run, protein) with distinct peptides and spectral counts. `join` writes the per-protein evidence table for the
Af293 proteome and the enrichment test. Both refuse an existing output and stop (FATAL:, exit 1) on malformed input. The PRIDE files are public but
untrusted: read with `python3 -I`; XML is parsed without external entities.
"""

import argparse
import re
import sys
from pathlib import Path
from xml.etree import ElementTree as ET  # noqa: N817

import pandas as pd

NS = "{http://regis-web.systemsbiology.net/protXML}"
SPECIES = {  # reference database name in the protXML header -> species, runs
    "AfumFDB": "A. fumigatus Af293",
    "AoerJGI": "A. oerlinghausenensis",
    "AfisFDB": "A. fischeri",
    "AlenNCBI": "A. lentulus",
}
AFU = re.compile(r"^Afu(\d)g(\d+)", re.IGNORECASE)
RUN_COLS = [
    "run", "run_suffix", "replicate", "species", "protein_id", "n_peptides", "spectral_count", "percent_coverage", "gene_product",
]  # fmt: skip


class PrideError(Exception):
    pass


def afu_to_afua(protein_id):
    m = AFU.match(protein_id)
    return f"AFUA_{m.group(1)}G{m.group(2)}" if m else ""


def parse_protxml(path):
    root = ET.parse(path).getroot()
    hdr = root.find(f"{NS}protein_summary_header")
    if hdr is None:
        raise PrideError(f"{path}: no protein_summary_header")
    db = hdr.get("reference_database", "").split(";")[0].split("_")[0]
    if db not in SPECIES:
        raise PrideError(f"{path}: unknown reference database {db!r}")
    rows = []
    for grp in root.iter(f"{NS}protein_group"):
        for prot in grp.findall(f"{NS}protein"):
            peps = prot.findall(f"{NS}peptide")
            ann = prot.find(f"{NS}annotation")
            desc = ann.get("protein_description", "") if ann is not None else ""
            m = re.search(r"gene_product=([^|]*)", desc)
            rows.append(
                {
                    "protein_id": prot.get("protein_name"),
                    "n_peptides": len({p.get("peptide_sequence") for p in peps}),
                    "spectral_count": sum(int(p.get("n_instances", "0")) for p in peps),
                    "percent_coverage": prot.get("percent_coverage"),
                    "gene_product": (m.group(1).strip() if m else ""),
                }
            )
    return SPECIES[db], rows


def detections(prot_dir):
    files = sorted(
        Path(prot_dir).glob("*.prot.xml"),
        key=lambda p: (int(re.match(r"\d+", p.name).group()), p.name),
    )
    if not files:
        raise PrideError(f"no *.prot.xml in {prot_dir}")
    out = []
    for f in files:
        stem = f.name.split(".")[0]
        num = int(re.match(r"\d+", stem).group())
        species, rows = parse_protxml(f)
        for r in rows:
            out.append(
                {
                    "run": stem,
                    "run_suffix": "s" if stem.endswith("s") else "plain",
                    "replicate": (num - 1) % 3 + 1,
                    "species": species,
                    **r,
                }
            )
    return pd.DataFrame(out)[RUN_COLS]


def join(runs, names, flags, mq=None):
    af = runs[runs["species"] == "A. fumigatus Af293"].copy()
    af["afua"] = af["protein_id"].map(afu_to_afua)
    unmapped = int((af["afua"] == "").sum())
    names = names.rename(
        columns={"Entry": "accession", "Gene Names (ORF)": "orf", "Protein names": "uniprot_name"}
    )
    names["afua"] = names["orf"].fillna("").str.extract(r"(AFUA_\w+)", expand=False).fillna("")
    names["accessions_for_locus"] = (
        names.groupby("afua")["accession"].transform("count").where(names["afua"] != "", 0)
    )
    fl = flags[flags["species"] == "Afum_Af293"].copy()
    ev = fl.merge(
        names[["accession", "afua", "uniprot_name", "accessions_for_locus"]],
        on="accession",
        how="left",
        validate="1:1",
    )
    per = (
        af.groupby("afua")
        .agg(
            runs_detected=("run", "nunique"),
            plain_runs=("run", lambda s: len({x for x in s if not x.endswith("s")})),
            s_runs=("run", lambda s: len({x for x in s if x.endswith("s")})),
            max_peptides=("n_peptides", "max"),
            total_spectra=("spectral_count", "sum"),
            pride_gene_product=("gene_product", "first"),
        )
        .reset_index()
    )
    ev = ev.merge(per, on="afua", how="left")
    for c in ("runs_detected", "plain_runs", "s_runs", "max_peptides", "total_spectra"):
        ev[c] = ev[c].fillna(0).astype(int)
    if mq is not None:
        mq = mq.copy()
        for c in ("samples_with_intensity", "intensity_sum"):
            mq[c] = pd.to_numeric(mq[c])
        m = mq[mq["af293_accession"] != ""].groupby("af293_accession").agg(
            mq_groups=("group_id", "nunique"), mq_samples=("samples_with_intensity", "max"), mq_intensity_sum=("intensity_sum", "sum")
        ).reset_index().rename(columns={"af293_accession": "accession"})  # fmt: skip
        ev = ev.merge(m, on="accession", how="left")
        for c in ("mq_groups", "mq_samples"):
            ev[c] = ev[c].fillna(0).astype(int)
        ev["mq_intensity_sum"] = ev["mq_intensity_sum"].fillna(0)
    ev["detected_any"] = (ev["runs_detected"] >= 1).astype(int)
    ev["detected_2plus_runs"] = (ev["runs_detected"] >= 2).astype(int)
    return ev, unmapped, af


KMER = 10


def _fasta(path):
    out, name, seq = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(seq)
            m = re.match(r">(?:sp|tr)\|([^|]+)\|", line)
            name, seq = (m.group(1) if m else line[1:].split()[0]), []
        else:
            seq.append(line.strip())
    if name:
        out[name] = "".join(seq)
    return out


def kmer_map(wanted, a1163, af293, k=KMER, min_share=0.8):
    """accession in `wanted` -> (best Af293 accession, share of its k-mers found in it); empty if below min_share or tied."""
    index = {}
    for acc, seq in af293.items():
        for i in range(len(seq) - k + 1):
            index.setdefault(seq[i : i + k], set()).add(acc)
    out = {}
    for acc in wanted:
        seq = a1163[acc]
        kmers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
        counts = {}
        for km in kmers:
            for tgt in index.get(km, ()):
                counts[tgt] = counts.get(tgt, 0) + 1
        if not counts:
            out[acc] = ("", 0.0)
            continue
        ranked = sorted(counts.items(), key=lambda kv: -kv[1])
        best, n = ranked[0]
        share = n / max(1, len(kmers))
        tied = len(ranked) > 1 and ranked[1][1] == n
        out[acc] = ("" if (share < min_share or tied) else best, round(share, 3))
    return out


def maxquant(groups, a1163_names, a1163_seqs, af293_seqs):
    g = groups[(groups["Reverse"] != "+") & (groups["Potential contaminant"] != "+")]
    g = g[g["Protein IDs"].str.contains("AFUB_", na=False)].copy()
    names = a1163_names.rename(columns={"Entry": "a1163_accession", "Gene Names (ORF)": "orf"})
    orf2acc = {
        o.split()[0]: a
        for a, o in zip(names["a1163_accession"], names["orf"].fillna(""), strict=True)
        if o
    }
    sample_cols = [c for c in g.columns if c.startswith("Intensity ") and c != "Intensity"]
    rows = []
    for _, r in g.iterrows():
        ids = [
            re.sub(r"-T$", "", x)
            for x in r["Majority protein IDs"].split(";")
            if x.startswith("AFUB_")
        ]
        accs = [orf2acc[i] for i in ids if i in orf2acc]
        for a in accs:
            rows.append({"a1163_accession": a, "group_id": r["id"], "n_in_group": len(accs), "peptides": r["Peptides"],
                         **{c: r[c] for c in sample_cols}})  # fmt: skip
    if not rows:
        raise PrideError(
            "no A. fumigatus protein groups could be mapped to A1163 UniProt accessions"
        )
    mq = pd.DataFrame(rows)
    mp = kmer_map(sorted(set(mq["a1163_accession"])), a1163_seqs, af293_seqs)
    mq["af293_accession"] = mq["a1163_accession"].map(lambda a: mp[a][0])
    mq["kmer_share"] = mq["a1163_accession"].map(lambda a: mp[a][1])
    for c in sample_cols:
        mq[c] = pd.to_numeric(mq[c], errors="coerce").fillna(0)
    mq["samples_with_intensity"] = (mq[sample_cols] > 0).sum(axis=1)
    mq["intensity_sum"] = mq[sample_cols].sum(axis=1)
    keep = [
        "a1163_accession",
        "af293_accession",
        "kmer_share",
        "group_id",
        "n_in_group",
        "peptides",
        "samples_with_intensity",
        "intensity_sum",
    ]
    return mq[keep + sample_cols]


def enrichment(ev):
    from scipy.stats import fisher_exact

    sets = {
        "signal peptide": ev["sp_yes"] == 1,
        "signal peptide and GPI": (ev["sp_yes"] == 1) & (ev["gpi_yes"] == 1),
        "broad rule": ev["broad_flag"] == 1,
        "narrow rule": ev["narrow_flag"] == 1,
        "2 of 3 models and signal peptide": (ev["n_models_call"] >= 2) & (ev["sp_yes"] == 1),
        "broad rule or 2 of 3 models with signal peptide": (ev["broad_flag"] == 1)
        | ((ev["n_models_call"] >= 2) & (ev["sp_yes"] == 1)),
    }
    rows = []
    for name, m in sets.items():
        for col in ("detected_any", "detected_2plus_runs"):
            a, b = int(ev.loc[m, col].sum()), int((~ev.loc[m, col].astype(bool) & m).sum())
            c, d = int(ev.loc[~m, col].sum()), int((~ev.loc[~m, col].astype(bool) & ~m).sum())
            odds, p = fisher_exact([[a, b], [c, d]])
            rows.append({"protein_set": name, "detection": col, "n_in_set": a + b, "detected_in_set": a, "rate_in_set": round(a / max(1, a + b), 3),
                         "n_outside": c + d, "detected_outside": c, "rate_outside": round(c / max(1, c + d), 4), "odds_ratio": round(odds, 2), "fisher_p": p})  # fmt: skip
    return pd.DataFrame(rows)


def _rd(path):
    return pd.read_csv(
        path, sep="\t", keep_default_na=False, na_values=[""], dtype=str, low_memory=False
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("detections")
    d.add_argument("--prot-dir", required=True)
    d.add_argument("--out", required=True)
    q = sub.add_parser("maxquant")
    q.add_argument("--groups", required=True)
    q.add_argument("--a1163-names", required=True)
    q.add_argument("--a1163-fasta", required=True)
    q.add_argument("--af293-fasta", required=True)
    q.add_argument("--out", required=True)
    j = sub.add_parser("join")
    j.add_argument("--runs", required=True)
    j.add_argument("--uniprot-names", required=True)
    j.add_argument("--flags", required=True)
    j.add_argument("--maxquant")
    j.add_argument("--out-dir", required=True)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "detections":
            if Path(a.out).exists():
                raise PrideError(f"{a.out} exists; nothing is overwritten")
            df = detections(a.prot_dir)
            df.to_csv(a.out, sep="\t", index=False)
            print(f"wrote {a.out}: {len(df)} detections in {df['run'].nunique()} runs")
        elif a.cmd == "maxquant":
            if Path(a.out).exists():
                raise PrideError(f"{a.out} exists; nothing is overwritten")
            mq = maxquant(
                _rd(a.groups), _rd(a.a1163_names), _fasta(a.a1163_fasta), _fasta(a.af293_fasta)
            )
            mq.to_csv(a.out, sep="\t", index=False)
            print(
                f"wrote {a.out}: {len(mq)} A1163 proteins, {int((mq['af293_accession'] != '').sum())} mapped to Af293"
            )
        else:
            out = Path(a.out_dir)
            if out.exists():
                raise PrideError(f"{out} exists; nothing is overwritten")
            runs = _rd(a.runs)
            for c in ("n_peptides", "spectral_count", "replicate"):
                runs[c] = runs[c].astype(int)
            flags = _rd(a.flags)
            for c in ("sp_yes", "gpi_yes", "broad_flag", "narrow_flag", "n_models_call"):
                flags[c] = (
                    flags[c].astype(float).fillna(0).astype(int)
                )  # a missing GPI call (protein of 40 aa or fewer) counts as no
            ev, unmapped, af = join(
                runs, _rd(a.uniprot_names), flags, _rd(a.maxquant) if a.maxquant else None
            )
            out.mkdir(parents=True)
            ev.to_csv(out / "afum_pride_surface_evidence.tsv", sep="\t", index=False)
            enr = enrichment(ev)
            enr.to_csv(out / "enrichment.tsv", sep="\t", index=False)
            det = ev[ev["detected_any"] == 1]
            print(f"Af293 runs: {af['run'].nunique()}; protein ids in the runs: {af['protein_id'].nunique()}; unmapped ids: {unmapped}; "
                  f"Af293 proteome rows: {len(ev)}; detected in any run: {len(det)}; in 2+ runs: {int(ev['detected_2plus_runs'].sum())}")  # fmt: skip
            print(enr.to_string(index=False))
    except (PrideError, OSError, ET.ParseError, KeyError, ValueError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
