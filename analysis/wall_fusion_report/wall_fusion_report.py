#!/usr/bin/env python3.12
"""Cell wall remodeling and cell fusion family report (report only, option c-i). NO call.

Same pattern as analysis/lectin_class_report/lectin_report.py. For each input proteome FASTA
(plain, .gz or .zst) the script runs `hmmsearch --cut_ga` with the models listed in a model
TSV (wall_fusion_models.tsv) against Pfam, and writes three files in --out-dir:

  wall_fusion_hits.tsv    one row per protein and model (family-only label, scores)
  wall_fusion_counts.tsv  per proteome and per model: proteins with a hit, best score
  wall_fusion_report.md   the report text

The script writes no status file and makes no call. The tables have no call or status
column. A hit is a Pfam domain hit. It is NOT evidence of adhesion, of cell fusion or of a
remodeling function. Wall enzymes are known false-positive sources for surface-glycoprotein
detection (docs/reports/cell-wall-fusion/07-model-hit-check.md: the ESM classifier calls 6 of
13 S288C wall enzymes above 0.5, and none of the 13 hits an active family-table row).

The two tables start with lines that begin with "#". They hold the no-call notice. Skip them
when you read the table (pandas: comment="#").

Usage:
  /usr/bin/python3.12 wall_fusion_report.py --out-dir OUT --models-tsv wall_fusion_models.tsv \
      --pfam-hmm Pfam-A.hmm [--set wall_remodeling|cell_fusion|all] [--tmpdir DIR] [--cpu 2] \
      [--fasta-list LIST] [FASTA ...]

Temp files go to --tmpdir, or to $SCRATCH when --tmpdir is not given. The script stops if
neither is set. Paths come from arguments only. No file location is derived from the
script path.
"""
# ruff: noqa: UP031

import argparse
import datetime
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict

DEFAULT_PFAM = "/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm"

SETS = ("wall_remodeling", "cell_fusion")
MODEL_COLUMNS = ("id", "name", "family", "set", "supporting_report", "verification_status")

NO_CALL_TEXT = (
    "This report makes no call. It writes no status file. It has no call or status column. "
    "A hit is a hit to a Pfam family model. It is NOT evidence of adhesion, of cell fusion or "
    "of a cell wall remodeling function."
)

FALSE_POSITIVE_NOTE = (
    "Cell wall enzymes are known false-positive sources for surface-glycoprotein detection. "
    "In the S288C check of docs/reports/cell-wall-fusion/07-model-hit-check.md the ESM "
    "classifier called 6 of 13 wall enzymes above 0.5 (CRH1, UTR2, SCW11, CTS1, SUN4, KRE9), "
    "and none of the 13 hit an active family-table row. Do not read a family hit here as "
    "an adhesin."
)

FUSION_NOTE = (
    "Fusion proteins such as Prm1 and Fig1 are membrane proteins. They lie outside the "
    "signal-peptide population of step 1 of the sorting hat. A hit to the Prm1 or Fig1 model "
    "does not show a role in cell fusion. The Fus1 SH3 hit (SH3_1) is not a fusion family "
    "and is not scanned."
)

FASTA_SUFFIXES = (".gz", ".zst", ".faa", ".fasta", ".fa", ".pep", ".proteins")


class ReportError(Exception):
    """A fatal input error."""


def table_header_lines():
    """Comment lines at the top of the hits and counts tables."""
    return [
        "# wall_fusion_report: family-only context table. " + NO_CALL_TEXT,
        "# " + FALSE_POSITIVE_NOTE,
    ]


# --------------------------------------------------------------------------- input


def load_models(path, which="all"):
    """Read the model TSV. Return a list of dicts (file order) for the chosen set."""
    if which != "all" and which not in SETS:
        raise ReportError("--set must be all or one of: " + ", ".join(SETS))
    rows = []
    with open(path) as fh:
        head = fh.readline().rstrip("\r\n").split("\t")
        if tuple(head) != MODEL_COLUMNS:
            raise ReportError(
                "%s: header must be %s, got %s" % (path, "\t".join(MODEL_COLUMNS), "\t".join(head))
            )
        for num, raw in enumerate(fh, 2):
            line = raw.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != len(MODEL_COLUMNS):
                raise ReportError(
                    "%s line %d: expected %d columns" % (path, num, len(MODEL_COLUMNS))
                )
            row = dict(zip(MODEL_COLUMNS, cols, strict=True))
            if not re.fullmatch(r"PF\d{5}", row["id"]):
                raise ReportError("%s line %d: bad Pfam ID %r" % (path, num, row["id"]))
            if row["set"] not in SETS:
                raise ReportError("%s line %d: unknown set %r" % (path, num, row["set"]))
            rows.append(row)
    ids = [r["id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ReportError(path + ": a Pfam ID occurs twice")
    if which != "all":
        rows = [r for r in rows if r["set"] == which]
    if not rows:
        raise ReportError("no models for set " + which)
    return rows


def proteome_name(path):
    """File name without the compression and FASTA suffixes."""
    name = os.path.basename(path)
    changed = True
    while changed:
        changed = False
        for suf in FASTA_SUFFIXES:
            if name.endswith(suf) and len(name) > len(suf):
                name = name[: -len(suf)]
                changed = True
    return name


def open_binary(path):
    """Return a binary stream and a closer for a plain, .gz or .zst file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rb"), (lambda: None)
    if path.endswith(".zst"):
        if shutil.which("zstd") is None:
            raise ReportError("zstd is not on PATH; it is needed to read " + path)
        proc = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE)

        def close():
            proc.stdout.close()
            if proc.wait() != 0:
                raise ReportError("zstd failed on " + path)

        return proc.stdout, close
    return open(path, "rb"), (lambda: None)


# --------------------------------------------------------------------------- models


def extract_models(pfam_hmm, accessions, out_path):
    """Write the models with the given accessions to out_path.

    One pass over the Pfam file. Return {acc: {"name", "ga_seq", "ga_dom"}}.
    """
    wanted = set(accessions)
    found = OrderedDict()
    block = []
    acc = None
    with open(pfam_hmm) as src, open(out_path, "w") as dst:
        for line in src:
            block.append(line)
            if line.startswith("ACC "):
                acc = line.split()[1].split(".")[0]
            if line.startswith("//"):
                if acc in wanted:
                    dst.writelines(block)
                    info = {}
                    for b in block:
                        if b.startswith("NAME "):
                            info["name"] = b.split()[1]
                        elif b.startswith("GA "):
                            v = b.split()[1:3]
                            info["ga_seq"] = float(v[0].rstrip(";"))
                            info["ga_dom"] = float(v[1].rstrip(";"))
                    found[acc] = info
                block = []
                acc = None
    missing = [a for a in accessions if a not in found]
    if missing:
        raise ReportError("models not found in " + pfam_hmm + ": " + ", ".join(missing))
    for a, info in found.items():
        if "ga_seq" not in info:
            raise ReportError("model %s has no GA line; --cut_ga cannot be used" % a)
    return found


def pfam_release(pfam_hmm):
    m = re.search(r"Pfam(\d+\.\d+)", os.path.realpath(pfam_hmm))
    return m.group(1) if m else "unknown"


def hmmsearch_version(hmmsearch):
    try:
        res = subprocess.run([hmmsearch, "-h"], capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ReportError("cannot run %s: %s" % (hmmsearch, exc)) from exc
    for line in res.stdout.splitlines():
        if line.startswith("# HMMER"):
            return line.lstrip("# ").strip()
    return "HMMER (version line not found)"


# --------------------------------------------------------------------------- search


def run_hmmsearch(hmmsearch, models_hmm, fasta, tbl_path, cpu, tmpdir):
    """Run hmmsearch --cut_ga on one FASTA. The table goes to tbl_path.

    A .gz or .zst file is decompressed into tmpdir and removed after.
    """
    plain = fasta
    made = None
    if fasta.endswith((".gz", ".zst")):
        made = os.path.join(tmpdir, "target.fa")
        stream, close = open_binary(fasta)
        try:
            with open(made, "wb") as out:
                shutil.copyfileobj(stream, out, 1 << 20)
        finally:
            close()
        plain = made
    err_path = os.path.join(tmpdir, "hmmsearch.err")
    try:
        with open(err_path, "wb") as err:
            rc = subprocess.run(
                [
                    hmmsearch,
                    "--cut_ga",
                    "--noali",
                    "--cpu",
                    str(cpu),
                    "--tblout",
                    tbl_path,
                    models_hmm,
                    plain,
                ],
                stdout=subprocess.DEVNULL,
                stderr=err,
            ).returncode
    finally:
        if made and os.path.exists(made):
            os.remove(made)
    if rc != 0:
        with open(err_path) as fh:
            msg = fh.read()[-500:]
        raise ReportError("hmmsearch failed on %s (exit %d): %s" % (fasta, rc, msg))


def count_records(fasta):
    """Return the number of FASTA records."""
    stream, close = open_binary(fasta)
    n = 0
    try:
        for line in stream:
            if line.startswith(b">"):
                n += 1
    finally:
        close()
    return n


def parse_tblout(path):
    """Yield (protein, acc, name, evalue, score, domain_score) for each hit line."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            # tblout: target, t_acc, query, q_acc, E, score, bias, dom E, dom score ...
            yield f[0], f[3].split(".")[0], f[2], float(f[4]), float(f[5]), float(f[8])


# --------------------------------------------------------------------------- output


def model_label(model):
    """Family-only label. It names the family and nothing else."""
    return "%s_family_hit: %s" % (model["set"], model["family"])


def fmt_score(x):
    return "NA" if x is None else ("%.1f" % x)


def write_atomic(path, text):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


def build_outputs(args, proteomes, models, info, hits, n_prot, versions, warnings=None):
    """Make the three output texts. hits is a list of dict rows."""
    meta = {m["id"]: m for m in models}
    accs = [m["id"] for m in models]
    head = table_header_lines()

    count_cols = [
        "proteome",
        "n_proteins",
        "set",
        "model_acc",
        "model_name",
        "family",
        "n_hit_proteins",
        "best_score",
        "n_hit_domain_ge_ga",
    ]
    count_lines = head + ["\t".join(count_cols)]
    stats = {}
    for p in proteomes:
        for a in accs:
            sel = [h for h in hits if h["proteome"] == p and h["model_acc"] == a]
            best = max((h["seq_score"] for h in sel), default=None)
            stats[(p, a)] = (len(sel), best)
            ndom = sum(1 for h in sel if h["domain_score"] >= info[a]["ga_dom"])
            count_lines.append(
                "\t".join(
                    [
                        p,
                        str(n_prot[p]),
                        meta[a]["set"],
                        a,
                        info[a]["name"],
                        meta[a]["family"],
                        str(len(sel)),
                        fmt_score(best),
                        str(ndom),
                    ]
                )
            )

    hit_cols = [
        "proteome",
        "protein",
        "set",
        "model_acc",
        "model_name",
        "family",
        "label",
        "seq_evalue",
        "seq_score",
        "domain_score",
        "domain_ge_ga",
        "supporting_report",
    ]
    hit_lines = head + ["\t".join(hit_cols)]
    for h in sorted(
        hits, key=lambda r: (r["proteome"], r["model_acc"], -r["seq_score"], r["protein"])
    ):
        m = meta[h["model_acc"]]
        hit_lines.append(
            "\t".join(
                [
                    h["proteome"],
                    h["protein"],
                    m["set"],
                    h["model_acc"],
                    info[h["model_acc"]]["name"],
                    m["family"],
                    model_label(m),
                    "%g" % h["seq_evalue"],
                    "%.1f" % h["seq_score"],
                    "%.1f" % h["domain_score"],
                    "yes" if h["domain_score"] >= info[h["model_acc"]]["ga_dom"] else "no",
                    m["supporting_report"],
                ]
            )
        )

    md = []
    md.append("# Cell wall remodeling and cell fusion family report (no call)\n")
    md.append(
        "Date: %s. Pattern: option c-i of the lectin report. Plan: "
        "docs/reports/cell-wall-fusion/04-repo-coverage-and-gaps.md section 4.\n"
        % datetime.date.today().isoformat()
    )
    md.append(NO_CALL_TEXT + "\n")
    md.append(FALSE_POSITIVE_NOTE + "\n")
    md.append("## Method\n")
    md.append(
        "- Search: `hmmsearch --cut_ga` (%s) with the models of `%s` against Pfam %s (`%s`). "
        "Set: %s." % (versions["hmmer"], args.models_tsv, versions["pfam"], args.pfam_hmm, args.set)
    )
    md.append(
        "- The cutoff is the Pfam GA cutoff of each model. A hit is a protein with a "
        "sequence score at or above the sequence GA. Scores in the tables are sequence scores."
    )
    md.append("- Column `domain_ge_ga` shows whether the best domain also reaches the domain GA.")
    md.append(
        "- Proteomes: %d. No signal peptide, TM or other condition is applied to the hits.\n"
        % len(proteomes)
    )

    md.append("## Models (all are family-only; none is a call)\n")
    md.append("| set | model | name | family | hit proteins | proteomes with hit | best score |")
    md.append("|---|---|---|---|---|---|---|")
    for a in accs:
        n = sum(stats[(p, a)][0] for p in proteomes)
        k = sum(1 for p in proteomes if stats[(p, a)][0] > 0)
        best = max(
            (stats[(p, a)][1] for p in proteomes if stats[(p, a)][1] is not None), default=None
        )
        md.append(
            "| %s | %s | %s | %s | %d | %d | %s |"
            % (meta[a]["set"], a, info[a]["name"], meta[a]["family"], n, k, fmt_score(best))
        )
    md.append(
        "\nSupporting report and verification status of every model are in `%s`.\n"
        % args.models_tsv
    )

    md.append("## Per proteome and per set\n")
    md.append("Distinct proteins with a hit to at least one model of the set.\n")
    md.append(
        "| proteome | proteins | "
        + " | ".join(s for s in SETS if any(meta[a]["set"] == s for a in accs))
        + " |"
    )
    present = [s for s in SETS if any(meta[a]["set"] == s for a in accs)]
    md.append("|---|---|" + "---|" * len(present))
    for p in proteomes:
        cells = []
        for s in present:
            prots = {
                h["protein"]
                for h in hits
                if h["proteome"] == p and meta[h["model_acc"]]["set"] == s
            }
            cells.append(str(len(prots)))
        md.append("| %s | %d | %s |" % (p, n_prot[p], " | ".join(cells)))
    md.append("")

    md.append("## Notes\n")
    if "cell_fusion" in present:
        md.append("- " + FUSION_NOTE)
    md.append(
        "- Large families (GH16, GH17, GH18, GH5) have intracellular and non-wall members. "
        "A hit says only that the protein carries the domain."
    )
    md.append(
        "- Models left out on purpose: PF01522 (chitin deacetylase, CE4 link not verified), "
        "PF00026 (generic aspartyl protease), the GPI-anchor biosynthesis models "
        "(ER enzymes), WSC (a sensor), the chitin-binding models (in the lectin report), and "
        "any ID that docs/reports/cell-wall-fusion/05-pfam-id-verification.md marks "
        "UNVERIFIABLE."
    )
    for w in warnings or []:
        md.append("- **Warning:** " + w)
    md.append("")

    return "\n".join(count_lines) + "\n", "\n".join(hit_lines) + "\n", "\n".join(md)


# --------------------------------------------------------------------------- main


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("fasta", nargs="*", help="proteome FASTA files (plain, .gz or .zst)")
    ap.add_argument("--fasta-list", default=None, help="file with one FASTA path per line")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--models-tsv", required=True, help="wall_fusion_models.tsv")
    ap.add_argument(
        "--set",
        default="all",
        help="gene set to scan: wall_remodeling, cell_fusion or all (default)",
    )
    ap.add_argument("--pfam-hmm", default=DEFAULT_PFAM)
    ap.add_argument("--hmmsearch", default="hmmsearch", help="hmmsearch executable")
    ap.add_argument("--tmpdir", default=None, help="temp dir (default: $SCRATCH)")
    ap.add_argument("--cpu", type=int, default=2)
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.fasta_list:
        with open(args.fasta_list) as fh:
            args.fasta = list(args.fasta) + [ln.strip() for ln in fh if ln.strip()]
    if not args.fasta:
        raise ReportError("no input FASTA given")
    models = load_models(args.models_tsv, args.set)
    base = args.tmpdir or os.environ.get("SCRATCH")
    if not base:
        raise ReportError("set --tmpdir or the SCRATCH variable (node-local temp space)")
    if shutil.which(args.hmmsearch) is None and not os.access(args.hmmsearch, os.X_OK):
        raise ReportError("hmmsearch not found: " + args.hmmsearch)
    if not os.path.isfile(args.pfam_hmm):
        raise ReportError("Pfam file not found: " + args.pfam_hmm)

    names = [proteome_name(f) for f in args.fasta]
    if len(set(names)) != len(names):
        raise ReportError("two inputs give the same proteome name: " + ", ".join(names))
    for f in args.fasta:
        if not os.path.isfile(f):
            raise ReportError("input not found: " + f)

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(base, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="wall_fusion_report.", dir=base)
    try:
        models_hmm = os.path.join(tmp, "models.hmm")
        info = extract_models(args.pfam_hmm, [m["id"] for m in models], models_hmm)
        versions = {"hmmer": hmmsearch_version(args.hmmsearch), "pfam": pfam_release(args.pfam_hmm)}
        hits, n_prot, warnings = [], {}, []
        for name, fasta in zip(names, args.fasta, strict=True):
            n_prot[name] = count_records(fasta)
            if n_prot[name] == 0:
                msg = "%s: the FASTA file has no sequences. No search run." % name
                warnings.append(msg)
                print("WARNING: " + msg, file=sys.stderr)
                continue
            tbl = os.path.join(tmp, name + ".tbl")
            run_hmmsearch(args.hmmsearch, models_hmm, fasta, tbl, args.cpu, tmp)
            best = {}
            for prot, acc, _name, ev, score, dscore in parse_tblout(tbl):
                if acc not in info:
                    continue
                key = (prot, acc)
                if key not in best or score > best[key]["seq_score"]:
                    best[key] = {
                        "proteome": name,
                        "protein": prot,
                        "model_acc": acc,
                        "seq_evalue": ev,
                        "seq_score": score,
                        "domain_score": dscore,
                    }
            hits.extend(best.values())
            os.remove(tbl)
        counts, hit_tsv, md = build_outputs(
            args, names, models, info, hits, n_prot, versions, warnings
        )
        write_atomic(os.path.join(args.out_dir, "wall_fusion_counts.tsv"), counts)
        write_atomic(os.path.join(args.out_dir, "wall_fusion_hits.tsv"), hit_tsv)
        write_atomic(os.path.join(args.out_dir, "wall_fusion_report.md"), md)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("wrote report for %d proteomes to %s" % (len(names), args.out_dir))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ReportError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        sys.exit(2)
