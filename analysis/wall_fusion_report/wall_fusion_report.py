#!/usr/bin/env python3.12
"""Cell wall remodeling and cell fusion family report (report only, option c-i). NO call.

Same pattern as analysis/lectin_class_report/lectin_report.py. For each input proteome FASTA
(plain, .gz or .zst) the script runs `hmmsearch --cut_ga` with the models listed in a model
TSV (wall_fusion_models.tsv) against Pfam, and writes three files in --out-dir:

  wall_fusion_hits.tsv    one row per protein and model (model-only label, scores)
  wall_fusion_counts.tsv  per proteome and per model: proteins with a hit, best score
  wall_fusion_report.md   the report text
  <table>.NOTICE          a sidecar file with the no-call notice, next to each table

The script writes no status file and makes no call. The tables have no call or status
column. A hit is a Pfam domain hit. It is NOT evidence of adhesion, of cell fusion or of a
remodeling function. Wall enzymes are known false-positive sources for surface-glycoprotein
detection (docs/reports/cell-wall-fusion/07-model-hit-check.md: the ESM classifier calls 6 of
13 S288C wall enzymes above 0.5, and none of the 13 hits an active family-table row).

The two tables start with lines that begin with "#". They hold the no-call notice. Skip them
when you read the table (pandas: comment="#"). A reader that does so loses the notice, so it is
ALSO written to wall_fusion_hits.tsv.NOTICE and wall_fusion_counts.tsv.NOTICE.

Inclusion rule for a model: it is verified in 05-pfam-id-verification.md (verification_level
05_interpro) OR it was measured in 07-model-hit-check.md (07_hmmscan_only). The script checks
the NAME of each model in the Pfam file against the name in the model TSV and stops on a
mismatch.

Usage:
  /usr/bin/python3.12 wall_fusion_report.py --out-dir OUT --models-tsv wall_fusion_models.tsv \
      --pfam-hmm Pfam-A.hmm [--set all|SET[,SET...]] [--tmpdir DIR] [--cpu 2] \
      [--fasta-list LIST] [FASTA ...]

Temp files go to --tmpdir, or to $SCRATCH when --tmpdir is not given. The script stops if
neither is set. Paths come from arguments only. No file location is derived from the
script path.
"""
# ruff: noqa: UP031

import argparse
import datetime
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict

DEFAULT_PFAM = "/srv/projects/db/pfam/2026-01-27-Pfam38.2/Pfam-A.hmm"

# wall_remodeling: enzymes that remodel the wall. wall_synthesis: chitin and glucan synthases.
# wall_structural_other: wall proteins that are not an enzyme family by themselves.
# cell_fusion_candidate_models: Pfam models of proteins named in the fusion literature.
SETS = (
    "wall_remodeling",
    "wall_synthesis",
    "wall_structural_other",
    "cell_fusion_candidate_models",
)
LEVELS = ("05_interpro", "07_hmmscan_only")
MODEL_COLUMNS = (
    "id",
    "name",
    "family",
    "set",
    "supporting_report",
    "verification_level",
    "verification_status",
)

# Custom (non-Pfam) models: a separate, optional table. They have no GA line; the threshold is
# the bit-score cutoff of the table (hmmsearch -T). Their hits are never merged with Pfam GA hits.
CUSTOM_COLUMNS = (
    "id",
    "name",
    "model_path",
    "cutoff_bits",
    "sha256",
    "kind",
    "use",
    "validation",
    "supporting_record",
)
CUSTOM_KIND = "custom"
CUSTOM_USE = "report_only"
CUSTOM_VALIDATION = "unvalidated_beyond_circular_holdout"

CUSTOM_NOTICE = (
    "CUSTOM MODEL SECTION. These models are not Pfam. Each is a custom, versioned model with "
    "NO GA line: the threshold is the bit-score cutoff stated per model (hmmsearch -T), set "
    "from the false-hit side only. They are report-only and UNVALIDATED beyond a circular "
    "held-out test (the held-out members came from the same search that supplied training). "
    "Their hits are kept apart from the Pfam GA counts and are never pooled with them. A hit "
    "is not a call, not a tier and not evidence of adhesion or of a remodeling function."
)

NO_CALL_TEXT = (
    "This report makes no call. It writes no status file. It has no call or status column. "
    "A hit is a hit to a Pfam family model. It is NOT evidence of adhesion, of cell fusion or "
    "of a cell wall remodeling function."
)

NOTICE_FILE_TEXT = (
    NO_CALL_TEXT + "\nThis notice is also in the first lines of the table and in the report.\n"
    "Pandas with comment='#' drops those lines. This file is written next to each table so the "
    "notice is not lost.\n"
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
    chosen = None
    if which != "all":
        chosen = [x.strip() for x in which.split(",") if x.strip()]
        bad = [x for x in chosen if x not in SETS]
        if bad or not chosen:
            raise ReportError("--set must be all or a comma list of: " + ", ".join(SETS))
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
            if row["verification_level"] not in LEVELS:
                raise ReportError(
                    "%s line %d: verification_level must be one of %s, got %r"
                    % (path, num, ", ".join(LEVELS), row["verification_level"])
                )
            rows.append(row)
    ids = [r["id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ReportError(path + ": a Pfam ID occurs twice")
    if chosen is not None:
        rows = [r for r in rows if r["set"] in chosen]
    if not rows:
        raise ReportError("no models for set " + which)
    return rows


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_custom_models(path):
    """Read the custom model TSV. Return a list of dicts (file order).

    model_path is resolved against the folder of the table when it is relative. The model
    file is NOT opened here; see verify_custom_model.
    """
    rows = []
    with open(path) as fh:
        head = fh.readline().rstrip("\r\n").split("\t")
        if tuple(head) != CUSTOM_COLUMNS:
            raise ReportError(
                "%s: header must be %s, got %s" % (path, "\t".join(CUSTOM_COLUMNS), "\t".join(head))
            )
        for num, raw in enumerate(fh, 2):
            line = raw.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != len(CUSTOM_COLUMNS):
                raise ReportError(
                    "%s line %d: expected %d columns" % (path, num, len(CUSTOM_COLUMNS))
                )
            row = dict(zip(CUSTOM_COLUMNS, cols, strict=True))
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", row["id"]):
                raise ReportError("%s line %d: bad custom model id %r" % (path, num, row["id"]))
            if re.fullmatch(r"PF\d{5}", row["id"]):
                raise ReportError(
                    "%s line %d: a custom id must not look like a Pfam ID" % (path, num)
                )
            try:
                row["cutoff"] = float(row["cutoff_bits"])
            except ValueError as exc:
                raise ReportError("%s line %d: cutoff_bits is not a number" % (path, num)) from exc
            if not re.fullmatch(r"[0-9a-f]{64}", row["sha256"]):
                raise ReportError(
                    "%s line %d: sha256 must be 64 lowercase hex digits" % (path, num)
                )
            for col, want in (
                ("kind", CUSTOM_KIND),
                ("use", CUSTOM_USE),
                ("validation", CUSTOM_VALIDATION),
            ):
                if row[col] != want:
                    raise ReportError("%s line %d: %s must be %r" % (path, num, col, want))
            if not os.path.isabs(row["model_path"]):
                row["model_path"] = os.path.normpath(
                    os.path.join(os.path.dirname(os.path.abspath(path)), row["model_path"])
                )
            rows.append(row)
    ids = [r["id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ReportError(path + ": a custom model id occurs twice")
    if not rows:
        raise ReportError(path + ": no custom models")
    return rows


def verify_custom_model(row):
    """Stop unless the model file exists, its sha256 is the pinned one and its NAME matches."""
    if not os.path.isfile(row["model_path"]):
        raise ReportError("custom model file not found: " + row["model_path"])
    got = sha256_file(row["model_path"])
    if got != row["sha256"]:
        raise ReportError(
            "custom model %s: sha256 of %s is %s, the table pins %s. Refusing to search."
            % (row["id"], row["model_path"], got, row["sha256"])
        )
    names = []
    with open(row["model_path"]) as fh:
        for line in fh:
            if line.startswith("NAME "):
                names.append(line.split()[1])
    if names != [row["name"]]:
        raise ReportError(
            "custom model %s: file has NAME %s, the table says %s"
            % (row["id"], ",".join(names) or "none", row["name"])
        )


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


def extract_models(pfam_hmm, accessions, out_path, expected_names=None):
    """Write the models with the given accessions to out_path.

    One pass over the Pfam file. Return {acc: {"name", "ga_seq", "ga_dom"}}.
    expected_names ({acc: name}) is checked against the NAME line of each model. A missing
    accession, a missing GA line or a different name stops the run.
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
    if expected_names:
        bad = [
            "%s: Pfam file says %s, model TSV says %s" % (a, info.get("name"), expected_names[a])
            for a, info in found.items()
            if info.get("name") != expected_names.get(a)
        ]
        if bad:
            raise ReportError(
                "model name differs between %s and the model TSV (renamed or renumbered "
                "accession?): %s" % (pfam_hmm, "; ".join(bad))
            )
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


def hmmsearch_command(hmmsearch, models_hmm, plain, tbl_path, cpu, cutoff_bits=None):
    """The hmmsearch command. Pfam: --cut_ga. A custom model (cutoff_bits given): -T, no --cut_ga."""
    thr = ["--cut_ga"] if cutoff_bits is None else ["-T", "%g" % cutoff_bits]
    return [hmmsearch, *thr, "--noali", "--cpu", str(cpu), "--tblout", tbl_path, models_hmm, plain]


def run_hmmsearch(hmmsearch, models_hmm, fasta, tbl_path, cpu, tmpdir, cutoff_bits=None):
    """Run hmmsearch (--cut_ga, or -T when cutoff_bits is given) on one FASTA.

    The table goes to tbl_path.

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
                hmmsearch_command(hmmsearch, models_hmm, plain, tbl_path, cpu, cutoff_bits),
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


def scan_fasta(fasta):
    """Return (number of records, set of IDs, number of repeated IDs inside the file).

    An ID is the first word of the header.
    """
    stream, close = open_binary(fasta)
    n = 0
    ids = set()
    try:
        for line in stream:
            if line.startswith(b">"):
                n += 1
                words = line[1:].split()
                ids.add(words[0].decode("utf-8", "replace") if words else "")
    finally:
        close()
    return n, ids, n - len(ids)


def parse_tblout(path):
    """Yield (protein, acc, name, evalue, score, domain_score) for each hit line."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            # tblout: target, t_acc, query, q_acc, E, score, bias, dom E, dom score ...
            yield f[0], f[3].split(".")[0], f[2], float(f[4]), float(f[5]), float(f[8])


def parse_tblout_by_name(path):
    """Yield (protein, model_name, evalue, score, domain_score). For models without ACC."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            yield f[0], f[2], float(f[4]), float(f[5]), float(f[8])


# --------------------------------------------------------------------------- output


def model_label(model):
    """Model-only label, as in the lectin report. It names the model and nothing else."""
    return "%s (%s) domain hit" % (model["family"], model["id"])


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
        "verification_level",
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
                    m["verification_level"],
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
        "- Proteomes: %d. No signal peptide, TM or other condition is applied to the hits."
        % len(proteomes)
    )
    md.append(
        "- `n_proteins` counts unique protein IDs (first word of the header). Records that "
        "repeat an ID are merged: one row per protein and model, with the best sequence "
        "score. Each repeat is reported as a warning below.\n"
    )

    md.append("## Models (all are family-only; none is a call)\n")
    md.append(
        "| set | model | name | family | verification | hit proteins | proteomes with hit "
        "| best score |"
    )
    md.append("|---|---|---|---|---|---|---|---|")
    for a in accs:
        n = sum(stats[(p, a)][0] for p in proteomes)
        k = sum(1 for p in proteomes if stats[(p, a)][0] > 0)
        best = max(
            (stats[(p, a)][1] for p in proteomes if stats[(p, a)][1] is not None), default=None
        )
        md.append(
            "| %s | %s | %s | %s | %s | %d | %d | %s |"
            % (
                meta[a]["set"],
                a,
                info[a]["name"],
                meta[a]["family"],
                meta[a]["verification_level"],
                n,
                k,
                fmt_score(best),
            )
        )
    md.append(
        "\nSupporting report and verification status of every model are in `%s`.\n"
        % args.models_tsv
    )

    md.append("## Per proteome and per set\n")
    md.append(
        "Distinct proteins with a hit to at least one model of the set. The sets are not "
        "added up: enzymes (wall_remodeling, wall_synthesis) and non-enzyme models "
        "(wall_structural_other, cell_fusion_candidate_models) are never in one total.\n"
    )
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
    if "cell_fusion_candidate_models" in present:
        md.append("- " + FUSION_NOTE)
    if "PF12454" in accs:
        md.append(
            "- Ecm33 (PF12454): S. cerevisiae Ecm33 itself has no hit to this model "
            "(07 row 22). A hit is not an Ecm33 call."
        )
    if "wall_structural_other" in present:
        md.append(
            "- wall_structural_other holds models that are not an enzyme family by themselves "
            "(Kre9/Knh1, Ecm33, Tos1, X8). SUN is kept here because its role is unclear: the "
            "InterPro name says beta-glucosidase, the reports show no enzyme role."
        )
    md.append(
        "- Large families (GH16, GH17, GH18, GH5) have intracellular and non-wall members. "
        "A hit says only that the protein carries the domain."
    )
    md.append(
        "- Inclusion rule: a model is listed if it is verified in "
        "docs/reports/cell-wall-fusion/05-pfam-id-verification.md (05_interpro) or measured in "
        "07-model-hit-check.md (07_hmmscan_only, name from Pfam-A, scores are hmmscan "
        "best-domain scores there). The script stops if a Pfam model NAME differs from "
        "the name in the model TSV."
    )
    md.append(
        "- Models left out on purpose: PF00026 (generic aspartyl protease), SH3_1 (PF00018, "
        "the Fus1 hit, generic domain), the GPI-anchor biosynthesis models (ER enzymes), "
        "WSC (a sensor), PF10333 (Pga1, an ER GPI enzyme), and the chitin-binding models "
        "(PF00187, PF01607, LysM; in the lectin report)."
    )
    md.append(
        "- Scores: 07 lists the hmmscan best-domain score; this report stores the sequence "
        "score (`seq_score`) and the best-domain score (`domain_score`). For PRM1 the "
        "values are 475.9 and 475.5. Same HMMER 3.4 and Pfam 38.2."
    )
    for w in warnings or []:
        md.append("- **Warning:** " + w)
    md.append("")

    return "\n".join(count_lines) + "\n", "\n".join(hit_lines) + "\n", "\n".join(md)


def build_custom_outputs(args, proteomes, custom, chits, n_prot):
    """Texts of the separate custom-model tables and the report section.

    Returns (counts_text, hits_text, md_section). Nothing here is added to the Pfam outputs.
    """
    head = ["# wall_fusion_report CUSTOM models. " + CUSTOM_NOTICE, "# " + NO_CALL_TEXT]
    ccols = [
        "proteome",
        "n_proteins",
        "custom_model",
        "model_sha256",
        "cutoff_bits",
        "n_hit_proteins",
        "best_score",
        "kind",
        "use",
        "validation",
    ]
    count_lines = head + ["\t".join(ccols)]
    for p in proteomes:
        for m in custom:
            sel = [h for h in chits if h["proteome"] == p and h["model_id"] == m["id"]]
            best = max((h["seq_score"] for h in sel), default=None)
            count_lines.append(
                "\t".join(
                    [
                        p,
                        str(n_prot[p]),
                        m["id"],
                        m["sha256"],
                        "%g" % m["cutoff"],
                        str(len(sel)),
                        fmt_score(best),
                        m["kind"],
                        m["use"],
                        m["validation"],
                    ]
                )
            )
    hcols = [
        "proteome",
        "protein",
        "custom_model",
        "model_sha256",
        "cutoff_bits",
        "seq_evalue",
        "seq_score",
        "domain_score",
        "kind",
        "use",
        "validation",
        "supporting_record",
    ]
    hit_lines = head + ["\t".join(hcols)]
    meta = {m["id"]: m for m in custom}
    for h in sorted(
        chits, key=lambda r: (r["proteome"], r["model_id"], -r["seq_score"], r["protein"])
    ):
        m = meta[h["model_id"]]
        hit_lines.append(
            "\t".join(
                [
                    h["proteome"],
                    h["protein"],
                    m["id"],
                    m["sha256"],
                    "%g" % m["cutoff"],
                    "%g" % h["seq_evalue"],
                    "%.1f" % h["seq_score"],
                    "%.1f" % h["domain_score"],
                    m["kind"],
                    m["use"],
                    m["validation"],
                    m["supporting_record"],
                ]
            )
        )
    md = ["## Custom models (separate from Pfam; report only; unvalidated)\n", CUSTOM_NOTICE + "\n"]
    md.append(
        "- Search: `hmmsearch -T <cutoff_bits>` (%s), no `--cut_ga`, models of `%s`. A hit is a "
        "protein with a sequence score at or above the cutoff. These counts are NOT part of any "
        "Pfam count above and are not added to it.\n" % (args.hmmer_version, args.custom_models_tsv)
    )
    md.append(
        "| custom model | sha256 | cutoff (bits) | hit proteins | proteomes with hit | best score |"
    )
    md.append("|---|---|---|---|---|---|")
    for m in custom:
        sel = [h for h in chits if h["model_id"] == m["id"]]
        k = len({h["proteome"] for h in sel})
        best = max((h["seq_score"] for h in sel), default=None)
        md.append(
            "| %s | `%s` | %g | %d | %d | %s |"
            % (
                m["id"],
                m["sha256"],
                m["cutoff"],
                len({(h["proteome"], h["protein"]) for h in sel}),
                k,
                fmt_score(best),
            )
        )
    md.append("")
    md.append("| proteome | proteins | " + " | ".join(m["id"] for m in custom) + " |")
    md.append("|---|---|" + "---|" * len(custom))
    for p in proteomes:
        cells = [
            str(len([h for h in chits if h["proteome"] == p and h["model_id"] == m["id"]]))
            for m in custom
        ]
        md.append("| %s | %d | %s |" % (p, n_prot[p], " | ".join(cells)))
    md.append("")
    for m in custom:
        md.append("- %s: supporting record `%s`." % (m["id"], m["supporting_record"]))
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
        help="gene set(s) to scan: all (default) or a comma list of " + ", ".join(SETS),
    )
    ap.add_argument(
        "--custom-models-tsv",
        default=None,
        help="optional wall_fusion_custom_models.tsv: custom models searched with -T, reported "
        "in a separate section and separate tables",
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
    custom = []
    if args.custom_models_tsv:
        custom = load_custom_models(args.custom_models_tsv)
        for row in custom:
            verify_custom_model(row)
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
        info = extract_models(
            args.pfam_hmm,
            [m["id"] for m in models],
            models_hmm,
            {m["id"]: m["name"] for m in models},
        )
        versions = {"hmmer": hmmsearch_version(args.hmmsearch), "pfam": pfam_release(args.pfam_hmm)}
        args.hmmer_version = versions["hmmer"]
        hits, chits, n_prot, warnings = [], [], {}, []
        seen_ids, dup_ids = set(), {}
        for name, fasta in zip(names, args.fasta, strict=True):
            n_records, ids, n_rep = scan_fasta(fasta)
            n_prot[name] = len(ids)  # unique IDs; n_records - len(ids) == n_rep
            for pid in ids & seen_ids:
                dup_ids.setdefault(pid, []).append(name)
            seen_ids |= ids
            if n_records == 0:
                msg = "%s: the FASTA file has no sequences. No search run." % name
                warnings.append(msg)
                print("WARNING: " + msg, file=sys.stderr)
                continue
            if n_rep:
                msg = (
                    "%s: %d records repeat an ID inside the file. Hits are counted once per "
                    "ID (best score kept); n_proteins counts unique IDs." % (name, n_rep)
                )
                warnings.append(msg)
                print("WARNING: " + msg, file=sys.stderr)
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
            for row in custom:
                ctbl = os.path.join(tmp, name + "." + row["id"] + ".tbl")
                run_hmmsearch(
                    args.hmmsearch, row["model_path"], fasta, ctbl, args.cpu, tmp, row["cutoff"]
                )
                cbest = {}
                for prot, _mname, ev, score, dscore in parse_tblout_by_name(ctbl):
                    if score < row["cutoff"]:
                        continue
                    if prot not in cbest or score > cbest[prot]["seq_score"]:
                        cbest[prot] = {
                            "proteome": name,
                            "protein": prot,
                            "model_id": row["id"],
                            "seq_evalue": ev,
                            "seq_score": score,
                            "domain_score": dscore,
                        }
                chits.extend(cbest.values())
                os.remove(ctbl)
        if dup_ids:
            ex = ", ".join(sorted(dup_ids)[:5])
            msg = (
                "%d protein IDs occur in more than one input proteome (examples: %s). Hits are "
                "kept per proteome and are not merged." % (len(dup_ids), ex)
            )
            warnings.append(msg)
            print("WARNING: " + msg, file=sys.stderr)
        counts, hit_tsv, md = build_outputs(
            args, names, models, info, hits, n_prot, versions, warnings
        )
        write_atomic(os.path.join(args.out_dir, "wall_fusion_counts.tsv"), counts)
        write_atomic(os.path.join(args.out_dir, "wall_fusion_hits.tsv"), hit_tsv)
        tables = ["wall_fusion_counts.tsv", "wall_fusion_hits.tsv"]
        if custom:
            ccounts, chit_tsv, cmd = build_custom_outputs(args, names, custom, chits, n_prot)
            md = md + "\n" + cmd
            write_atomic(os.path.join(args.out_dir, "wall_fusion_custom_counts.tsv"), ccounts)
            write_atomic(os.path.join(args.out_dir, "wall_fusion_custom_hits.tsv"), chit_tsv)
            tables += ["wall_fusion_custom_counts.tsv", "wall_fusion_custom_hits.tsv"]
        write_atomic(os.path.join(args.out_dir, "wall_fusion_report.md"), md)
        for table in tables:
            text = NOTICE_FILE_TEXT if "custom" not in table else CUSTOM_NOTICE + "\n"
            write_atomic(os.path.join(args.out_dir, table + ".NOTICE"), text)
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
