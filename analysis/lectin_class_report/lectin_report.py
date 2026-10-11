#!/usr/bin/env python3.12
"""Lectin-like family report (spec option c-i). It makes NO call.

Spec: docs/superpowers/specs/2026-10-10-lectin-class-design.md (sections 3.1, 3.4, 6.2).

For each input proteome FASTA (plain, .gz or .zst) the script runs
`hmmsearch --cut_ga` with the lectin-like models of the spec (M1 and M2) against Pfam,
and writes three files in --out-dir:

  lectin_report_counts.tsv  per proteome and per model: proteins with a hit, best score
  lectin_report_hits.tsv    one row per protein and model (score, E-value, label)
  lectin_report.md          the report text (fixed texts of spec 6.2, notes)

The script writes no status file and makes no call. A hit is a domain hit. It does not
show location, secretion or host-glycan binding.

The optional --signalp file gives signal peptide data. Only then the report has the
"(a2) group" section (PF07938 score >= 100 and a signal peptide).

Usage:
  /usr/bin/python3.12 lectin_report.py --out-dir OUT --pfam-hmm Pfam-A.hmm \
      [--signalp FILE] [--tmpdir DIR] [--cpu 2] [--fasta-list LIST] [FASTA ...]

Optional --name-list FILE: proteins to name in the report (ID, then a free note per line).
With --signalp the table goes in the (a2) section. Without it, a "Named proteins" section.

ID assumption (owner decision 2026-10-10): signal peptide data and the name list are matched
on the protein ID alone (the first word of the FASTA header; for sp|ACC|NAME also ACC). This is
safe for the ID sources in use (Fungi_5k, UniProt, NCBI). The script warns on stderr, and in the
report, if the same ID occurs in two input proteomes. It does not stop. --no-dup-check skips the
check. The check keeps all IDs of the run in memory (about 100 bytes per protein).

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

# The label of spec decision D7, option A. Do not change the text.
PF07938_LABEL = "PF07938 (AAL-like / FleA-family) domain hit"

# Fixed text of spec 6.2 item 3.
FIXED_TEXT = (
    "This is a hit to the PF07938 (AAL-like / FleA-family) domain model. "
    "It does not show that the protein is on the cell surface or secreted. "
    "It does not show that the protein binds a host glycan. "
    "Host-glycan binding is shown in papers for FleA and SapL1 only; "
    "the role of FleA is disputed."
)

A2_FLOOR = 100.0

# (accession, set, H level, U on whole proteomes, 653 genera, from spec 3.1 table).
# Names are read from the Pfam file. H and U are copied from the spec. They are not
# measured by this script.
MODELS = [
    ("PF07938", "M1", "H1a", "392 genera (60.0%); 104 (15.9%) at score 100 or more"),
    ("PF01476", "M2", "H1b", "not scanned; at least 511 (78%), pool S value"),
    ("PF08881", "M2", "H2", "not scanned; at least 233 (36%), pool S value"),
    ("PF07367", "M2", "H2", "68 genera (10.4%)"),
    ("PF21691", "M2", "H0", "36 genera (5.5%)"),
    ("PF00337", "M2", "H0", "31 genera (4.7%)"),
    ("PF21595", "M2", "H0", "2 genera (0.3%)"),
    ("PF18647", "M2", "H0", "not scanned; at least 170 (26%), pool S value"),
    ("PF22803", "M2", "H0", "not scanned; at least 51 (8%), pool S value"),
    ("PF01419", "M2", "H0", "not scanned; at least 256 (39%), pool S value"),
    ("PF00652", "M2", "H0", "not scanned; at least 218 (33%), pool S value"),
    ("PF14200", "M2", "H0", "not scanned; at least 122 (19%), pool S value"),
    ("PF13385", "M2", "H0", "not scanned; at least 54 (8%), pool S value"),
    ("PF01607", "M2", "H1b", "31 genera (4.7%)"),
    ("PF00187", "M2", "H1b-like", "572 genera (87.6%)"),
]

# M3 house-keeping models. They are listed as a note and are not scanned.
M3_NOTE = [
    ("PF03388", "Lectin_leg-like", "ER-to-Golgi cargo receptors (Emp46p, Emp47p)"),
    ("PF00262", "Calreticulin", "calnexin and calreticulin"),
    ("PF07915 and PF02157", "PRKCSH and Man-6-P_recep", "P-type lectin (MRH) models"),
    ("PF11721", "Malectin", "ER glycoprotein quality control"),
]

PARACOCCIN_NOTE = (
    "Known gap (spec decision D9, no custom model): paracoccin of Paracoccidioides "
    "(UniProt C1G842) is a published GlcNAc-binding, laminin-binding protein "
    "(Coltri 2006, PMID 16476564). Its only Pfam hit is a partial GH18 chitinase "
    "domain (Glyco_hydro_18 PF00704). Its predicted lectin sites lie inside that domain. "
    "No lectin model of this report hits it, so this report does not show it."
)

ENDEMIC_NOTE = (
    "Coccidioides, Histoplasma and Blastomyces: the literature searches of spec F17 "
    "found no fungal lectin literature. The hits in this report are domain hits only."
)

F27_NOTE = (
    "Spec F27 names six full-length PF07938 proteins (score 288 to 331) in six Onygenales "
    "genera. No paper on the function of these proteins was found: a PubMed search for a "
    "dermatophyte AAL-like or fucose-binding lectin gave 0 hits (spec F27, review 3). "
    '"Not found" is not "does not exist". Kerr 2016 (PMID 27058347) reports putative '
    "FleA orthologs in dermatophytes, found by BLAST of AAL against NCBI sequences "
    "(S1 Fig legend). Orthology of these proteins to FleA is not checked by this report."
)

NO_CALL_TEXT = (
    "This report makes no call. It writes no status file. It does not claim a function "
    "for any protein. Every model here is reported only (spec decision D1, option c)."
)

FASTA_SUFFIXES = (".gz", ".zst", ".faa", ".fasta", ".fa", ".pep", ".proteins")


class ReportError(Exception):
    """A fatal input error."""


# --------------------------------------------------------------------------- input


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


def read_signalp(path):
    """Read signal peptide data. Return {protein id: True or False}.

    Accepted lines (tab or space separated). Lines that start with # are skipped.
      ID                 the protein has a signal peptide
      ID  PREDICTION     SignalP 6 style; OTHER or NO_SP means no signal peptide, any
                         other word (SP, LIPO, TAT, TATLIPO, PILIN) means a signal peptide
    An ID is the first word of the FASTA header. For a UniProt style header
    (sp|ACC|NAME) the ACC also matches.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            if "\t" in line:
                # SignalP 6 file: the ID column can hold the whole FASTA header
                cols = line.split("\t")
                parts = [cols[0].split()[0] if cols[0].split() else "", *cols[1:]]
            else:
                parts = line.split()
            pid = parts[0]
            if len(parts) == 1:
                has = True
            else:
                pred = parts[1].strip().upper()
                if pred == "ID" or pid == "ID":
                    continue
                has = pred not in ("OTHER", "NO_SP")
            out[pid] = has
    return out


def lookup_sp(sp, protein):
    """Return True, False or None (no data) for a protein id."""
    if protein in sp:
        return sp[protein]
    if "|" in protein:
        bits = protein.split("|")
        if len(bits) >= 2 and bits[1] in sp:
            return sp[bits[1]]
    return None


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
    """Run hmmsearch on one FASTA. The table goes to tbl_path.

    hmmsearch cannot read a stream with several models (the target must be rewindable).
    A .gz or .zst file is decompressed into tmpdir (node-local scratch) and removed after.
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


def scan_fasta(fasta):
    """Return (number of records, set of IDs). An ID is the first word of the header."""
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
    return n, ids


def read_name_list(path):
    """Read the list of proteins to name. Return {protein id: note}.

    One protein per line. The first tab- or space-separated word is the ID. The rest of the
    line is a free note (for example the genus). Lines that start with # are skipped.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            out[parts[0]] = parts[1].strip() if len(parts) > 1 else ""
    return out


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


def model_label(acc, name):
    if acc == "PF07938":
        return PF07938_LABEL
    return "%s (%s) domain hit" % (name, acc)


def fmt_score(x):
    return "NA" if x is None else ("%.1f" % x)


def write_atomic(path, text):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


def named_section(hits, sp, named, heading):
    """Lines for the table of named proteins (the --name-list option)."""
    md = [heading, ""]
    md.append(
        "These proteins were named by the person who ran the report (--name-list). "
        "Each row shows its PF07938 result in the input proteomes. No call is made. "
        "No function is claimed.\n"
    )
    md.append("| protein | note | proteome | PF07938 score | signal peptide | label |")
    md.append("|---|---|---|---|---|---|")
    by_prot = {}
    for h in hits:
        if h["model_acc"] == "PF07938":
            by_prot.setdefault(h["protein"], []).append(h)
    for pid, note in named.items():
        rows = by_prot.get(pid)
        sp_val = lookup_sp(sp, pid) if sp is not None else None
        sp_txt = "NA" if sp_val is None else ("yes" if sp_val else "no")
        if not rows:
            md.append(
                "| %s | %s | not in the input proteomes, or no PF07938 hit | NA | %s | NA |"
                % (pid, note, sp_txt)
            )
            continue
        for h in rows:
            md.append(
                "| %s | %s | %s | %.1f | %s | %s |"
                % (pid, note, h["proteome"], h["seq_score"], sp_txt, PF07938_LABEL)
            )
    md.append("")
    return md


def build_outputs(args, proteomes, info, hits, n_prot, sp, versions, named=None, dup_ids=None):
    """Make the three output texts. hits is a list of dict rows."""
    accs = [m[0] for m in MODELS]
    meta = {m[0]: m for m in MODELS}

    # counts table
    count_lines = [
        "proteome\tn_proteins\tmodel_acc\tmodel_name\tset\tn_hit_proteins\tbest_score\tn_hit_domain_ge_ga"
    ]
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
                        a,
                        info[a]["name"],
                        meta[a][1],
                        str(len(sel)),
                        fmt_score(best),
                        str(ndom),
                    ]
                )
            )

    # hits table
    hit_cols = [
        "proteome",
        "protein",
        "model_acc",
        "model_name",
        "set",
        "label",
        "seq_evalue",
        "seq_score",
        "domain_score",
        "domain_ge_ga",
    ]
    if sp is not None:
        hit_cols.append("signal_peptide")
    hit_lines = ["\t".join(hit_cols)]
    for h in sorted(
        hits, key=lambda r: (r["proteome"], r["model_acc"], -r["seq_score"], r["protein"])
    ):
        row = [
            h["proteome"],
            h["protein"],
            h["model_acc"],
            info[h["model_acc"]]["name"],
            meta[h["model_acc"]][1],
            model_label(h["model_acc"], info[h["model_acc"]]["name"]),
            "%g" % h["seq_evalue"],
            "%.1f" % h["seq_score"],
            "%.1f" % h["domain_score"],
            "yes" if h["domain_score"] >= info[h["model_acc"]]["ga_dom"] else "no",
        ]
        if sp is not None:
            v = lookup_sp(sp, h["protein"])
            row.append("NA" if v is None else ("yes" if v else "no"))
        hit_lines.append("\t".join(row))

    # report
    md = []
    md.append("# Lectin-like family report (no call)\n")
    md.append(
        "Date: %s. Spec: docs/superpowers/specs/2026-10-10-lectin-class-design.md "
        "(option c-i; sections 3.1 and 6.2).\n" % datetime.date.today().isoformat()
    )
    md.append(NO_CALL_TEXT + "\n")
    md.append("## Method\n")
    md.append(
        "- Search: `hmmsearch --cut_ga` (%s) with the models of the table below against "
        "Pfam %s (`%s`)." % (versions["hmmer"], versions["pfam"], args.pfam_hmm)
    )
    md.append(
        "- The cutoff is the Pfam GA cutoff of each model. A hit is a protein with a "
        "sequence score at or above the sequence GA (the `--tblout` table). Scores in the "
        "tables are sequence scores."
    )
    md.append(
        "- Column `domain_ge_ga` (hits table) and `n_hit_domain_ge_ga` (counts table) show "
        "whether the best domain also reaches the domain GA. The tool's own Pfam table "
        "comes from the domain table, so the tool can count fewer proteins for a "
        "model than this report (for example RicinB_lectin_2 in C. neoformans H99)."
    )
    md.append(
        "- Proteomes: %d. No signal peptide or other condition is applied to the hits."
        % len(proteomes)
    )
    md.append(
        "- Signal peptide data: %s.\n"
        % ("supplied (`%s`)" % args.signalp if sp is not None else "not supplied")
    )

    md.append("## PF07938 hits per proteome\n")
    md.append("Label of every PF07938 hit: **%s**.\n" % PF07938_LABEL)
    md.append('Fixed text for every PF07938 hit: "%s"\n' % FIXED_TEXT)
    md.append(
        "| proteome | proteins | PF07938 hit proteins | best score | proteins with score >= %g |"
        % A2_FLOOR
    )
    md.append("|---|---|---|---|---|")
    for p in proteomes:
        n, best = stats[(p, "PF07938")]
        hi = sum(
            1
            for h in hits
            if h["proteome"] == p and h["model_acc"] == "PF07938" and h["seq_score"] >= A2_FLOOR
        )
        md.append("| %s | %d | %d | %s | %d |" % (p, n_prot[p], n, fmt_score(best), hi))
    tot = sum(stats[(p, "PF07938")][0] for p in proteomes)
    md.append(
        "\nTotal PF07938 hit proteins: %d in %d of %d proteomes.\n"
        % (tot, sum(1 for p in proteomes if stats[(p, "PF07938")][0] > 0), len(proteomes))
    )

    md.append("## Lectin-like families (not called)\n")
    md.append(
        "Counts per model over all input proteomes. H is the host-interaction evidence "
        "level of spec 3.1. U is the whole-proteome ubiquity of spec 3.1 (653 genera of "
        "the discovery pass). H and U are copied from the spec. This script does not "
        'measure them. The column "proteomes with hit" is for the proteomes of this run.\n'
    )
    md.append(
        "| set | model | name | H | U (spec 3.1) | hit proteins | proteomes with hit | best score |"
    )
    md.append("|---|---|---|---|---|---|---|---|")
    for a, st, hl, u in MODELS:
        n = sum(stats[(p, a)][0] for p in proteomes)
        k = sum(1 for p in proteomes if stats[(p, a)][0] > 0)
        best = max(
            (stats[(p, a)][1] for p in proteomes if stats[(p, a)][1] is not None), default=None
        )
        md.append(
            "| %s | %s | %s | %s | %s | %d | %d | %s |"
            % (st, a, info[a]["name"], hl, u, n, k, fmt_score(best))
        )
    md.append(
        "\nM1 here means the models with H1a evidence (PF07938 only, spec 3.1). M1 is "
        "reported, not called. M2 models are reported, not called.\n"
    )

    md.append("## Per proteome and per model\n")
    md.append(
        "Hit proteins (best score) per proteome. A dash means no hit. The full table is "
        "`lectin_report_counts.tsv`.\n"
    )
    md.append("| proteome | " + " | ".join(info[a]["name"] for a in accs) + " |")
    md.append("|---|" + "---|" * len(accs))
    for p in proteomes:
        cells = []
        for a in accs:
            n, best = stats[(p, a)]
            cells.append("-" if n == 0 else "%d (%s)" % (n, fmt_score(best)))
        md.append("| %s | %s |" % (p, " | ".join(cells)))
    md.append("")

    if sp is not None:
        md.append("## (a2) group: PF07938 score >= %g and a signal peptide\n" % A2_FLOOR)
        md.append(
            "This group is named for the record only. It is the group of spec option (a2). "
            "No call is made. No function is claimed. Orthology to FleA is not checked. "
            "The spec states that the floor of 100 was fitted to the same data that "
            "found the group.\n"
        )
        grp, n_no, n_missing = [], 0, 0
        for h in hits:
            if h["model_acc"] != "PF07938" or h["seq_score"] < A2_FLOOR:
                continue
            v = lookup_sp(sp, h["protein"])
            if v is None:
                n_missing += 1
            elif v:
                grp.append(h)
            else:
                n_no += 1
        md.append(
            "PF07938 hits at score %g or more: %d with a signal peptide, %d without, "
            "%d with no signal peptide data.\n" % (A2_FLOOR, len(grp), n_no, n_missing)
        )
        if grp:
            md.append("| proteome | protein | score | label |")
            md.append("|---|---|---|---|")
            for h in sorted(grp, key=lambda r: (r["proteome"], -r["seq_score"], r["protein"])):
                md.append(
                    "| %s | %s | %.1f | %s |"
                    % (h["proteome"], h["protein"], h["seq_score"], PF07938_LABEL)
                )
            md.append('\nFixed text for every row: "%s"\n' % FIXED_TEXT)
        else:
            md.append("No protein is in the group.\n")
        md.append(F27_NOTE + "\n")
        if named:
            md.extend(named_section(hits, sp, named, "### Named proteins"))
    elif named:
        md.extend(named_section(hits, sp, named, "## Named proteins"))
        md.append(F27_NOTE + "\n")

    md.append("## Notes\n")
    md.append(
        "- M3 house-keeping models are not scanned. They are listed so the reader sees what "
        "is left out:"
    )
    for acc, name, what in M3_NOTE:
        md.append("  - %s %s: %s." % (acc, name, what))
    md.append(
        "- Models owned by the narrow tier (PA14 PF07691, GLEYA PF10528, Candida_ALS_N "
        "PF11766) are not part of this report."
    )
    md.append("- " + PARACOCCIN_NOTE)
    md.append("- " + ENDEMIC_NOTE)
    md.append(
        "- The list of models comes from a name search for lectin-like words (spec F4). "
        "A carbohydrate-binding protein whose model has no such word can be missed."
    )
    md.append(
        "- The report items 4 and 5 of spec 6.2 (UniLectin3D sequence table, T4 outcomes) "
        "are not made by this script."
    )
    if sp is not None or named:
        md.append(
            "- Signal peptide data and the name list are matched on the protein ID alone, not "
            "on the proteome name. Equal IDs in two proteomes would share one value."
        )
        if dup_ids:
            md.append(
                "- **Warning:** %d protein IDs occur in more than one input proteome. "
                "Matching can mix proteomes." % len(dup_ids)
            )
    md.append("")

    return "\n".join(count_lines) + "\n", "\n".join(hit_lines) + "\n", "\n".join(md)


# --------------------------------------------------------------------------- main


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("fasta", nargs="*", help="proteome FASTA files (plain, .gz or .zst)")
    ap.add_argument("--fasta-list", default=None, help="file with one FASTA path per line")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--pfam-hmm", default=DEFAULT_PFAM)
    ap.add_argument("--hmmsearch", default="hmmsearch", help="hmmsearch executable")
    ap.add_argument(
        "--signalp", default=None, help="optional signal peptide table; enables the (a2) section"
    )
    ap.add_argument(
        "--name-list",
        default=None,
        help="optional list of proteins to name in the report (ID, then a free note)",
    )
    ap.add_argument(
        "--no-dup-check",
        action="store_true",
        help="skip the check for equal protein IDs in two proteomes (saves memory)",
    )
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
    sp = read_signalp(args.signalp) if args.signalp else None
    named = read_name_list(args.name_list) if args.name_list else None

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(base, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="lectin_report.", dir=base)
    try:
        models_hmm = os.path.join(tmp, "models.hmm")
        info = extract_models(args.pfam_hmm, [m[0] for m in MODELS], models_hmm)
        versions = {"hmmer": hmmsearch_version(args.hmmsearch), "pfam": pfam_release(args.pfam_hmm)}
        hits, n_prot = [], {}
        check_dups = not args.no_dup_check
        seen_ids, dup_ids = set(), {}
        for name, fasta in zip(names, args.fasta, strict=False):
            tbl = os.path.join(tmp, name + ".tbl")
            run_hmmsearch(args.hmmsearch, models_hmm, fasta, tbl, args.cpu, tmp)
            n_prot[name], ids = scan_fasta(fasta)
            if check_dups:
                for pid in ids & seen_ids:
                    dup_ids.setdefault(pid, []).append(name)
                seen_ids |= ids
            for prot, acc, _name, ev, score, dscore in parse_tblout(tbl):
                hits.append(
                    {
                        "proteome": name,
                        "protein": prot,
                        "model_acc": acc,
                        "seq_evalue": ev,
                        "seq_score": score,
                        "domain_score": dscore,
                    }
                )
            os.remove(tbl)
        if dup_ids:
            ex = ", ".join(sorted(dup_ids)[:5])
            print(
                "WARNING: %d protein IDs occur in more than one input proteome (examples: %s). "
                "Signal peptide and name-list matching uses the protein ID alone, so equal IDs "
                "share one value." % (len(dup_ids), ex),
                file=sys.stderr,
            )
        counts, hit_tsv, md = build_outputs(
            args, names, info, hits, n_prot, sp, versions, named, dup_ids
        )
        write_atomic(os.path.join(args.out_dir, "lectin_report_counts.tsv"), counts)
        write_atomic(os.path.join(args.out_dir, "lectin_report_hits.tsv"), hit_tsv)
        write_atomic(os.path.join(args.out_dir, "lectin_report.md"), md)
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
