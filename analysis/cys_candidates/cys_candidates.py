#!/usr/bin/env python3.12
"""Candidate list of short, cysteine-rich, secreted proteins (read-only analysis).

Python 3.12, standard library only. No model, no training, no accuracy claim.

Inputs per proteome: a protein FASTA (plain or .gz), a SignalP 6 prediction_results.txt
(plain or .gz), and an optional hmmsearch --domtblout file (plain or .gz) for the known
families. Output: one row per protein with sequence features and a tier.

Tiers (thresholds are command-line options, see --help):
  cys_rich_sp_unassigned    SP called, mature length <= L, max Cys in a W-window >= K,
                            no cfem / bys1 / hydrophobin hit. PF28404 (pra3_like_family)
                            hits stay here and are flagged.
  cys_rich_sp_known_family  same, with a cfem / bys1 / hydrophobin hit.
  cys_rich_no_sp            no SP call; length <= L and max Cys in a W-window >= K on the
                            full sequence. Family flags are reported, not used.
  other                     everything else.

STOP contract: any input problem prints "STOP: <reason>" to stderr and exits 2. No output
file is left behind: outputs are written to temporary names and moved with os.replace
only after every proteome has been processed.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

DEFAULT_MAX_MATURE_LEN = 300  # L
DEFAULT_MIN_CYS_WINDOW = 8  # K
DEFAULT_WINDOW = 60  # W

# Pfam accession (no version) -> family flag.
FAMILY_OF_PFAM = {
    "PF05730": "cfem",
    "PF04681": "bys1",
    "PF01185": "hydrophobin",
    "PF06766": "hydrophobin",
    "PF28987": "hydrophobin",  # DewD, "Unclassified hydrophobin dewD"
    "PF28404": "pra3_like_family",
}
FAMILIES = ("cfem", "bys1", "hydrophobin", "pra3_like_family")
EXCLUDING_FAMILIES = ("cfem", "bys1", "hydrophobin")

TIERS = (
    "cys_rich_sp_unassigned",
    "cys_rich_sp_known_family",
    "cys_rich_no_sp",
    "other",
)

COLUMNS = (
    "proteome",
    "protein_id",
    "header",
    "length",
    "sp_call",
    "sp_prob",
    "cs_end",
    "mature_length",
    "cys_count",
    "cys_frac",
    "max_cys_window",
    "window_start",
    "cc_pairs",
    "cxc",
    "pest_frac",
    "cfem",
    "cfem_evalue",
    "bys1",
    "bys1_evalue",
    "hydrophobin",
    "hydrophobin_evalue",
    "pra3_like_family",
    "pra3_like_family_evalue",
    "tier",
)


class StopArgParser(argparse.ArgumentParser):
    """argparse that reports 'STOP: ...' and exits 2 (the STOP contract)."""

    def error(self, message):
        print(f"STOP: {message}", file=sys.stderr)
        sys.exit(2)


class Stop(Exception):  # noqa: N818
    """Input problem. Reported as 'STOP: <reason>' with exit code 2."""


# ----------------------------------------------------------------------------- I/O


def open_text(path):
    """Open a plain or gzip text file (detected by magic bytes, not by suffix)."""
    path = Path(path)
    try:
        with open(path, "rb") as fh:
            magic = fh.read(2)
    except OSError as exc:
        raise Stop(f"cannot read {path}: {exc}") from exc
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, encoding="utf-8", newline="")


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def short_id(header):
    """First whitespace-delimited token: the id hmmsearch reports as the target name."""
    return header.split(None, 1)[0] if header.split() else ""


def read_fasta(path):
    """Return a list of (header, sequence). STOP on duplicate headers or short ids."""
    if not Path(path).is_file():
        raise Stop(f"FASTA not found: {path}")
    records = []
    seen_full = set()
    seen_short = set()
    header = None
    parts = []

    def flush():
        if header is None:
            return
        sid = short_id(header)
        if not sid:
            raise Stop(f"empty sequence id in {path}")
        if header in seen_full or sid in seen_short:
            raise Stop(f"duplicate id {sid!r} in {path}")
        seen_full.add(header)
        seen_short.add(sid)
        # Rule: lowercase is uppercased; ONE trailing '*' (stop symbol) is removed; any other
        # non-letter (internal '*', '-', digit, space) or an empty sequence is a STOP.
        seq = "".join(parts).upper()
        if seq.endswith("*"):
            seq = seq[:-1]
        if not seq:
            raise Stop(f"{path}: empty sequence for {sid!r}")
        bad = re.search(r"[^A-Z]", seq)
        if bad:
            raise Stop(
                f"{path}: non-letter {bad.group()!r} at position {bad.start() + 1} in {sid!r}"
            )
        records.append((header, seq))

    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                flush()
                header = line[1:]
                parts = []
            elif header is None:
                if line.strip():
                    raise Stop(f"{path}: sequence before the first header")
            else:
                parts.append(line.strip())
    flush()
    if not records:
        raise Stop(f"no sequences in {path}")
    return records


_CS_RE = re.compile(r"CS pos:\s*(\d+)-(\d+)")


def parse_cs(text):
    """Return the last residue of the signal peptide (the N in 'CS pos: N-M'), or None."""
    m = _CS_RE.search(text or "")
    return int(m.group(1)) if m else None


def read_signalp(path):
    """Return {full_header: (call, sp_prob, cs_end or None)} from prediction_results.txt."""
    if not Path(path).is_file():
        raise Stop(f"SignalP file not found: {path}")
    out = {}
    header_seen = False
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                cols = line.lstrip("# ").split("\t")
                if len(cols) >= 4 and cols[0] == "ID":
                    if cols[1] != "Prediction" or not cols[3].startswith("SP"):
                        raise Stop(f"{path}: unexpected SignalP column header: {line!r}")
                    header_seen = True
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                raise Stop(f"{path}: row with fewer than 4 tab-separated columns: {line[:80]!r}")
            pid, call = cols[0], cols[1]
            if call not in ("SP", "OTHER"):
                raise Stop(f"{path}: unknown SignalP prediction {call!r} for {short_id(pid)}")
            try:
                prob = float(cols[3])
            except ValueError as exc:
                raise Stop(f"{path}: bad SP probability for {short_id(pid)}") from exc
            cs_end = parse_cs(cols[4]) if len(cols) > 4 else None
            if call == "SP" and cs_end is None:
                raise Stop(f"{path}: SP call without a CS position for {short_id(pid)}")
            if pid in out:
                raise Stop(f"duplicate id {short_id(pid)!r} in {path}")
            out[pid] = (call, prob, cs_end)
    if not header_seen:
        raise Stop(f"{path}: no '# ID<TAB>Prediction...' header line; not a SignalP 6 file")
    if not out:
        raise Stop(f"no rows in {path}")
    return out


def read_domtbl(path):
    """Return {short_id: {family: best i-Evalue}} from an hmmsearch --domtblout file."""
    if not Path(path).is_file():
        raise Stop(f"domtbl not found: {path}")
    hits = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split(None, 22)
            if len(f) < 22:
                raise Stop(f"{path}: domtbl row with fewer than 22 columns")
            target, qname, qacc = f[0], f[3], f[4]
            acc = (qacc if qacc != "-" else qname).split(".")[0]
            family = FAMILY_OF_PFAM.get(acc)
            if family is None:
                raise Stop(f"{path}: unexpected model {qname!r} ({qacc}); not one of the five")
            try:
                ev = float(f[12])
            except ValueError as exc:
                raise Stop(f"{path}: bad i-Evalue for {target}") from exc
            fam = hits.setdefault(target, {})
            if family not in fam or ev < fam[family]:
                fam[family] = ev
    return hits


# ------------------------------------------------------------------------ features


def max_cys_window(seq, window):
    """Return (max Cys count in any window of `window` residues, 1-based start of the first
    such window). A sequence shorter than the window counts as one window."""
    n = len(seq)
    if n == 0:
        return 0, 1
    w = min(window, n)
    cur = seq[:w].count("C")
    best, best_start = cur, 0
    for i in range(1, n - w + 1):
        cur += (seq[i + w - 1] == "C") - (seq[i - 1] == "C")
        if cur > best:
            best, best_start = cur, i
    return best, best_start + 1


def adjacent_cc(seq):
    """Number of positions i with seq[i] == seq[i+1] == 'C' (overlapping: CCC gives 2)."""
    return sum(1 for i in range(len(seq) - 1) if seq[i] == "C" and seq[i + 1] == "C")


def c_x_c(seq):
    """Number of positions i with seq[i] == seq[i+2] == 'C' (x is any residue)."""
    return sum(1 for i in range(len(seq) - 2) if seq[i] == "C" and seq[i + 2] == "C")


def pest_fraction(seq):
    if not seq:
        return 0.0
    return sum(seq.count(a) for a in "PSTE") / len(seq)


def assign_tier(has_sp, mature_len, max_cys, families, max_mature_len, min_cys_window):
    """Pure tier rule. `families` is a set of flags with a hit."""
    passes = mature_len <= max_mature_len and max_cys >= min_cys_window
    if not passes:
        return "other"
    if not has_sp:
        return "cys_rich_no_sp"
    if any(f in families for f in EXCLUDING_FAMILIES):
        return "cys_rich_sp_known_family"
    return "cys_rich_sp_unassigned"


def features(seq, cs_end, window):
    """Sequence features on the mature sequence (seq[cs_end:]) or the full sequence."""
    mature = seq[cs_end:] if cs_end else seq
    mx, start = max_cys_window(mature, window)
    n_c = mature.count("C")
    return {
        "mature": mature,
        "mature_length": len(mature),
        "cys_count": n_c,
        "cys_frac": (n_c / len(mature)) if mature else 0.0,
        "max_cys_window": mx,
        "window_start": start,
        "cc_pairs": adjacent_cc(mature),
        "cxc": c_x_c(mature),
        "pest_frac": pest_fraction(mature),
    }


def build_rows(proteome, records, sp, domtbl, params):
    """One row dict per protein. STOP when FASTA, SignalP and domtbl ids disagree."""
    fasta_headers = {h for h, _ in records}
    missing_in_fasta = [h for h in sp if h not in fasta_headers]
    if missing_in_fasta:
        raise Stop(
            f"{proteome}: {len(missing_in_fasta)} SignalP id(s) not in the FASTA, "
            f"first {short_id(missing_in_fasta[0])!r}"
        )
    missing_in_sp = [h for h, _ in records if h not in sp]
    if missing_in_sp:
        raise Stop(
            f"{proteome}: {len(missing_in_sp)} FASTA id(s) without a SignalP row, "
            f"first {short_id(missing_in_sp[0])!r}"
        )
    if domtbl is not None:
        shorts = {short_id(h) for h in fasta_headers}
        bad = [t for t in domtbl if t not in shorts]
        if bad:
            raise Stop(
                f"{proteome}: {len(bad)} domtbl target(s) not in the FASTA, first {bad[0]!r}"
            )
    rows = []
    for header, seq in records:
        call, prob, cs_end = sp[header]
        has_sp = call == "SP"
        if has_sp and cs_end > len(seq):
            raise Stop(f"{proteome}: CS position {cs_end} beyond length {len(seq)}")
        ft = features(seq, cs_end if has_sp else None, params["window"])
        sid = short_id(header)
        fam_hits = domtbl.get(sid, {}) if domtbl is not None else {}
        row = {
            "proteome": proteome,
            "protein_id": sid,
            "header": header,
            "length": len(seq),
            "sp_call": call,
            "sp_prob": f"{prob:.6g}",
            "cs_end": cs_end if has_sp else "NA",
            "mature_length": ft["mature_length"],
            "cys_count": ft["cys_count"],
            "cys_frac": f"{ft['cys_frac']:.4f}",
            "max_cys_window": ft["max_cys_window"],
            "window_start": ft["window_start"],
            "cc_pairs": ft["cc_pairs"],
            "cxc": ft["cxc"],
            "pest_frac": f"{ft['pest_frac']:.4f}",
        }
        for fam in FAMILIES:
            if domtbl is None or fam in params.get("missing_models", ()):
                # No search was run for this family: report NA, not 0.
                row[fam] = "NA"
                row[f"{fam}_evalue"] = "NA"
                continue
            row[fam] = 1 if fam in fam_hits else 0
            row[f"{fam}_evalue"] = f"{fam_hits[fam]:.3g}" if fam in fam_hits else "NA"
        row["tier"] = assign_tier(
            has_sp,
            ft["mature_length"],
            ft["max_cys_window"],
            {f for f in FAMILIES if fam_hits.get(f) is not None},
            params["max_mature_len"],
            params["min_cys_window"],
        )
        rows.append(row)
    return rows


# -------------------------------------------------------------------------- output


def write_gz(path, text):
    """Write gzip without file name or time stamp (gzip -n equivalent)."""
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            gz.write(text.encode("utf-8"))


def rows_to_text(rows):
    lines = ["\t".join(COLUMNS)]
    for r in rows:
        lines.append("\t".join(str(r[c]) for c in COLUMNS))
    return "\n".join(lines) + "\n"


def git_commit(repo):
    try:
        res = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return res.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "NA"


def read_manifest(path):
    """Manifest TSV: name<TAB>fasta<TAB>signalp<TAB>domtbl ('-' or empty for none)."""
    jobs = []
    names = set()
    if not Path(path).is_file():
        raise Stop(f"manifest not found: {path}")
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 3:
                raise Stop(f"manifest row needs name, fasta, signalp[, domtbl]: {line!r}")
            name = c[0]
            if name in names:
                raise Stop(f"duplicate proteome name {name!r} in manifest")
            names.add(name)
            dom = c[3] if len(c) > 3 and c[3] not in ("", "-") else None
            jobs.append((name, c[1], c[2], dom))
    if not jobs:
        raise Stop(f"no proteomes in manifest {path}")
    return jobs


def read_provenance(path, missing_models):
    """Read the pfam_provenance.json written by 01_known_family_hmm.sh. STOP if absent or
    if a model needed for a searched family is not listed."""
    if not path or not Path(path).is_file():
        raise Stop(
            "a domtbl was given but pfam_provenance.json is missing; "
            "pass --pfam-provenance (written by 01_known_family_hmm.sh)"
        )
    try:
        prov = json.loads(Path(path).read_text(encoding="utf-8"))
        accs = {m["acc"].split(".")[0] for m in prov["models"]}
        for key in ("pfam_release_dir", "pfam_hmm_resolved", "known_families_sha256"):
            if not prov[key]:
                raise KeyError(key)
    except (ValueError, KeyError, TypeError) as exc:
        raise Stop(f"bad pfam provenance file {path}: {exc!r}") from exc
    for acc, fam in FAMILY_OF_PFAM.items():
        if acc not in accs and fam not in missing_models:
            raise Stop(f"{path}: model {acc} ({fam}) is not in the provenance model list")
    return prov


def run(jobs, out_dir, params, repo, provenance=None):
    out_dir = Path(out_dir)
    all_rows = []
    per = {}
    inputs = {}
    for name, fasta, sp_path, dom_path in jobs:
        records = read_fasta(fasta)
        sp = read_signalp(sp_path)
        dom = read_domtbl(dom_path) if dom_path else None
        rows = build_rows(name, records, sp, dom, params)
        per[name] = rows
        all_rows.extend(rows)
        inputs[name] = {
            "fasta": {"path": str(fasta), "sha256": sha256_of(fasta)},
            "signalp": {"path": str(sp_path), "sha256": sha256_of(sp_path)},
            "domtbl": (
                {"path": str(dom_path), "sha256": sha256_of(dom_path)} if dom_path else None
            ),
        }
    cands = [r for r in all_rows if r["tier"] != "other"]
    summary = ["proteome\tn_proteins\tn_sp\t" + "\t".join(TIERS)]
    for name, rows in per.items():
        counts = {t: 0 for t in TIERS}
        for r in rows:
            counts[r["tier"]] += 1
        n_sp = sum(1 for r in rows if r["sp_call"] == "SP")
        summary.append(f"{name}\t{len(rows)}\t{n_sp}\t" + "\t".join(str(counts[t]) for t in TIERS))
    run_json = {
        "pfam": provenance,
        "inputs": inputs,
        "parameters": params,
        "git_commit": git_commit(repo),
        "python": sys.version.split()[0],
        "library": "none",
    }
    finals = {}
    for name, rows in per.items():
        finals[out_dir / f"{name}.tsv.gz"] = ("gz", rows_to_text(rows))
    finals[out_dir / "candidates.tsv.gz"] = ("gz", rows_to_text(cands))
    if out_dir.is_dir():
        stale = sorted(f.name for f in out_dir.glob("*.tsv.gz") if f not in finals)
        if stale:
            raise Stop(
                f"stale per-proteome files from an earlier run in {out_dir}: "
                + ", ".join(stale)
                + "; remove them or use a new --out-dir"
            )
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        probe = out_dir / ".write_test.tmp"
        probe.write_text("")
        probe.unlink()
    except OSError as exc:
        raise Stop(f"cannot write to --out-dir {out_dir}: {exc}") from exc
    finals[out_dir / "summary.tsv"] = ("txt", "\n".join(summary) + "\n")
    finals[out_dir / "run.json"] = (
        "txt",
        json.dumps(run_json, indent=2, sort_keys=True) + "\n",
    )
    tmps = {}
    try:
        for final, (kind, text) in finals.items():
            tmp = final.with_name(final.name + ".tmp")
            if kind == "gz":
                write_gz(tmp, text)
            else:
                tmp.write_text(text, encoding="utf-8")
            tmps[tmp] = final
        for tmp, final in tmps.items():
            os.replace(tmp, final)
    except BaseException:
        for tmp in tmps:
            tmp.unlink(missing_ok=True)
        raise
    return per, cands


def parse_args(argv):
    ap = StopArgParser(
        description="Tier short, Cys-rich, secreted proteins. Read-only; no model.",
    )
    ap.add_argument("--manifest", help="TSV: name, fasta, signalp, domtbl ('-' for none)")
    ap.add_argument("--name", help="proteome name (single-proteome mode)")
    ap.add_argument("--fasta", help="protein FASTA, plain or .gz (single-proteome mode)")
    ap.add_argument("--signalp", help="SignalP 6 prediction_results.txt (single mode)")
    ap.add_argument("--domtbl", help="hmmsearch --domtblout file, optional (single mode)")
    ap.add_argument("--out-dir", required=True, help="output directory")
    ap.add_argument(
        "--max-mature-len",
        type=int,
        default=DEFAULT_MAX_MATURE_LEN,
        help="L: keep proteins with mature length <= L (default %(default)s)",
    )
    ap.add_argument(
        "--min-cys-window",
        type=int,
        default=DEFAULT_MIN_CYS_WINDOW,
        help="K: keep proteins with >= K Cys in some window (default %(default)s)",
    )
    ap.add_argument(
        "--window",
        type=int,
        default=DEFAULT_WINDOW,
        help="W: window size in residues (default %(default)s)",
    )
    ap.add_argument(
        "--missing-models",
        default="",
        help="comma list of families with no HMM search (for example pra3_like_family); "
        "their flag columns are NA",
    )
    ap.add_argument(
        "--pfam-provenance",
        help="pfam_provenance.json from 01_known_family_hmm.sh (required with a domtbl)",
    )
    ap.add_argument("--repo", default=".", help="git checkout recorded in run.json")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        if args.manifest:
            jobs = read_manifest(args.manifest)
        elif args.name and args.fasta and args.signalp:
            jobs = [(args.name, args.fasta, args.signalp, args.domtbl)]
        else:
            raise Stop("give --manifest, or --name with --fasta and --signalp")
        if args.window < 1 or args.max_mature_len < 1 or args.min_cys_window < 1:
            raise Stop("--window, --max-mature-len and --min-cys-window must be >= 1")
        missing = sorted(m for m in args.missing_models.split(",") if m)
        unknown = [m for m in missing if m not in FAMILIES]
        if unknown:
            raise Stop(f"--missing-models: unknown family {unknown[0]!r}; use {FAMILIES}")
        params = {
            "max_mature_len": args.max_mature_len,
            "min_cys_window": args.min_cys_window,
            "window": args.window,
            "missing_models": missing,
        }
        prov = None
        if any(j[3] for j in jobs):
            prov = read_provenance(args.pfam_provenance, missing)
        run(jobs, args.out_dir, params, args.repo, prov)
    except Stop as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
