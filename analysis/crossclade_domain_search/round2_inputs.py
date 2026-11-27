"""Build the inputs of round 2 of the cross-clade domain search: AlphaFold DB models for structure search with real models.

Runs once on the HPCC login node (about 1,200 small HTTPS downloads from the AlphaFold DB; no heavy work) and writes a NEW
directory outside the repo (built as <out>.partial and renamed only when complete):

  queries/<query_id>.pdb   the AFDB model of each round 1 query protein, cut to the round 1 crop (crop_start..crop_end); the crop's
                           sequence must equal the round 1 query sequence (queries.faa)
  targets/<acc>.pdb        AFDB target models with residues of pLDDT < --min-plddt removed (long low-confidence disordered regions give
                           spurious structure similarity); a model with fewer than --min-residues residues left is not used. Targets:
                             - candidate:  every target hit by any query at E <= 0.1 in round 1 (hits.tsv, either method)
                             - pf21671:    the PF21671 reference members among the targets (minus query accessions)
                             - background: per species, proteins drawn in a seeded random order from the WHOLE proteome (minus queries
                               and PF21671 members; a candidate can be drawn too) until 2 x --background-per-species models are usable;
                               they alternate between half A (the null of the structure methods) and half B (held out to measure how
                               often a protein with no known relation passes the rule)
  models.tsv               id, role, bg_half, species, afdb entry and version, residues, residues kept, mean pLDDT, status
  provenance.txt, SHA256SUMS

  python3 -I analysis/crossclade_domain_search/round2_inputs.py --round1-inputs /path/round1_in \
      --round1-hits /path/round1_summary/hits.tsv --background-per-species 250 --out /path/NEW_DIR

A model is used only if its sequence equals the round 1 sequence (status ok); others are listed as no_model, sequence_mismatch,
numbering or too_short. Stdlib only. Fails loudly (FATAL:, exit 1) on an existing --out or --out.partial, an --out inside the repo,
missing round 1 files, a network failure after retries, or a species with too few usable background models.
"""

import argparse
import csv
import hashlib
import json
import random
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

AFDB = "https://alphafold.ebi.ac.uk/api/prediction"
UA = {"User-Agent": "crossclade-domain-search/2 (python urllib)"}
REPO = Path(__file__).resolve().parents[2]
RETRIES = 4


class BuildError(Exception):
    pass


class NotFoundError(Exception):
    pass


def http_get(url, binary=False, retries=RETRIES, backoff=2.0):
    """GET with retries on network errors and 5xx/429; 404/422 raise NotFoundError; other failures raise BuildError."""
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            return data if binary else data.decode()
        except urllib.error.HTTPError as e:
            if e.code in (404, 422):
                raise NotFoundError(url) from e
            if e.code not in (429, 500, 502, 503, 504) or attempt == retries:
                raise BuildError(f"{url}: HTTP {e.code}") from e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            if attempt == retries:
                raise BuildError(f"{url}: {e}") from e
        time.sleep(backoff * 2**attempt)
    raise BuildError(f"{url}: no response")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path):
    out, name, seq = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(seq)
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        out[name] = "".join(seq)
    return out


def afdb_entry(acc, fetch=http_get):
    """The F1 entry of the AFDB prediction API for acc (legacy or current field names), or None when AFDB has no model."""
    try:
        entries = json.loads(fetch(f"{AFDB}/{acc}"))
    except NotFoundError:
        return None
    for e in entries:
        if (e.get("entryId") or e.get("modelEntityId")) == f"AF-{acc}-F1":
            return e
    return None


def entry_fields(e):
    return {
        "entry": e.get("entryId") or e.get("modelEntityId") or "NA",
        "version": e.get("latestVersion", "NA"),
        "seq": e.get("uniprotSequence") or e.get("sequence") or "",
        "start": int(e.get("uniprotStart") or e.get("sequenceStart") or 1),
        "pdb_url": e.get("pdbUrl"),
    }


def _ca(pdb_text):
    return [
        ln for ln in pdb_text.splitlines() if ln.startswith("ATOM") and ln[12:16].strip() == "CA"
    ]


def pdb_residues(pdb_text):
    """Residue numbers that have a CA atom, in file order."""
    return [int(ln[22:26]) for ln in _ca(pdb_text)]


def mean_plddt(pdb_text):
    """Mean B-factor (pLDDT in AFDB models) over CA atoms; None for no CA."""
    vals = [float(ln[60:66]) for ln in _ca(pdb_text)]
    return round(sum(vals) / len(vals), 1) if vals else None


def crop_pdb(pdb_text, start, end):
    """Keep ATOM records of residues start..end (1-based, inclusive) and close the chain."""
    keep = [
        ln
        for ln in pdb_text.splitlines()
        if ln.startswith("ATOM") and start <= int(ln[22:26]) <= end
    ]
    if not keep:
        raise BuildError(f"no atoms in residues {start}-{end}")
    return "\n".join(keep + ["TER", "END"]) + "\n"


def trim_low_plddt(pdb_text, min_plddt):
    """Drop residues whose CA pLDDT is below min_plddt; returns (pdb_text, residues kept)."""
    keep_res = {int(ln[22:26]) for ln in _ca(pdb_text) if float(ln[60:66]) >= min_plddt}
    keep = [
        ln for ln in pdb_text.splitlines() if ln.startswith("ATOM") and int(ln[22:26]) in keep_res
    ]
    return "\n".join(keep + ["TER", "END"]) + "\n", len(keep_res)


def background_order(targets, exclude, seed):
    """{species: [accessions in seeded random order]} over the whole proteome minus exclude."""
    by_sp = {}
    for r in targets:
        if r["accession"] not in exclude:
            by_sp.setdefault(r["species"], []).append(r["accession"])
    out = {}
    for sp, pool in sorted(by_sp.items()):
        pool = sorted(pool)
        random.Random(f"{seed}:{sp}").shuffle(pool)
        out[sp] = pool
    return out


def plan_scored(hits, targets, truth, query_accs):
    """{accession: role} for the scored targets: round 1 candidates and PF21671 reference members."""
    sp = {r["accession"] for r in targets}
    pf = {r["accession"] for r in truth if r["in_targets"] == "yes"} - query_accs
    cand = {h["target"] for h in hits} - query_accs
    for a in cand | pf:
        if a not in sp:
            raise BuildError(f"{a} is not a round 1 target")
    return {a: ("pf21671" if a in pf else "candidate") for a in sorted(cand | pf)}


def fetch_model(acc, expected_seq, fetch=http_get):
    """(status, fields, pdb_text). status: ok, no_model, sequence_mismatch, numbering."""
    e = afdb_entry(acc, fetch)
    if e is None:
        return "no_model", None, None
    f = entry_fields(e)
    if f["seq"] != expected_seq:
        return "sequence_mismatch", f, None
    if f["start"] != 1:
        return "numbering", f, None
    pdb = fetch(f["pdb_url"], binary=True).decode()
    if pdb_residues(pdb) != list(range(1, len(f["seq"]) + 1)):
        return "sequence_mismatch", f, None
    return "ok", f, pdb


def main(argv=None, fetch=http_get, pause=0.1):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--round1-inputs", required=True)
    ap.add_argument("--round1-hits", required=True)
    ap.add_argument(
        "--background-per-species",
        type=int,
        default=250,
        help="per half; 2x this many per species in total",
    )
    ap.add_argument("--min-plddt", type=float, default=50.0)
    ap.add_argument("--min-residues", type=int, default=40)
    ap.add_argument("--seed", default="round2")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    try:
        out = Path(a.out).resolve()
        part = out.with_name(out.name + ".partial")
        for p in (out, part):
            if p.exists():
                raise BuildError(f"{p} exists; nothing is overwritten")
        if out == REPO or REPO in out.parents:
            raise BuildError(f"{out} is inside the repo; write outside it")
        r1 = Path(a.round1_inputs)
        domains = read_tsv(r1 / "query_domains.tsv")
        targets = read_tsv(r1 / "target_proteins.tsv")
        truth = read_tsv(r1 / "truth_pf21671.tsv")
        seqs = read_fasta(r1 / "targets.faa")
        qseqs = read_fasta(r1 / "queries.faa")
        hits = read_tsv(a.round1_hits)
        species = {r["accession"]: r["species"] for r in targets}
        for r in targets:
            if r["accession"] not in seqs:
                raise BuildError(f"{r['accession']} has no sequence in targets.faa")
        qaccs = {d["accession"] for d in domains}
        scored = plan_scored(hits, targets, truth, qaccs)
        pf = {r["accession"] for r in truth if r["in_targets"] == "yes"}
        (part / "queries").mkdir(parents=True)
        (part / "targets").mkdir()
        rows, done = [], {}

        def target(acc):
            """Fetch, trim and write one target model once; returns its status."""
            if acc in done:
                return done[acc]
            st, f, pdb = fetch_model(acc, seqs[acc], fetch)
            kept, plddt = "NA", "NA"
            if st == "ok":
                plddt = mean_plddt(pdb)
                trimmed, kept = trim_low_plddt(pdb, a.min_plddt)
                if kept < a.min_residues:
                    st = "too_short"
                else:
                    (part / "targets" / f"{acc}.pdb").write_text(trimmed)
            done[acc] = st
            rows.append(
                {
                    "id": acc,
                    "role": scored.get(acc, "background"),
                    "bg_half": "",
                    "species": species[acc],
                    "afdb_entry": f["entry"] if f else "NA",
                    "afdb_version": f["version"] if f else "NA",
                    "residues": len(seqs[acc]),
                    "residues_kept": kept,
                    "mean_plddt": plddt,
                    "status": st,
                }
            )
            time.sleep(pause)
            return st

        for d in domains:
            f = None
            e = afdb_entry(d["accession"], fetch)
            cs, ce = int(d["crop_start"]), int(d["crop_end"])
            plddt = "NA"
            if e is None:
                st = "no_model"
            else:
                f = entry_fields(e)
                if f["seq"][cs - 1 : ce] != qseqs.get(d["query_id"]):
                    st = "sequence_mismatch"
                elif f["start"] != 1:
                    st = "numbering"
                else:
                    pdb = fetch(f["pdb_url"], binary=True).decode()
                    if pdb_residues(pdb) != list(range(1, len(f["seq"]) + 1)):
                        st = "sequence_mismatch"
                    else:
                        st = "ok"
                        cropped = crop_pdb(pdb, cs, ce)
                        plddt = mean_plddt(cropped)
                        (part / "queries" / f"{d['query_id']}.pdb").write_text(cropped)
            rows.append(
                {
                    "id": d["query_id"],
                    "role": "query",
                    "bg_half": "",
                    "species": d["species"],
                    "afdb_entry": f["entry"] if f else "NA",
                    "afdb_version": f["version"] if f else "NA",
                    "residues": ce - cs + 1,
                    "residues_kept": ce - cs + 1 if st == "ok" else "NA",
                    "mean_plddt": plddt,
                    "status": st,
                }
            )
            time.sleep(pause)
        for acc in scored:
            target(acc)
        halves = {}
        for sp, order in background_order(targets, qaccs | pf, a.seed).items():
            n_ok = 0
            for acc in order:
                if n_ok == 2 * a.background_per_species:
                    break
                if target(acc) == "ok":
                    halves[acc] = "A" if n_ok % 2 == 0 else "B"
                    n_ok += 1
            if n_ok < 2 * a.background_per_species:
                raise BuildError(
                    f"{sp}: only {n_ok} usable background models for {2 * a.background_per_species}"
                )
        for r in rows:
            r["bg_half"] = halves.get(r["id"], "")
        with open(part / "models.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        counts = {}
        for r in rows:
            k = (r["role"] + (f"/{r['bg_half']}" if r["bg_half"] else ""), r["status"])
            counts[k] = counts.get(k, 0) + 1
        (part / "provenance.txt").write_text(
            "\n".join(
                [
                    "script: analysis/crossclade_domain_search/round2_inputs.py",
                    f"round1_inputs: {r1.resolve()}",
                    f"round1_hits: {Path(a.round1_hits).resolve()} sha256 "
                    + hashlib.sha256(Path(a.round1_hits).read_bytes()).hexdigest(),
                    f"background_per_species_per_half: {a.background_per_species}; seed: {a.seed}",
                    f"min_plddt: {a.min_plddt}; min_residues: {a.min_residues}",
                    f"afdb_api: {AFDB}",
                    *(f"{role} {st}: {n}" for (role, st), n in sorted(counts.items())),
                ]
            )
            + "\n"
        )
        files = sorted(p for p in part.rglob("*") if p.is_file())
        (part / "SHA256SUMS").write_text(
            "".join(
                f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(part)}\n"
                for p in files
            )
        )
        part.rename(out)
        ok = sum(1 for r in rows if r["status"] == "ok")
        print(f"wrote {out}: {ok}/{len(rows)} models ok")
    except BuildError as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
