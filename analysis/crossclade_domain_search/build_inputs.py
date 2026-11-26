"""Build the inputs of the cross-clade domain search (strategy A of docs/reports/2026-11-25-plm-survey-rhodotorula-adhesion.md).

Runs once on the HPCC login node (a few dozen small HTTPS requests, no heavy work) and writes into a NEW directory outside the repo:

  queries.faa         the Pfam domain of each query in queries.tsv, cut at the InterPro envelope +- --pad residues
  decoys.faa          per query: the crop reversed (`dec00_<id>`) and --n-shuffles local-window shuffles (`decNN_<id>`, seeded by the id);
                      searched exactly like the queries, they give the empirical false-hit rate (shuffles keep local composition and length)
  query_domains.tsv   query id, panel, accession, gene, pfam, domain start/end, crop start/end, number of Pfam locations, protein length
  targets.faa         the target proteomes, taken from an existing staged proteome directory (proteins.tsv + chunks/), never re-downloaded
  target_proteins.tsv accession, species, length
  truth_pf21671.tsv   UniProt accessions of each target species with Pfam PF21671 (CPL1-like); the reference for the basidio panel
  provenance.txt, SHA256SUMS

  python3 -I analysis/crossclade_domain_search/build_inputs.py --queries analysis/crossclade_domain_search/queries.tsv \
      --staged /path/nonyeast_test/in --species Rtor_NP11=1130832 --species Cneo_H99=235443 --out /path/NEW_DIR

Stdlib only. Fails loudly (FATAL:, exit 1): existing --out, --out inside the repo, a query without its Pfam domain, a staged species
missing, a sequence that differs in length from InterPro's record.
"""

import argparse
import csv
import hashlib
import json
import random
import sys
import urllib.parse
import urllib.request
from pathlib import Path

UNIPROT = "https://rest.uniprot.org/uniprotkb"
INTERPRO = "https://www.ebi.ac.uk/interpro/api"
UA = {"User-Agent": "crossclade-domain-search/1 (python urllib)"}
REPO = Path(__file__).resolve().parents[2]


class BuildError(Exception):
    pass


def http_get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode(), dict(r.headers)


def parse_fasta_text(text):
    """Return [(id, seq)]; id is the first word after '>' (UniProt sp|ACC|NAME is reduced to ACC)."""
    out, name, seq = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(seq)))
            word = line[1:].split()[0] if line[1:].strip() else ""
            parts = word.split("|")
            name = parts[1] if len(parts) >= 3 and parts[0] in ("sp", "tr") else word
            seq = []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def pick_location(locations):
    """Pick the Pfam location to crop: the best (lowest) score; ties go to the N-terminal one. Returns (start, end, n_locations)."""
    if not locations:
        raise BuildError("no Pfam location")
    frags = []
    for loc in locations:
        fs = loc["fragments"]
        frags.append(
            (
                float(loc["score"]) if loc.get("score") is not None else float("inf"),
                min(f["start"] for f in fs),
                max(f["end"] for f in fs),
            )
        )
    score, start, end = sorted(frags, key=lambda t: (t[0], t[1]))[0]
    return start, end, len(frags)


def crop(seq, start, end, pad):
    """1-based inclusive domain [start, end] widened by pad on each side and clipped to the sequence. Returns (sub, cstart, cend)."""
    if not (1 <= start <= end <= len(seq)):
        raise BuildError(f"domain {start}-{end} outside a sequence of {len(seq)}")
    cs, ce = max(1, start - pad), min(len(seq), end + pad)
    return seq[cs - 1 : ce], cs, ce


def window_shuffle(seq, window, seed):
    """Shuffle residues inside consecutive windows of `window` residues (keeps local composition, destroys order). Deterministic."""
    rng = random.Random(seed)
    out = []
    for i in range(0, len(seq), window):
        w = list(seq[i : i + window])
        rng.shuffle(w)
        out.extend(w)
    return "".join(out)


def decoys_for(qid, sub, n_shuffles, window):
    """[(decoy_id, seq)]: dec00 = reversed, dec01..decNN = window shuffles."""
    out = [(f"dec00_{qid}", sub[::-1])]
    out += [
        (f"dec{k:02d}_{qid}", window_shuffle(sub, window, f"{qid}:{k}"))
        for k in range(1, n_shuffles + 1)
    ]
    return out


def query_id(acc, pfam, cs, ce):
    return f"{acc}_{pfam}_{cs}_{ce}"


def wrap(seq, width=60):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def staged_targets(staged, species):
    """Read proteins.tsv + chunks/*.faa of an aga1_class staged directory; return {acc: (species, seq)} for the given species."""
    staged = Path(staged)
    rows = [r for r in read_tsv(staged / "proteins.tsv") if r["species"] in species]
    have = {r["species"] for r in rows}
    missing = set(species) - have
    if missing:
        raise BuildError(f"species not in {staged}/proteins.tsv: {sorted(missing)}")
    want = {r["accession"]: r for r in rows}
    out = {}
    for chunk in sorted({r["chunk"] for r in rows}):
        for acc, seq in parse_fasta_text((staged / "chunks" / f"{chunk}.faa").read_text()):
            if acc in want:
                if len(seq) != int(want[acc]["length"]):
                    raise BuildError(
                        f"{acc}: length {len(seq)} differs from proteins.tsv {want[acc]['length']}"
                    )
                out[acc] = (want[acc]["species"], seq)
    if set(out) != set(want):
        raise BuildError(
            f"{len(set(want) - set(out))} accessions of proteins.tsv not found in the chunks"
        )
    return out


def uniprot_pfam_members(taxon, pfam):
    q = urllib.parse.urlencode(
        {
            "query": f"organism_id:{taxon} AND xref:pfam-{pfam}",
            "fields": "accession",
            "format": "tsv",
            "size": 500,
        }
    )
    text, headers = http_get(f"{UNIPROT}/search?{q}")
    accs = [line.strip() for line in text.splitlines()[1:] if line.strip()]
    if headers.get("Link") or headers.get("link"):
        raise BuildError(f"more than 500 {pfam} members for taxon {taxon}; paging not implemented")
    return accs, headers.get("X-UniProt-Release") or headers.get("x-uniprot-release", "?")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--queries", required=True)
    ap.add_argument(
        "--staged", required=True, help="staged proteome directory (proteins.tsv, chunks/)"
    )
    ap.add_argument(
        "--species", action="append", required=True, help="NAME=TAXID, as named in proteins.tsv"
    )
    ap.add_argument("--out", required=True)
    ap.add_argument("--pad", type=int, default=10)
    ap.add_argument(
        "--n-shuffles",
        type=int,
        default=50,
        help="window-shuffled decoys per query (plus one reversal); at most 99",
    )
    ap.add_argument("--window", type=int, default=10, help="shuffle window in residues")
    a = ap.parse_args(argv)
    out = Path(a.out).resolve()
    try:
        if not 0 <= a.n_shuffles <= 99 or a.window < 2:
            raise BuildError("--n-shuffles must be 0..99 and --window at least 2")
        if out.exists():
            raise BuildError(f"{out} exists; nothing is overwritten")
        if out == REPO or REPO in out.parents:
            raise BuildError("--out is inside the repo; no sequence is written there")
        species = dict(s.split("=", 1) for s in a.species)

        qrows, dom_rows, qfa, dfa = read_tsv(a.queries), [], [], []
        releases = set()
        for r in qrows:
            acc, pfam = r["accession"], r["pfam"]
            text, h = http_get(f"{UNIPROT}/{acc}.fasta")
            releases.add(h.get("X-UniProt-Release") or h.get("x-uniprot-release", "?"))
            recs = parse_fasta_text(text)
            if len(recs) != 1:
                raise BuildError(f"{acc}: expected one FASTA record, got {len(recs)}")
            seq = recs[0][1]
            js = json.loads(http_get(f"{INTERPRO}/entry/pfam/{pfam}/protein/uniprot/{acc}")[0])
            prot = js["proteins"][0]
            if int(prot["protein_length"]) != len(seq):
                raise BuildError(
                    f"{acc}: UniProt length {len(seq)} != InterPro length {prot['protein_length']}"
                )
            start, end, n = pick_location(prot["entry_protein_locations"])
            sub, cs, ce = crop(seq, start, end, a.pad)
            qid = query_id(acc, pfam, cs, ce)
            qfa.append(f">{qid} {r['gene']} {r['domain']}\n{wrap(sub)}")
            dfa.extend(
                f">{did} decoy\n{wrap(dseq)}"
                for did, dseq in decoys_for(qid, sub, a.n_shuffles, a.window)
            )
            dom_rows.append(
                [
                    qid,
                    r["panel"],
                    acc,
                    r["species"],
                    r["gene"],
                    pfam,
                    start,
                    end,
                    cs,
                    ce,
                    n,
                    len(seq),
                ]
            )

        targets = staged_targets(a.staged, set(species))
        truth = []
        for sp, tax in sorted(species.items()):
            accs, rel = uniprot_pfam_members(tax, "PF21671")
            releases.add(rel)
            for acc in sorted(accs):
                truth.append(
                    [acc, sp, "PF21671", "yes" if acc in targets else "not_in_staged_proteome"]
                )

        out.mkdir(parents=True)
        (out / "queries.faa").write_text("\n".join(qfa) + "\n")
        (out / "decoys.faa").write_text("\n".join(dfa) + "\n")
        with open(out / "query_domains.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(
                [
                    "query_id",
                    "panel",
                    "accession",
                    "species",
                    "gene",
                    "pfam",
                    "dom_start",
                    "dom_end",
                    "crop_start",
                    "crop_end",
                    "n_pfam_locations",
                    "protein_length",
                ]
            )
            w.writerows(dom_rows)
        with (
            open(out / "targets.faa", "w") as fh,
            open(out / "target_proteins.tsv", "w", newline="") as th,
        ):
            w = csv.writer(th, delimiter="\t", lineterminator="\n")
            w.writerow(["accession", "species", "length"])
            for acc in sorted(targets):
                sp, seq = targets[acc]
                fh.write(f">{acc}\n{wrap(seq)}\n")
                w.writerow([acc, sp, len(seq)])
        with open(out / "truth_pf21671.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["accession", "species", "pfam", "in_targets"])
            w.writerows(truth)
        files = [
            "queries.faa",
            "decoys.faa",
            "query_domains.tsv",
            "targets.faa",
            "target_proteins.tsv",
            "truth_pf21671.tsv",
        ]
        sums = {f: hashlib.sha256((out / f).read_bytes()).hexdigest() for f in files}
        (out / "SHA256SUMS").write_text("".join(f"{v}  {k}\n" for k, v in sums.items()))
        staged_sha = Path(a.staged) / "SHA256SUMS"
        (out / "provenance.txt").write_text(
            f"script: analysis/crossclade_domain_search/build_inputs.py\n"
            f"queries: {a.queries} sha256 {hashlib.sha256(Path(a.queries).read_bytes()).hexdigest()}\n"
            f"uniprot_release: {','.join(sorted(releases))}\n"
            f"pad: {a.pad}\ndecoys: 1 reversal + {a.n_shuffles} window shuffles (window {a.window}) per query\nstaged: {Path(a.staged).resolve()} (SHA256SUMS sha256 "
            f"{hashlib.sha256(staged_sha.read_bytes()).hexdigest() if staged_sha.exists() else 'absent'})\n"
            f"species: {', '.join(f'{k}={v}' for k, v in sorted(species.items()))}\n"
            f"n_queries: {len(qfa)}\nn_targets: {len(targets)}\nn_truth_pf21671: {len(truth)}\n"
        )
        print(
            f"wrote {out}: {len(qfa)} queries, {len(targets)} target proteins, {len(truth)} PF21671 members"
        )
    except BuildError as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
