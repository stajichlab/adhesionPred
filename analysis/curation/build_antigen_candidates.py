#!/usr/bin/env python
"""Rank Coccidioides surface proteins as candidate antigens for recombinant production
(data/curated/antigens/coccidioides_candidates.tsv).

The target application is immunodiagnostics and serology: proteins the host immune system
is likely to see and make antibodies against, which can then be produced recombinantly and
tested against patient sera. IEDB curates only a handful of Coccidioides antigens, so
candidates are ranked by transferable evidence instead:

  +3  homologous to a protein with curated immune-assay evidence in IEDB (any fungus)
  +2  known Coccidioides antigen already in IEDB (positive control for the ranking)
  +2  conserved between C. immitis and C. posadasii (a diagnostic must detect both)
  +1  GPI-anchored or cell-wall (displayed on the spherule/hypha surface)
  +1  has a signal peptide (secreted: reaches host fluids, so detectable in serology)
  +1  no close human homolog (reduces antibody cross-reactivity, and self-tolerance
      would otherwise blunt the response)
  -2  close human homolog (conserved housekeeping proteins such as enolase or HSP70 are
      common false leads: immunogenic but cross-reactive and non-specific)

This ranking is a hypothesis generator for wet-lab prioritization, not a claim that any
protein is protective or a vaccine target; those require immunological testing.

Requires MMseqs2 on PATH (or --mmseqs).

Usage:
    python analysis/curation/build_antigen_candidates.py [--mmseqs /path/to/mmseqs]
"""

import argparse
import csv
import subprocess
import tempfile
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

from curation_lib import CURATED, get, log, write_table

CUR = CURATED / "antigens"
COCCI = {"Coccidioides immitis RS": "Cimm", "Coccidioides posadasii C735 delta SOWgp": "Cpos"}
HUMAN_PROTEOME = "UP000005640"


def fetch_fasta(accs, out_path):
    accs = sorted({a for a in accs if a})
    with open(out_path, "w") as f:
        for i in range(0, len(accs), 90):
            q = urllib.parse.urlencode(
                {"accessions": ",".join(accs[i : i + 90]), "format": "fasta"}
            )
            f.write(get(f"https://rest.uniprot.org/uniprotkb/accessions?{q}", "text/plain"))
    return out_path


def fetch_human(out_path):
    log("downloading human reference proteome (cross-reactivity screen)...")
    url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=" + urllib.parse.quote(
        f"proteome:{HUMAN_PROTEOME}"
    )
    with open(out_path, "wb") as f:
        f.write(urllib.request.urlopen(url, timeout=600).read())
    return out_path


def search(mmseqs, query, target, out, tmp, min_id, cov):
    subprocess.run(
        [
            mmseqs,
            "easy-search",
            str(query),
            str(target),
            str(out),
            str(tmp),
            "--min-seq-id",
            str(min_id),
            "-c",
            str(cov),
            "--cov-mode",
            "0",
            "-v",
            "1",
            "--format-output",
            "query,target,pident,evalue,bits",
        ],
        check=True,
    )
    best = {}
    for line in open(out):
        q, t, pid, ev, bits = line.rstrip("\n").split("\t")
        acc = q.split("|")[1] if "|" in q else q
        tacc = t.split("|")[1] if "|" in t else t
        if acc not in best or float(bits) > float(best[acc][2]):
            best[acc] = (tacc, float(pid), float(bits), float(ev))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mmseqs", default="mmseqs")
    ap.add_argument(
        "--skip-human", action="store_true", help="skip the human cross-reactivity screen"
    )
    args = ap.parse_args()

    surface = list(csv.DictReader(open(CURATED / "surface" / "surface.tsv"), delimiter="\t"))
    cocci = [r for r in surface if r["genome"] in COCCI]
    antigens = {
        r["accession"]: r for r in csv.DictReader(open(CUR / "antigens.tsv"), delimiter="\t")
    }
    log(f"{len(cocci)} Coccidioides surface proteins; {len(antigens)} curated fungal antigens")

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        q = fetch_fasta([r["accession"] for r in cocci], td / "cocci.fa")
        a = fetch_fasta(antigens, td / "antigens.fa")

        log("searching against curated fungal antigens...")
        antigen_hit = search(args.mmseqs, q, a, td / "ag.m8", td / "t1", 0.3, 0.5)

        log("cross-species conservation (C. immitis <-> C. posadasii)...")
        by_sp = defaultdict(list)
        for r in cocci:
            by_sp[COCCI[r["genome"]]].append(r["accession"])
        fetch_fasta(by_sp["Cimm"], td / "cimm.fa")
        fetch_fasta(by_sp["Cpos"], td / "cpos.fa")
        # Each species is searched against the OTHER one only; searching the combined set
        # against either species would match every protein to itself.
        cons = search(
            args.mmseqs, td / "cimm.fa", td / "cpos.fa", td / "c1.m8", td / "t2", 0.8, 0.8
        )
        cons.update(
            search(args.mmseqs, td / "cpos.fa", td / "cimm.fa", td / "c2.m8", td / "t3", 0.8, 0.8)
        )

        human = {}
        if not args.skip_human:
            fetch_human(td / "human.fa")
            log("human cross-reactivity screen...")
            human = search(args.mmseqs, q, td / "human.fa", td / "hs.m8", td / "t4", 0.4, 0.5)

    rows = []
    for r in cocci:
        acc = r["accession"]
        known = acc in antigens
        hom = antigen_hit.get(acc)
        homolog_of_antigen = bool(hom) and not known
        conserved = acc in cons
        hs = human.get(acc)
        score = (
            (2 if known else 0)
            + (3 if homolog_of_antigen else 0)
            + (2 if conserved else 0)
            + (1 if r["gpi_anchor"] == "yes" else 0)
            + (1 if r["signal_peptide"] == "yes" else 0)
            + (-2 if hs else (0 if args.skip_human else 1))
        )
        rows.append(
            {
                "accession": acc,
                "gene": r["gene"],
                "protein_name": r["protein_name"],
                "genome": r["genome"],
                "candidate_score": score,
                "known_iedb_antigen": "yes" if known else "no",
                "iedb_assays": r["antigen_assays"],
                "antigen_homolog": hom[0] if homolog_of_antigen else "",
                "antigen_homolog_pident": f"{hom[1]:.1f}" if homolog_of_antigen else "",
                "conserved_both_species": "yes" if conserved else "no",
                "human_homolog": hs[0] if hs else "",
                "human_homolog_pident": f"{hs[1]:.1f}" if hs else "",
                "adhesion_status": r["adhesion_status"],
                "length": r["length"],
                "signal_peptide": r["signal_peptide"],
                "gpi_anchor": r["gpi_anchor"],
                "ec": r["ec"],
                "pfam": r["pfam"],
                "needs_review": "yes",
            }
        )

    cols = list(rows[0].keys())
    write_table(
        rows,
        CUR / "coccidioides_candidates.tsv",
        cols,
        lambda r: (-r["candidate_score"], r["accession"]),
    )


if __name__ == "__main__":
    main()
