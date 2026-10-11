#!/usr/bin/env python3.12
"""Build the lectin validation truth table and the sequence file.

Inputs are files fetched by fetch_unilectin.py, fetch_uniprot.py and select_assumed.py, the
hand-curated table curated_positives.tsv, and the shared hard negatives of the repo.
Output (--out-dir): truth_all.tsv, sequences.faa, rejected.tsv (appended by this script).
Leakage is set to `unknown` here and filled by finalize_truth.py (after seed identity).

Classes (column `label_class`):
  positive classes: structure_backed, host_binding_h1a, chitin_masking_h1b, literature_binding,
                    known_gap, pa14_owned_elsewhere
  negative classes: er_housekeeping, enzyme_lectin_like_domain, n1_shared, assumed
Rules are in LOG.md. No hmmsearch result is used here.
"""

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict

DATE = "2026-10-11"
CURATOR = "claude-sonnet-5-5"
NAMES = "/srv/projects/db/taxonomy/names.dmp"
NODES = "/srv/projects/db/taxonomy/nodes.dmp"
SUGAR = re.compile(r"(pyranos|furanos|neuraminic)", re.I)
BIND = re.compile(
    r"\b(binds?|binding|specific(ally)? (for|to)|recogni[sz]es?|agglutinat\w+|hemagglutin\w+)\b",
    re.I,
)
SUGARW = re.compile(
    r"(chitin|fucos|lactos|galactos|mannos|glcnac|acetylglucosamine|oligosaccharide|glycan|sialic|carbohydrate|sugar)",
    re.I,
)
M3_PF = {"PF03388", "PF00262", "PF07915", "PF02157", "PF11721"}

COLS = [
    "accession",
    "gene",
    "species",
    "taxon_id",
    "label",
    "stratum",
    "evidence_level",
    "source",
    "quote",
    "tuned",
    "leakage",
    "curator",
    "date",
    "notes",
    "label_class",
    "tier",
    "h_flag",
    "protein_name",
    "length",
    "strain_taxid",
    "unilectin_class",
    "pdb_ids",
    "pdb_sugar_comps",
    "uniprot_status",
    "negative_evidence",
    "phylum",
]


def load_taxonomy(taxids, names_wanted):
    """Return (species_of[taxid] -> (species_taxid, name), sciname_to_taxid)."""
    parent, rank = {}, {}
    for line in open(NODES):
        f = line.split("\t|\t")
        parent[int(f[0])] = int(f[1])
        rank[int(f[0])] = f[2]
    sci = {}
    byname = {}
    for line in open(NAMES):
        f = [x.strip() for x in line.rstrip("\n").rstrip("|").split("\t|\t")]
        if f[3].startswith("scientific name"):
            sci[int(f[0])] = f[1]
        if f[1] in names_wanted and f[3].startswith("scientific name"):
            byname[f[1]] = int(f[0])

    def phylum(t):
        x, seen = t, 0
        while x in parent and rank.get(x) != "phylum" and x != parent[x] and seen < 60:
            x = parent[x]
            seen += 1
        return sci.get(x, "") if rank.get(x) == "phylum" else ""

    out = {}
    for t in taxids:
        x = t
        seen = 0
        while x in parent and rank.get(x) != "species" and x != parent[x] and seen < 40:
            x = parent[x]
            seen += 1
        out[t] = (x, sci.get(x, "")) if rank.get(x) == "species" else (t, sci.get(t, ""))
    return out, byname, phylum


def read_fasta(path):
    seqs, cur = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line.strip()
    return seqs


def uniprot_fields(e):
    fn = [c for c in e.get("comments", []) if c["commentType"] == "FUNCTION"]
    text, codes, pm = "", set(), set()
    for c in fn:
        for t in c.get("texts", []):
            text += (" " if text else "") + t["value"]
            for ev in t.get("evidences", []):
                codes.add(ev["evidenceCode"])
                if ev.get("source") == "PubMed":
                    pm.add(ev["id"])
    pd = e.get("proteinDescription", {})
    name = pd.get("recommendedName", {}).get("fullName", {}).get("value") or (
        pd.get("submissionNames") or [{}]
    )[0].get("fullName", {}).get("value", "")
    genes = e.get("genes", [])
    gene = ""
    if genes:
        g = genes[0]
        gene = g.get("geneName", {}).get("value") or (g.get("orfNames") or [{}])[0].get("value", "")
    pf = sorted({x["id"] for x in e.get("uniProtKBCrossReferences", []) if x["database"] == "Pfam"})
    return {
        "text": text,
        "codes": codes,
        "pmids": pm,
        "name": name,
        "gene": gene,
        "pfam": pf,
        "taxon": e["organism"]["taxonId"],
        "length": e["sequence"]["length"],
        "reviewed": e["entryType"].startswith("UniProtKB reviewed"),
        "seq": e["sequence"]["value"],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="analysis/lectin_validation/data")
    ap.add_argument("--curated", required=True)
    ap.add_argument("--unilectin-tsv", required=True)
    ap.add_argument("--unilectin-faa", required=True)
    ap.add_argument("--shared", required=True, help="data/controls/hard-negatives-shared")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    ul_rows = {r["uniprot"]: r for r in csv.DictReader(open(a.unilectin_tsv), delimiter="\t")}
    docs = json.load(gzip.open(a.data + "/unilectin_docs.json.gz", "rt"))
    rcsb = json.load(gzip.open(a.data + "/rcsb_ligands.json.gz", "rt"))
    pdbs = defaultdict(list)
    ul_pmids = defaultdict(set)
    for d in docs:
        pdbs[d["uniprot"]].append(d["pdb"])
        if d.get("pubmed") and d["pubmed"] != "NULL":
            ul_pmids[d["uniprot"]].add(d["pubmed"])
    uni = {
        e["primaryAccession"]: uniprot_fields(e)
        for e in json.load(gzip.open(a.data + "/raw/uniprot_candidates.json.gz", "rt"))
    }
    uni_assumed = {
        e["primaryAccession"]: uniprot_fields(e)
        for e in json.load(gzip.open(a.data + "/raw/uniprot_assumed.json.gz", "rt"))
    }
    ul_seq = read_fasta(a.unilectin_faa)
    curated = {r["accession"]: r for r in csv.DictReader(open(a.curated), delimiter="\t")}
    shared = list(csv.DictReader(open(a.shared + "/controls.tsv"), delimiter="\t"))
    shared_seq = read_fasta(a.shared + "/sequences.faa")

    def sugar_comps(acc):
        names = set()
        for p in pdbs.get(acc, []):
            for n in (rcsb.get(p) or {}).get("nonpolymer_entities") or []:
                nm = n["pdbx_entity_nonpoly"]["name"]
                if SUGAR.search(nm):
                    names.add(n["pdbx_entity_nonpoly"]["comp_id"])
        return sorted(names)

    rows, seqs, rejected = [], {}, []

    def reject(acc, reason, where=""):
        rejected.append({"accession": acc, "candidate_class": where, "reason": reason})

    # taxonomy: collect taxids first
    taxids = set()
    for d in list(uni.values()) + list(uni_assumed.values()):
        taxids.add(d["taxon"])
    for r in shared:
        taxids.add(int(r["taxon_id"]))
        if r["strain_taxid"]:
            taxids.add(int(r["strain_taxid"]))
    sp_of, byname, phylum_of = load_taxonomy(taxids | {5271}, {"Sclerotium rolfsii"})
    sp_of_extra = {}
    if "Sclerotium rolfsii" in byname:
        t = byname["Sclerotium rolfsii"]
        sp_of_extra[t] = (t, "Sclerotium rolfsii")

    def sp(taxid):
        if taxid in sp_of:
            return sp_of[taxid]
        return sp_of_extra[taxid]

    # ---------------------------------------------------------------- positives from UniLectin3D
    for acc, r in ul_rows.items():
        cls = r["unilectin_class"]
        if acc in ("Q12396", "P43555"):
            continue  # handled as er_housekeeping negatives below
        if cls.startswith("AAL-like") and r["species"].startswith("Ralstonia"):
            reject(
                acc,
                "Bacterial protein (Ralstonia solanacearum) labelled fungal in UniLectin3D (spec F9). Not a fungal lectin.",
                "unilectin3d",
            )
            continue
        if acc.startswith("UPI"):
            seq = ul_seq.get(acc, "")
            if not seq:
                reject(
                    acc,
                    "UniParc ID without a sequence in the committed sequence file.",
                    "unilectin3d",
                )
                continue
            taxon = byname["Sclerotium rolfsii"]
            species, st = "Sclerotium rolfsii", taxon
            length, name, gene, status = len(seq), r["protein_name"], "", "UniParc only"
            uni_d = None
        else:
            uni_d = uni[acc]
            seq = uni_d["seq"]
            st, species = sp(uni_d["taxon"])
            length, name, gene = uni_d["length"], uni_d["name"], uni_d["gene"]
            status = "reviewed" if uni_d["reviewed"] else "unreviewed"
        sc = sugar_comps(acc)
        c = curated.get(acc)
        if cls.startswith("PA14"):
            stratum = "pa14_owned_elsewhere"
        elif c:
            stratum = c["stratum"]
        else:
            stratum = "structure_backed"
        ev_level = (
            c["evidence_level"]
            if c
            else (
                "unilectin3d_structure_ligand_not_retrieved"
                if not pdbs[acc]
                else (
                    "unilectin3d_structure_sugar_comp_in_pdb"
                    if sc
                    else "unilectin3d_structure_no_sugar_comp_in_pdb"
                )
            )
        )
        tier = c["tier"] if c else "T2"
        pm = sorted(ul_pmids.get(acc, []))
        src = "UniLectin3D:" + ",".join(sorted(pdbs[acc])) + "; UniProtKB:" + acc
        no_pdb_note = ""
        if not pdbs[acc]:
            src = (
                f"UniLectin3D ({r['n_pdb_structures']} structures in the committed table analysis/lectin_class_spec/unilectin_fungal_sequences.tsv; PDB ids not returned by the re-fetch of this work); UniProtKB:"
                + acc
            )
            no_pdb_note = "PDB ids and ligands not retrieved: the API re-fetch of this work did not return this accession (cap of 100 documents per call). "
        if pm:
            src += (
                "; "
                + "; ".join("PMID:" + x for x in pm)
                + " (from the UniLectin3D record, abstract not read)"
            )
        if c:
            src = c["source"] + "; UniLectin3D:" + ",".join(sorted(pdbs[acc]))
        quote = (
            f"UniLectin3D class {cls}, family {r['unilectin_family']}; {len(pdbs[acc]) or r['n_pdb_structures']} PDB entries; "
            f"sugar components in PDB non-polymer entities: {(','.join(sc) if sc else 'none') if pdbs[acc] else 'not retrieved'}"
        )
        notes = no_pdb_note
        if c:
            quote = c["quote"] + " | " + quote
            notes = (notes + c["notes"]).strip()
        if (
            uni_d
            and uni_d["reviewed"]
            and "ECO:0000269" in uni_d["codes"]
            and BIND.search(uni_d["text"])
            and SUGARW.search(uni_d["text"])
        ):
            notes = (
                notes + " " if notes else ""
            ) + "Swiss-Prot FUNCTION has ECO:0000269 binding text."
        if sc:
            notes = (
                (notes + " " if notes else "")
                + "A sugar component in a PDB entry can be a glycosylation site and not a ligand; ligand binding was not checked per structure."
            )
        rows.append(
            {
                "accession": acc,
                "gene": gene,
                "species": species,
                "taxon_id": st,
                "label": "positive",
                "stratum": stratum,
                "evidence_level": ev_level,
                "source": src,
                "quote": quote,
                "notes": notes,
                "label_class": stratum,
                "tier": tier,
                "h_flag": c["h_flag"] if c else "H0",
                "protein_name": name,
                "length": length,
                "strain_taxid": uni_d["taxon"] if uni_d else st,
                "unilectin_class": cls,
                "pdb_ids": ",".join(sorted(pdbs[acc])),
                "pdb_sugar_comps": ",".join(sc),
                "uniprot_status": status,
                "negative_evidence": "",
            }
        )
        seqs[acc] = seq

    # ---------------------------------------------------------------- positives not in UniLectin3D
    for acc, c in curated.items():
        if acc in ul_rows:
            continue
        d = uni[acc]
        st, species = sp(d["taxon"])
        rows.append(
            {
                "accession": acc,
                "gene": d["gene"],
                "species": species,
                "taxon_id": st,
                "label": "positive",
                "stratum": c["stratum"],
                "evidence_level": c["evidence_level"],
                "source": c["source"],
                "quote": c["quote"],
                "notes": c["notes"],
                "label_class": c["stratum"],
                "tier": c["tier"],
                "h_flag": c["h_flag"],
                "protein_name": d["name"],
                "length": d["length"],
                "strain_taxid": d["taxon"],
                "unilectin_class": "",
                "pdb_ids": "",
                "pdb_sugar_comps": "",
                "uniprot_status": "reviewed" if d["reviewed"] else "unreviewed",
                "negative_evidence": "",
            }
        )
        seqs[acc] = d["seq"]

    # ---------------------------------------------------------------- negatives: er_housekeeping
    pos_acc = {r["accession"] for r in rows}
    for acc, d in uni.items():
        if acc in pos_acc:
            continue
        in_m3 = bool(set(d["pfam"]) & M3_PF)
        if acc in ("Q12396", "P43555"):
            in_m3 = True
        if in_m3:
            if not d["reviewed"]:
                reject(
                    acc, "Not reviewed; no function annotation to give a reason.", "er_housekeeping"
                )
                continue
            if not d["text"]:
                reject(
                    acc,
                    "Swiss-Prot entry without FUNCTION text: no reason for a negative label (uncharacterised).",
                    "er_housekeeping",
                )
                continue
            ev = "experimental_eco269" if "ECO:0000269" in d["codes"] else "similarity_or_unstated"
            st, species = sp(d["taxon"])
            src = "UniProtKB:" + acc
            if d["pmids"]:
                src += "; " + "; ".join("PMID:" + x for x in sorted(d["pmids"]))
            if acc in ("Q12396", "P43555"):
                src += "; PMID:12134087"
            rows.append(
                {
                    "accession": acc,
                    "gene": d["gene"],
                    "species": species,
                    "taxon_id": st,
                    "label": "negative",
                    "stratum": "er_housekeeping",
                    "evidence_level": ev,
                    "source": src,
                    "quote": "UniProt FUNCTION: " + d["text"][:300],
                    "notes": "Reason: ER or Golgi house-keeping lectin-fold protein (spec 3.1 row M3; Lebreton 2021, PMID 34200153). The fold binds glycans; the label means 'not a host-interaction lectin', not 'no glycan binding'."
                    + (
                        " In UniLectin3D as 'Yeast Emp L-type' (structure)."
                        if acc in ul_rows
                        else ""
                    ),
                    "label_class": "er_housekeeping",
                    "tier": "",
                    "h_flag": "",
                    "protein_name": d["name"],
                    "length": d["length"],
                    "strain_taxid": d["taxon"],
                    "unilectin_class": ul_rows.get(acc, {}).get("unilectin_class", ""),
                    "pdb_ids": ",".join(sorted(pdbs.get(acc, []))),
                    "pdb_sugar_comps": "",
                    "uniprot_status": "reviewed",
                    "negative_evidence": ev,
                }
            )
            seqs[acc] = d["seq"]
    # ---------------------------------------------------------------- negatives: enzymes with a lectin-like domain
    for acc, d in uni.items():
        if acc in pos_acc or acc in seqs:
            continue
        if d["reviewed"] and re.search(r"alpha-galactosidase", d["name"], re.I):
            st, species = sp(d["taxon"])
            ev = "experimental_eco269" if "ECO:0000269" in d["codes"] else "similarity_or_unstated"
            rows.append(
                {
                    "accession": acc,
                    "gene": d["gene"],
                    "species": species,
                    "taxon_id": st,
                    "label": "negative",
                    "stratum": "enzyme_lectin_like_domain",
                    "evidence_level": ev,
                    "source": "UniProtKB:" + acc,
                    "quote": "UniProt FUNCTION: " + d["text"][:200],
                    "notes": "Reason: glycoside hydrolase (EC 3.2.1.22) with a ricin-B-like (PF00652) module. UniProt carries the keyword Lectin for it. The reason is protein-level function. The ricin-B-like module may bind sugar. Report apart from the other negatives.",
                    "label_class": "enzyme_lectin_like_domain",
                    "tier": "",
                    "h_flag": "",
                    "protein_name": d["name"],
                    "length": d["length"],
                    "strain_taxid": d["taxon"],
                    "unilectin_class": "",
                    "pdb_ids": "",
                    "pdb_sugar_comps": "",
                    "uniprot_status": "reviewed",
                    "negative_evidence": ev,
                }
            )
            seqs[acc] = d["seq"]
    # PA14 glucosidases from the shared set (N2 rows with an enzyme function)
    for r in shared:
        if r["accession"] in ("Q4WL79", "Q4WU49", "Q4WLY1", "Q4WA69"):
            st, species = sp(int(r["taxon_id"]))
            rows.append(
                {
                    "accession": r["accession"],
                    "gene": r["gene"],
                    "species": species,
                    "taxon_id": st,
                    "label": "negative",
                    "stratum": "enzyme_lectin_like_domain",
                    "evidence_level": "shared_set_N2_" + r["evidence_level"],
                    "source": "data/controls/hard-negatives-shared/controls.tsv:"
                    + r["accession"]
                    + "; "
                    + r["source"],
                    "quote": r["quote"][:300],
                    "notes": "Reason: beta-glucosidase (glycoside hydrolase) with a PA14 module (PF07691). PA14 is not one of the 15 models. Row is N2 in the shared set (adhesin-like by domain only). "
                    + r["notes"][:200],
                    "label_class": "enzyme_lectin_like_domain",
                    "tier": "",
                    "h_flag": "",
                    "protein_name": r["protein_name"],
                    "length": len(shared_seq[r["accession"]]),
                    "strain_taxid": r["taxon_id"],
                    "unilectin_class": "",
                    "pdb_ids": "",
                    "pdb_sugar_comps": "",
                    "uniprot_status": r["uniprot_status"][:40],
                    "negative_evidence": "shared_set_N2",
                }
            )
            seqs[r["accession"]] = shared_seq[r["accession"]]
    # ---------------------------------------------------------------- negatives: shared N1
    used = set(seqs)
    for r in shared:
        if r["evidence_level"] != "N1" or r["accession"] in used:
            continue
        st, species = sp(int(r["taxon_id"]))
        rows.append(
            {
                "accession": r["accession"],
                "gene": r["gene"],
                "species": species,
                "taxon_id": st,
                "label": "negative",
                "stratum": "n1_shared:" + r["stratum"],
                "evidence_level": "shared_set_N1",
                "source": "data/controls/hard-negatives-shared/controls.tsv:"
                + r["accession"]
                + "; "
                + r["source"],
                "quote": r["quote"][:300],
                "notes": "Reason: characterised non-adhesive function with a source (shared set N1). "
                + r["notes"][:200],
                "label_class": "n1_shared",
                "tier": "",
                "h_flag": "",
                "protein_name": r["protein_name"],
                "length": len(shared_seq[r["accession"]]),
                "strain_taxid": r["taxon_id"],
                "unilectin_class": "",
                "pdb_ids": "",
                "pdb_sugar_comps": "",
                "uniprot_status": r["uniprot_status"][:40],
                "negative_evidence": "shared_set_N1",
            }
        )
        seqs[r["accession"]] = shared_seq[r["accession"]]
    # ---------------------------------------------------------------- assumed
    taken = set(seqs)
    for acc, d in uni_assumed.items():
        if acc in taken:
            continue
        st, species = sp(d["taxon"])
        rows.append(
            {
                "accession": acc,
                "gene": d["gene"],
                "species": species,
                "taxon_id": st,
                "label": "negative",
                "stratum": "assumed",
                "evidence_level": "assumed_no_reason",
                "source": "UniProtKB:" + acc + " (release 2026_03); select_assumed.py",
                "quote": "No lectin keyword, no lectin-like Pfam cross-reference in UniProt, no lectin word in name or function text.",
                "notes": "ASSUMED. Not a reasoned negative. Absence of an annotation is not a reason (COMMON-RULES).",
                "label_class": "assumed",
                "tier": "",
                "h_flag": "",
                "protein_name": d["name"],
                "length": d["length"],
                "strain_taxid": d["taxon"],
                "unilectin_class": "",
                "pdb_ids": "",
                "pdb_sugar_comps": "",
                "uniprot_status": "reviewed",
                "negative_evidence": "assumed",
            }
        )
        seqs[acc] = d["seq"]

    # ---------------------------------------------------------------- candidates not used
    manual = {
        "Q8NJT4": "Same protein as Q4WW81 (FleA) under another accession (taxon 746128). Duplicate; identity checked in finalize step.",
        "Q8TGE0": "Same protein as Q2UNX8 (AOL) under another accession (taxon 5062). Duplicate.",
        "A0A1X7S147": "Mg1LysM of another Zymoseptoria tritici strain (ST99CH). Same gene as F9XHX3. Duplicate.",
        "A0A1X7S1V5": "MgxLysM of another strain; the abstract of PMID 21467214 does not show binding for MgxLysM (only Mg1 and Mg3). Not used.",
        "P9WEM2": "MgxLysM; the abstract of PMID 21467214 does not show binding for MgxLysM (only Mg1 and Mg3). UniProt cites PMID 33797163, which I did not read. Not used.",
        "A0A1X7S690": "Mg3LysM of another strain (ST99CH). Same gene as F9XMT4. Duplicate.",
        "G1X5Y2": "Arthrobotrys oligospora AOL_s00054g276 'fucose-specific lectin'. The abstract of PMID 35142828 reports a deletion phenotype (trap formation). It does not report a binding experiment. UniProt cites the paper for fucose binding. Tier T3 (name or similarity), not used as truth.",
        "D4AL26": "Arthroderma benhamiae 'Fucose-specific lectin': UniProt function 'Probable L-fucose-binding lectin', evidence ECO:0000250 (by similarity to P18891). Tier T3. Not used as truth. It is an Onygenales PF07938-family protein; see report.",
        "Q4WER9": "LysM effector ldpA of A. fumigatus: UniProt FUNCTION 'might have a role in sequestration of chitin' (ECO:0000305). Spec 4.3 lists it as a look-alike candidate. Probably a LysM effector, so it is not a negative. Not used.",
        "Q4WRR0": "LysM effector ldpB of A. fumigatus: same as Q4WER9. Not used.",
        "I1C083": "Mucoricin of Rhizopus delemar: ricin-like toxin and N-glycosylase (PMID 33462434). A ricin-like toxin may bind host glycans through its ricin-B-like module. Ambiguous. Not used.",
        "P83973": "Rhizopus stolonifer fucose-specific lectin: UniProt sequence is a 28 aa fragment. Fragment; no gene model.",
        "P86922": "Ganoderma lucidum lectin: UniProt sequence is an 8 aa fragment. Fragment.",
        "P84526": "Xerocomus spadiceus inulin-specific lectin: UniProt sequence is a 13 aa fragment. Fragment.",
        "Q9P403": "Colletotrichum lindemuthianum intracellular hyphae protein 1: function 'may have roles' (ECO:0000303). No binding evidence.",
    }
    for acc, d in uni.items():
        if acc in seqs or acc in {x["accession"] for x in rejected}:
            continue
        if acc in manual:
            reject(acc, manual[acc], "candidate")
        elif d["length"] < 50:
            reject(acc, "Sequence fragment under 50 aa.", "candidate")
        elif "PF01476" in d["pfam"]:
            reject(
                acc,
                "LysM protein without an experimental binding statement in UniProt FUNCTION (prediction, 'might have a role', or no text). Not a positive. Not a negative: LysM proteins bind chitin.",
                "candidate",
            )
        elif d["reviewed"] and not d["text"]:
            reject(
                acc,
                "Swiss-Prot entry without FUNCTION text. No reason for either label.",
                "candidate",
            )
        else:
            reject(
                acc,
                "No inclusion rule applies (not experimental sugar binding; not house-keeping; not an enzyme with a lectin-like domain).",
                "candidate",
            )

    # ---------------------------------------------------------------- finish
    seen = set()
    for r in rows:
        if r["accession"] in seen:
            raise SystemExit("duplicate accession " + r["accession"])
        seen.add(r["accession"])
        r.setdefault("tuned", "no")
        r["tuned"] = "no"
        r["leakage"] = "unknown"
        r["phylum"] = phylum_of(int(r["taxon_id"]))
        r["curator"] = CURATOR
        r["date"] = DATE
    rows.sort(key=lambda r: (r["label"], r["label_class"], r["accession"]))
    with open(a.out_dir + "/truth_all.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(a.out_dir + "/sequences.faa", "w") as fh:
        for r in rows:
            s = seqs[r["accession"]]
            fh.write(">" + r["accession"] + "\n")
            for i in range(0, len(s), 80):
                fh.write(s[i : i + 80] + "\n")
    with open(a.out_dir + "/rejected_build.tsv", "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["accession", "candidate_class", "reason"],
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rejected)
    from collections import Counter

    print(len(rows), "rows;", Counter((r["label"], r["label_class"]) for r in rows))


main()
