#!/usr/bin/env python3
"""Build the per-protein table (ids and labels only, no sequences) of every protein usable for the AGA1-class study.

Sources (all committed): truth_v3 annotated tables (S288C, C. albicans), repeat-mechanism consensus + curation table,
hard-negatives-shared controls, repeat-adhesins-other-species tables, wall-fusion-kiss2019 controls.
Clades come from the UniProt lineage (cached outside the repo) with the rule of analysis/wall_fusion_truth/map_kiss.py,
or from the Kiss table's own clade columns.

CLASS POSITIVE (pre-declared, operational): curated class `adhesin` AND repeat-bearing, where repeat-bearing means
the curators' mechanism label is 2a (repeat/avidity) OR the truth_v3 repeat label is 1. `pos_basis` records which.
Tier A = independent repeat evidence (truth_v3 repeat label 1, or mechanism basis stated in a paper / database annotation);
tier B = mechanism 2a by family inference only. GPI/signal peptide are UniProt annotations (feature lines), reported but NOT part of the positive rule.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from common import (  # noqa: E402
    CTRL,
    OUT,
    ROOT,
    all_sequences,
    clades_from_lineage_text,
    family_from_lineage_text,
    read_tsv,
    truth_tables,
    uniprot_cache,
    write_tsv,
)

# Candidate narrow classes, by gene name first, then by the curated Pfam list (data/curated/adhesins/adhesins.tsv).
# First match wins. These are CANDIDATES for the report, not decisions; the order is an assumption stated in the README.
GENE_RULES = [
    ("AGA1", r"^AGA1$"),
    ("Als_like", r"^ALS\d"),
    ("Flo11_type", r"^(FLO11|MUC1)\b"),
    ("Flo1_type", r"^FLO(1|5|9|10)\b"),
    ("Hwp1_Eap1_like", r"^(HWP|EAP|RBT)\d"),
    ("Hyr_Iff_like", r"^(HYR|IFF|RBR)\d"),
    ("Epa_like", r"^EPA\d"),
    ("SOWgp", r"^SOWGP"),
    ("Pfl_like", r"^PFL\d"),
]
PFAM_RULES = [  # (class, Pfam accessions that place a protein there)
    ("Als_like", {"PF05792", "PF11766"}),
    ("Flo11_type", {"PF10182"}),
    ("GLEYA_type", {"PF10528"}),
    ("Flocculin_repeat_type", {"PF00624", "PF13928", "PF07691"}),
    ("Hyphal_reg_CWP_type", {"PF11765", "PF15789"}),
]
COLS = [
    "accession", "gene", "organism", "taxon_id", "clade_coarse", "clade_fine", "clade_family", "role", "pos_basis", "pos_tier", "curated_class",
    "evidence", "mechanism_label", "mechanism_basis", "truth_proteome", "truth_repeat_label", "hn_stratum",
    "narrow_class", "pfam_curated", "cluster_source", "cluster_id", "uniprot_gpi_feature", "uniprot_signal_feature",
    "uniprot_repeat_n", "uniprot_reviewed", "sequence_source", "sources",
]  # fmt: skip


# two hard negatives are not UniProt accessions (GenBank/RefSeq ids); clade set from the species named in the table
MANUAL_CLADE = {
    "34373": ("Pezizomycotina", "Leotiomycetes"),
    "5039": ("Pezizomycotina", "Eurotiomycetes"),
}


def narrow_of(gene, pfams):
    g = re.sub(r"\s*\(.*", "", gene or "").strip().upper()
    for name, rx in GENE_RULES:
        if re.search(rx, g):
            return name
    for name, accs in PFAM_RULES:
        if accs & set(pfams):
            return name
    return ""


def main():
    cache = uniprot_cache()
    _, seqsrc = all_sequences()
    pf = {
        r["accession"]: r["pfam"]
        for r in read_tsv(ROOT / "data/curated/adhesins/adhesins.tsv")
        if r["accession"]
    }
    for r in read_tsv(CTRL / "hard-negatives-shared/controls.tsv"):
        pf.setdefault(r["accession"], r["family_domain_pfam"])
    P = {}

    def rec(acc):
        return P.setdefault(acc, {"accession": acc, "sources": set()})

    for prot, rows in truth_tables().items():
        for r in rows:
            d = rec(r["accession"])
            d["sources"].add(f"truth_v3.{prot}")
            d.update(gene=r["gene"], truth_proteome=prot, truth_repeat_label=r["label"], curated_class=r["curated_class"],
                     evidence=r["evidence"], cluster_source="truth_v3", cluster_id=r["cluster"])  # fmt: skip
    cons = {
        r["accession"]: r
        for r in read_tsv(CTRL / "repeat-mechanism/adjudication/curation_table.consensus.tsv")
    }
    for r in read_tsv(CTRL / "repeat-mechanism/curation_table.tsv"):
        d = rec(r["accession"])
        d["sources"].add("repeat-mechanism")
        d.setdefault("gene", r["gene"])
        d["organism"], d["taxon_id"] = r["organism"], r["taxon_id"]
        d["curated_class"], d["evidence"] = r["cls"], r["evidence_level"]
        d["cluster_source"], d["cluster_id"] = "repeat-mechanism", r["cluster_id"]
        c = cons.get(r["accession"])
        if c:
            d["mechanism_label"], d["mechanism_basis"] = c["final_label"], c["final_tier"]
    for r in read_tsv(CTRL / "repeat-adhesins-other-species/curation_table.other_species.tsv"):
        d = rec(r["accession"])
        d["sources"].add("repeat-adhesins-other-species")
        d.update(gene=r["gene"], organism=r["organism"], taxon_id=r["taxon_id"], curated_class=r["cls"],
                 evidence=r["evidence_level"], mechanism_label=r["mechanism_label"], mechanism_basis=r["evidence_basis"])  # fmt: skip
    for r in read_tsv(CTRL / "repeat-adhesins-other-species/hard_negatives_other_species.tsv"):
        d = rec(r["accession"])
        d["sources"].add("repeat-adhesins-other-species")
        d.update(gene=r["gene"], organism=r["organism"], taxon_id=r["taxon_id"], curated_class="hard_negative",
                 evidence=r["evidence_level"], hn_stratum=r["stratum"])  # fmt: skip
    for r in read_tsv(CTRL / "hard-negatives-shared/controls.tsv"):
        d = rec(r["accession"])
        d["sources"].add("hard-negatives-shared")
        d.setdefault("gene", r["gene"])
        d["organism"], d["taxon_id"] = r["species"], r["taxon_id"]
        d["curated_class"] = "hard_negative"
        d["evidence"], d["hn_stratum"] = r["evidence_level"], r["stratum"]
        if not d.get("cluster_id"):
            d["cluster_source"], d["cluster_id"] = "hard-negatives-shared", r["cluster_id"]
    for r in read_tsv(CTRL / "wall-fusion-kiss2019/controls.tsv"):
        d = rec(r["accession"])
        d["sources"].add("wall-fusion-kiss2019")
        d.setdefault("gene", r["gene"])
        d.setdefault("organism", r["species"])
        d.setdefault("taxon_id", r["taxon_id"])
        d["kiss"] = (r["clade_coarse"], r["clade_fine"], r["stratum"], r["label"], r["cluster_id"])
        d.setdefault("cluster_source", "kiss")
        d.setdefault("cluster_id", r["cluster_id"])

    rows = []
    for acc, d in sorted(P.items()):
        u = cache.get(acc)
        if u:
            d["clade_coarse"], d["clade_fine"] = clades_from_lineage_text(u["Taxonomic lineage"])
            d["clade_family"] = family_from_lineage_text(u["Taxonomic lineage"])
            d["uniprot_gpi_feature"] = "yes" if "GPI-anchor" in u["Lipidation"] else "no"
            d["uniprot_signal_feature"] = "yes" if u["Signal peptide"].strip() else "no"
            d["uniprot_repeat_n"] = str(u["Repeat"].count("REPEAT "))
            d["uniprot_reviewed"] = u["Reviewed"]
            d.setdefault("organism", u["Organism"])
            d.setdefault("taxon_id", u["Organism (ID)"])
            d.setdefault("gene", u["Gene Names (primary)"])
        elif d.get("kiss"):
            d["clade_coarse"], d["clade_fine"] = d["kiss"][0], d["kiss"][1]
        elif d.get("taxon_id") in MANUAL_CLADE:
            d["clade_coarse"], d["clade_fine"] = MANUAL_CLADE[d["taxon_id"]]
        d["sequence_source"] = seqsrc.get(acc, "none")
        d["pfam_curated"] = pf.get(acc, "")
        d["narrow_class"] = narrow_of(
            d.get("gene", ""), [x for x in re.split(r"[;,]", d["pfam_curated"]) if x]
        )
        cc = d.get("curated_class", "")
        tl = d.get("truth_repeat_label", "")
        mech = d.get("mechanism_label", "")
        d["pos_basis"] = d["pos_tier"] = ""
        if cc == "adhesin" and (mech == "2a" or tl == "1"):
            d["role"] = "class_positive"
            d["pos_basis"] = "+".join(x for x in (("mechanism_2a" if mech == "2a" else ""), ("truth_repeat_1" if tl == "1" else "")) if x)  # fmt: skip
            strong = tl == "1" or d.get("mechanism_basis", "") in (
                "stated_in_paper",
                "database_annotation",
            )
            d["pos_tier"] = "A" if strong else "B"
        elif cc == "adhesin":
            d["role"] = "adhesin_other_mechanism"
        elif cc == "hard_negative":
            d["role"] = "hard_negative"
        elif cc in ("surface_other_adhesion_phenotype", "indirect_regulator"):
            d["role"] = "excluded_" + cc
        elif cc == "uniprot":
            d["role"] = "uncurated_repeat_bearing" if tl == "1" else "background_assumed_negative"
        elif d.get("kiss"):
            d["role"] = "kiss_row_" + d["kiss"][3]
        else:
            d["role"] = "unclassified"
        d["sources"] = ";".join(sorted(d["sources"]))
        rows.append(d)
    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT / "proteins.tsv", rows, COLS)
    from collections import Counter

    print(len(rows), Counter(r["role"] for r in rows))


if __name__ == "__main__":
    main()
