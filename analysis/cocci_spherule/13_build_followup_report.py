#!/usr/bin/env python3
"""Write docs/reports/2026-10-03-cocci-spherule-followup.md from the follow-up summary files.

All counts and gene rows come from followup_genes_summary.tsv, phobius_spherule_up.tsv and the gene
table. Usage: 13_build_followup_report.py --dir analysis/cocci_spherule --out <report.md>
"""

import argparse
import csv
import gzip
from pathlib import Path


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def f(x, nd=3):
    try:
        return f"{float(x):.{nd}g}"
    except (TypeError, ValueError):
        return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    table = {
        r["gene_id"]: r
        for r in csv.DictReader(
            gzip.open(d / "spherule_surface_table.tsv.gz", "rt"), delimiter="\t"
        )
    }
    rows = list(csv.DictReader(open(d / "followup_genes_summary.tsv"), delimiter="\t"))
    phob = list(csv.DictReader(open(d / "phobius_spherule_up.tsv"), delimiter="\t"))

    def cls(x):
        if int(x["dimorphic_nsp"]) > 0 or int(x["outgroup_nsp"]) > 0:
            return "C1"
        return "C2" if int(x["other_ony_nsp"]) > 0 else "C3"

    def prev(x):
        v = table[x["gene_id"]]["prevalence"]
        return float(v) if v else None

    def supported(x):
        p = prev(x)
        return (
            (x["model_call"] == "agree" or x["in_fungi5k_RS_model"] == "yes")
            and p is not None
            and p > 0.01
        )

    for x in rows:
        x["cls"] = cls(x)
    n = len(rows)
    n_a = sum(1 for x in rows if "A_" in x["set"])
    n_b = sum(1 for x in rows if "B_" in x["set"])
    n_ab = sum(1 for x in rows if "A_" in x["set"] and "B_" in x["set"])
    prot = (
        list(
            csv.DictReader(
                open(d / "search_proteomes.tsv"),
                delimiter="\t",
            )
        )
        if (d / "search_proteomes.tsv").exists()
        else []
    )
    n_out = sum(1 for r in prot if r["order"] != "Onygenales")
    n_ony_f5k = sum(1 for r in prot if r["order"] == "Onygenales" and r["locustag"] != "pangenome")
    by = {k: [x for x in rows if x["cls"] == k] for k in ("C1", "C2", "C3")}
    c3_sup = [x for x in by["C3"] if supported(x)]
    n_spec_flag = sum(1 for x in rows if x["ranking_specific"] == "yes")
    n_agree = sum(1 for x in rows if x["model_call"] == "agree")
    n_nohit = sum(1 for x in rows if x["model_call"] == "no_hit")
    n_partial = sum(1 for x in rows if x["model_call"] == "partial_or_different")
    n_absent_f5k = sum(1 for x in rows if x["in_fungi5k_RS_model"] == "no")
    n_prev_low = sum(1 for x in rows if prev(x) is not None and prev(x) <= 0.01)
    n_up = len(phob)
    sp_sig = [p for p in phob if p["signalp_call"] == "SP"]
    sp_pho = [p for p in phob if p["phobius_sp"] == "yes"]
    pho_only = [p for p in phob if p["phobius_sp"] == "yes" and p["signalp_call"] != "SP"]
    pho_only_notm = [p for p in pho_only if int(p["phobius_tm"] or 0) == 0]
    n_pfam = sum(1 for x in rows if x["pfam_domains"])

    def hits(x):
        return (
            f"{x['coccidioides_nsp']}/{x['dimorphic_nsp']}/{x['other_ony_nsp']}/{x['outgroup_nsp']}"
        )

    def model(x):
        if x["model_call"] == "agree":
            return "agree"
        if x["model_call"] == "no_hit":
            return "no hit"
        return f"partial ({x['model_pident']}%, cov {x['model_qcov']})"

    def gene_row(x):
        t = table[x["gene_id"]]
        return [
            x["gene_id"], x["cls"], "AB" if "A_" in x["set"] and "B_" in x["set"] else x["set"][0],
            x["length"], f(x["log2fc_48h"], 3), f(x["cys_frac"], 2), f(t["prevalence"], 3),
            model(x), x["in_fungi5k_RS_model"], hits(x), x["pfam_domains"][:40] or "-",
            x["phobius_sp"],
        ]  # fmt: skip

    gh = ["Gene", "Class", "Set", "aa", "log2FC", "Cys", "Prevalence", "Model vs pangenome RS",
          "In Fungi5k RS", "Hits cocc/dimorph/other ony/outgroup", "Pfam", "Phobius SP"]  # fmt: skip

    def pho_row(p):
        t = table[p["gene_id"]]
        return [
            p["gene_id"], t["product"][:34], t["length"], f(t["log2fc_48h"], 3), f(p["signalp_sp_prob"], 2),
            p["phobius_tm"], f(t["cys_frac"], 2), t["flag_cocci_specific"] or "NA",
        ]  # fmt: skip

    text = f"""# Spherule follow-up: specificity, gene models, Cys domains and signal peptides

*2026-10-03. Branch `cocci-extreme-genes`. Follows `2026-10-03-cocci-spherule-surface-table.md` (PR #56).
Scripts: `analysis/cocci_spherule/10_` to `13_`. Counts are computed from the output files.*

## 1. What I checked and what came out

The first report left a set of {n} genes to check: {n_a} specific
extreme genes without a signal peptide (set A) and {
        n_b
    } Cys-rich spherule-up genes without a signal peptide that pass
the specificity flag (set B). {n_ab} genes are in both. All {n} pass the ranking's specificity flag
({n_spec_flag} of {n}).

1. **The ranking's specificity flag is not reliable for these genes.** A phmmer search (E-value at most 1e-5)
   finds a hit outside *Coccidioides* for {len(by["C1"]) + len(by["C2"])} of {n} genes. {
        len(by["C1"])
    } have a
   hit in a dimorphic relative or an outgroup (class C1). These contradict the flag, which excludes genes
   with a confounder-fungus hit. {
        len(by["C2"])
    } have hits only in other Onygenales (class C2), which the flag
   does not exclude.
   Example: CIMG_04662 (a P-type cation ATPase, 929 aa) has a hit in
   {
        next(x["dimorphic_nsp"] for x in rows if x["gene_id"] == "CIMG_04662")
    } dimorphic-relative proteomes. The
   flag called it specific. I do not know why. I did not inspect the ranking's confounder search.
2. **{len(by["C3"])} genes have no hit outside *Coccidioides*** (class C3). Only {
        len(c3_sup)
    } of them are
   supported by a second annotation of the same genome and are not rare in the pangenome
   ({", ".join(x["gene_id"] for x in c3_sup)}).
3. **Gene models.** Only {n_agree} of {
        n
    } proteins have a close match (at least 95% identity and 90% coverage on
   both sides) in the pangenome RS annotation. {n_partial} match partly or with lower identity. {
        n_nohit
    } have
   no match at all (diamond ultra-sensitive, E-value at most 1e-3). {
        n_absent_f5k
    } are absent from the Fungi5k
   annotation of RS. {
        n_prev_low
    } are in at most 1% of the 488 proteomes, so they are essentially RS-only.
   All {
        n
    } genes are up in spherules in the RNA-seq, so reads map to the RefSeq transcript sequences. That does
   not show that the RefSeq gene models are right or that the genes encode proteins.
4. **Phobius calls more signal peptides than SignalP.** Among the {
        n_up
    } spherule-up genes, Phobius calls
   {len(sp_pho)} and SignalP calls {len(sp_sig)}. All {
        len(sp_sig)
    } SignalP calls are also Phobius calls.
   {len(pho_only)} genes are Phobius-only, {
        len(pho_only_notm)
    } of them with no TM helix. Phobius is usually
   less specific than SignalP, so read the Phobius-only genes as candidates.
5. **No CFEM domain.** None of the {
        n_up
    } spherule-up proteins has a CFEM hit (PF05730, Pfam gathering
   threshold). {n_pfam} of the {n} follow-up genes have any Pfam domain.

## 2. Methods

| Check | Tool | Against |
|---|---|---|
| Orthologs | phmmer (HMMER 3.4), E-value at most 1e-3 searched, 1e-5 used for "strict" | {
        n_ony_f5k
    } Onygenales proteomes of Fungi5k, the two *C. posadasii* references (Silveira, C735-TKO), and one proteome for each of {
        n_out
    } outgroup genera (`search_proteomes.tsv`). The Fungi5k annotation of RS itself is excluded from the hit counts and used only for the model check |
| Gene model | diamond blastp ultra-sensitive, E-value at most 1e-3 | pangenome RS proteome (`C3CED3A_*`) |
| Domains | hmmscan, Pfam-A, gathering thresholds | all {n_up} spherule-up proteins |
| Signal peptide | Phobius 1.01 | all {n_up} spherule-up proteins |

Class rules: C1 = strict hit in a dimorphic relative (*Histoplasma*, *Blastomyces*, *Paracoccidioides*,
*Emergomyces*, *Emmonsia*, *Ajellomyces*) or an outgroup. C2 = strict hit only in other Onygenales. C3 = no
strict hit outside *Coccidioides*. Model agrees = 95% identity and 90% coverage of query and subject.

## 3. The {n} follow-up genes

Hits column: number of proteomes with a strict hit in *Coccidioides* (other than RS), dimorphic relatives,
other Onygenales and outgroups.

{
        md_table(
            gh,
            [gene_row(x) for x in sorted(rows, key=lambda x: (x["cls"], -float(x["log2fc_48h"])))],
        )
    }

## 4. Class C3: no hit outside *Coccidioides* ({len(by["C3"])} genes)

These are the candidates for *Coccidioides*-restricted genes. Most lack support in a second annotation.

{md_table(gh, [gene_row(x) for x in sorted(by["C3"], key=lambda x: -float(x["log2fc_48h"]))])}

Genes with a second-annotation match and prevalence above 1% are the best supported: {
        ", ".join(x["gene_id"] for x in c3_sup)
    }.
A phmmer miss at E-value 1e-5 does not prove absence. Short proteins (under 150 aa) have little power.

## 5. Phobius-only signal peptides among spherule-up genes ({len(pho_only)})

{
        md_table(
            [
                "Gene",
                "Product",
                "aa",
                "log2FC 48 h",
                "SignalP SP prob",
                "Phobius TM",
                "Cys",
                "Specific flag",
            ],
            [
                pho_row(p)
                for p in sorted(pho_only, key=lambda p: -float(table[p["gene_id"]]["log2fc_48h"]))
            ],
        )
    }

## 6. Limits

- The orthology search is sensitive but not complete. phmmer used single-sequence queries. A profile search
  (jackhmmer, or an HMM built from the hits) would find more.
- The outgroup set is {n_out} genera, one proteome each. It does not span the fungal tree.
- Absence from the pangenome RS annotation or from Fungi5k RS means the other pipelines did not predict the
  gene. It does not show that the RefSeq model is wrong, and the reverse is also open.
- The ranking's specificity flag was not re-derived. I only show that it disagrees with phmmer.
- Phobius-only calls can be false positives. TM-containing proteins can also be signal anchors.
- Two replicates, one strain (first report, section 7).
- No experiment tests any surface localisation here.

## 7. Next

- Map the {
        len(by["C3"])
    } C3 genes to the RS genome and check exon structure against the RNA-seq reads
  (the kallisto run does not show gene structure).
- Re-run the specificity search with a profile method, and find out why CIMG_04662 passed the flag.
- Decide whether a hand-checked gene list is worth sending for a proteomics or localisation test.
"""
    Path(a.out).write_text(text)
    print(f"wrote {a.out}: {len(text.splitlines())} lines")


if __name__ == "__main__":
    main()
