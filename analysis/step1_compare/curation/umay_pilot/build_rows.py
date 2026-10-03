import csv
import re

from qc import norm


def seg(pmc, start, end):
    t = open(f"texts/{pmc}.txt").read()
    t = re.sub(r"\s+", " ", t)
    i = t.find(start)
    if i < 0:
        raise SystemExit(f"START not found {pmc}: {start}")
    j = t.find(end, i)
    if j < 0:
        raise SystemExit(f"END not found {pmc}: {end}")
    return t[i : j + len(end)]


SP = "Umay_MYCMD"
SPN = "Ustilago maydis 521"
R = []  # dict rows


def add(sym, acc, go, code, lab, sel, pmid, pmc, segs, note, idnote, full="full text"):
    q = " ... ".join(seg(pmc, s, e) for s, e in segs)
    R.append(
        dict(  # noqa: C408
            sym=sym,
            acc=acc,
            go=go,
            code=code,
            lab=lab,
            sel=sel,
            pmid=pmid,
            pmc=pmc,
            quote=q,
            note=note,
            idnote=idnote,
            full=full,
        )
    )


# ---- Rsp3 PMC5923269
RS = "Rsp3 (UMAG_03274)"
add(
    "Rsp3",
    "A0A0D1DYI3",
    "GO:0005618",
    "IDA",
    "P-ext",
    "unknown",
    "29703884",
    "PMC5923269",
    [
        (
            "To investigate whether Rsp3 when expressed from its native promoter attaches",
            "could not be detected (Fig. 4c, lower panel).",
        ),
        ("Secreted Rsp3-HA binds to the U. maydis cell wall.", "cell wall."),
    ],
    "strain: SG200 derivatives named in the quote (SG200Δrsp3-rsp3-HA, native promoter, non-permeabilised anti-HA immunostain of infected leaves); wall term rests on the paper's figure title; C-terminal HA on a protein the paper says has no GPI anchor",
    "accession by UMAG_03274 stated in paper; checked in UniProt REST and proteome FASTA",
)
add(
    "Rsp3",
    "A0A0D1DYI3",
    "GO:0005576",
    "IDA",
    "P-ext",
    "unknown",
    "29703884",
    "PMC5923269",
    [("To visualize secretion, we generated SG200Δrsp3 strains", "anti-HA antibody (Fig. 3b).")],
    "strain: SG200Δrsp3 derivatives named; Rsp3-HA from otef promoter (overexpressed); supernatant versus pellet western with tubulin as cytosolic lysis control; non-GPI protein",
    "accession by UMAG_03274",
)
# ---- Stp (PMC8159752)
STPQ = (
    "To test whether Stp complex members are secreted",
    "Stp4 proteins could be detected by western blot",
)
STPL = (
    "a-d, HA-specific western blot to reveal secretion of Stp proteins",
    "SG200Δkex2 (b) or SG200 (c, d).",
)
stp = [
    ("Stp1", "A0A0D1CU07", "UMAG_02475"),
    ("Stp2", "A0A0D1E167", "UMAG_10067"),
    ("Stp3", "A0A0D1EAN2", "UMAG_00715"),
    ("Stp4", "A0A0D1C5E3", "UMAG_12197"),
]
for s, a, u in stp:
    add(
        s,
        a,
        "GO:0005576",
        "IDA",
        "P-ext",
        "yes",
        "33941900",
        "PMC8159752",
        [STPQ, STPL],
        "strain: AB33, SG200Δkex2 or SG200 depending on panel (panel for this protein not checked; Stp1 needed the kex2 deletion strain); HA-tagged, otef promoter (overexpressed); selected as predicted effector genes; supernatant western with tubulin lysis control",
        f"accession by {u} from the paper's data availability statement",
    )
for s, a, u in stp[::2]:
    add(
        s,
        a,
        "GO:0005576",
        "IDA",
        "P-ext",
        "yes",
        "33941900",
        "PMC8159752",
        [
            ("To localize Stp1 and Stp3 at higher resolution", "for both Stp1-HA and Stp3-HA."),
            (
                "In an untagged SG200 strain, protrusions were seen, but no specific labelling was detected (Fig. 2d).",
                "(Fig. 2d).",
            ),
        ],
        "strain: SG200 derivatives (HA fusions from native promoters, untagged SG200 control named); immuno-EM, label outside the hyphal wall; the paper also describes the Stp complex as membrane anchored",
        "accession by " + u,
    )
for s, a, u in [("Stp5", "A0A0D1CK59", "UMAG_04342"), ("Stp6", "A0A0D1E805", "UMAG_01695")]:
    add(
        s,
        a,
        "GO:0005886",
        "IDA",
        "N-sec",
        "yes",
        "33941900",
        "PMC8159752",
        [
            (
                "Instead, HA-Stp5 and Stp6-HA localized to the fungal plasma membrane fraction",
                "(Extended Data Fig. 4e).",
            ),
            (
                "Total extract (T), plasma membrane fraction (M) and soluble fraction (Sol) from SG200-derived strains",
                "were loaded.",
            ),
        ],
        "strain: SG200-derived; HA fusion from the otef promoter (overexpressed); fractionation with mCherry cytosol, Cmu1 secreted and Pit1 plasma-membrane controls; the paper also finds Stp6-mCherry in speckles outside hyphae from the native promoter, which is not counted here",
        "accession by " + u,
    )
# ---- Pep1 PMC2631132
add(
    "Pep1",
    "G0X7E8",
    "GO:0005576",
    "IDA",
    "P-ext",
    "unknown",
    "19197359",
    "PMC2631132",
    [
        (
            "To overcome this problem, pep1 under control of its own promoter was fused",
            "strain SG200Δpep1.",
        ),
        ("SG200Δpep1-pep1M was used to infect maize lines", "(Figure 5A)."),
    ],
    "strain: SG200Δpep1 named; Pep1-mCherry from native promoter, fusion fully functional; signal around intracellular hyphae and in the apoplast, partial co-localisation with the maize plasma membrane marker PIN1-YFP",
    "accession by UMAG_01987 from the paper's data availability statement",
)
# ---- Pit2 and Erc1 (PMC9556619)
ERQ = (
    "Confocal microscopy confirmed that both Erc1-mCherry and Pit2-mCherry",
    "inside the fungal cell (Fig. 3a).",
)
ERS = (
    "The SG200 strain expressing Pit2-mCherry, an effector that was previously shown",
    "positive control for secretion.",
)
add(
    "Pit2",
    "A0A0D1EAR7",
    "GO:0005576",
    "IDA",
    "P-ext",
    "no",
    "36224193",
    "PMC9556619",
    [ERS, ERQ],
    "strain: SG200 named; Pit2-mCherry used as the positive control for secretion in this paper; cytosolic mCherry as negative control",
    "mapped by exact gene name PIT2 (paper text does not give the UMAG ID); proteome header GN=PIT2 UMAG_01375",
)
add(
    "Pit2",
    "A0A0D1EAR7",
    "GO:0005576",
    "IDA",
    "P-ext",
    "no",
    "33941900",
    "PMC8159752",
    [
        (
            "For microscopic analysis, a strain expressing the apoplastic effector Pit2-mCherryHA",
            "(Fig. 2a)",
        ),
        (
            "a, SG200 derived strains expressing the indicated fusion proteins growing in the epidermal layer",
            "mCherry signal, red.",
        ),
    ],
    "strain: SG200 derived (figure legend); control strain, native promoter; the sentence reads 'uniformly distributed around biotrophic hyphae' which is taken as extracellular",
    "mapped by exact gene name PIT2; proteome header GN=PIT2 UMAG_01375",
)
add(
    "Erc1",
    "A0A0D1E3H5",
    "GO:0005576",
    "IDA",
    "P-ext",
    "yes",
    "36224193",
    "PMC9556619",
    [
        (
            "Erc1 with a fused C-terminal mCherry tag was expressed in the SG200Δerc1 mutant under the control of its native promoter.",
            "native promoter.",
        ),
        ERQ,
        ("After plasmolysis with 1 M sodium chloride solution", "(Fig. 3a)."),
    ],
    "strain: SG200Δerc1; native promoter; confocal with plasmolysis; Erc1 is not described as GPI anchored; chosen as predicted effector",
    "accession by UMAG_01829",
)
add(
    "Erc1",
    "A0A0D1E3H5",
    "GO:0005576",
    "IDA",
    "P-ext",
    "yes",
    "36224193",
    "PMC9556619",
    [
        (
            "TEM micrographs obtained from immunogold labeling with a monoclonal antibody for HA depicted secretion of Erc1-HA to the biotrophic interface.",
            "biotrophic interface.",
        ),
        (
            "we performed transmission electron microscopy (TEM) of immunogold labeled maize leaf sections infected by an Erc1-HA expressing SG200 strain",
            "SG200 strain (Fig. 3b).",
        ),
    ],
    "strain: SG200 (named for the Erc1-HA strain in the paper's text: 'Erc1-HA expressing SG200 strain'); immunogold TEM; the same paper reports no specific accumulation in fungal wall, so only extracellular region is claimed",
    "accession by UMAG_01829",
)
# ---- Cmu1 PMC5923269
add(
    "Cmu1",
    "A0A0D1DWQ2",
    "GO:0005576",
    "IDA",
    "P-ext",
    "no",
    "29703884",
    "PMC5923269",
    [
        (
            "To determine whether Rsp3 might also reside as non-attached form in the biotrophic interface",
            "(Supplementary Fig. 4d).",
        )
    ],
    "strain: SG200-cmu1-AvitagHA named; apoplastic fluid from infected maize, western; Cmu1 is used as the secreted control here",
    "mapped by exact gene name CMU1; proteome header GN=CMU1 UMAG_05731",
)
# ---- Sts2 PMC10593772
add(
    "Sts2",
    "A0A0D1BWA6",
    "GO:0005576",
    "IDA",
    "P-ext",
    "unknown",
    "37872143",
    "PMC10593772",
    [
        ("To test for secretion, we expressed Sts2-mCherry in the CR-Sts2 mutant", "(Fig. 2a, b)."),
        (
            "We generated an open reading frame shift knockout of UMAG_05318 in U. maydis strain SG200",
            "CRISPR-Cas9 mutagenesis.",
        ),
    ],
    "strain: SG200 (CR-Sts2 is the SG200 knockout); C-terminal mCherry from the pit2 promoter (strong); signal on the hyphal edge, Pit2-mCherry as the comparison; not a GPI protein",
    "accession by UMAG_05318",
)
# ---- Llp1 PMC11857070
add(
    "Llp1",
    "A0A0D1E6R1",
    "GO:0005618",
    "IDA",
    "P-ext",
    "unknown",
    "39997458",
    "PMC11857070",
    [
        (
            "To assess whether Llp1, expressed under control of its native promoter",
            "was confirmed in this experiment (Figure 3b).",
        )
    ],
    "strain: SG200 Llp1-HA named; HA tag at the native locus (C-terminal); GPI status not stated in the text read; partially macerated infected leaves, anti-HA immunostain",
    "accession by UMAG_00027",
)
# ---- Xyn PMC8706147
XQ = ("As expected, given the presence of signal peptides", "(Figure 1B,C).")
XS = ("(C) Secretion of Xyn2, Xyn11A, and Xyn3 (tagged with GFP)", "(SG200-GFP).")
add(
    "Xyn2",
    "Q4P902",
    "GO:0005576",
    "IDA",
    "P-ext",
    "no",
    "34947062",
    "PMC8706147",
    [XQ, XS],
    "strain: SG200 (named in the legend for the control; the tagged strains are not named in the quote, SG200 background inferred); GFP-tagged colony secretion assay with cytoplasmic GFP lysis control",
    "accession by UMAG_03411",
)
add(
    "Xyn11A",
    "Q4P0L3",
    "GO:0005576",
    "IDA",
    "P-ext",
    "no",
    "34947062",
    "PMC8706147",
    [XQ, XS],
    "strain: SG200 (as for Xyn2); GFP-tagged colony secretion assay",
    "accession by UMAG_06350",
)
add(
    "Xyn11A",
    "Q4P0L3",
    "GO:0005576",
    "IDA",
    "P-ext",
    "no",
    "34947062",
    "PMC8706147",
    [
        (
            "we infected plants with SG200 strains harboring Xyn1-, Xyn2-, or Xyn11A:mCherry-HA fusion proteins",
            "is the first time this enzyme has been identified in the apoplast of infected plants.",
        )
    ],
    "strain: SG200 named; mCherry-HA fusion from the pit2 promoter; apoplastic fluid western; a free mCherry band is also present, so the fusion is partly cleaved",
    "accession by UMAG_06350",
)
# ---- Pdi1 PMC6881057
add(
    "Pdi1",
    "A0A0D1E191",
    "GO:0005783",
    "IDA",
    "N-sec",
    "no",
    "31730668",
    "PMC6881057",
    [
        (
            "Microscopy co-localization analysis of Pdi1:GFP and the ER marker mRFP:HDEL confirmed",
            "(Fig 5B)",
        ),
        (
            "To corroborate ER Pdi1 localization, experimental SG200pdi1:gfp and ER localization control SG200cals:mrfp:HDEL cells",
            "cells",
        ),
    ],
    "strain: SG200pdi1:gfp named; C-terminal GFP on a protein with an HDEL motif (tag may mask the motif, but the signal still co-localised with the ER marker); Pdi1 chosen from a glycoproteomic screen, not by prediction",
    "accession by UMAG_10156",
)
# ---- Sirtuins PMC10126416
SQN = (
    "Consistent with their localization motifs, we observed that Sir2 and Hst4 displayed nuclear localization (Figure 1B).",
    "(Figure 1B).",
)
SQM = (
    "By contrast, Hst5 and 6 were localized in the mitochondria, as verified by Mito-tracker colocalization",
    "Supplementary Figure S2A).",
)
SQL = (
    "(B) Subcellular localization of the indicated U. maydis sirtuins tagged with eGFP in its endogenous loci.",
    "endogenous loci.",
)
for s, a, u, go, q in [
    ("Sir2", "A0A0D1CD62", "UMAG_00963", "GO:0005634", SQN),
    ("Hst4", "A0A0D1DSY4", "UMAG_05758", "GO:0005634", SQN),
    ("Hst5", "A0A0D1E6Y0", "UMAG_05239", "GO:0005739", SQM),
    ("Hst6", "A0A0D1CIL7", "UMAG_12006", "GO:0005739", SQM),
]:
    add(
        s,
        a,
        go,
        "IDA",
        "N-int",
        "no",
        "37113216",
        "PMC10126416",
        [q, SQL],
        "strain: NOT NAMED for the eGFP-tagged strains in the text read; the deletions were made in SG200 and the controls are called wild-type, so SG200 is probable but not stated (strain_unstated); eGFP at the endogenous locus; the sentence names two sirtuins and is shared with the other protein of the pair",
        "accession by " + u,
    )

# ---- Pleiades PMC8224859 (late addition, see candidates_addendum.tsv)
add(
    "Mer1",
    "A0A0D1CND6",
    "GO:0005576",
    "IDA",
    "P-ext",
    "yes",
    "34166468",
    "PMC8224859",
    [
        (
            "a. Secretion of Mer1 during maize infection. Left: Plants infected with the U. maydis strain SG200ΔC",
            "apoplastic space.",
        )
    ],
    "strain: SG200ΔC (SG200 derivative) named; C-terminal mCherry-3xHA behind the Mer1 signal peptide from the tay1 promoter, plasmolysis; control without signal peptide stays diffuse in hyphae; selected as SignalP-predicted secreted Pleiades",
    "accession by UMAG_03753 (Tay1 UMAG_03752 and Mer1 UMAG_03753 stated in the results)",
)
add(
    "Tay1",
    "A0A0D1C393",
    "GO:0005576",
    "IDA",
    "P-ext",
    "yes",
    "34166468",
    "PMC8224859",
    [
        (
            "c. Secretion of Tay1 during maize infection. Left: Plants infected with the U. maydis strain SG200ΔC",
            "apoplastic space.",
        )
    ],
    "strain: SG200ΔC (SG200 derivative) named; N-terminal mCherry behind the Tay1 signal peptide from the tay1 promoter, plasmolysis; selected as SignalP-predicted secreted Pleiades",
    "accession by UMAG_03752",
)
# ---- Cce1 = Stp4 (PMC6638113) and Cpl1 (PMC10257043), added after the candidate hash
add(
    "Stp4",
    "A0A0D1C5E3",
    "GO:0005576",
    "IDA",
    "P-ext",
    "yes",
    "29745456",
    "PMC6638113",
    [
        (
            "The fusion proteins containing a secretion signal were tested using the same SG200Δcce1 complementation strains",
            "complementation strains as in the previous virulence assay.",
        ),
        ("The fusion proteins Cce1", "indicating that these proteins are secreted."),
    ],
    "strain: SG200Δcce1 complementation strains named; cce1 is UMAG_12197, which the Stp paper (33941900) calls Stp4 (Cce1); mCherry-HA from the native promoter; signal at the hyphal surroundings in the biotrophic interface; second paper for this gene",
    "accession by UMAG_12197 (stated in this paper's results)",
)
add(
    "Cpl1",
    "A0A0D1E4Q7",
    "GO:0005618",
    "IDA",
    "P-ext",
    "unknown",
    "37171083",
    "PMC10257043",
    [
        ("To consolidate our findings", "subsequently lysed."),
        ("To analyse if Cpl1", "fungal cell wall."),
    ],
    "strain: SG200Δcpl1 complemented with cpl1-HA named; otef promoter (overexpressed), C-terminal HA, anti-HA immunostain of filaments, no signal around hyphae; the paper also finds only faint Cpl1-HA in the supernatant; Cpl1 not described as GPI",
    "accession by UMAG_01820 (stated in the paper)",
)
# write
hdr = "source_id species gene_id uniprot_accession symbol go_term evidence_code expected_label evidence_level selected_by_predictor pmid evidence_note reviewer review_date".split()
with open("draft_rows.tsv", "w", newline="") as f, open("row_provenance.tsv", "w", newline="") as g:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(hdr)
    p = csv.writer(g, delimiter="\t", lineterminator="\n")
    p.writerow("row symbol uniprot_accession pmid text_file text_type quote_found".split())
    for n, r in enumerate(R, 1):
        note = f"{r['quote']} | strain/experiment: {r['note']}; accession: {r['idnote']}; text checked: {r['pmc']} EuropePMC XML; retrieved 2026-10-02"
        w.writerow(
            [
                SP,
                SPN,
                r["acc"],
                r["acc"],
                r["sym"],
                r["go"],
                r["code"],
                r["lab"],
                "direct",
                r["sel"],
                r["pmid"],
                note,
                "",
                "",
            ]
        )
        ok = all(
            norm(s) in norm(open(f"texts/{r['pmc']}.txt").read()) for s in r["quote"].split(" ... ")
        )
        p.writerow(
            [n, r["sym"], r["acc"], r["pmid"], r["pmc"] + ".txt", r["full"], "yes" if ok else "NO"]
        )
print(len(R), "rows")
