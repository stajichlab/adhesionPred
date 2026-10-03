import csv

import qc

SP = "Cryptococcus neoformans var. grubii H99"
DATE = "retrieved 2026-10-02"
FILES = {
    "22354955": "PMC3280450.txt",
    "30459196": "PMC6247093.txt",
    "27212659": "PMC5186401.txt",
    "17947228": "PMID17947228_abstract.txt",
    "16524904": "PMC1398056.txt",
    "25227465": "PMC4172073.txt",
    "17101662": "PMC1828480.txt",
    "26453029": "PMC4600298.txt",
    "27977806": "PMC5158083.txt",
    "23251520": "PMC3520850.txt",
    "33567338": "PMC7961099.txt",
    "41313167": "PMC12797982.txt",
    "31932719": "PMC7036007.txt",
    "37099613": "PMC10166503.txt",
    "40424315": "PMC12113280.txt",
    "42085564": "PMC13220255.txt",
    "11500433": "PMC98673.txt",
}
SRC = {
    "22354955": "PMC3280450 EuropePMC XML",
    "30459196": "PMC6247093 EuropePMC XML",
    "27212659": "PMC5186401 EuropePMC XML",
    "17947228": "abstract only (EuropePMC abstract; no PMC copy)",
    "16524904": "PMC1398056 PMC HTML page",
    "25227465": "PMC4172073 EuropePMC XML",
    "17101662": "PMC1828480 PMC HTML page",
    "26453029": "PMC4600298 EuropePMC XML",
    "27977806": "PMC5158083 EuropePMC XML",
    "23251520": "PMC3520850 EuropePMC XML",
    "33567338": "PMC7961099 EuropePMC XML",
    "41313167": "PMC12797982 EuropePMC XML",
    "31932719": "PMC7036007 EuropePMC XML",
    "37099613": "PMC10166503 EuropePMC XML",
    "40424315": "PMC12113280 EuropePMC XML",
    "42085564": "PMC13220255 NCBI efetch XML",
    "11500433": "PMC98673 PMC HTML page",
}
GENE = {
    "CDA2": "J9VND2",
    "CDA1": "J9VPD7",
    "QSP1": "J9VPP4",
    "PQP1": "J9VDC2",
    "PLB1": "Q9P8P2",
    "LAC1": "J9VY90",
    "APH1": "J9VHR6",
    "CNAG_05312": "J9VPN6",
    "cnap1": "J9VS02",
    "CRZ1": "J9VE33",
    "SOD2": "J9VWW9",
    "SOD1": "J9VLJ9",
    "CIG1": "J9VUB8",
    "BIM1": "J9VHN6",
    "CEL1": "J9VH79",
    "UGG1": "J9VG93",
    "MNS1": "J9VN97",
    "MNS101": "J9VSV0",
    "NOP1": "J9VR32",
    "VCX1": "J9VDQ4",
}
R = []  # (sym,go,code,label,sel,pmid,quote,strain_note,origin)


def add(sym, go, code, label, sel, pmid, quote, note, origin):
    R.append((sym, go, code, label, sel, pmid, quote, note, origin))


# ---- carried / corrected
Q_CDA2a = "Digestion with PI-PLC released Cda2, which ran at ~76 kDa by SDS-PAGE, from membranes into the supernatant"
# quotes carried from pass 1 are loaded from the pass-1 file by symbol/go/pmid to avoid retyping
old = list(csv.DictReader(open("../h99_pilot/draft_rows.tsv"), delimiter="\t"))


def oq(sym, go, pmid):
    for r in old:
        if r["symbol"] == sym and r["go_term"] == go and r["pmid"] == pmid:
            return r["evidence_note"].split(" | ")[0].replace("[abstract] ", "")
    raise KeyError((sym, go, pmid))


KN = "strain: KN99 (serotype A, H99-congenic); Results text: 'wild-type serotype A strain KN99'"
add(
    "CDA2",
    "GO:0005886",
    "IDA",
    "P-ext",
    "no",
    "22354955",
    oq("CDA2", "GO:0005886", "22354955"),
    KN + "; crude membranes, PI-PLC release and Triton X-114 partitioning, anti-Cda2 mAb western",
    "carried",
)
add(
    "CDA2",
    "GO:0005618",
    "IDA",
    "P-ext",
    "no",
    "22354955",
    oq("CDA2", "GO:0005618", "22354955"),
    KN + "; cell wall fraction released by beta-1,3-glucanase, cda2 deletion control",
    "carried",
)
add(
    "CDA2",
    "GO:0005576",
    "IDA",
    "P-ext",
    "no",
    "22354955",
    oq("CDA2", "GO:0005576", "22354955"),
    KN
    + "; culture secretions; Cda2 also seen in a cytosolic fraction, which the authors attribute to shearing from the wall, so no cytosol row",
    "carried",
)
add(
    "CDA1",
    "GO:0016020",
    "IDA",
    "P-ext",
    "no",
    "30459196",
    "Mutation of the catalytic residues of the Cda1 protein (Cda1CS) had no measurable effect on its localization in the membrane fraction compared to its wild-type Cda1 counterpart as can be seen on lanes marked M for membrane and CW for cell wall fractions in Fig. 5.",
    "strain: KN99 (H99-congenic); crude 'cell membrane fraction' western, evidence mainly for the catalytic mutant; text does not say whether Cda1 is in the CW lane; term GO:0016020 chosen because plasma membrane vs internal membrane is not separated. Gene-level label P-ext from H99 GOA ISS wall/extracellular terms per reviewer (GOA not re-opened by curator)",
    "corrected",
)
add(
    "CDA1",
    "GO:0016020",
    "IDA",
    "P-ext",
    "no",
    "40424315",
    "Cda1 was primarily detected in the insoluble cellular protein fraction, encompassing the cell wall and membrane compartments (Figure 7G, left), consistent with its cell surface localization via a GPI anchor.",
    "strain: paper builds mutants in H99 (Methods: 'introduced into the C. neoformans serotype A strain H99'); the WT strain used for this blot is not named as H99 in the sentence; anti-Cda1 western of total, soluble, insoluble fractions; wall and membrane not separated; secreted fraction also blotted (Fig 7G right) but WT secreted level not stated in text",
    "new",
)
add(
    "QSP1",
    "GO:0005576",
    "IDA",
    "P-ext",
    "no",
    "27212659",
    oq("QSP1", "GO:0005576", "27212659"),
    "strain: KN99alpha; Methods: 'The rest of the strains reported are derived from the KN99α wild-type since that isolate is more closely related to the H99 clinical isolate'; culture supernatant immunoblot and ELISA",
    "carried",
)
add(
    "PQP1",
    "GO:0009986",
    "IDA",
    "P-ext",
    "no",
    "27212659",
    oq("PQP1", "GO:0009986", "27212659"),
    "strain: KN99alpha (same sentence as QSP1); intact cells vs supernatant, fluorogenic substrate; GO:0009986 has no wall or extracellular ancestor; gene is P-ext through extracellular-region evidence (reviewer)",
    "carried",
)
add(
    "PLB1",
    "GO:0005618",
    "IDA",
    "P-ext",
    "no",
    "17947228",
    oq("PLB1", "GO:0005618", "17947228"),
    "strain: H99 (stated in the quote); H99 cell walls, beta-1,3-glucanase release; abstract only, no PMC full text exists",
    "carried",
)
add(
    "PLB1",
    "GO:0005886",
    "IDA",
    "P-ext",
    "unknown",
    "16524904",
    "Treatment of viable cryptococcal cells in suspension with M\u03b2CD also released PLB1 protein and enzyme activity, consistent with localization of PLB1 in plasma membrane rafts prior to secretion.",
    "strain: H99, now resolved from full text PMC1398056: 'A highly virulent, high PLB1-producing strain (H99) of C. neoformans var. grubii'; plasma membrane is the authors' inference from raft fractions plus MbetaCD release",
    "carried; strain resolved",
)
PL = "strain: H99, Methods: 'Wild-type C. neoformans var. grubii strain H99 (serotype A, MATa) and the SEC14-1 deletion mutant (28) were used in this study.'; Fig. 6A legend 'WT cells expressing Plb1-GFP'"
add(
    "PLB1",
    "GO:0005618",
    "IDA",
    "P-ext",
    "no",
    "25227465",
    oq("PLB1", "GO:0005618", "25227465"),
    PL + "; Plb1-GFP microscopy with calcofluor white",
    "carried; strain fixed",
)
add(
    "PLB1",
    "GO:0005886",
    "IDA",
    "P-ext",
    "no",
    "25227465",
    oq("PLB1", "GO:0005886", "25227465"),
    PL + "; same experiment, FM 4-64 membrane co-staining",
    "carried; strain fixed",
)
add(
    "LAC1",
    "GO:0005618",
    "IDA",
    "P-ext",
    "no",
    "17101662",
    oq("LAC1", "GO:0005618", "17101662"),
    "strain: H99, resolved from full text PMC1828480: 'Strain H99 lac1Δ ura5 (41) was employed as a recipient strain for expression studies.'; N-terminal GFP-laccase fusion; location depends on pH",
    "carried; strain resolved",
)
add(
    "LAC1",
    "GO:0005618",
    "IDA",
    "P-ext",
    "no",
    "11500433",
    "A monoclonal antibody to the C. neoformans laccase was generated and used to show localization in the cell walls of representative serotype A (H99) and serotype D (B-3501) strains by immunoelectron microscopy.",
    "strain: H99 stated in the quote (plus B-3501); text names the construct CNLAC1-GFP; the antibody sentence says 'the C. neoformans laccase', mapped to LAC1 by that CNLAC1 mention, not by an ID in the sentence",
    "new",
)
AP = "strain: H99, from the quote"
add(
    "APH1",
    "GO:0005576",
    "HDA",
    "P-ext",
    "yes",
    "25227465",
    "We performed a proteomic analysis on secretions obtained from the more-virulent encapsulated serotype A strain of C. neoformans, H99, and identified CNAG_02944 (designated Aph1) as the acid phosphatase responsible for the extracellular acid phosphatase activity previously observed by others",
    AP
    + "; H99 culture secretions by mass spectrometry; selected_by_predictor yes: paper says only Aph1 is predicted secreted among candidate acid phosphatases",
    "corrected",
)
add(
    "APH1",
    "GO:0005773",
    "IDA",
    "P-ext",
    "yes",
    "25227465",
    oq("APH1", "GO:0005773", "25227465"),
    "strain: H99 (Methods sentence as for PLB1 row, same paper); Aph1-GFP and Aph1-DsRed microscopy, vacuole and cell periphery; vacuole is secretory in labels.py, so gene label is P-ext with the extracellular row",
    "corrected",
)
add(
    "CNAG_05312",
    "GO:0005576",
    "HDA",
    "P-ext",
    "no",
    "26453029",
    oq("CNAG_05312", "GO:0005576", "26453029"),
    "strain: H99, Methods: 'The C. neoformans var. grubii wild-type strain H99 (WT) and the PGAL7::PKA1 strain with galactose-inducible/glucose repressible expression of PKA1 were used for this study'; culture supernatant LC-MS/MS",
    "carried; strain fixed",
)
add(
    "cnap1",
    "GO:0005576",
    "IDA",
    "P-ext",
    "yes",
    "27977806",
    oq("cnap1", "GO:0005576", "27977806"),
    "strain: H99-derived, Methods: 'strains used in this study were derived from strain H99'; ID check: paper lists CNAG_05872 / MAY1, UniProt J9VS02 ORF CNAG_05872 (reviewer); selected_by_predictor yes: deletion candidates were MS-identified peptidases with predicted signal sequences (reviewer)",
    "corrected",
)
add(
    "CRZ1",
    "GO:0005634",
    "IDA",
    "N-int",
    "no",
    "23251520",
    "Exposure to elevated temperature (30\u201337\u00b0C vs 25\u00b0C) and extracellular calcium caused calcineurin-dependent nuclear accumulation of Crz1-GFP.",
    "strain: H99, Methods: 'The strains Crz1-GFP, Δcrz1, Crz1-GFP:Pab-dsRed were created from H99'; location is condition dependent (cytosol and nuclei at 25 C; puncta under salt and heat shock)",
    "carried; strain fixed",
)
add(
    "SOD2",
    "GO:0005739",
    "IDA",
    "N-int",
    "no",
    "33567338",
    "As expected, Cu-sufficient cells expressed predominantly mitochondrially localized Sod2 (Fig. 6A).",
    "strain: H99 per full text: 'original C. neoformans H99 clinically isolated strain ... a lab evolved “wimp” strain emerged ... compared to the original “stud” strain' and 'the current studies presented here used the fully virulent “stud” strain'; H99 is inferred from these two sentences, no sentence says 'we used H99' for this assay; crude mitochondria vs cytosol fractionation (Fig. 6)",
    "corrected: quote changed from abstract to full-text sentence, code IDA",
)
add(
    "SOD2",
    "GO:0005829",
    "IDA",
    "N-int",
    "no",
    "33567338",
    "However, cells grown under Cu-limiting conditions demonstrated significant accumulation of Sod2 in the cytosol, with levels of mitochondrial Sod2 approximately equivalent to that found in cells grown under Cu-replete conditions.",
    "strain: as SOD2 mitochondrion row (H99 inferred via 'stud' sentences); cytosolic isoform under copper limitation",
    "corrected: quote changed from abstract to full-text sentence, code IDA",
)
# ---- new
add(
    "SOD1",
    "GO:0005829",
    "IDA",
    "N-int",
    "no",
    "33567338",
    "High levels of cytosolic Sod1 were detected during Cu sufficiency, with a small fraction associated with the mitochondria, presumably Sod1 that localizes to the mitochondrial IMS.",
    "strain: as SOD2 rows (H99 inferred via 'stud' sentences); fractionation into crude mitochondria and cytosol. Caveat: PMID 16524904 (abstract of raft paper, H99) reports SOD1 enriched in raft membrane fractions; if added as a membrane row the gene label would no longer be N-int",
    "new",
)
add(
    "CIG1",
    "GO:0005576",
    "HDA",
    "P-ext",
    "no",
    "26453029",
    "These proteins included three enzymes (α-amylase, acid phosphatase, and glyoxal oxidase), the Cig1 protein (cytokine-inducing glycoprotein) associated with virulence and heme uptake, and a novel protein containing a carbohydrate-binding domain (CNAG_05312).",
    "strain: H99 (same Methods sentence as CNAG_05312 row); secretome LC-MS/MS; ID: the paper's tables list 'CNAG_01653 Cytokine-inducing glycoprotein' and UniProt J9VUB8 has GN CNAG_01653; the sentence itself does not give the ID",
    "new",
)
add(
    "CIG1",
    "GO:0005886",
    "IDA",
    "P-ext",
    "no",
    "41313167",
    "We found that the mCherry signal was primarily localized to the plasma membrane, consistent with its function as a GPI-anchored mannoprotein (Fig. 3C).",
    "strain: H99, Methods: 'Cryptococcus neoformans clinical strain H99 and its mutants were cultured on the yeast extract-peptone-dextrose (YPD) medium'; Cig1:mCherry expressed from the ACT1 promoter in a strain also expressing Fbp1:FLAG (overexpression); ID given in the same paper: 'Cig1, CNAG_01653'",
    "new",
)
BM = "strain: wild-type control H99 (DTY758) named in Fig. 3b legend; the bim1Δ parent of the Bim1-HA strains is not named in the text read (Supplementary Dataset 1 not read), so H99 background for these strains is unverified; ID in paper: BIM1 (CNAG_02775, Uniprot ID J9VHN6)"
bq = "Bim1-HA fractionated with the plasma membrane, was released from the cell wall after treatment with fungal cell wall lytic enzymes and also partially fractionated in the culture supernatant (Figure 3a)."
add(
    "BIM1",
    "GO:0005886",
    "IDA",
    "P-ext",
    "no",
    "31932719",
    bq,
    BM + "; subcellular fractionation, Triton X-114",
    "new",
)
add(
    "BIM1",
    "GO:0005618",
    "IDA",
    "P-ext",
    "no",
    "31932719",
    bq,
    BM
    + "; release by lytic enzymes; indirect immunofluorescence also showed cell surface association",
    "new",
)
add(
    "BIM1",
    "GO:0005576",
    "IDA",
    "P-ext",
    "no",
    "31932719",
    bq,
    BM + "; only partial fractionation into supernatant; GPI-deleted form mostly in supernatant",
    "new",
)
add(
    "CEL1",
    "GO:0005618",
    "IDA",
    "P-ext",
    "no",
    "37099613",
    "In contrast, a dose dependent decrease of the Cel1-4xFLAG signal was observed when cells were treated with increasing concentration of the Zymolyase yeast lytic enzyme mix, indicating that Cel1-4xFLAG is likely associated with the Cn cell wall.",
    "strain: H99, Methods: 'All strains were generated in the C. neoformans var. grubii H99 background.'; ID in paper: CNAG_00601/CEL1; UniProt J9VH79 has GN=CEL1 (name match, ORF not shown in header); protein not found in supernatant (Fig 2F)",
    "new",
)
ER = "strain: H99; Methods of the same paper build deletions in 'the C. neoformans serotype A strain H99'; the GFP strains were made in 'the WT strain', not named H99 in that sentence; GFP fusion colocalised with ER-Tracker; "
eq = "Furthermore, green fluorescence protein (GFP)-tagged Ugg1, Mns1, and Mns101 proteins colocalized with the ER marker, supporting their role as confirming them to be functional ERQC components based on their subcellular localization in the ER (Figure 3—figure supplement 2)."
add(
    "UGG1",
    "GO:0005783",
    "IDA",
    "N-sec",
    "no",
    "40424315",
    eq,
    ER + "ID: CNAG_03648 (Methods)",
    "new",
)
add(
    "MNS1",
    "GO:0005783",
    "IDA",
    "N-sec",
    "no",
    "40424315",
    eq,
    ER + "ID: CNAG_02081 (Methods)",
    "new",
)
add(
    "MNS101",
    "GO:0005783",
    "IDA",
    "N-sec",
    "no",
    "40424315",
    eq,
    ER + "ID: CNAG_03240 (Methods)",
    "new",
)
MK = "strain: Methods: 'constructs were integrated into the MATα H99 strain' (mCherry markers; GFP markers in KN99a); marker expressed from the constitutive histone H3 promoter at the safe haven locus, so this is an overexpressed fusion; "
add(
    "NOP1",
    "GO:0005730",
    "IDA",
    "N-int",
    "no",
    "42085564",
    "Localization of the nucleolar marker Nop1 (Lee and Heitman 2012; Tollervey et al. 1991) was validated by DAPI staining, which showed colocalization of the fluorescent signal with the nucleolar region (Fig. 2c).",
    MK
    + "accession by protein name only (the single fibrillarin entry J9VR32, CNAG_06919); the paper gives no ID in the text read",
    "new",
)
add(
    "VCX1",
    "GO:0005773",
    "IDA",
    "N-sec",
    "no",
    "42085564",
    "We evaluated the vacuolar marker Vcx1 (Kmetzsch et al. 2010) under hypotonic and hypertonic conditions, as well as growth in RPMI medium. As expected, hypotonic growth produced enlarged vacuoles, whereas hypertonic conditions led to fragmented punctate structures (Fig. 2g).",
    MK
    + "weak: the sentence reports vacuole morphology, not an explicit Vcx1 localisation; accession J9VDQ4 has GN=VCX1 (a second vacuolar calcium transporter J9VIR4 exists)",
    "new",
)
# ---- write
hdr = "source_id species gene_id uniprot_accession symbol go_term evidence_code expected_label evidence_level selected_by_predictor pmid evidence_note reviewer review_date".split()
fail = 0
out = []
for sym, go, code, label, sel, pmid, quote, note, origin in R:
    ok = qc.check(quote, FILES[pmid])
    if not ok:
        fail += 1
        print("QUOTE NOT FOUND", sym, go, pmid, quote[:70])
    acc = GENE[sym]
    en = f"{quote} | strain/experiment: {note}; text checked: {SRC[pmid]}; {DATE}"
    out.append(
        ["Cneo_H99_GOA", SP, acc, acc, sym, go, code, label, "direct", sel, pmid, en, "", ""]
    )
    R[R.index((sym, go, code, label, sel, pmid, quote, note, origin))] = (
        sym,
        go,
        code,
        label,
        sel,
        pmid,
        quote,
        note,
        origin,
    )
with open("draft_rows.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(hdr)
    w.writerows(out)
with open("row_provenance.tsv", "w") as f:
    f.write("row\tsymbol\tgo_term\tpmid\torigin\ttext_file\tquote_machine_check\n")
    for i, (sym, go, _code, _label, _sel, pmid, quote, _note, origin) in enumerate(R, 2):
        f.write(
            f"{i}\t{sym}\t{go}\t{pmid}\t{origin}\t{FILES[pmid]}\t{'pass' if qc.check(quote,FILES[pmid]) else 'FAIL'}\n"
        )
print(len(R), "rows; quote failures", fail)
