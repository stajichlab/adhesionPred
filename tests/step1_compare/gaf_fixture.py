"""Golden 50-row GAF fixture for the D1 extractor (spec section 7, test_gaf_extract_golden).

Each tuple is (db, gene_id, symbol, qualifier, term, evidence, aspect, object_type, taxon).
The comment after each gene block gives its expected label under the non-IEA policy.
"""

from pathlib import Path

TAXON = "taxon:559292"
P, F = "GO:0008150", "GO:0003674"  # biological_process, molecular_function roots
WALL, EXT, PM = "GO:0009277", "GO:0005576", "GO:0005886"
CYTO, NUC, MITO, MITO_IM = "GO:0005829", "GO:0005634", "GO:0005739", "GO:0005743"
ER, VAC, MEM, PERIPH, OBSOLETE = (
    "GO:0005783",
    "GO:0000324",
    "GO:0016020",
    "GO:0071944",
    "GO:0031225",
)

ROWS = [
    # 1-3 G01 wall IDA -> P-ext (wall); F and P aspect rows are dropped from labels
    ("SGD", "G01", "WALL1", "located_in", WALL, "IDA", "C", "protein", TAXON),
    ("SGD", "G01", "WALL1", "enables", F, "IDA", "F", "protein", TAXON),
    ("SGD", "G01", "WALL1", "involved_in", P, "IEA", "P", "protein", TAXON),
    # 4-5 G02 extracellular IDA (+ IEA copy) -> P-ext (extracellular-only)
    ("SGD", "G02", "SECR1", "located_in", EXT, "IDA", "C", "protein", TAXON),
    ("SGD", "G02", "SECR1", "located_in", EXT, "IEA", "C", "protein", TAXON),
    # 6-7 G03 wall + plasma membrane IDA -> P-ext (wall), pm_candidate
    ("SGD", "G03", "GPIPM", "located_in", WALL, "IDA", "C", "protein", TAXON),
    ("SGD", "G03", "GPIPM", "located_in", PM, "IDA", "C", "protein", TAXON),
    # 8-9 G04 extracellular IDA + cytosol IDA -> ambiguous
    ("SGD", "G04", "MOON1", "located_in", EXT, "IDA", "C", "protein", TAXON),
    ("SGD", "G04", "MOON1", "located_in", CYTO, "IDA", "C", "protein", TAXON),
    # 10-11 G05 cytosol IDA, nucleus IEA -> N-int
    ("SGD", "G05", "CYTO1", "located_in", CYTO, "IDA", "C", "protein", TAXON),
    ("SGD", "G05", "CYTO1", "located_in", NUC, "IEA", "C", "protein", TAXON),
    # 12-13 G06 nucleus IDA, membrane IEA -> unlabelled (membrane at any level blocks N-int)
    ("SGD", "G06", "NUCMEM", "located_in", NUC, "IDA", "C", "protein", TAXON),
    ("SGD", "G06", "NUCMEM", "located_in", MEM, "IEA", "C", "protein", TAXON),
    # 14 G07 ER IDA -> N-sec (ER is part_of endomembrane system)
    ("SGD", "G07", "ER1", "located_in", ER, "IDA", "C", "protein", TAXON),
    # 15-16 G08 plasma membrane IDA + nucleus IDA -> N-sec (internal terms allowed in N-sec)
    ("SGD", "G08", "PMNUC", "located_in", PM, "IDA", "C", "protein", TAXON),
    ("SGD", "G08", "PMNUC", "located_in", NUC, "IDA", "C", "protein", TAXON),
    # 17-19 G09 wall IEA + cytosol IEA + PM IDA -> unlabelled (IEA wall blocks N-sec)
    ("SGD", "G09", "IEAWALL", "located_in", WALL, "IEA", "C", "protein", TAXON),
    ("SGD", "G09", "IEAWALL", "located_in", CYTO, "IEA", "C", "protein", TAXON),
    ("SGD", "G09", "IEAWALL", "located_in", PM, "IDA", "C", "protein", TAXON),
    # 20-22 G10 wall IDA + wall IBA + cytosol IBA -> ambiguous; P-ext without homology codes
    ("SGD", "G10", "IBAWALL", "located_in", WALL, "IDA", "C", "protein", TAXON),
    ("SGD", "G10", "IBAWALL", "is_active_in", WALL, "IBA", "C", "protein", TAXON),
    ("SGD", "G10", "IBAWALL", "located_in", CYTO, "IBA", "C", "protein", TAXON),
    # 23 G11 extracellular IBA only -> P-ext (extracellular-only); unlabelled without homology
    ("SGD", "G11", "IBAONLY", "located_in", EXT, "IBA", "C", "protein", TAXON),
    # 24-25 G12 NOT wall IDA + vacuole IDA -> N-sec (NOT row skipped)
    ("SGD", "G12", "NOTWALL", "NOT|located_in", WALL, "IDA", "C", "protein", TAXON),
    ("SGD", "G12", "NOTWALL", "located_in", VAC, "IDA", "C", "protein", TAXON),
    # 26-27 G13 obsolete term IDA + ER IDA -> N-sec; obsolete row counted
    ("SGD", "G13", "OBSOL", "located_in", OBSOLETE, "IDA", "C", "protein", TAXON),
    ("SGD", "G13", "OBSOL", "located_in", ER, "IDA", "C", "protein", TAXON),
    # 28-29 G14 wall IDA + mitochondrion HDA -> ambiguous (as GAS1 in sgd.gaf.gz)
    ("SGD", "G14", "HDAMITO", "located_in", WALL, "IDA", "C", "protein", TAXON),
    ("SGD", "G14", "HDAMITO", "located_in", MITO, "HDA", "C", "protein", TAXON),
    # 30 G15 vacuole IEA only -> unlabelled; counted in genes_cc, not in genes_noniea_cc
    ("SGD", "G15", "VAC1", "located_in", VAC, "IEA", "C", "protein", TAXON),
    # 31-32 G16 wall ISS + ER IDA -> P-ext (wall); unlabelled without homology codes
    ("SGD", "G16", "ISSWALL", "located_in", WALL, "ISS", "C", "protein", TAXON),
    ("SGD", "G16", "ISSWALL", "located_in", ER, "IDA", "C", "protein", TAXON),
    # 33 G17 mitochondrial inner membrane IDA -> unlabelled (membrane ancestor blocks N-int)
    ("SGD", "G17", "MITOMEM", "located_in", MITO_IM, "IDA", "C", "protein", TAXON),
    # 34 G18 cell periphery IDA only -> unlabelled
    ("SGD", "G18", "PERIPH", "located_in", PERIPH, "IDA", "C", "protein", TAXON),
    # 35-36 complex rows from ComplexPortal -> dropped (other_db)
    ("ComplexPortal", "CPX-1", "cpx1", "part_of", CYTO, "IDA", "C", "protein_complex", TAXON),
    ("ComplexPortal", "CPX-1", "cpx1", "part_of", NUC, "IDA", "C", "protein_complex", TAXON),
    # 37-38 ncRNA rows from the primary database -> dropped (object_type)
    ("SGD", "R01", "SNR1", "located_in", NUC, "IDA", "C", "ncRNA", TAXON),
    ("SGD", "R01", "SNR1", "located_in", CYTO, "IDA", "C", "ncRNA", TAXON),
    # 39-41 duplicate-database rows (UniProtKB in a MOD file) -> dropped (other_db)
    ("UniProtKB", "P99999", "DUP1", "located_in", WALL, "IBA", "C", "protein", TAXON),
    ("UniProtKB", "P99999", "DUP1", "involved_in", P, "IBA", "P", "protein", TAXON),
    ("UniProtKB", "P99999", "DUP1", "enables", F, "IBA", "F", "protein", TAXON),
    # 42-43 G21 other strain -> P-ext without a taxon filter, dropped with taxon 559292
    ("SGD", "G21", "STRAIN2", "located_in", WALL, "IDA", "C", "protein", "taxon:999999"),
    ("SGD", "G21", "STRAIN2", "involved_in", P, "IDA", "P", "protein", "taxon:999999"),
    # 44 G22 NCBITaxon prefix for the same taxon -> kept, P-ext (extracellular-only)
    ("SGD", "G22", "NCBI1", "located_in", EXT, "IDA", "C", "gene_product", "NCBITaxon:559292"),
    # 45-50 more rows for genes above
    ("SGD", "G05", "CYTO1", "involved_in", P, "IDA", "P", "protein", TAXON),
    ("SGD", "G07", "ER1", "enables", F, "IEA", "F", "protein", TAXON),
    ("SGD", "G03", "GPIPM", "located_in", PM, "IEA", "C", "protein", TAXON),
    ("SGD", "G08", "PMNUC", "located_in", NUC, "HDA", "C", "protein", TAXON),
    ("SGD", "G11", "IBAONLY", "enables", F, "IBA", "F", "protein", TAXON),
    ("SGD", "G04", "MOON1", "located_in", EXT, "EXP", "C", "protein", TAXON),
]

EXPECTED_LABELS = {
    "G01": ("P-ext", "wall"),
    "G02": ("P-ext", "extracellular-only"),
    "G03": ("P-ext", "wall"),
    "G04": ("ambiguous", ""),
    "G05": ("N-int", ""),
    "G06": ("unlabelled", ""),
    "G07": ("N-sec", ""),
    "G08": ("N-sec", ""),
    "G09": ("unlabelled", ""),
    "G10": ("ambiguous", ""),
    "G11": ("P-ext", "extracellular-only"),
    "G12": ("N-sec", ""),
    "G13": ("N-sec", ""),
    "G14": ("ambiguous", ""),
    "G15": ("unlabelled", ""),
    "G16": ("P-ext", "wall"),
    "G17": ("unlabelled", ""),
    "G18": ("unlabelled", ""),
    "G21": ("P-ext", "wall"),
    "G22": ("P-ext", "extracellular-only"),
}


def write_golden_gaf(path: Path) -> Path:
    lines = ["!gaf-version: 2.2", "!date-generated: 2026-05-21T09:00", ""]
    for db, gene, symbol, qual, term, ev, aspect, otype, taxon in ROWS:
        cols = [
            db,
            gene,
            symbol,
            qual,
            term,
            "PMID:1",
            ev,
            "",
            aspect,
            symbol.lower(),
            f"{symbol}|{gene}-syn",
            otype,
            taxon,
            "20260101",
            "SGD",
            "",
            "",
        ]
        lines.append("\t".join(cols))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
