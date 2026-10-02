"""Score lookup for whole proteomes, agreement counts, named panel (Phase C spec 4).

score_source of a sequence (one value per hash):
- oof: the hash is in an S1 test fold (GO rows of the training sources and T-c rows); its
  score comes from the S1 model that did not train on it or its cluster.
- in_sample: the hash is in a FULL training table (V-go or V-kw) but in no S1 fold.
- final: the hash is in no training table; its score comes from the FULL model.
"""

PANEL = (
    # (name, lookup id, why) from the parent spec section 6; the lookup id is a UniProt
    # accession of the keyword set, a truth gene_id, or a proteome gene_id
    ("FLO1", "P32768", "1,537 aa GPI adhesin (C-terminal variant)"),
    ("SAG1", "P20840", "GPI wall protein"),
    ("CWP1", "P28319", "GPI wall protein"),
    ("CCW12", "Q12127", "GPI wall protein, 133 aa"),
    ("GAS1", "P22146", "GPI enzyme; ambiguous under D1"),
    ("PIR1", "Q03178", "wall protein without GPI; ambiguous under D1"),
    ("MSB2", "P32334", "PM-TM negative"),
    ("HKR1", "P41809", "PM-TM, 1,802 aa; surface evidence IBA only"),
    ("SUC2", "P00724", "secreted enzyme; ambiguous"),
    ("PHO5", "P00635", "secreted enzyme; P-ext wall"),
    ("ALS3", "Q59L12", "C. albicans GPI protein"),
    ("HWP1", "P46593", "C. albicans GPI protein"),
    ("SAP9", "Q59SU1", "C. albicans GPI protein"),
    ("ECM33 (C. albicans)", "A0A1D8PCY4", "P-ext wall; gpi_anchor=no in keywords"),
    ("ENO1 (C. albicans)", "CAL0000185645", "ambiguous stratum"),
    ("TDH3 (C. albicans)", "CAL0000197744", "ambiguous; internal term IBA only"),
    ("SOWgp58", "Q8NK60", "Onygenales Pro-rich surface protein"),
    ("SOWgp (RS)", "CIMG_04613-t26_1-p1", "C. immitis RS proteome gene CIMG_04613"),
    ("CTS1", "Q1E3R8", "Eurotiomycetes literature row"),
    ("CspA", "Q4WXC4", "literature row; P-ext wall in ASPFU"),
    ("cfmA", "Q4WLB9", "literature row; P-ext wall in ASPFU"),
    ("HSP60", "P50142", "moonlighting; never positive"),
)
SOURCE_ORDER = ("oof", "in_sample", "final")


def score_source(h: str, s1_hashes, full_train) -> str:
    """s1_hashes and full_train: containers of hashes (a dict of S1 folds works too)."""
    if h in s1_hashes:
        return "oof"
    if h in full_train:
        return "in_sample"
    return "final"


def find_panel_hashes(panel, members) -> dict[str, dict]:
    """lookup id -> {seq_sha256, set_ids}. A keyword-set accession wins over other sets."""
    found: dict[str, dict] = {}
    wanted = {p[1] for p in panel}
    for m in members:
        gid = m["gene_id"]
        if gid not in wanted:
            continue
        hit = found.setdefault(gid, {"seq_sha256": m["seq_sha256"], "set_ids": []})
        hit["set_ids"].append(m["set_id"])
        if m["set_id"] == "uniprot_kw":
            hit["seq_sha256"] = m["seq_sha256"]
    return found


def agreement_counts(rule_calls, ml_calls) -> dict:
    """Counts of (rule call, ML call) pairs."""
    out = {"rule1_ml1": 0, "rule1_ml0": 0, "rule0_ml1": 0, "rule0_ml0": 0}
    for r, m in zip(rule_calls, ml_calls, strict=True):
        out[f"rule{int(bool(r))}_ml{int(bool(m))}"] += 1
    return out
