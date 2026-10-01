"""A small, realistic Phase A + Phase B work directory for the Phase C tests.

`make_work(root)` writes every file that 08_build_eval_tables.py reads, with the real column
sets, consistent run JSONs (the M3 hash chain holds) and small embedding matrices (8M: 6
columns, 35M: 8 columns). Sources: Scer_SGD and Calb_CGD (train), Spom_PomBase (test_species),
Spom_SCHPO-mod (alternate_file), Afum_ASPFU (test_clade, Eurotiomycetes), Cneo_H99_GOA
(test_clade, Basidiomycota), Umay_MYCMD (undecided, Basidiomycota). Features and embeddings
carry a class signal plus noise, so the models can learn.

Planted cases (the tests rely on them):
- `SHARED_NSEC`, `SHARED_UNRES`, `SHARED_AMBIG`, `SHARED_POS`: T-c rows with the hash of an
  N-sec, a pm-unresolved, an ambiguous and a P-ext gene of Scer_SGD (ruling C-6).
- `TC_DUP`: two T-c accessions with one sequence.
- the alternate source Spom_SCHPO-mod has one gene (`ALT_ONLY`) whose sequence is in no other
  source; it must not enter the table.
- `LONG_1022` and `LONG_1023`: Calb_CGD proteins of 1,022 and 1,023 aa.
- families: groups of genes that share the first 6 residues (the stub MMseqs2 clusters them).
- proteome sets Scer_proteome and Cimm_RS_proteome (members of truth and new sequences).
"""

import gzip
import hashlib
import json
import random
from pathlib import Path

import numpy as np
import seqhash
import truth_table

AA = "ACDEFGHIKLMNPQRSTVWY"
MODELS = {"esm2_t6_8M_UR50D": 6, "esm2_t12_35M_UR50D": 8}
TRUTH_SHA = "5" * 64
SPECIES_COLUMNS = (
    "source_id", "species", "taxon_id", "taxon_filter", "in_clade", "role", "role_note",
    "gaf_file", "fasta_file", "id_mapping",
)  # fmt: skip
SOURCES = [
    # source_id, species, taxon, clade, role, (pos, nint, nsec, pmtm, unres, ambig)
    ("Scer_SGD", "Saccharomyces cerevisiae S288C", "559292", "Saccharomycotina", "train",
     (16, 16, 12, 2, 1, 2)),
    ("Calb_CGD", "Candida albicans SC5314", "237561", "Saccharomycotina", "train",
     (16, 16, 12, 2, 1, 2)),
    ("Spom_PomBase", "Schizosaccharomyces pombe 972h-", "284812", "Taphrinomycotina",
     "test_species", (6, 8, 6, 0, 0, 1)),
    ("Spom_SCHPO-mod", "Schizosaccharomyces pombe 972h-", "284812", "Taphrinomycotina",
     "alternate_file", (0, 0, 0, 0, 0, 0)),
    ("Afum_ASPFU", "Aspergillus fumigatus Af293", "330879", "Eurotiomycetes", "test_clade",
     (6, 6, 6, 0, 1, 0)),
    ("Cneo_H99_GOA", "Cryptococcus neoformans H99", "235443", "Basidiomycota", "test_clade",
     (4, 4, 4, 0, 0, 0)),
    ("Umay_MYCMD", "Ustilago maydis 521", "5270", "Basidiomycota", "undecided",
     (4, 4, 4, 0, 0, 0)),
]  # fmt: skip
TC_TAXA = [("330879", "Aspergillus fumigatus Af293", 5), ("246410", "Coccidioides immitis RS", 3),
           ("237561", "Candida albicans SC5314", 3), ("559292", "Saccharomyces cerevisiae S288C", 3),
           ("498019", "Candidozyma auris B8441", 2)]  # fmt: skip
TC_CLADES = {"330879": "Eurotiomycetes", "246410": "Eurotiomycetes", "237561": "Saccharomycotina",
             "559292": "Saccharomycotina", "498019": "Saccharomycotina"}  # fmt: skip


def _seq(rng: random.Random, n: int, st: float, prefix: str = "") -> str:
    body = []
    for _ in range(n - len(prefix)):
        body.append(rng.choice("ST") if rng.random() < st else rng.choice(AA))
    return prefix + "".join(body)


def _truth_row(src, gene_id, label, subset, d8, homology_only, htp, symbol=""):
    source_id, species, taxon, clade, role, _ = src
    stratum = d8 or subset or label
    return {
        "source_id": source_id, "species": species, "taxon_id": taxon, "in_clade": clade,
        "role": role, "gene_id": gene_id, "symbol": symbol, "synonym1": "", "label": label,
        "subset": subset, "stratum": stratum, "tier": "T-a", "label_no_homology": label,
        "label_experimental": label, "homology_only": homology_only,
        "pm_candidate": "yes" if d8 else "no", "evidence_codes": "IDA,IEA",
        "surface_evidence": "IDA" if label in ("P-ext", "ambiguous") else "",
        "internal_evidence": "HDA" if htp else ("IDA" if label in ("N-int", "ambiguous") else ""),
        "internal_evidence_htp_only": "yes" if htp else "no",
        "secretory_evidence": "IDA" if label == "N-sec" else "", "source_file": "f.gaf",
        "source_sha256": "0" * 64, "source_date": "2026-05-21", "obo_sha256": "1" * 64,
        "d8_class": d8, "d8_reason": "",
    }  # fmt: skip


def _features(rng: random.Random, cls: str, seq: str) -> dict:
    if cls == "pos":
        sp = rng.random() < 0.85
        gpi = rng.choices(["highly_probable", "probable", "weakly", "none"], [4, 1, 1, 4])[0]
    elif cls == "nsec":
        sp = rng.random() < 0.5
        gpi = rng.choices(["highly_probable", "none"], [1, 9])[0]
    else:
        sp = rng.random() < 0.05
        gpi = "none"
    if len(seq) <= 40:
        gpi = "too_short"
    sp_prob = rng.uniform(0.7, 1.0) if sp else rng.uniform(0.0, 0.2)
    gpi_prob = {"highly_probable": 1.0, "probable": 0.7, "weakly": 0.55}.get(gpi, 0.0)
    st = (seq.count("S") + seq.count("T")) / len(seq)
    return {
        "ser_thr_frac": f"{st:.6f}", "sp_prediction": "SP" if sp else "OTHER",
        "sp_prob": f"{sp_prob:.4f}", "sp_other_prob": f"{1 - sp_prob:.4f}",
        "sp_cs_end": "22" if sp else "", "sp_cs_prob": "0.9" if sp else "", "gpi_call": gpi,
        "gpi_prob": f"{gpi_prob}", "gpi_omega": "", "gpi_fpr": "0.5", "gpi_svm": "-1.0",
    }  # fmt: skip


def _gz_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_work(root: Path, seed: int = 7) -> dict:
    """Write the work directory under root/work and the side files under root. Return names."""
    rng = random.Random(seed)
    root = Path(root)
    work = root / "work"
    (work / "phaseb" / "emb").mkdir(parents=True)
    seqs: dict[str, str] = {}  # hash -> sequence
    cls_of: dict[str, str] = {}  # hash -> signal class for features and embeddings
    truth, members = [], []
    named: dict[str, str] = {}
    families = [_seq(rng, 6, 0.0) for _ in range(6)]

    def add_member(set_id, source_id, gene_id, seq, cls):
        h = seqhash.seq_sha256(seq)
        seqs[h] = seq
        cls_of.setdefault(h, cls)
        members.append({"set_id": set_id, "source_id": source_id, "gene_id": gene_id,
                        "seq_sha256": h, "length": str(len(seq))})  # fmt: skip
        return h

    for src in SOURCES:
        source_id = src[0]
        npos, nint, nsec, pmtm, unres, ambig = src[5]
        k = 0
        plan = (
            [("P-ext", "wall" if i % 3 else "extracellular-only", "", "pos") for i in range(npos)]
            + [("N-int", "", "", "nint")] * nint
            + [("N-sec", "", "", "nsec")] * nsec
            + [("P-ext", "wall", "PM-TM", "nsec")] * pmtm
            + [("P-ext", "wall", "pm-unresolved", "pos")] * unres
            + [("ambiguous", "", "", "pos")] * ambig
        )
        for label, subset, d8, cls in plan:
            k += 1
            gene_id = f"{source_id[:4].upper()}{k:04d}"
            st = 0.30 if cls == "pos" else 0.12
            n = rng.randint(80, 400)
            prefix = families[k % 6] if k % 7 == 0 else ""
            if source_id == "Calb_CGD" and k == 1:
                n = 1022
            if source_id == "Calb_CGD" and k == 2:
                n = 1023
            seq = _seq(rng, n, st, prefix)
            hom = "yes" if k % 5 == 0 else "no"
            htp = label == "ambiguous" and k % 2 == 0
            truth.append(_truth_row(src, gene_id, label, subset, d8, hom, htp))
            h = add_member("truth", source_id, gene_id, seq, cls)
            if source_id == "Calb_CGD" and k in (1, 2):
                named[f"LONG_{n}"] = h
            if source_id == "Scer_SGD":
                key = {
                    "N-sec": "SHARED_NSEC",
                    "pm-unresolved": "SHARED_UNRES",
                    "ambiguous": "SHARED_AMBIG",
                    "P-ext": "SHARED_POS",
                }.get(d8 or label)
                if key and key not in named:
                    named[key] = h
            if source_id == "Spom_PomBase":
                truth.append(_truth_row(SOURCES[3], gene_id, label, subset, d8, hom, htp))
                add_member("truth", "Spom_SCHPO-mod", gene_id, seq, cls)
        # one unlabelled gene per source: never in the table
        k += 1
        seq = _seq(rng, 120, 0.1)
        truth.append(_truth_row(src, f"{source_id[:4].upper()}{k:04d}", "unlabelled", "", "",
                                "no", False))  # fmt: skip
        add_member("truth", source_id, truth[-1]["gene_id"], seq, "nint")
    alt_seq = _seq(rng, 150, 0.1)
    truth.append(_truth_row(SOURCES[3], "ALT0001", "N-int", "", "", "no", False))
    named["ALT_ONLY"] = add_member("truth", "Spom_SCHPO-mod", "ALT0001", alt_seq, "nint")

    # T-c rows (keyword tier) and the uniprot_kw members
    kw = []
    j = 0
    for taxon, genome, n in TC_TAXA:
        for _ in range(n):
            j += 1
            seq = _seq(rng, rng.randint(100, 300), 0.30)
            acc = f"Q{j:05d}"
            kw.append({"accession": acc, "gene": "", "genome": genome, "taxon_id": taxon,
                       "length": str(len(seq)), "seq_sha256": add_member(
                           "uniprot_kw", "uniprot_kw", acc, seq, "pos"), "tier": "T-c"})  # fmt: skip
    for key in ("SHARED_NSEC", "SHARED_UNRES", "SHARED_AMBIG", "SHARED_POS"):
        j += 1
        acc = f"Q{j:05d}"
        h = named[key]
        kw.append({"accession": acc, "gene": "", "genome": "Saccharomyces cerevisiae S288C",
                   "taxon_id": "559292", "length": str(len(seqs[h])), "seq_sha256": h,
                   "tier": "T-c"})  # fmt: skip
        add_member("uniprot_kw", "uniprot_kw", acc, seqs[h], cls_of[h])
    dup = _seq(rng, 210, 0.30)
    for acc in ("QDUP01", "QDUP02"):
        h = add_member("uniprot_kw", "uniprot_kw", acc, dup, "pos")
        kw.append({"accession": acc, "gene": "", "genome": "Candidozyma auris B8441",
                   "taxon_id": "498019", "length": str(len(dup)), "seq_sha256": h,
                   "tier": "T-c"})  # fmt: skip
    named["TC_DUP"] = h

    # literature seeds: an adhesin, a hard_negative, a moonlighting row, a row without accession
    lit_seqs = {"P90001": _seq(rng, 300, 0.35), "P90002": _seq(rng, 250, 0.30),
                "P90003": _seq(rng, 500, 0.08)}  # fmt: skip
    for acc, seq in lit_seqs.items():
        add_member("uniprot_kw", "uniprot_kw", acc, seq, "pos" if acc != "P90003" else "nint")
    seeds = root / "seeds.tsv"
    seeds.write_text(
        "# test seeds\n"
        "gene\tuniprot_query\tspecies\torder\tfamily\tclass\tevidence_level\tmoonlighting\t"
        "pmids\tevidence_summary\n"
        "LIT1\taccession:P90001\tCoccidioides immitis\tOnygenales\tf\tadhesin\tE1\tno\t1\ts\n"
        "LIT2\taccession:P90002\tAspergillus fumigatus\tEurotiales\tf\thard_negative\tN1\tno\t-\ts\n"
        "LIT3\taccession:P90003\tHistoplasma capsulatum\tOnygenales\tf\tadhesin\tE1\tYES\t1\ts\n"
        "LIT4\t\tHistoplasma capsulatum\tOnygenales\tf\tadhesin\tE2\tno\t1\ts\n"
    )

    # proteome sets: two truth sequences, one T-c sequence and new sequences
    scer_truth = [m for m in members if m["source_id"] == "Scer_SGD"][:2]
    for i, m in enumerate(scer_truth):
        add_member("Scer_proteome", "Scer_proteome", f"YPR{i:03d}W", seqs[m["seq_sha256"]], "pos")
    add_member("Scer_proteome", "Scer_proteome", "YPR900W", _seq(rng, 200, 0.1), "nint")
    for i in range(3):
        add_member("Cimm_RS_proteome", "Cimm_RS_proteome", f"CIMG_{i:05d}-t26_1-p1",
                   _seq(rng, 220, 0.25 if i == 0 else 0.1), "pos" if i == 0 else "nint")  # fmt: skip
    add_member(
        "Cimm_RS_proteome",
        "Cimm_RS_proteome",
        "CIMG_09999-t26_1-p1",
        seqs[kw[5]["seq_sha256"]],
        "pos",
    )

    # unique sequences, features, embeddings
    hashes = sorted(seqs)
    unique, fu, crow_of = [], [], {}
    crow = 0
    for i, h in enumerate(hashes):
        seq = seqs[h]
        long = len(seq) > 1022
        crow_of[h] = str(crow) if long else ""
        unique.append({"row": str(i), "seq_sha256": h, "length": str(len(seq)),
                       "cterm_row": crow_of[h], "sequence": seq})  # fmt: skip
        fu.append({"row": str(i), "seq_sha256": h, "length": str(len(seq)),
                   "cterm_row": crow_of[h], **_features(rng, cls_of[h], seq)})  # fmt: skip
        crow += long
    row_of = {h: str(i) for i, h in enumerate(hashes)}
    fu_by = {r["seq_sha256"]: r for r in fu}
    truth_by = {(t["source_id"], t["gene_id"]): t for t in truth}
    frows = []
    for m in members:
        t = truth_by.get((m["source_id"], m["gene_id"])) if m["set_id"] == "truth" else None
        tcols = {c: (t[c] if t else "") for c in
                 ("label", "subset", "stratum", "d8_class", "homology_only", "role")}  # fmt: skip
        f = fu_by[m["seq_sha256"]]
        frows.append({**m, **tcols, **{c: f[c] for c in FEATURE_COLUMNS},
                      "emb_row": row_of[m["seq_sha256"]],
                      "emb_cterm_row": crow_of[m["seq_sha256"]]})  # fmt: skip
    pb = work / "phaseb"
    truth_table.write_tsv(pb / "unique_sequences.tsv.gz", UNIQUE_COLUMNS, unique)
    truth_table.write_tsv(pb / "features_unique.tsv.gz", UNIQUE_FEATURE_COLUMNS, fu)
    truth_table.write_tsv(pb / "features.tsv.gz", MEMBER_FEATURE_COLUMNS, frows)
    truth_table.write_tsv(pb / "sequence_members.tsv.gz", MEMBER_COLUMNS, members)
    truth_table.write_tsv(work / "truth_set_triaged.tsv.gz", list(truth[0]), truth)
    truth_table.write_tsv(work / "keyword_tier.tsv.gz", KEYWORD_COLUMNS, kw)
    nprng = np.random.default_rng(seed)
    signal = {"pos": 1.5, "nsec": -0.5, "nint": -1.5}
    emb_models = {}
    n_long = sum(1 for h in hashes if crow_of[h] != "")
    for model, dim in MODELS.items():
        nterm = nprng.normal(size=(len(hashes), dim)).astype(np.float32)
        nterm[:, 0] += np.array([signal[cls_of[h]] for h in hashes], dtype=np.float32)
        cterm = nprng.normal(size=(n_long, dim)).astype(np.float32)
        np.save(pb / "emb" / f"{model}.nterm.npy", nterm)
        np.save(pb / "emb" / f"{model}.cterm.npy", cterm)
        emb_models[model] = {
            w: {"array_sha256": hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest(),
                "dtype": "float32", "shape": list(a.shape)}
            for w, a in (("nterm", nterm), ("cterm", cterm))
        }  # fmt: skip
    unique_sha = _sha(pb / "unique_sequences.tsv.gz")
    _gz_json(pb / "emb" / "embedding_run.json", {"models": emb_models,
             "unique_sequences": len(hashes), "unique_sequences_sha256": unique_sha})  # fmt: skip
    _gz_json(pb / "features_run.json", {
        "all_sources": True, "truth_set_sha256": TRUTH_SHA,
        "input_sha256": {"truth_set_triaged.tsv.gz": _sha(work / "truth_set_triaged.tsv.gz"),
                         "sequence_members.tsv.gz": _sha(pb / "sequence_members.tsv.gz"),
                         "unique_sequences.tsv.gz": unique_sha}})  # fmt: skip
    _gz_json(work / "d8_run.json", {"all_sources": True, "truth_set_sha256": TRUTH_SHA})
    _gz_json(work / "keyword_tier_run.json", {"all_sources": True, "truth_set_sha256": TRUTH_SHA})
    species = root / "species.tsv"
    truth_table.write_tsv(species, SPECIES_COLUMNS, [
        {"source_id": s[0], "species": s[1], "taxon_id": s[2], "taxon_filter": "",
         "in_clade": s[3], "role": s[4], "role_note": "", "gaf_file": "", "fasta_file": "",
         "id_mapping": ""} for s in SOURCES])  # fmt: skip
    clades = root / "tc_taxon_clades.tsv"
    truth_table.write_tsv(clades, ("taxon_id", "organism", "clade"), [
        {"taxon_id": t, "organism": g, "clade": TC_CLADES[t]} for t, g, _ in TC_TAXA])  # fmt: skip
    sets = root / "sequence_sets.tsv"
    sets.write_text("set_id\tkind\tlocation\tnote\ntruth\ttruth\ttruth_sequences.tsv.gz\t\n"
                    "uniprot_kw\tkeyword\tkeyword_sequences.fasta.gz\t\n"
                    "Scer_proteome\tdownload\tx.fasta.gz\t\n"
                    "Cimm_RS_proteome\tsite\tcocci:x.fasta\t\n")  # fmt: skip
    return {"work": work, "species": species, "seeds": seeds, "clades": clades, "sets": sets,
            "named": named, "seqs": seqs}  # fmt: skip


def build_argv(fx: dict) -> list[str]:
    return ["--work-dir", str(fx["work"]), "--species", str(fx["species"]),
            "--seeds", str(fx["seeds"]), "--tc-clades", str(fx["clades"])]  # fmt: skip


def rewrite_json(path: Path, **changes) -> None:
    obj = json.loads(Path(path).read_text())
    for key, value in changes.items():
        obj[key] = value
    Path(path).write_text(json.dumps(obj))


def gz_lines(path: Path) -> list[str]:
    return gzip.decompress(Path(path).read_bytes()).decode().splitlines()


FEATURE_COLUMNS = (
    "ser_thr_frac", "sp_prediction", "sp_prob", "sp_other_prob", "sp_cs_end", "sp_cs_prob",
    "gpi_call", "gpi_prob", "gpi_omega", "gpi_fpr", "gpi_svm",
)  # fmt: skip
UNIQUE_COLUMNS = ("row", "seq_sha256", "length", "cterm_row", "sequence")
UNIQUE_FEATURE_COLUMNS = ("row", "seq_sha256", "length", "cterm_row", *FEATURE_COLUMNS)
MEMBER_COLUMNS = ("set_id", "source_id", "gene_id", "seq_sha256", "length")
MEMBER_FEATURE_COLUMNS = (
    "set_id", "source_id", "gene_id", "seq_sha256", "length", "label", "subset", "stratum",
    "d8_class", "homology_only", "role", *FEATURE_COLUMNS, "emb_row", "emb_cterm_row",
)  # fmt: skip
KEYWORD_COLUMNS = ("accession", "gene", "genome", "taxon_id", "length", "seq_sha256", "tier")
