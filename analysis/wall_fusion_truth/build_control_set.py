#!/usr/bin/env python3
"""Build the control set data/controls/wall-fusion-kiss2019/ (plan Task 5.2; spec D8, D19, D20).

DRAFT, not run on HPCC (there is nothing to run there: this script is offline and stdlib only).
Source: analysis/wall_fusion_truth/results/row_mapping.tsv (Task 5.1), the Kiss table
docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv, results/organism_resolution.tsv and the
UniProt answers in results/query_cache.json. Nothing is read from the JGI files and no JGI
sequence can reach the output: the only sequences written are UniProt sequences from the cache.

Owner decisions of 2026-10-22 (docs/reports/2026-10-18-kiss-feasibility.md section 10):
  * all 133 strictly mapped rows are accepted; the flag jgi_accept (row_verified, locus_confirmed,
    portal_validated) stays in the table so a later analysis can filter on it;
  * the three name_conflict rows (hyr1, bhp2, bhp3) are not controls (rejected.tsv);
  * fine clade is the primary level, coarse the second; cluster_stratum_fine and
    cluster_stratum_coarse say whether the Kiss cluster of the row also has members in another clade;
  * fusion rows are descriptive only (stratum fusion);
  * family representatives (35 Afu rows) are never test controls: they go to family_reps.tsv and
    family_rep_sequences.faa, not to controls.tsv, and never count toward a floor.

Files written to the output directory:
  controls.tsv                one row per protein with a public (UniProt) sequence; the L1 validator checks it
  controls_local_sequence.tsv rows whose sequence is not public (no UniProt entry): same columns, key
                              `jgi:<portal>:<proteinId>`; the sequence stays in the JGI file on the owner's disk
  sequences.faa               UniProt sequences of controls.tsv
  clusters.tsv                accession, cluster_id (Kiss cluster per slice: a PROXY until Task 5.3)
  family_reps.tsv, family_rep_sequences.faa   the 35 family representatives (public UniProt sequences)
  rejected.tsv, counts.md, unverified_notes.md, README.md

Task 5.3 gives MMseqs2 clusters and seed classes. Run this script again with
--clusters-tsv and --split-members (written by seed_split.py) to replace the proxy.

Usage:  python3 -I analysis/wall_fusion_truth/build_control_set.py [--out DIR] [--date YYYY-MM-DD]
            [--clusters-tsv clusters.tsv] [--split-members members_split.tsv]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
KISS_TSV = ROOT / "docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv"
OUT_DEFAULT = ROOT / "data/controls/wall-fusion-kiss2019"
DATE_DEFAULT = "2026-10-22"
CURATOR = "claude-sonnet-5-5 (build_control_set.py; owner decisions of 2026-10-22)"
PMID = "PMID:31501435"
STRICT = ("mapped_jgi_file", "mapped_locus_tag")
EVIDENCE_ORDER = ("row_verified", "locus_tag", "locus_confirmed", "portal_validated")
SERVES = {
    "adhesion": "adhesion_family_recall",
    "wall": "wall_model_recall,narrow_false_hit_panel",
    "fusion": "fusion_descriptive",
}
FAMILY_REP_CODE = "Aspfu1"  # the family representatives are A. fumigatus loci
AA_RUN = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{40,}")
QUOTE_MAX = 220

BASE_COLUMNS = [
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
]
EXTRA_COLUMNS = [
    "control_role",
    "serves",
    "clade_coarse",
    "clade_fine",
    "strain_taxid",
    "kiss_cluster_number",
    "cluster_id",
    "cluster_stratum_fine",
    "cluster_stratum_coarse",
    "jgi_accept",
    "map_status",
    "protein_id_jgi",
    "jgi_seq_sha256",
    "sequence_sha256",
    "sequence_source",
    "sequence_matches_jgi",
    "sequence_in_repo",
    "seed_class_slice",
    "seed_class_all39",
]
COLUMNS = BASE_COLUMNS + EXTRA_COLUMNS
FAMILY_COLUMNS = [
    "accession",
    "gene_or_family_id",
    "kiss_cluster_number",
    "slice",
    "afu_locus",
    "species",
    "taxon_id",
    "strain_taxid",
    "control_role",
    "may_be_test_control",
    "counts_toward_floor",
    "accession_also_in_controls",
    "sequence_sha256",
    "source",
    "quote",
]
REJECTED_COLUMNS = [
    "gene_or_family_id",
    "kiss_cluster_number",
    "slice",
    "protein_id_jgi",
    "organism_code",
    "map_status",
    "control_role",
    "reason",
    "see",
]


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(c, "") for c in header])


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def uniprot_sequences(cache_path: Path) -> dict[str, dict]:
    """{accession: {sequence, entry_name, genes}} from every UniProt answer in the cache.

    The cache holds the public UniProt TSV answers (columns of map_kiss.UNIPROT_FIELDS:
    accession, reviewed, id, gene_names, organism_name, organism_id, length, sequence).
    An accession with two different sequences in the cache is an error.
    """
    cache = json.loads(cache_path.read_text())
    out: dict[str, dict] = {}
    for url in sorted(cache):
        text = cache[url]
        if "rest.uniprot.org" not in url or not text.startswith("Entry"):
            continue
        for line in text.strip("\n").split("\n")[1:]:
            f = line.split("\t")
            if len(f) < 8 or not f[7]:
                continue
            seq = f[7].strip().upper()
            if f[0] in out and out[f[0]]["sequence"] != seq:
                raise ValueError(f"{f[0]}: two different sequences in the cache")
            out[f[0]] = {"sequence": seq, "entry_name": f[2], "genes": f[3]}
    return out


def kiss_quotes(path: Path) -> dict[tuple, dict]:
    """{(gene, cluster, jgi id, category): row} of the Kiss table; the first row of a key wins."""
    out: dict[tuple, dict] = {}
    for r in read_tsv(path):
        key = (r["gene_or_family_id"], r["cluster_number"], r["protein_id_jgi"], r["category"])
        out.setdefault(key, r)
    return out


def make_quote(row: dict, kiss: dict) -> tuple[str, str]:
    """(source, quote): the field of the Kiss table that supports the row."""
    k = kiss.get(
        (row["gene_or_family_id"], row["cluster_number"], row["protein_id_jgi"], row["category"])
    )
    sheet = k["source_sheet"] if k else "Kiss 2019 supplementary table"
    ann = (k["annotation"] if k else "").strip()
    if len(ann) > QUOTE_MAX:
        ann = ann[: QUOTE_MAX - 3] + "..."
    quote = (
        f"{sheet}; category {row['category']}; cluster {row['cluster_number']}; "
        f"annotation: {ann or '(none)'}"
    )
    source = f"{PMID}; {sheet}; docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv"
    return source, quote


def evidence_of(row: dict) -> str:
    return row["jgi_accept"] or "locus_tag"


def clade_strata(rows: list[dict]) -> dict[tuple, dict]:
    """{(slice, cluster): {fine: shared|clade_only, coarse: ...}} over the strictly mapped rows."""
    seen: dict[tuple, dict] = defaultdict(lambda: {"fine": set(), "coarse": set()})
    for r in rows:
        if r["map_status"] in STRICT and r["slice"]:
            seen[(r["slice"], r["cluster_number"])]["fine"].add(r["clade_fine"])
            seen[(r["slice"], r["cluster_number"])]["coarse"].add(r["clade_coarse"])
    return {
        k: {lv: ("shared" if len(v[lv]) > 1 else "clade_only") for lv in ("fine", "coarse")}
        for k, v in seen.items()
    }


def jgi_key(row: dict) -> str:
    m = re.match(r"^(\d+)_(.+)$", row["protein_id_jgi"].strip())
    if not m:
        raise ValueError(
            f"{row['gene_or_family_id']}: no JGI protein ID in {row['protein_id_jgi']!r}"
        )
    return f"jgi:{row['jgi_portal'] or m.group(2)}:{m.group(1)}"


def build(rows: list[dict], orgs: dict[str, dict], seqs: dict[str, dict], kiss: dict, date: str):
    """Return (controls, local, family_reps, rejected, family_sequences) as lists of dicts."""
    strata = clade_strata(rows)
    by_sha: dict[str, list[str]] = defaultdict(list)
    for acc, rec in seqs.items():
        by_sha[sha256_text(rec["sequence"])].append(acc)
    found_by_sequence: set[str] = set()
    groups: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        if r["map_status"] in STRICT and r["control_role"] == "kiss_row":
            acc = r["accession"]
            if not acc and len(by_sha.get(r["jgi_seq_sha256"], [])) == 1:
                # no accession by name, but one cached UniProt entry has this very sequence
                acc = by_sha[r["jgi_seq_sha256"]][0]
                found_by_sequence.add(acc)
            groups[acc if acc else jgi_key(r)].append(r)

    controls, local = [], []
    for key in sorted(groups, key=lambda k: _sort_key(groups[k][0], k)):
        g = groups[key]
        slices = {r["slice"] for r in g}
        if len(slices) != 1:
            raise ValueError(f"{key}: the same protein is in two slices {sorted(slices)}")
        first = g[0]
        org = orgs[first["organism_code"]]
        genes = "/".join(dict.fromkeys(r["gene_or_family_id"] for r in g))
        clusters = sorted({r["cluster_number"] for r in g}, key=int)
        evidence = min((evidence_of(r) for r in g), key=EVIDENCE_ORDER.index)
        accept = next((e for e in EVIDENCE_ORDER if any(r["jgi_accept"] == e for r in g)), "")
        public = key in seqs
        sha_jgi = next((r["jgi_seq_sha256"] for r in g if r["jgi_seq_sha256"]), "")
        if public:
            seq = seqs[key]["sequence"]
            sha = sha256_text(seq)
            matches = ("yes" if sha == sha_jgi else "no") if sha_jgi else "not_checked"
            source_of = "uniprot"
        else:
            sha, matches, source_of = sha_jgi, "yes", "jgi_file_local"
            if not sha:
                raise ValueError(f"{key}: no public sequence and no JGI sequence sha256")
        source, quote = make_quote(first, kiss)
        notes = []
        if len(g) > 1:
            notes.append(
                "merged Kiss rows of one protein: "
                + ", ".join(f"{r['gene_or_family_id']} (cluster {r['cluster_number']})" for r in g)
            )
        if first["accession"] and not public:
            notes.append("UniProt accession without a sequence in the cache")
        if key in found_by_sequence and not first["accession"]:
            notes.append(
                "accession found by an identical sequence (sha256) in the UniProt answers cached "
                "in Task 5.1, not by name"
            )
        if matches == "no":
            notes.append(
                "the UniProt sequence differs from the JGI protein (map_kiss jgi_vs_uniprot "
                "similar); the sequence of this row is the UniProt one"
            )
        if accept == "portal_validated":
            notes.append(
                "no row-level name evidence; accepted because its portal is validated by other rows"
            )
        if accept == "locus_confirmed":
            notes.append(
                "the UniProt entry of the JGI header locus tag carries only that tag as gene name"
            )
        rec = {
            "accession": key,
            "gene": genes,
            "species": org["ncbi_name"],
            "taxon_id": org["ncbi_taxid"],
            "label": "positive",
            "stratum": first["slice"],
            "evidence_level": evidence,
            "source": source,
            "quote": quote,
            "tuned": "no",
            "leakage": "unknown",
            "curator": CURATOR,
            "date": date,
            "notes": "; ".join(notes),
            "control_role": "kiss_row",
            "serves": SERVES[first["slice"]],
            "clade_coarse": first["clade_coarse"],
            "clade_fine": first["clade_fine"],
            "strain_taxid": org["scope_taxid"],
            "kiss_cluster_number": "+".join(clusters),
            "cluster_id": f"KISS-{first['slice']}-{clusters[0]}",
            "cluster_stratum_fine": strata[(first["slice"], clusters[0])]["fine"],
            "cluster_stratum_coarse": strata[(first["slice"], clusters[0])]["coarse"],
            "jgi_accept": accept,
            "map_status": first["map_status"],
            "protein_id_jgi": "+".join(dict.fromkeys(r["protein_id_jgi"] for r in g)),
            "jgi_seq_sha256": sha_jgi,
            "sequence_sha256": sha,
            "sequence_source": source_of,
            "sequence_matches_jgi": matches,
            "sequence_in_repo": "yes" if public else "no",
            "seed_class_slice": "",
            "seed_class_all39": "",
        }
        (controls if public else local).append(rec)

    controls_acc = {c["accession"] for c in controls}
    reps, fam_seqs = [], {}
    for r in sorted(
        (x for x in rows if x["control_role"] == "family_rep"),
        key=lambda x: (x["slice"], int(x["cluster_number"]), x["gene_or_family_id"]),
    ):
        if r["accession"] not in seqs:
            raise ValueError(f"family representative {r['accession']} has no sequence in the cache")
        org = orgs[FAMILY_REP_CODE]
        source, quote = make_quote(r, kiss)
        fam_seqs[r["accession"]] = seqs[r["accession"]]["sequence"]
        reps.append(
            {
                "accession": r["accession"],
                "gene_or_family_id": r["gene_or_family_id"],
                "kiss_cluster_number": r["cluster_number"],
                "slice": r["slice"],
                "afu_locus": r["match_term"],
                "species": org["ncbi_name"],
                "taxon_id": org["ncbi_taxid"],
                "strain_taxid": r["uniprot_organism_id"],
                "control_role": "family_rep",
                "may_be_test_control": "no",
                "counts_toward_floor": "no",
                "accession_also_in_controls": "yes" if r["accession"] in controls_acc else "no",
                "sequence_sha256": sha256_text(seqs[r["accession"]]["sequence"]),
                "source": source,
                "quote": quote,
            }
        )

    rejected = []
    for r in sorted(
        (x for x in rows if x["control_role"] != "family_rep" and x["map_status"] not in STRICT),
        key=lambda x: (x["slice"], int(x["cluster_number"]), x["gene_or_family_id"]),
    ):
        conflict = "differs_from_name_candidate" in r["map_reason"]
        reason = r["map_reason"] or r["map_status"]
        rejected.append(
            {
                "gene_or_family_id": r["gene_or_family_id"],
                "kiss_cluster_number": r["cluster_number"],
                "slice": r["slice"],
                "protein_id_jgi": r["protein_id_jgi"],
                "organism_code": r["organism_code"],
                "map_status": r["map_status"],
                "control_role": r["control_role"],
                "reason": ("name_conflict: " if conflict else "") + reason,
                "see": "results/name_conflict_review.tsv"
                if conflict
                else "results/row_mapping.tsv",
            }
        )
    return controls, local, reps, rejected, fam_seqs


def _sort_key(row: dict, key: str):
    return (row["slice"], int(row["cluster_number"]), row["gene_or_family_id"], key)


def guard_tables_have_no_sequences(paths: list[Path]) -> None:
    """Stop when a table cell holds a run of 40 or more amino-acid letters (a sequence)."""
    for p in paths:
        for n, line in enumerate(p.read_text().splitlines(), 1):
            for cell in line.split("\t"):
                if AA_RUN.search(cell):
                    raise RuntimeError(f"{p.name} line {n}: a cell looks like a protein sequence")


def apply_clusters(controls: list[dict], local: list[dict], mapping: dict[str, str]) -> None:
    missing = [c["accession"] for c in controls + local if c["accession"] not in mapping]
    if missing:
        raise ValueError(
            f"{len(missing)} accessions are not in the cluster file, e.g. {missing[0]}"
        )
    for c in controls + local:
        c["cluster_id"] = mapping[c["accession"]]


def apply_split(controls: list[dict], local: list[dict], members: list[dict]) -> None:
    by = {m["accession"]: m for m in members}
    for c in controls + local:
        m = by.get(c["accession"])
        if m is None:
            raise ValueError(f"{c['accession']} is not in the seed split table")
        c["seed_class_slice"] = m["seed_class_slice"]
        c["seed_class_all39"] = m["seed_class_all39"]
        c["leakage"] = m["leakage_suggested"]


def read_fasta_text(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    cur = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif cur is not None:
            out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def fasta_text(seqs: dict[str, str]) -> str:
    return "".join(f">{a}\n{seqs[a]}\n" for a in sorted(seqs))


def count_lines(controls: list[dict], local: list[dict], reps: list[dict], rejected: list[dict]):
    def clusters(rows):
        return len({r["cluster_id"] for r in rows})

    return {
        "rows": len(controls),
        "clusters": clusters(controls),
        "local_rows": len(local),
        "local_clusters": clusters(local),
        "all_clusters": clusters(controls + local),
        "reps": len(reps),
        "rep_sequences": len({r["accession"] for r in reps}),
        "rejected": len(rejected),
    }


def counts_md(controls, local, reps, rejected, proxy: bool) -> str:
    n = count_lines(controls, local, reps, rejected)
    both = controls + local
    cl_label = (
        "Independent clusters (Kiss cluster per slice; a proxy until Task 5.3)"
        if proxy
        else "Independent clusters (MMseqs2, 30% identity, 50% coverage)"
    )
    by_slice = Counter(r["stratum"] for r in both)
    by_slice_pub = Counter(r["stratum"] for r in controls)
    by_ev = Counter(r["evidence_level"] for r in both)
    lines = [
        "# Counts: wall, fusion and adhesion controls from Kiss 2019",
        "",
        "*Written by `analysis/wall_fusion_truth/build_control_set.py` (DRAFT, Task 5.2). "
        "Every row is `label=positive` for the call named in `stratum` and `serves`. "
        "A wall row is never a negative and never a truth row of the adhesion call (D20). "
        "No number is pooled across slices.*",
        "",
        "## Summary (`controls.tsv`: rows with a public UniProt sequence; checked by the L1 validator)",
        "",
        "| Item | Value |",
        "|---|---|",
        f"| Rows (proteins) | {n['rows']} |",
        f"| Sequences in `sequences.faa` | {n['rows']} |",
        f"| {cl_label} | {n['clusters']} |",
        f"| Rows positive | {n['rows']} |",
        f"| Sequences positive | {n['rows']} |",
        f"| Clusters positive | {n['clusters']} |",
        "| Rows negative | 0 |",
        "| Sequences negative | 0 |",
        "| Clusters negative | 0 |",
        "",
        "This set has no negatives on purpose. The hard negatives are the shared set of 161 rows "
        "(`data/controls/hard-negatives-shared`; owner decision 2026-10-22, that set only). "
        "**Anything that reads this set must read both `controls.tsv` and `controls_local_sequence.tsv`.**",
        "",
        "## Rows by slice",
        "",
        "| Slice | In `controls.tsv` | In `controls_local_sequence.tsv` | Both |",
        "|---|---|---|---|",
    ]
    for sl in ("adhesion", "wall", "fusion"):
        lines.append(
            f"| {sl} | {by_slice_pub[sl]} | {by_slice[sl] - by_slice_pub[sl]} | {by_slice[sl]} |"
        )
    lines += [
        "",
        "## Rows kept out of `controls.tsv` because the sequence is not public",
        "",
        f"`controls_local_sequence.tsv` holds {n['local_rows']} rows in {n['local_clusters']} Kiss "
        "clusters. Their sequence is in a JGI file that stays on the owner's disk (key "
        "`jgi:<portal>:<proteinId>`; `sequence_sha256` is the check). The L1 validator "
        "does not read this file: it requires every accession to have a sequence in "
        "`sequences.faa`. `build_kiss_fasta.py` joins them back for Task 5.3.",
        f"Both files together: {len(both)} rows in {n['all_clusters']} clusters.",
        "",
        "## Evidence level (both files)",
        "",
        "| evidence_level | rows |",
        "|---|---|",
    ]
    for e in EVIDENCE_ORDER:
        lines.append(f"| {e} | {by_ev[e]} |")
    lines += [
        "",
        "## Other files",
        "",
        "| File | Rows |",
        "|---|---|",
        f"| `family_reps.tsv` (never test controls) | {n['reps']} rows, {n['rep_sequences']} sequences |",
        f"| `rejected.tsv` | {n['rejected']} |",
        "",
        "## Size targets (COMMON-RULES rule 7)",
        "",
        "The rule gives 20 clusters per class for `estimated`. The L1 validator reports the "
        "distance only. A clade floor of 5 clusters per stratum applies to the hold-out "
        "(owner decision 2026-10-22); see `docs/reports/2026-10-18-kiss-feasibility.md` "
        "section 10 and `analysis/wall_fusion_truth/results/clusters_per_clade.tsv`.",
        "",
    ]
    return "\n".join(lines)


UNVERIFIED = """# Notes that are not checked

*Written by `build_control_set.py` (DRAFT, Task 5.2).*

- **Clusters.** `clusters.tsv` uses the Kiss `cluster_number` per slice (`KISS-<slice>-<n>`). It
  is an orthology and MCL family, not a 30% identity cluster. COMMON-RULES rule 4 asks for
  MMseqs2. Task 5.3 makes them. Until then every count of clusters is a proxy.
- **Leakage.** `leakage` is `unknown` for every row: nobody has compared these proteins with the
  Pfam seed sequences yet. Task 5.3 (`seed_split.py`) sets `none` or `in_reference`,
  `seed_class_slice` and `seed_class_all39`.
- **portal_validated rows.** These rows have no name evidence of their own. They are accepted by
  the owner (2026-10-22) because their portal agrees on other rows. The flag stays; filter on it
  if a result looks sensitive.
- **locus_confirmed rows.** The UniProt entry found by the JGI header locus tag has that tag as its
  only gene name. The file protein equals the entry. This does not show that the Kiss row means
  this protein.
- **Sequence differs from JGI.** Where `sequence_matches_jgi` is `no`, the committed sequence is
  the UniProt one, not the Kiss (JGI) one.
- **Accession found by sequence.** For a row with no accession by name, the accession is taken from
  the one cached UniProt entry whose sequence has the same sha256 as the JGI protein (the `notes`
  column says so). The Kiss gene name was not used and may differ from the UniProt gene name.
- **Species names.** `species` and `taxon_id` come from NCBI Taxonomy answers cached by `map_kiss.py`
  (`results/organism_resolution.tsv`). They were not checked against `names.dmp`.
- **Wall rows.** The wall rows are a descriptive panel for the wall models and a false-hit panel for
  the narrow call. They are positives for the wall slice only.
- **Family representatives.** The 35 Afu rows stand for a Kiss family. The Afu protein is the
  ortholog named in the Kiss table, not the protein that the authors measured.
- **Bcine1.** Two Botrytis rows (bhp2, bhp3) are not controls: the 2025 annotation renumbered the
  JGI IDs (see `rejected.tsv`).
"""


def readme_md(controls, local, reps, rejected, proxy: bool) -> str:
    n = count_lines(controls, local, reps, rejected)
    return f"""# Control set: wall, fusion and adhesion genes of Kiss et al. 2019

*Written by `analysis/wall_fusion_truth/build_control_set.py` (DRAFT, plan Task 5.2). Do not edit
the files by hand; run the script again.*

Source: Kiss et al. 2019, Nat Commun ({PMID}), Supplementary Data 3, 4 and 8, as read in
`docs/reports/cell-wall-fusion/06-kiss2019-hypha-genes.tsv` and mapped to sequences by
`analysis/wall_fusion_truth/map_kiss.py` (Task 5.1; report
`docs/reports/2026-10-18-kiss-feasibility.md`). The owner decisions of 2026-10-22 (section 10 of
that report) bind this set.

## Files

| File | Content |
|---|---|
| `controls.tsv` | {n['rows']} positive controls with a public UniProt sequence. Checked by the L1 validator. |
| `controls_local_sequence.tsv` | {n['local_rows']} positive controls whose sequence is not public. Same columns. The validator does not read it. **These two files are one set: anything that reads the set must read both.** |
| `sequences.faa` | UniProt sequences of `controls.tsv`. Header = accession. |
| `clusters.tsv` | `accession<TAB>cluster_id`. {"Kiss cluster per slice: a proxy until Task 5.3." if proxy else "MMseqs2 clusters (Task 5.3)."} |
| `family_reps.tsv`, `family_rep_sequences.faa` | {n['reps']} family representatives (Afu orthologs). Training seeds and baseline report only; not clustered; each keeps its Kiss `cluster_number`; a flagged group in the report. |
| `rejected.tsv` | {n['rejected']} Kiss rows that are not controls, with the reason. |
| `counts.md`, `unverified_notes.md` | counts and unchecked claims. |

## Keys

The key (`accession`) is the UniProt accession wherever the row has one (owner preference: NCBI or
UniProt accessions, not JGI IDs). A row with no UniProt accession gets the local key
`jgi:<portal>:<proteinId>` (for example `jgi:Neucr2:4224`). A local key is not a public identifier.
Its sequence is not in this repository (JGI data are never committed); `sequence_sha256` is the
check, and `build_kiss_fasta.py` reads the sequence from the downloaded JGI file.

## Columns

The required columns of `docs/agent-tasks/COMMON-RULES.md`, then: `control_role` (`kiss_row`),
`serves` (the calls the row may be used for), `clade_coarse`, `clade_fine`, `strain_taxid`,
`kiss_cluster_number`, `cluster_id`, `cluster_stratum_fine` and `cluster_stratum_coarse`
(`shared`: the Kiss cluster of this slice has strictly mapped members in another clade at that
level; `clade_only`: it has not), `jgi_accept` (`row_verified`, `locus_confirmed`,
`portal_validated`, or empty for a row mapped by a locus tag of its own), `map_status`,
`protein_id_jgi`, `jgi_seq_sha256`, `sequence_sha256`, `sequence_source`, `sequence_matches_jgi`,
`sequence_in_repo`, `seed_class_slice` and `seed_class_all39` (empty until Task 5.3; `seed-like` when the
Kiss cluster has any member at 30% identity or more to a Pfam seed: of the row's own slice models, or of
all 39 models).

`label` is `positive` for the call named in `stratum` (adhesion, wall or fusion). Read the rules
before you use a row:

- adhesion rows: positives of the adhesion-family recall;
- wall rows: positives for the wall models only, and a descriptive false-hit panel for the
  narrow call. A wall row is never a negative and never a truth row of the adhesion call;
- fusion rows: descriptive only (at most 2 clusters per clade: no floor is met);
- no number is pooled across slices (D20).

## Rules of this set

- All strictly mapped rows are accepted (owner, 2026-10-22). The three `jgi_accept` flags stay;
  filter later if a result looks sensitive.
- hyr1, bhp2 and bhp3 are not controls (`name_conflict`, listed in `rejected.tsv`).
- The 35 family representatives are not in `controls.tsv`. They may be training seeds and appear
  in the baseline recall report; they are never test controls and never count toward a floor.
  They are not clustered (the `--include-family-reps` option of `build_kiss_fasta.py` exists and is not used).
- The hard negatives for this set are the shared set of 161 rows, and only that set (owner decision
  2026-10-22; `data/curated/adhesins/hard_negative_seeds.tsv` is excluded and may be added later).
  Eight proteins are both a control here and a shared negative; they keep both roles (the negative
  copy is `<accession>#neg` in the clustering input) and their clusters are flagged `mixed_pos_neg`.
- This set is two files. Anything that reads it must read `controls.tsv` and `controls_local_sequence.tsv`.
- Two Kiss rows that point to the same protein are one row here (`notes` says which).
- `leakage` is `unknown` until `seed_split.py` has compared each row with the Pfam seeds.

## How it was built

```
python3 -I analysis/wall_fusion_truth/build_control_set.py
```

The script is offline and deterministic (two runs give the same bytes). It reads
`analysis/wall_fusion_truth/results/{{row_mapping,organism_resolution}}.tsv`,
`results/query_cache.json` (the UniProt sequences) and the Kiss table. It never reads the JGI
files. After Task 5.3 run it again with `--clusters-tsv` and `--split-members` (written by
`seed_split.py`) to replace the proxy clusters and fill `leakage`, `seed_class_slice` and `seed_class_all39`.

Check: `python3 -m cellsurface_sorting_hat.controls data/controls/wall-fusion-kiss2019`
(with `PYTHONPATH=src`).
"""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--date", default=DATE_DEFAULT)
    ap.add_argument("--row-mapping", type=Path, default=RESULTS / "row_mapping.tsv")
    ap.add_argument("--organisms", type=Path, default=RESULTS / "organism_resolution.tsv")
    ap.add_argument("--cache", type=Path, default=RESULTS / "query_cache.json")
    ap.add_argument("--kiss-table", type=Path, default=KISS_TSV)
    ap.add_argument("--clusters-tsv", type=Path, help="accession<TAB>cluster_id from seed_split.py")
    ap.add_argument("--split-members", type=Path, help="members_split.tsv from seed_split.py")
    a = ap.parse_args(argv)

    rows = read_tsv(a.row_mapping)
    orgs = {r["code"]: r for r in read_tsv(a.organisms)}
    seqs = uniprot_sequences(a.cache)
    controls, local, reps, rejected, fam_seqs = build(
        rows, orgs, seqs, kiss_quotes(a.kiss_table), a.date
    )
    proxy = a.clusters_tsv is None
    if a.clusters_tsv:
        mapping = {}
        for line in a.clusters_tsv.read_text().splitlines():
            f = line.split("\t")
            if len(f) >= 2 and f[0] != "accession":
                mapping[f[0]] = f[1]
        apply_clusters(controls, local, mapping)
    if a.split_members:
        apply_split(controls, local, read_tsv(a.split_members))

    out = a.out
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "controls.tsv", COLUMNS, controls)
    write_tsv(out / "controls_local_sequence.tsv", COLUMNS, local)
    (out / "sequences.faa").write_text(
        fasta_text({c["accession"]: seqs[c["accession"]]["sequence"] for c in controls})
    )
    (out / "clusters.tsv").write_text(
        "accession\tcluster_id\n"
        + "".join(f"{c['accession']}\t{c['cluster_id']}\n" for c in controls)
    )
    write_tsv(out / "family_reps.tsv", FAMILY_COLUMNS, reps)
    (out / "family_rep_sequences.faa").write_text(fasta_text(fam_seqs))
    write_tsv(out / "rejected.tsv", REJECTED_COLUMNS, rejected)
    (out / "counts.md").write_text(counts_md(controls, local, reps, rejected, proxy))
    (out / "unverified_notes.md").write_text(UNVERIFIED)
    (out / "README.md").write_text(readme_md(controls, local, reps, rejected, proxy))
    guard_tables_have_no_sequences(
        [
            out / n
            for n in (
                "controls.tsv",
                "controls_local_sequence.tsv",
                "family_reps.tsv",
                "rejected.tsv",
            )
        ]
    )
    n = count_lines(controls, local, reps, rejected)
    print(
        f"controls.tsv {n['rows']} rows in {n['clusters']} clusters; local-sequence rows {n['local_rows']}; "
        f"family reps {n['reps']}; rejected {n['rejected']}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
