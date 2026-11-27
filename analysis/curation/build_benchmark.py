"""Build the frozen adhesion benchmark (`data/benchmark/adhesion_<version>/`) from the curated sources.

Owner decisions of 2026-11-26 (plan `docs/superpowers/plans/2026-11-26-adhesion-round2-benchmark-cruft-and-layout-study.md`):

* **Tiers, not one label.** A protein lands in exactly one tier:
  - ``positive`` - direct experimental adhesion evidence AND a secretory architecture to predict from;
  - ``positive_moonlighting`` - direct evidence but a cytoplasmic protein found on the surface (no signal peptide, no
    anchor). Kept on record, left out of both contrasts and of the layout study, because no architecture method can
    predict it;
  - ``hard_negative`` - a surface or wall protein with a demonstrated non-adhesion function. Includes proteins whose
    deletion reduces adhesion by an indirect route (``indirect=yes``): they are the trap a predictor must not fall into;
  - ``excluded`` - everything that fails a rule, with ``exclude_reason`` saying which. Nothing is dropped silently.
* **Two contrasts.** ``C1`` (positives vs hard negatives) decides whether a component is kept; ``C2`` (positives vs all
  other secreted or GPI proteins of the same proteomes) is descriptive and is assembled at scoring time, not here.
* **Surface filter.** A hard negative must look like a surface protein (signal peptide, GPI anchor, or an explicit
  ``surface_evidence`` note). This is what keeps transcription factors out of a "surface protein" negative set, where
  they would make the contrast easier than it is.
* **Frozen.** The output carries a VERSION and SHA256SUMS; a later change is a new version directory, never an edit.

Inputs:

* ``data/curated/adhesins/adhesins.tsv`` - the existing curated table (``cls``, ``evidence_level``, ``signal_peptide``,
  ``gpi_anchor``, ``surface`` columns).
* ``data/benchmark/<version>/sources/*.tsv`` - literature curation passes in the 15-column schema that
  ``REQUIRED_SOURCE_COLUMNS`` lists, one file per clade group, each row carrying its PMIDs and evidence code.
* ``data/benchmark/clade_map.tsv`` - reviewed taxon-to-class map for the leave-one-clade-out splits. A scored row
  whose taxon is not mapped is excluded rather than guessed, so a new organism must be classified on purpose.

Usage::

    python3 -I analysis/curation/build_benchmark.py --version adhesion_v1 --out data/benchmark/adhesion_v1
    python3 -I analysis/curation/build_benchmark.py --version adhesion_v1 --out /tmp/check --dry-run
"""

import argparse
import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

REQUIRED_SOURCE_COLUMNS = (
    "accession",
    "accession_db",
    "gene",
    "protein_name",
    "organism",
    "taxon_id",
    "clade_class",
    "label",
    "mechanism",
    "evidence_code",
    "assay_summary",
    "pmids",
    "uniprot_length",
    "verified",
    "notes",
)

OUTPUT_COLUMNS = (
    "accession",
    "accession_db",
    "gene",
    "protein_name",
    "organism",
    "taxon_id",
    "clade_class",
    "tier",
    "mechanism",
    "evidence_code",
    "indirect",
    "assay_summary",
    "pmids",
    "length",
    "signal_peptide",
    "gpi_anchor",
    "source",
    "exclude_reason",
    "notes",
)

DIRECT_EVIDENCE = frozenset({"ko_assay", "het_expr", "binding"})
NEGATIVE_EVIDENCE = frozenset({"function_shown"})
LABELS = frozenset({"adhesin", "hard_negative"})
TIERS = ("positive", "positive_moonlighting", "hard_negative", "excluded")

# Classes of the existing curated table and how they enter the benchmark.
CURATED_POSITIVE_LEVELS = frozenset({"E1"})  # E2/E3 are family or weak evidence, not direct
CURATED_NEGATIVE_CLASSES = frozenset({"hard_negative", "surface_other_adhesion_phenotype"})
INDIRECT_CLASS = "surface_other_adhesion_phenotype"

# A fold must hold enough of both sides to give a usable number. Below either threshold the class is
# scored inside its higher group; a group that is still too thin is reported as "not scored" with its
# counts, never silently folded into another clade.
MIN_FOLD_POSITIVES = 5
MIN_FOLD_NEGATIVES = 5


class BuildError(Exception):
    pass


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_clade_map(path):
    """{taxon_id: clade_class} from the reviewed map. Comment lines start with '#'."""
    out = {}
    with open(path, newline="") as fh:
        lines = [ln for ln in fh.read().splitlines() if not ln.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        taxon = (row.get("taxon_id") or "").strip()
        clade = (row.get("clade_class") or "").strip()
        if not taxon or not clade:
            raise BuildError(f"{path}: a row has no taxon_id or no clade_class")
        if taxon in out and out[taxon] != clade:
            raise BuildError(f"{path}: taxon {taxon} is mapped to two classes")
        out[taxon] = clade
    if not out:
        raise BuildError(f"{path}: no clade rows")
    return out


def read_clade_groups(path):
    """{clade_class: group} from the reviewed map. Comment lines start with '#'."""
    out = {}
    with open(path, newline="") as fh:
        lines = [ln for ln in fh.read().splitlines() if not ln.startswith("#")]
    for row in csv.DictReader(lines, delimiter="\t"):
        clade = (row.get("clade_class") or "").strip()
        group = (row.get("group") or "").strip()
        if not clade or not group:
            raise BuildError(f"{path}: a row has no clade_class or no group")
        out[clade] = group
    if not out:
        raise BuildError(f"{path}: no group rows")
    return out


def plan_folds(rows, groups, min_pos=MIN_FOLD_POSITIVES, min_neg=MIN_FOLD_NEGATIVES):
    """Decide what each clade is scored as. Returns a list of fold rows, most positives first.

    A class that reaches both thresholds is its own fold. Otherwise its proteins count towards its
    group, and the group is a fold if the pooled counts reach the thresholds. Anything left is
    reported with ``scored = no`` and its counts, so a thin clade is visible rather than assumed.

    ``n_positive`` counts proteins here. Clusters (MMseqs2, 30% identity) are the real unit and are
    filled in before the set is frozen; until then a fold list built on proteins is provisional.
    """
    scored = [r for r in rows if r["tier"] in ("positive", "hard_negative")]
    by_class = {}
    for r in scored:
        c = r["clade_class"]
        key = by_class.setdefault(c, {"positive": 0, "hard_negative": 0})
        key[r["tier"]] += 1

    unmapped = sorted(c for c in by_class if c not in groups)
    if unmapped:
        raise BuildError(
            f"classes with no group: {unmapped}; add them to data/benchmark/clade_groups.tsv"
        )

    folds, pooled = [], {}
    for clade, n in sorted(by_class.items()):
        if n["positive"] >= min_pos and n["hard_negative"] >= min_neg:
            folds.append(
                {
                    "fold": clade,
                    "level": "class",
                    "members": clade,
                    "n_positive": n["positive"],
                    "n_hard_negative": n["hard_negative"],
                    "scored": "yes",
                }
            )
        else:
            g = pooled.setdefault(groups[clade], {"positive": 0, "hard_negative": 0, "members": []})
            g["positive"] += n["positive"]
            g["hard_negative"] += n["hard_negative"]
            g["members"].append(clade)

    for group, n in sorted(pooled.items()):
        ok = n["positive"] >= min_pos and n["hard_negative"] >= min_neg
        folds.append(
            {
                "fold": group,
                "level": "group",
                "members": ",".join(sorted(n["members"])),
                "n_positive": n["positive"],
                "n_hard_negative": n["hard_negative"],
                "scored": "yes" if ok else "no",
            }
        )
    folds.sort(key=lambda f: (-f["n_positive"], f["fold"]))
    return folds


def apply_clades(rows, clade_map):
    """Fill clade_class from the map; exclude a scored row whose taxon is unmapped.

    A row with no class cannot take part in a leave-one-clade-out split, and guessing its class
    would put it in the wrong fold, so it is excluded with a reason instead.
    """
    for row in rows:
        if not row["clade_class"]:
            row["clade_class"] = clade_map.get(row["taxon_id"], "")
        if row["tier"] != "excluded" and not row["clade_class"]:
            row["tier"] = "excluded"
            row["exclude_reason"] = (
                f"taxon {row['taxon_id'] or 'none'} is not in the clade map; "
                "add it to data/benchmark/clade_map.tsv"
            )
    return rows


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def yes(value):
    return str(value).strip().lower() in {"yes", "y", "true", "1"}


def looks_secreted(signal_peptide, gpi_anchor):
    """A protein an architecture method could predict from: it enters the secretory pathway."""
    return yes(signal_peptide) or yes(gpi_anchor)


def check_source_schema(path, rows):
    """Every source file must carry the agreed columns; a renamed column is a silent data loss."""
    if not rows:
        raise BuildError(f"{path}: no data rows")
    missing = [c for c in REQUIRED_SOURCE_COLUMNS if c not in rows[0]]
    if missing:
        raise BuildError(f"{path}: missing columns {missing}")


def classify_source_row(row, path):
    """(tier, indirect, exclude_reason) for one literature-curation row."""
    label = (row["label"] or "").strip()
    if label not in LABELS:
        return "excluded", "", f"label {label!r} is not one of {sorted(LABELS)}"
    if not (row["accession"] or "").strip():
        return "excluded", "", "no accession"
    if not (row["pmids"] or "").strip():
        return "excluded", "", "no PMID"

    codes = {c.strip() for c in (row["evidence_code"] or "").split(";") if c.strip()}
    indirect = "yes" if "indirect" in (row["notes"] or "").lower() else ""

    if label == "hard_negative":
        if not codes & NEGATIVE_EVIDENCE:
            return "excluded", indirect, f"hard negative without {sorted(NEGATIVE_EVIDENCE)}"
        return "hard_negative", indirect, ""

    if not codes & DIRECT_EVIDENCE:
        return "excluded", "", f"adhesin without direct evidence (has {sorted(codes)})"
    if (row["mechanism"] or "").strip() == "moonlighting":
        return "positive_moonlighting", "", ""
    return "positive", "", ""


def from_sources(source_dir):
    """Rows of every `*.tsv` in the sources directory, classified into tiers."""
    out = []
    paths = sorted(Path(source_dir).glob("*.tsv"))
    if not paths:
        raise BuildError(f"{source_dir}: no source *.tsv files")
    for path in paths:
        rows = read_tsv(path)
        check_source_schema(path, rows)
        for row in rows:
            tier, indirect, reason = classify_source_row(row, path)
            out.append(
                {
                    "accession": row["accession"].strip(),
                    "accession_db": row["accession_db"].strip(),
                    "gene": row["gene"].strip(),
                    "protein_name": row["protein_name"].strip(),
                    "organism": row["organism"].strip(),
                    "taxon_id": row["taxon_id"].strip(),
                    "clade_class": row["clade_class"].strip(),
                    "tier": tier,
                    "mechanism": row["mechanism"].strip(),
                    "evidence_code": row["evidence_code"].strip(),
                    "indirect": indirect,
                    "assay_summary": row["assay_summary"].strip(),
                    "pmids": row["pmids"].strip(),
                    "length": row["uniprot_length"].strip(),
                    "signal_peptide": "",
                    "gpi_anchor": "",
                    "source": f"curation:{path.stem}",
                    "exclude_reason": reason,
                    "notes": row["notes"].strip(),
                }
            )
    return out


def classify_curated_row(row):
    """(tier, indirect, exclude_reason) for one row of the existing curated adhesins table."""
    cls = (row.get("cls") or "").strip()
    level = (row.get("evidence_level") or "").strip()
    sp, gpi = row.get("signal_peptide", ""), row.get("gpi_anchor", "")

    if cls == "adhesin":
        if level not in CURATED_POSITIVE_LEVELS:
            return "excluded", "", f"adhesin evidence level {level or 'none'} is not direct"
        if yes(row.get("moonlighting")):
            return "positive_moonlighting", "", ""
        if not looks_secreted(sp, gpi):
            # No signal peptide and no anchor: either a moonlighting protein or a broken gene model.
            return "positive_moonlighting", "", ""
        return "positive", "", ""

    if cls in CURATED_NEGATIVE_CLASSES:
        indirect = "yes" if cls == INDIRECT_CLASS else ""
        if not looks_secreted(sp, gpi):
            return (
                "excluded",
                indirect,
                "negative candidate is not secreted (no signal peptide, no GPI anchor)",
            )
        return "hard_negative", indirect, ""

    return "excluded", "", f"class {cls or 'none'} is not used by the benchmark"


def from_curated(path):
    out = []
    for row in read_tsv(path):
        tier, indirect, reason = classify_curated_row(row)
        out.append(
            {
                "accession": (row.get("accession") or "").strip(),
                "accession_db": "UniProtKB",
                "gene": (row.get("gene") or "").strip(),
                "protein_name": (row.get("protein_name") or "").strip(),
                "organism": (row.get("organism") or "").strip(),
                "taxon_id": (row.get("taxon_id") or "").strip(),
                "clade_class": "",
                "tier": tier,
                "mechanism": (row.get("family") or "").strip(),
                "evidence_code": (row.get("evidence_codes") or "").strip(),
                "indirect": indirect,
                "assay_summary": (row.get("evidence_summary") or "").strip(),
                "pmids": (row.get("pmids") or "").strip(),
                "length": (row.get("length") or "").strip(),
                "signal_peptide": (row.get("signal_peptide") or "").strip(),
                "gpi_anchor": (row.get("gpi_anchor") or "").strip(),
                "source": "curated:adhesins.tsv",
                "exclude_reason": reason,
                "notes": "",
            }
        )
    return out


def merge(rows):
    """One row per accession. A later source wins only when the earlier row was excluded.

    Raises on a real conflict: the same accession in two tiers from two sources is a curation error
    that must be seen, not resolved silently.
    """
    by_acc = {}
    for row in rows:
        acc = row["accession"]
        if not acc:
            continue
        seen = by_acc.get(acc)
        if seen is None:
            by_acc[acc] = row
            continue
        if seen["tier"] == row["tier"]:
            continue
        if seen["tier"] == "excluded":
            by_acc[acc] = row
        elif row["tier"] != "excluded":
            raise BuildError(
                f"{acc}: tier conflict between {seen['source']} ({seen['tier']}) "
                f"and {row['source']} ({row['tier']})"
            )
    return [by_acc[a] for a in sorted(by_acc)]


def summarize(rows):
    counts = Counter(r["tier"] for r in rows)
    per_clade = Counter(
        (r["clade_class"] or "unassigned", r["tier"])
        for r in rows
        if r["tier"] in ("positive", "hard_negative")
    )
    return counts, per_clade


def write_outputs(out_dir, version, rows, counts, per_clade, inputs, folds):
    out_dir.mkdir(parents=True, exist_ok=True)
    proteins = out_dir / "proteins.tsv"
    with open(proteins, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(OUTPUT_COLUMNS), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    (out_dir / "VERSION").write_text(
        f"{version}\n"
        + "".join(f"input {p}\tsha256 {h}\n" for p, h in inputs)
        + "".join(f"tier {t}\t{counts.get(t, 0)}\n" for t in TIERS)
    )

    lines = ["\t".join(("clade_class", "positive", "hard_negative"))]
    clades = sorted({c for c, _ in per_clade})
    for clade in clades:
        lines.append(
            f"{clade}\t{per_clade.get((clade, 'positive'), 0)}"
            f"\t{per_clade.get((clade, 'hard_negative'), 0)}"
        )
    (out_dir / "counts.tsv").write_text("\n".join(lines) + "\n")

    fold_cols = ("fold", "level", "members", "n_positive", "n_hard_negative", "scored")
    with open(out_dir / "folds.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fold_cols), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(folds)

    sums = []
    for path in sorted(out_dir.iterdir()):
        if path.name != "SHA256SUMS" and path.is_file():
            sums.append(f"{sha256_file(path)}  {path.name}")
    (out_dir / "SHA256SUMS").write_text("\n".join(sums) + "\n")
    return proteins


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--version", required=True, help="version name, e.g. adhesion_v1")
    ap.add_argument("--out", required=True, help="output directory (must not exist unless --force)")
    ap.add_argument("--curated", default=str(ROOT / "data/curated/adhesins/adhesins.tsv"))
    ap.add_argument(
        "--sources", default=None, help="directory of curation *.tsv (default <out>/sources)"
    )
    ap.add_argument("--clade-map", default=str(ROOT / "data/benchmark/clade_map.tsv"))
    ap.add_argument("--clade-groups", default=str(ROOT / "data/benchmark/clade_groups.tsv"))
    ap.add_argument("--dry-run", action="store_true", help="print the summary, write nothing")
    ap.add_argument("--force", action="store_true", help="overwrite an existing output directory")
    a = ap.parse_args(argv)

    try:
        out_dir = Path(a.out)
        sources = Path(a.sources) if a.sources else out_dir / "sources"
        if out_dir.exists() and not (a.force or a.dry_run):
            raise BuildError(f"{out_dir} exists; pass --force to overwrite or use a new version")

        clade_map = read_clade_map(a.clade_map)
        groups = read_clade_groups(a.clade_groups)
        rows = apply_clades(merge(from_curated(a.curated) + from_sources(sources)), clade_map)
        folds = plan_folds(rows, groups)
        counts, per_clade = summarize(rows)

        if not counts.get("positive"):
            raise BuildError(
                "no positives survived the rules; refusing to write an empty benchmark"
            )

        inputs = [
            (a.curated, sha256_file(a.curated)),
            (a.clade_map, sha256_file(a.clade_map)),
            (a.clade_groups, sha256_file(a.clade_groups)),
        ] + [(str(p), sha256_file(p)) for p in sorted(sources.glob("*.tsv"))]

        for tier in TIERS:
            print(f"{tier}: {counts.get(tier, 0)}")
        print("\nclade\tpositive\thard_negative")
        for clade in sorted({c for c, _ in per_clade}):
            print(
                f"{clade}\t{per_clade.get((clade, 'positive'), 0)}"
                f"\t{per_clade.get((clade, 'hard_negative'), 0)}"
            )

        print("\nfold\tlevel\tpositive\thard_negative\tscored")
        for f in folds:
            print(
                f"{f['fold']}\t{f['level']}\t{f['n_positive']}"
                f"\t{f['n_hard_negative']}\t{f['scored']}"
            )

        if a.dry_run:
            print("\ndry run: nothing written")
            return 0

        proteins = write_outputs(out_dir, a.version, rows, counts, per_clade, inputs, folds)
        print(
            f"\nwrote {proteins} ({len(rows)} rows) and VERSION, counts.tsv, folds.tsv, SHA256SUMS"
        )
    except (BuildError, OSError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
