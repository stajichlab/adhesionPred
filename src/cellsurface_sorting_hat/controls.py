"""Check a committed control set (L1, D15) against ``docs/agent-tasks/COMMON-RULES.md``.

A control set is a directory with ``controls.tsv``, ``sequences.faa``, ``clusters.tsv`` and
``counts.md``. A defect is a broken rule: it fails the set. The size floor of rule 7 is only
reported, because a small set is not a defect. The taxon name check needs ``names.dmp`` and is
skipped when the file is not given.

    python -m cellsurface_sorting_hat.controls data/controls/<dir> [--names-dmp names.dmp]
"""

import argparse
import csv
import datetime
import gzip
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

from cellsurface_sorting_hat.status import LEAKAGE_VALUES

REQUIRED_COLUMNS = (
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
)
LABELS = ("positive", "negative")
FLOOR = 20  # rule 7: rows and clusters per class
MIXED_COLUMN = "mixed_cluster_with_positive"
_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
_COUNT_KEY = re.compile(
    r"^(rows|sequences|independent clusters|clusters)"
    r"(?: \([^)]*\))?(?: in sequences\.faa)?(?: (positive|negative))?$"
)


@dataclass(frozen=True)
class Defect:
    file: str
    row: str  # accession, or "" for a defect of the whole file
    rule: str
    message: str

    def __str__(self):
        where = f"row {self.row}: " if self.row else ""
        return f"{self.file}: {where}rule {self.rule}: {self.message}"


@dataclass
class Report:
    path: str
    defects: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    size: dict = field(default_factory=dict)
    counts: dict = field(default_factory=dict)

    @property
    def ok(self):
        return not self.defects

    def add(self, file, row, rule, message):
        self.defects.append(Defect(file, row, rule, message))


def _open(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return open(path, encoding="utf-8-sig", newline="")


def _find(directory, name):
    for candidate in (directory / name, directory / f"{name}.gz"):
        if candidate.is_file():
            return candidate
    return None


def _read_controls(path, report):
    with _open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = [dict(r) for r in reader]
    missing = [c for c in REQUIRED_COLUMNS if c not in header]
    if missing:
        report.add(path.name, "", "columns", f"missing required column(s): {', '.join(missing)}")
    return header, rows


def _read_fasta_ids(path):
    ids = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                parts = line[1:].split()
                ids.append(parts[0] if parts else "")
    return ids


def _read_clusters(path):
    """``accession -> cluster_id``; a header line is skipped; a first duplicate wins."""
    mapping = {}
    with _open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or (f[0], f[1]) == ("accession", "cluster_id"):
                continue
            mapping.setdefault(f[0], f[1])
    return mapping


def _check_rows(name, rows, report):
    seen = set()
    for r in rows:
        acc = (r.get("accession") or "").strip()
        if not acc:
            report.add(name, "", "accession", "a row has an empty accession")
            continue
        if acc in seen:
            report.add(name, acc, "duplicate", "accession occurs more than once (one row each)")
        seen.add(acc)

        def val(col, r=r):
            return (r.get(col) or "").strip() if col in r else None

        if val("label") is not None and val("label") not in LABELS:
            report.add(name, acc, "label", f"label {val('label')!r} is not positive or negative")
        if val("tuned") is not None and val("tuned") not in ("yes", "no"):
            report.add(name, acc, "tuned", f"tuned {val('tuned')!r} is not yes or no")
        if val("leakage") is not None and val("leakage") not in LEAKAGE_VALUES:
            report.add(
                name,
                acc,
                "leakage (rule 6)",
                f"leakage {val('leakage')!r} is not one of {', '.join(LEAKAGE_VALUES)}",
            )
        for col, rule in (("source", "source (rule 1)"), ("quote", "quote (rule 1)")):
            if val(col) == "":
                report.add(name, acc, rule, f"{col} is empty; every row needs a source and a quote")
        if val("species") == "":
            report.add(name, acc, "species (rule 5)", "species is empty")
        if val("taxon_id") is not None and not re.fullmatch(r"\d+", val("taxon_id")):
            report.add(
                name,
                acc,
                "taxon_id (rule 5)",
                f"taxon_id {val('taxon_id')!r} is not an integer",
            )
        d = val("date")
        if d is not None:
            try:
                if not _DATE.fullmatch(d):
                    raise ValueError(d)
                datetime.date.fromisoformat(d)
            except ValueError:
                report.add(name, acc, "date", f"date {d!r} is not a valid YYYY-MM-DD date")
    return seen


def _read_names(path):
    names = {}
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            f = [x.strip() for x in line.rstrip("\n").split("|")]
            if len(f) >= 4 and f[3] == "scientific name" and f[0].isdigit():
                names[int(f[0])] = f[1]
    return names


def _check_names(name, rows, names_dmp, report):
    names = _read_names(names_dmp)
    for r in rows:
        tid, species = (r.get("taxon_id") or "").strip(), (r.get("species") or "").strip()
        if not tid.isdigit() or not species:
            continue
        acc = r.get("accession", "")
        if int(tid) not in names:
            report.add(name, acc, "names.dmp (rule 5)", f"taxon_id {tid} is not in names.dmp")
        elif names[int(tid)] != species:
            report.add(
                name,
                acc,
                "names.dmp (rule 5)",
                f"species {species!r} is not the scientific name {names[int(tid)]!r} of taxon {tid}",
            )


def _computed_counts(rows, clusters):
    """Counts from the files: rows, clusters, and the same per class."""
    out = {"rows": len(rows)}
    cl = {r["accession"]: clusters.get(r["accession"]) for r in rows if r.get("accession")}
    out["clusters"] = len({c for c in cl.values() if c})
    for label in LABELS:
        accs = [r["accession"] for r in rows if r.get("label") == label and r.get("accession")]
        out[f"rows {label}"] = len(accs)
        out[f"clusters {label}"] = len({cl[a] for a in accs if cl.get(a)})
    return out


def _parse_counts_md(path):
    """``(kind, class) -> (value, line number)`` from table rows such as ``| Rows (proteins) | 4 |``."""
    found = {}
    with open(path, encoding="utf-8-sig") as fh:
        for n, line in enumerate(fh, 1):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if len(cells) < 2:
                continue
            key = re.sub(r"[`*]", "", cells[0]).strip().lower()
            m = _COUNT_KEY.match(key)
            if not m or not re.fullmatch(r"\d+", cells[1]):
                continue
            kind = "clusters" if "clusters" in m.group(1) else m.group(1)
            found.setdefault((kind, m.group(2) or ""), (int(cells[1]), n))
    return found


def _check_counts(path, computed, n_sequences, report):
    found = _parse_counts_md(path)
    if not found:
        report.notes.append(
            f"{path.name}: no count lines in the checkable form; counts not checked"
        )
        return
    expected = {
        ("rows", ""): computed["rows"],
        ("sequences", ""): n_sequences,
        ("clusters", ""): computed["clusters"],
    }
    for label in LABELS:
        expected[("rows", label)] = computed[f"rows {label}"]
        expected[("clusters", label)] = computed[f"clusters {label}"]
        expected[("sequences", label)] = computed[f"rows {label}"]  # one sequence per row
    for key, (value, line) in found.items():
        want = expected[key]
        if value != want:
            what = f"{key[0]} {key[1]}".strip()
            report.add(
                path.name,
                "",
                "counts (rule 4)",
                f"line {line}: {what} is {value} in counts.md but {want} from the files",
            )


def _size_floor(rows, clusters, report):
    size = {}
    for label in LABELS:
        n_rows = sum(1 for r in rows if r.get("label") == label)
        n_clu = len(
            {clusters.get(r.get("accession")) for r in rows if r.get("label") == label} - {None}
        )
        size[f"{label}_rows"] = n_rows
        size[f"{label}_clusters"] = n_clu
        size[f"{label}_rows_to_floor"] = max(0, FLOOR - n_rows)
        size[f"{label}_clusters_to_floor"] = max(0, FLOOR - n_clu)
        report.notes.append(
            f"size floor (rule 7), {label}: {n_rows} rows ({max(0, FLOOR - n_rows)} to {FLOOR}), "
            f"{n_clu} clusters ({max(0, FLOOR - n_clu)} to {FLOOR})"
        )
    report.size = size


def validate_controls(path, names_dmp=None):
    """Check one control-set directory; return a ``Report`` (``report.ok`` is False on a defect)."""
    directory = Path(path)
    report = Report(str(directory))
    if not directory.is_dir():
        report.add(str(directory), "", "file", "the control-set directory does not exist")
        return report
    files = {n: _find(directory, n) for n in ("controls.tsv", "sequences.faa", "clusters.tsv")}
    counts_path = directory / "counts.md"
    for name, found in {
        **files,
        "counts.md": counts_path if counts_path.is_file() else None,
    }.items():
        if found is None:
            report.add(name, "", "file", f"required file is missing in {directory}")
    if files["controls.tsv"] is None:
        return report
    header, rows = _read_controls(files["controls.tsv"], report)
    accessions = _check_rows(files["controls.tsv"].name, rows, report)
    if names_dmp is not None:
        _check_names(files["controls.tsv"].name, rows, names_dmp, report)

    clusters = {}
    if files["clusters.tsv"] is not None:
        clusters = _read_clusters(files["clusters.tsv"])
        for acc in sorted(accessions - set(clusters)):
            report.add("clusters.tsv", acc, "cluster (rule 4)", "accession is not in clusters.tsv")
        _check_mixed(rows, clusters, header, report)
    n_sequences = 0
    if files["sequences.faa"] is not None:
        ids = _read_fasta_ids(files["sequences.faa"])
        n_sequences = len(set(ids))
        for acc in sorted(accessions - set(ids)):
            report.add("sequences.faa", acc, "sequence", "accession has no sequence")
    computed = _computed_counts(rows, clusters)
    report.counts = {**computed, "sequences": n_sequences}
    if counts_path.is_file():
        _check_counts(counts_path, computed, n_sequences, report)
    _size_floor(rows, clusters, report)
    return report


def _check_mixed(rows, clusters, header, report):
    """A cluster with a positive and a negative needs ``mixed_cluster_with_positive=yes`` on the negative."""
    positives = {
        clusters.get(r.get("accession"))
        for r in rows
        if r.get("label") == "positive" and clusters.get(r.get("accession"))
    }
    for r in rows:
        acc = r.get("accession", "")
        if r.get("label") != "negative" or clusters.get(acc) not in positives:
            continue
        flag = (r.get(MIXED_COLUMN) or "").strip().lower() if MIXED_COLUMN in header else ""
        if flag not in ("yes", "true", "1"):
            report.add(
                "clusters.tsv",
                acc,
                "mixed cluster (rule 4)",
                f"cluster {clusters[acc]} holds a positive and this negative; "
                f"the row must say {MIXED_COLUMN}=yes",
            )


def build_parser():
    ap = argparse.ArgumentParser(
        prog="python -m cellsurface_sorting_hat.controls",
        description="Check control-set directories made as docs/agent-tasks/COMMON-RULES.md says.",
    )
    ap.add_argument("sets", nargs="+", help="directory with controls.tsv, sequences.faa, ...")
    ap.add_argument(
        "--names-dmp", help="NCBI names.dmp; without it the species names are not checked"
    )
    return ap


def run(sets, names_dmp=None, out=None):
    """Print one report per set; return 1 when any set has a defect, else 0."""
    out = out or sys.stdout
    status = 0
    for s in sets:
        rep = validate_controls(s, names_dmp=names_dmp)
        print(f"{'PASS' if rep.ok else 'FAIL'}\t{rep.path}\t{len(rep.defects)} defect(s)", file=out)
        for d in rep.defects:
            print(f"  {d}", file=out)
        for n in rep.notes:
            print(f"  note: {n}", file=out)
        status = status or (0 if rep.ok else 1)
    return status


def main(argv=None):
    args = build_parser().parse_args(argv)
    return run(args.sets, names_dmp=args.names_dmp)


if __name__ == "__main__":
    sys.exit(main())
