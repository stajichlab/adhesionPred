#!/usr/bin/env python3
"""Compare the old and the new rule of cell_wall_adhesion_candidate on one finished run.

Reads only (never writes to) an `out_final` directory of the classifier:
`calls.long.tsv.gz` and/or `calls.wide.tsv.gz`, and `evidence.tsv.gz` for the family split.
Standard library only.

    python compare_old_new_rule.py --out counts.tsv OUT_FINAL_DIR [OUT_FINAL_DIR ...]

The label of a proteome is the name of the work directory (the parent of `out_final`), or
`--label` when one directory is given.

What it counts, per proteome and per signal-peptide variant:
  new      the narrow `cell_wall_adhesion_candidate` as the run wrote it.
  old      the rule before D4: signal_peptide_protein AND (tandem_repeat_protein OR
           wall_family_domain), read from the same run. Kleene logic as in logic.py:
           not_applicable acts as not_called inside and/or; not_assessable stays unknown.
  old_only proteins the old rule calls and the new rule does not call.
  new_only proteins the new rule calls and the old rule does not call (expected 0).
  surface  `surface_attachment_candidate`, for reference.
  rule_mismatch  proteins where the written narrow call differs from the narrow rule recomputed
           from tandem_repeat_protein, wall_family_narrow_domain and signal_peptide_protein
           (expected 0; a non-zero value means the run is not what this script assumes).
Values are counted in four separate columns: called, not_called, not_assessable and
not_applicable. not_assessable and not_applicable are never counted as not_called.

The family split (`<out>.by_family.tsv`) takes the old_only proteins and counts them by the
accession in the pfam_adhesion `families` evidence field. A protein with several families is
counted once in each; `n_proteins` is the number of distinct proteins. It is skipped, with a
message, when evidence.tsv.gz is missing. The evidence table has the families of the broad
flag only, so it shows through which families the broad flag fired, not which one the rule used.

This script makes no call and measures no sensitivity or specificity. It only counts rows of a run.
"""

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

CALLED, NOT_CALLED, NOT_ASSESSABLE, NOT_APPLICABLE = (
    "called",
    "not_called",
    "not_assessable",
    "not_applicable",
)
VALUES = (CALLED, NOT_CALLED, NOT_ASSESSABLE, NOT_APPLICABLE)
LONG_NEEDED = ["protein", "call", "variant", "value"]
LEAF_CALLS = [
    "signal_peptide_protein",
    "tandem_repeat_protein",
    "wall_family_domain",
    "wall_family_narrow_domain",
]
OUT_CALLS = ["cell_wall_adhesion_candidate", "surface_attachment_candidate"]
EVIDENCE_NEEDED = ["protein", "module", "field", "value"]


class InputError(Exception):
    pass


def _and(*values):
    values = [NOT_CALLED if v == NOT_APPLICABLE else v for v in values]
    if NOT_CALLED in values:
        return NOT_CALLED
    if NOT_ASSESSABLE in values:
        return NOT_ASSESSABLE
    return CALLED


def _or(*values):
    values = [NOT_CALLED if v == NOT_APPLICABLE else v for v in values]
    if CALLED in values:
        return CALLED
    if NOT_ASSESSABLE in values:
        return NOT_ASSESSABLE
    return NOT_CALLED


def _open(path):
    return gzip.open(path, "rt", newline="")


def _reader(path, needed):
    fh = _open(path)
    reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
    missing = [c for c in needed if c not in (reader.fieldnames or [])]
    if missing:
        fh.close()
        raise InputError(f"{path}: missing column(s) {', '.join(missing)}")
    return fh, reader


def _check_value(path, line, value):
    if value not in VALUES:
        raise InputError(f"{path}: row {line}: not a call value: {value!r}")


def read_long(path):
    """Return (proteins in file order, {(protein, call, variant): value})."""
    proteins, table = [], {}
    fh, reader = _reader(path, LONG_NEEDED)
    with fh:
        seen = set()
        for n, row in enumerate(reader, start=2):
            _check_value(path, n, row["value"])
            key = (row["protein"], row["call"], row["variant"])
            if key in table:
                raise InputError(f"{path}: row {n}: duplicate row for {key}")
            table[key] = row["value"]
            if row["protein"] not in seen:
                seen.add(row["protein"])
                proteins.append(row["protein"])
    return proteins, table


def read_wide(path):
    """Same result as read_long, from the wide table (columns `call` or `call[variant]`)."""
    fh, reader = _reader(path, ["protein"])
    table, proteins = {}, []
    with fh:
        names = []
        for col in reader.fieldnames:
            if col == "protein" or col.endswith("_status") or col.endswith("_basis"):
                continue
            if col in ("primary_label", "primary_basis"):
                continue
            call, _, rest = col.partition("[")
            names.append((col, call, rest.rstrip("]")))
        for n, row in enumerate(reader, start=2):
            proteins.append(row["protein"])
            for col, call, variant in names:
                if row[col] == "":
                    continue
                _check_value(path, n, row[col])
                table[(row["protein"], call, variant)] = row[col]
    return proteins, table


def load_calls(out_dir):
    long_path = out_dir / "calls.long.tsv.gz"
    wide_path = out_dir / "calls.wide.tsv.gz"
    if long_path.is_file():
        result = read_long(long_path)
        if wide_path.is_file():
            wide = read_wide(wide_path)
            if wide[1] != result[1]:
                raise InputError(f"{out_dir}: calls.long and calls.wide disagree")
        return result
    if wide_path.is_file():
        return read_wide(wide_path)
    raise InputError(f"{out_dir}: neither calls.long.tsv.gz nor calls.wide.tsv.gz found")


def families_by_protein(out_dir):
    path = out_dir / "evidence.tsv.gz"
    if not path.is_file():
        return None
    fh, reader = _reader(path, EVIDENCE_NEEDED)
    found = {}
    with fh:
        for row in reader:
            if row["module"] == "pfam_adhesion" and row["field"] == "families":
                found[row["protein"]] = [a for a in row["value"].split(",") if a]
    return found


def _get(table, protein, call, variant):
    """A leaf call has variant '' unless the call is per variant; None if the row is missing."""
    for v in (variant, ""):
        if (protein, call, v) in table:
            return table[(protein, call, v)]
    return None


def compare(out_dir):
    """Return (rows of the count table, rows of the family split or None)."""
    proteins, table = load_calls(out_dir)
    present = {call for (_, call, _) in table}
    absent = [c for c in LEAF_CALLS + OUT_CALLS if c not in present]
    if absent:
        raise InputError(f"{out_dir}: missing call(s) {', '.join(absent)}")
    variants = sorted({v for (_, c, v) in table if c == "cell_wall_adhesion_candidate"})
    families = families_by_protein(out_dir)
    count_rows, family_rows = [], defaultdict(set)
    for variant in variants:
        c = {k: Counter() for k in ("new", "old", "surface")}
        n_old_only = n_new_only = mismatch = 0
        for p in proteins:
            leaf = {}
            for call in LEAF_CALLS:
                leaf[call] = _get(
                    table, p, call, variant if call == "signal_peptide_protein" else ""
                )
            new = _get(table, p, "cell_wall_adhesion_candidate", variant)
            surface = _get(table, p, "surface_attachment_candidate", variant)
            if new is None or surface is None or None in leaf.values():
                raise InputError(
                    f"{out_dir}: protein {p}, variant {variant!r}: a needed row is missing"
                )
            old = _and(
                _or(leaf["tandem_repeat_protein"], leaf["wall_family_domain"]),
                leaf["signal_peptide_protein"],
            )
            recomputed = _and(
                _or(leaf["tandem_repeat_protein"], leaf["wall_family_narrow_domain"]),
                leaf["signal_peptide_protein"],
            )
            if recomputed != new and NOT_APPLICABLE not in (new, recomputed):
                mismatch += 1
            c["new"][new] += 1
            c["old"][old] += 1
            c["surface"][surface] += 1
            if old == CALLED and new != CALLED:
                n_old_only += 1
                if families is not None:
                    for acc in families.get(p, ["no_families_row"]) or ["no_families_row"]:
                        family_rows[(variant, acc)].add(p)
            if new == CALLED and old != CALLED:
                n_new_only += 1
        row = {"proteome": out_dir.parent.name, "variant": variant, "n_proteins": len(proteins)}
        for name in ("new", "old", "surface"):
            for value in VALUES:
                row[f"{name}_{value}"] = c[name][value]
        row.update(old_only=n_old_only, new_only=n_new_only, rule_mismatch=mismatch)
        count_rows.append(row)
    fam = None
    if families is not None:
        fam = [
            {"proteome": out_dir.parent.name, "variant": v, "family": a, "n_proteins": len(ps)}
            for (v, a), ps in sorted(family_rows.items())
        ]
    return count_rows, fam


HEADER = [
    "#compare_old_new_rule.py: counts rows of finished runs; it makes no call and measures nothing.",
    "#old = signal_peptide_protein AND (tandem_repeat_protein OR wall_family_domain); "
    "new = the narrow cell_wall_adhesion_candidate as written. Same run, no re-run.",
    "#not_assessable and not_applicable are counted apart from not_called.",
]


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as fh:
        for line in HEADER:
            fh.write(line + "\n")
        writer = csv.DictWriter(fh, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def print_table(rows, columns):
    widths = {c: max(len(c), *(len(str(r[c])) for r in rows)) for c in columns}
    print("  ".join(c.ljust(widths[c]) for c in columns))
    for r in rows:
        print("  ".join(str(r[c]).ljust(widths[c]) for c in columns))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out_final", nargs="+", type=Path, help="out_final directories (read only)")
    ap.add_argument("--out", required=True, type=Path, help="TSV to write")
    ap.add_argument("--label", help="proteome label (one directory only)")
    args = ap.parse_args(argv)
    if args.label and len(args.out_final) != 1:
        ap.error("--label needs exactly one directory")
    counts, fams = [], []
    try:
        for d in args.out_final:
            if not d.is_dir():
                raise InputError(f"{d}: not a directory")
            rows, fam = compare(d)
            if args.label:
                for r in rows:
                    r["proteome"] = args.label
                for r in fam or []:
                    r["proteome"] = args.label
            counts += rows
            if fam is None:
                print(f"{d}: no evidence.tsv.gz, family split skipped", file=sys.stderr)
            else:
                fams += fam
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    columns = ["proteome", "variant", "n_proteins"]
    columns += [f"{n}_{v}" for n in ("new", "old", "surface") for v in VALUES]
    columns += ["old_only", "new_only", "rule_mismatch"]
    write_tsv(args.out, counts, columns)
    print_table(counts, columns)
    if fams:
        fam_path = args.out.with_name(args.out.name + ".by_family.tsv")
        write_tsv(fam_path, fams, ["proteome", "variant", "family", "n_proteins"])
        print(f"\nold_only proteins by family of the broad flag: {fam_path}")
    if any(r["rule_mismatch"] for r in counts):
        print("warning: written narrow call differs from the recomputed rule", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
