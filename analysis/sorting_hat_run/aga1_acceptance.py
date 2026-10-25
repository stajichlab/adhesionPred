"""Acceptance run of the AGA1 modules and calls: build module tables for each proteome from the tool outputs already in hand, run the engine, and
compare its calls with the flags of `analysis/aga1_class/` (stdlib only; local, no SLURM, no sequence is written inside the repository).

For every --set (the staged chunks of one proteome set, its proteins.tsv, the merged raw tool outputs, the detector-14 table, the classifier score tables
and the committed per-protein flags table) the script splits the set by species, runs the wrappers `signalp`, `gpi`, `disorder`, `repeat14` and
`aga1_model` for each species, runs `cellsurface_sorting_hat` (modules without a table stay absent and give not_assessable), and compares, protein by protein:
  ser_thr_gpi_wall_candidate (default gate) called  <=>  broad_flag == 1
  flocculin_like_candidate   (default gate) called  <=>  narrow_flag == 1
  model_surface_candidate    (default gate) called  <=>  at least 2 of 3 models and a signal peptide (n_models_call >= 2 and sp_yes == 1)
Any difference is a bug in a wrapper or a call: the exit code is 1 and the differences are listed.

  python3 -I analysis/sorting_hat_run/aga1_acceptance.py --taxdump nodes.dmp --out-root NEW_DIR \\
     --set NAME CHUNKS_DIR PROTEINS_TSV MERGED_DIR REP14_TSV SCORES_DIR FLAGS_TSV   [--set ...]

Writes OUT_ROOT/<species>/{p.faa, wd/, out/} and OUT_ROOT/acceptance_counts.tsv. Refuses an existing --out-root.
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
TAXA = Path(__file__).resolve().parent / "aga1_acceptance_taxa.tsv"
MODELS = ("logreg", "rf", "gboost")
SIGNALP_VERSION = "6.0h"
PREDGPI_VERSION = "202001"
AIUPRED_VERSION = "2.x"


class AcceptanceError(Exception):
    pass


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def species_taxa():
    return {r["species"]: int(r["taxon_id"]) for r in read_tsv(TAXA)}


def _records(path):
    name, seq = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def split_fasta(chunks_dir, members):
    """{species: fasta text} for the proteins of `members` ({accession: species}); chunk files have bare-accession headers."""
    out, seen = {}, set()
    for f in sorted(Path(chunks_dir).glob("*.faa")):
        for name, seq in _records(f):
            sp = members.get(name)
            if sp is not None:
                out.setdefault(sp, []).append(f">{name}\n{seq}\n")
                seen.add(name)
    missing = set(members) - seen
    if missing:
        raise AcceptanceError(
            f"{len(missing)} proteins of proteins.tsv are not in {chunks_dir}, first {sorted(missing)[:3]}"
        )
    return {sp: "".join(parts) for sp, parts in out.items()}


def split_signalp(path, members):
    """Header lines go to every species; result rows by the first token of the ID."""
    heads, rows = [], {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            heads.append(line)
        elif line.strip():
            sp = members.get(line.split("\t")[0].split()[0])
            if sp is not None:
                rows.setdefault(sp, []).append(line)
    return {sp: "\n".join(heads + r) + "\n" for sp, r in rows.items()}


def split_gff3(path, members):
    rows = {}
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            sp = members.get(line.split("\t")[0])
            if sp is not None:
                rows.setdefault(sp, []).append(line)
    return {sp: "\n".join(r) + "\n" for sp, r in rows.items()}


def split_disorder_to_files(path, members, folder):
    """Stream the per-residue file and append each record to the file of its species (never held whole in memory)."""
    handles, cur = {}, None
    try:
        with open(path, encoding="utf-8-sig") as fh:
            for line in fh:
                if line.startswith("#>") or line.startswith(">"):
                    name = line.lstrip("#>").split()[0]
                    sp = members.get(name)
                    cur = None
                    if sp is not None:
                        if sp not in handles:
                            handles[sp] = open(Path(folder) / f"{sp}.disorder_raw.txt", "w")
                        cur = handles[sp]
                if cur is not None:
                    cur.write(line)
    finally:
        for h in handles.values():
            h.close()
    return {sp: Path(folder) / f"{sp}.disorder_raw.txt" for sp in handles}


def split_table(path, key, members):
    """{species: text} for a TSV with a header, split by the value of column `key`."""
    with open(path, newline="") as fh:
        lines = fh.read().splitlines()
    head = lines[0]
    col = head.split("\t").index(key)
    rows = {}
    for line in lines[1:]:
        sp = members.get(line.split("\t")[col])
        if sp is not None:
            rows.setdefault(sp, []).append(line)
    return {sp: "\n".join([head, *r]) + "\n" for sp, r in rows.items()}


def read_calls(
    out_dir,
    calls=("ser_thr_gpi_wall_candidate", "flocculin_like_candidate", "model_surface_candidate"),
):
    """{call: {protein: value}} at the default gate (the R0 variant for the per_variant calls)."""
    res = {c: {} for c in calls}
    with gzip.open(Path(out_dir) / "calls.long.tsv.gz", "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["call"] in res and r["variant"] in ("", "R0"):
                res[r["call"]][r["protein"]] = r["value"]
    return res


def compare(calls, flags):
    """List of differences between the engine's calls and the flags rows ({accession: row of the flags table})."""
    diffs = []
    for acc, f in flags.items():
        sp = f.get("sp_yes", "0") in ("1", "1.0")
        expect = {
            "ser_thr_gpi_wall_candidate": f["broad_flag"] == "1",
            "flocculin_like_candidate": f["narrow_flag"] == "1",
            "model_surface_candidate": int(f["n_models_call"]) >= 2 and sp,
        }
        for call, want in expect.items():
            got = calls[call].get(acc)
            if got is None:
                diffs.append((acc, call, "missing from the engine output"))
            elif (got == "called") != want:
                diffs.append(
                    (
                        acc,
                        call,
                        f"engine {got}, analysis flags say {'called' if want else 'not called'}",
                    )
                )
    return diffs


def run_species(sp, taxon, fasta_text, files, nodes, out_root, thresholds):
    from cellsurface_sorting_hat.cli import main as engine
    from cellsurface_sorting_hat.modules.cli import main as module

    base = Path(out_root) / sp
    base.mkdir(parents=True)
    fa = base / "p.faa"
    fa.write_text(fasta_text)
    wd, out = base / "wd", base / "out"
    common = ["--fasta", str(fa), "--workdir", str(wd)]
    steps = [
        ["signalp", *common, "--results", str(files["signalp"]), "--signalp-version", SIGNALP_VERSION],
        ["gpi", *common, "--gff3", str(files["gpi"]), "--predgpi-version", PREDGPI_VERSION],
        ["disorder", *common, "--raw", str(files["disorder"]), "--tool", "aiupred", "--tool-version", AIUPRED_VERSION],
        ["repeat14", *common, "--table", str(files["rep14"]), "--script", str(ROOT / "analysis/cocci_repeats/14_repeat_detect_general.py")],
        ["aga1_model", *common, "--taxon", str(taxon), "--taxdump", str(nodes), "--thresholds", str(thresholds),
         "--scores", *[f"{m}={files[m]}" for m in MODELS]],
    ]  # fmt: skip
    for argv in steps:
        if module(argv) != 0:
            raise AcceptanceError(f"{sp}: module {argv[0]} failed")
    if (
        engine(
            [
                "--fasta",
                str(fa),
                "--taxon",
                str(taxon),
                "--taxdump",
                str(nodes),
                "--workdir",
                str(wd),
                "--out",
                str(out),
            ]
        )
        != 0
    ):
        raise AcceptanceError(f"{sp}: the engine failed")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--taxdump", required=True, help="nodes.dmp")
    ap.add_argument("--out-root", required=True)
    ap.add_argument(
        "--thresholds", default=str(ROOT / "data/sorting_hat/aga1_model_thresholds.tsv")
    )
    ap.add_argument("--set", dest="sets", action="append", nargs=7, required=True,
                    metavar=("NAME", "CHUNKS_DIR", "PROTEINS_TSV", "MERGED_DIR", "REP14_TSV", "SCORES_DIR", "FLAGS_TSV"))  # fmt: skip
    a = ap.parse_args(argv)
    out_root = Path(a.out_root)
    if out_root.exists():
        print(f"FATAL: {out_root} exists; nothing is overwritten", file=sys.stderr)
        return 1
    out_root.mkdir(parents=True)
    taxa = species_taxa()
    counts, all_diffs = [], []
    try:
        for name, chunks, proteins_tsv, merged, rep14, scores, flags_tsv in a.sets:
            members = {r["accession"]: r["species"] for r in read_tsv(proteins_tsv)}
            unknown = sorted(set(members.values()) - set(taxa))
            if unknown:
                raise AcceptanceError(f"{name}: species without a taxon in {TAXA.name}: {unknown}")
            fastas = split_fasta(chunks, members)
            tmp = out_root / f"_split_{name}"
            tmp.mkdir()
            per = {
                "signalp": split_signalp(Path(merged) / "signalp/prediction_results.txt", members),
                "gpi": split_gff3(Path(merged) / "predgpi/predgpi_raw.gff3", members),
                "rep14": split_table(rep14, "protein", members),
            }
            for m in MODELS:
                per[m] = split_table(Path(scores) / f"scores_{m}.tsv", "accession", members)
            disorder = split_disorder_to_files(
                Path(merged) / "aiupred/disorder_raw.txt", members, tmp
            )
            flags = {r["accession"]: r for r in read_tsv(flags_tsv)}
            for sp in sorted(fastas):
                files = {}
                for key, texts in per.items():
                    if sp not in texts:
                        raise AcceptanceError(f"{name}/{sp}: no rows for {key} after the split")
                    p = tmp / f"{sp}.{key}.txt"
                    p.write_text(texts[sp])
                    files[key] = p
                files["disorder"] = disorder[sp]
                out = run_species(
                    sp, taxa[sp], fastas[sp], files, a.taxdump, out_root, a.thresholds
                )
                calls = read_calls(out)
                sp_flags = {acc: f for acc, f in flags.items() if members.get(acc) == sp}
                diffs = compare(calls, sp_flags)
                all_diffs += [(sp, *d) for d in diffs]
                n = {c: sum(v == "called" for v in calls[c].values()) for c in calls}
                counts.append({"species": sp, "taxon": taxa[sp], "proteins": len(sp_flags), "broad_engine": n["ser_thr_gpi_wall_candidate"],
                               "broad_flags": sum(f["broad_flag"] == "1" for f in sp_flags.values()),
                               "narrow_engine": n["flocculin_like_candidate"], "narrow_flags": sum(f["narrow_flag"] == "1" for f in sp_flags.values()),
                               "model_engine": n["model_surface_candidate"],
                               "model_flags": sum(int(f["n_models_call"]) >= 2 and f["sp_yes"] in ("1", "1.0") for f in sp_flags.values()),
                               "differences": len(diffs)})  # fmt: skip
                print(counts[-1])
    except (AcceptanceError, OSError, KeyError, ValueError) as e:
        print(f"FATAL: {e}", file=sys.stderr)
        return 1
    with open(out_root / "acceptance_counts.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(counts[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(counts)
    for d in all_diffs[:40]:
        print("DIFF", *d)
    print(f"{len(counts)} proteomes, {len(all_diffs)} differences")
    return 1 if all_diffs else 0


if __name__ == "__main__":
    sys.exit(main())
