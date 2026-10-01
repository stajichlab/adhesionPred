#!/usr/bin/env python3
"""Phase B step 4: merge SignalP 6, PredGPI and Ser+Thr into one feature table (D2).

Reads from $STEP1_WORKDIR/phaseb/: unique_sequences.tsv.gz, sequence_members.tsv.gz,
signalp/part_*/prediction_results.txt.gz and output.gff3.gz, predgpi/part_*.tsv.gz (all from
J1). Reads $STEP1_WORKDIR/truth_set_triaged.tsv.gz (03) for the truth columns. Writes to
$STEP1_WORKDIR/phaseb/:

  features_unique.tsv.gz  one row per unique sequence (row, seq_sha256, features)
  features.tsv.gz         one row per member (set_id, source_id, gene_id), truth columns for
                          truth members, the features, and the embedding rows emb_row and
                          emb_cterm_row (rows of <model>.nterm.npy and <model>.cterm.npy)
  feature_coverage.tsv    per set_id: members, unique sequences, no SignalP call,
                          no PredGPI call, longer than 1,022 aa, PredGPI too_short
  features_run.json       input hashes, counts, all_sources, git commit, arguments

STOP (exit 2, no output): a tool output id that is not a unique seq_sha256 (stale output); an id
in two parts; SignalP prediction_results and output.gff3 that disagree; a truth member with no
row in truth_set_triaged.tsv.gz; d8_run.json without all_sources: true (unless
--allow-partial-truth-set); a unique sequence without a SignalP or PredGPI call (unless
--allow-missing-calls, which writes empty feature fields and counts them in the coverage).
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import feature_parsers as fp
import manifest
import paths
import runinfo
import seqhash
import truth_table

OUTPUT_NAMES = (
    "features_unique.tsv.gz",
    "features.tsv.gz",
    "feature_coverage.tsv",
    "features_run.json",
)
FEATURE_COLUMNS = (
    "ser_thr_frac",
    "sp_prediction",
    "sp_prob",
    "sp_other_prob",
    "sp_cs_end",
    "sp_cs_prob",
    "gpi_call",
    "gpi_prob",
    "gpi_omega",
    "gpi_fpr",
    "gpi_svm",
)
UNIQUE_FEATURE_COLUMNS = ("row", "seq_sha256", "length", "cterm_row", *FEATURE_COLUMNS)
TRUTH_COLUMNS = ("label", "subset", "stratum", "d8_class", "homology_only", "role")
MEMBER_FEATURE_COLUMNS = (
    "set_id",
    "source_id",
    "gene_id",
    "seq_sha256",
    "length",
    *TRUTH_COLUMNS,
    *FEATURE_COLUMNS,
    "emb_row",
    "emb_cterm_row",
)
COVERAGE_COLUMNS = (
    "set_id",
    "members",
    "unique_sequences",
    "no_signalp",
    "no_predgpi",
    "over_1022",
    "gpi_too_short",
)


class FeatureError(ValueError):
    """The tool outputs do not cover the unique sequences as expected."""


def _fmt(value) -> str:
    return "" if value is None else str(value)


def merge_parts(parts: list[dict], what: str) -> dict:
    merged: dict = {}
    for part in parts:
        dup = merged.keys() & part.keys()
        if dup:
            raise FeatureError(f"{what}: id {sorted(dup)[0]} occurs in two parts")
        merged.update(part)
    return merged


def check_inputs(unique, members, truth_rows) -> dict:
    """Check the unique sequences, the members and the truth rows. Return truth rows by key."""
    seen = set()
    for u in unique:
        h = u["seq_sha256"]
        if h in seen:
            raise FeatureError(f"unique_sequences has seq_sha256 {h} twice")
        seen.add(h)
        seq = u["sequence"]
        if seq != seqhash.clean(seq):
            raise FeatureError(f"sequence {h} is not cleaned upper case; re-run 05")
        if seqhash.seq_sha256(seq) != h:
            raise FeatureError(f"sequence {h} does not match its seq_sha256; re-run 05")
        if str(len(seq)) != u["length"]:
            raise FeatureError(f"sequence {h} has a length that differs from its row")
    keys = set()
    for m in members:
        key = (m["set_id"], m["source_id"], m["gene_id"])
        if key in keys:
            raise FeatureError(f"member {':'.join(key)} occurs twice in sequence_members")
        keys.add(key)
        if m["seq_sha256"] not in seen:
            raise FeatureError(
                f"member {':'.join(key)} has seq_sha256 {m['seq_sha256']} that is not in "
                "unique_sequences"
            )
    truth_by_key: dict = {}
    for r in truth_rows:
        key = (r["source_id"], r["gene_id"])
        if key in truth_by_key:
            raise FeatureError(f"truth_set_triaged.tsv.gz has {key[0]}:{key[1]} twice")
        truth_by_key[key] = r
    return truth_by_key


def check_truth_set_current(work: Path) -> str:
    """The truth set that 03 triaged must be the one that 02 attached sequences to."""
    triaged = json.loads((work / "d8_run.json").read_text()).get("truth_set_sha256")
    attached = json.loads((work / "sequence_run.json").read_text()).get("truth_set_sha256")
    if not triaged or triaged != attached:
        raise FeatureError(
            f"d8_run.json truth_set_sha256 ({triaged}) differs from sequence_run.json "
            f"({attached}); re-run 02, 03 and 05 on the same truth set"
        )
    return triaged


def feature_rows(unique, sp_calls, gpi_calls, allow_missing: bool):
    known = {r["seq_sha256"] for r in unique}
    for what, calls in (("SignalP", sp_calls), ("PredGPI", gpi_calls)):
        stale = calls.keys() - known
        if stale:
            raise FeatureError(
                f"{what} output has {len(stale)} ids that are not unique sequences "
                f"(for example {sorted(stale)[0]}); re-run J1 on the current unique_sequences"
            )
    missing_sp = [r["seq_sha256"] for r in unique if r["seq_sha256"] not in sp_calls]
    missing_gpi = [r["seq_sha256"] for r in unique if r["seq_sha256"] not in gpi_calls]
    if (missing_sp or missing_gpi) and not allow_missing:
        raise FeatureError(
            f"{len(missing_sp)} sequences have no SignalP call and {len(missing_gpi)} no "
            "PredGPI call; use --allow-missing-calls to write the table anyway"
        )
    rows = []
    for u in unique:
        sp = sp_calls.get(u["seq_sha256"])
        gpi = gpi_calls.get(u["seq_sha256"])
        rows.append(
            {
                "row": u["row"],
                "seq_sha256": u["seq_sha256"],
                "length": u["length"],
                "cterm_row": u["cterm_row"],
                "ser_thr_frac": f"{fp.ser_thr_fraction(u['sequence']):.6f}",
                "sp_prediction": sp.prediction if sp else "",
                "sp_prob": _fmt(sp.sp_prob if sp else None),
                "sp_other_prob": _fmt(sp.other_prob if sp else None),
                "sp_cs_end": _fmt(sp.cs_end if sp else None),
                "sp_cs_prob": _fmt(sp.cs_prob if sp else None),
                "gpi_call": gpi.call if gpi else "",
                "gpi_prob": _fmt(gpi.prob if gpi else None),
                "gpi_omega": _fmt(gpi.omega if gpi else None),
                "gpi_fpr": _fmt(gpi.fpr if gpi else None),
                "gpi_svm": _fmt(gpi.svm if gpi else None),
            }
        )
    return rows


def member_rows(members, features_by_hash, truth_by_key):
    out = []
    for m in members:
        f = features_by_hash[m["seq_sha256"]]
        truth = {c: "" for c in TRUTH_COLUMNS}
        if m["set_id"] == "truth":
            t = truth_by_key.get((m["source_id"], m["gene_id"]))
            if t is None:
                raise FeatureError(
                    f"truth member {m['source_id']}:{m['gene_id']} has no row in "
                    "truth_set_triaged.tsv.gz; re-run 03_triage_pm.py and 05"
                )
            truth = {c: t.get(c, "") for c in TRUTH_COLUMNS}
        out.append(
            {
                **{c: m[c] for c in ("set_id", "source_id", "gene_id", "seq_sha256", "length")},
                **truth,
                **{c: f[c] for c in FEATURE_COLUMNS},
                "emb_row": f["row"],
                "emb_cterm_row": f["cterm_row"],
            }
        )
    return out


def coverage_rows(members, features_by_hash):
    by_set: dict[str, list] = {}
    for m in members:
        by_set.setdefault(m["set_id"], []).append(m)
    rows = []
    for set_id, ms in by_set.items():
        hashes = {m["seq_sha256"] for m in ms}
        feats = [features_by_hash[h] for h in hashes]
        rows.append(
            {
                "set_id": set_id,
                "members": str(len(ms)),
                "unique_sequences": str(len(hashes)),
                "no_signalp": str(sum(f["sp_prediction"] == "" for f in feats)),
                "no_predgpi": str(sum(f["gpi_call"] == "" for f in feats)),
                "over_1022": str(sum(f["cterm_row"] != "" for f in feats)),
                "gpi_too_short": str(sum(f["gpi_call"] == "too_short" for f in feats)),
            }
        )
    return rows


def read_tool_outputs(out: Path) -> tuple[dict, dict, dict]:
    sp_dirs = sorted((out / "signalp").glob("part_*"))
    gpi_files = sorted((out / "predgpi").glob("part_*.tsv.gz"))
    if not sp_dirs or not gpi_files:
        raise FeatureError(f"no J1 output under {out}/signalp or {out}/predgpi")
    hashes, sp_parts = {}, []
    for d in sp_dirs:
        pred, gff = d / "prediction_results.txt.gz", d / "output.gff3.gz"
        calls = fp.parse_signalp(pred)
        fp.check_signalp_consistency(calls, fp.parse_signalp_gff(gff))
        sp_parts.append(calls)
        hashes[str(pred.relative_to(out))] = manifest.sha256_file(pred)
        hashes[str(gff.relative_to(out))] = manifest.sha256_file(gff)
    gpi_parts = []
    for f in gpi_files:
        gpi_parts.append(fp.parse_predgpi_scores(f))
        hashes[str(f.relative_to(out))] = manifest.sha256_file(f)
    return merge_parts(sp_parts, "SignalP"), merge_parts(gpi_parts, "PredGPI"), hashes


def run(work: Path, allow_missing: bool, provenance: dict | None = None):
    out = Path(work) / "phaseb"
    unique = truth_table.read_tsv(out / "unique_sequences.tsv.gz")
    members = truth_table.read_tsv(out / "sequence_members.tsv.gz")
    truth_path = Path(work) / "truth_set_triaged.tsv.gz"
    truth_by_key = check_inputs(unique, members, truth_table.read_tsv(truth_path))
    truth_set_sha256 = check_truth_set_current(Path(work))
    sp_calls, gpi_calls, hashes = read_tool_outputs(out)
    feats = feature_rows(unique, sp_calls, gpi_calls, allow_missing)
    by_hash = {f["seq_sha256"]: f for f in feats}
    mrows = member_rows(members, by_hash, truth_by_key)
    cov = coverage_rows(members, by_hash)
    log = {
        "tool_outputs_sha256": hashes,
        "truth_set_sha256": truth_set_sha256,
        "input_sha256": {
            "truth_set_triaged.tsv.gz": manifest.sha256_file(truth_path),
            "sequence_members.tsv.gz": manifest.sha256_file(out / "sequence_members.tsv.gz"),
            "unique_sequences.tsv.gz": manifest.sha256_file(out / "unique_sequences.tsv.gz"),
        },
        "unique_sequences": len(feats),
        "members": len(mrows),
        "missing_signalp": sum(f["sp_prediction"] == "" for f in feats),
        "missing_predgpi": sum(f["gpi_call"] == "" for f in feats),
        "sp_predictions": dict(Counter(f["sp_prediction"] for f in feats)),
        "gpi_calls": dict(Counter(f["gpi_call"] for f in feats)),
        "all_sources": None,
        "git_commit": runinfo.git_commit(),
        "python": runinfo.python_version(),
        "arguments": [],
    }
    log.update(provenance or {})
    runinfo.atomic_write_all(
        out,
        {
            "features_unique.tsv.gz": lambda p: truth_table.write_tsv(
                p, UNIQUE_FEATURE_COLUMNS, feats
            ),
            "features.tsv.gz": lambda p: truth_table.write_tsv(p, MEMBER_FEATURE_COLUMNS, mrows),
            "feature_coverage.tsv": lambda p: truth_table.write_tsv(p, COVERAGE_COLUMNS, cov),
            "features_run.json": lambda p: runinfo.write_json(p, log),
        },
    )
    return feats, mrows, cov, log


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=None, help="default: $STEP1_WORKDIR")
    parser.add_argument("--allow-missing-calls", action="store_true")
    parser.add_argument("--allow-partial-truth-set", action="store_true")
    args = parser.parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else paths.workdir()
    try:
        runinfo.require_full(
            work / "d8_run.json", "truth_set_triaged.tsv.gz", args.allow_partial_truth_set
        )
        runinfo.require_full(
            work / "phaseb" / "prepare_run.json",
            "unique_sequences.tsv.gz",
            args.allow_partial_truth_set,
        )
        provenance = {
            "all_sources": runinfo.says_all_sources(work / "d8_run.json")
            and runinfo.says_all_sources(work / "phaseb" / "prepare_run.json"),
            "arguments": list(argv) if argv is not None else sys.argv[1:],
        }
        _, _, cov, log = run(work, args.allow_missing_calls, provenance)
    except (FeatureError, fp.OutputFormatError, ValueError, OSError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 2
    for c in cov:
        print("\t".join(f"{k}={c[k]}" for k in COVERAGE_COLUMNS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
