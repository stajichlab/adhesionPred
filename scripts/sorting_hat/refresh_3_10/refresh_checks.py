"""Small read-only checks for the Task 3.10 job scripts (standard library only).

Subcommands (each exits 0 when the check holds, 1 when it fails, 3 for a measured difference, 4 for a missing record):
  identity-equal MODULE_JSON STATUS_JSON   module identity of a rebuilt module equals the one in a status file
  json-equal A B                           two JSON files hold the same data
  notes STATUS_JSON                        print the --notes text of the first entry (before the counts)
  measure STATUS_JSON SENS SPEC            print and compare sensitivity and specificity (3 decimals)
  models-searched PROVENANCE_JSON FAMILY_TABLE
                                           every row of the family table is in the model list of the Pfam provenance
                                           (exit 1: a table row is not in the record; exit 4: the provenance has no list)
  r0-check STATUS_JSON                     every entry carries leakage_exception and negative_basis
"""

import csv
import json
import sys

IDENTITY_KEYS = ("module", "version", "params_hash", "artefact_hash")


def load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def identity_equal(module_json, status_json):
    new, old = load(module_json), load(status_json)
    bad = [k for k in IDENTITY_KEYS if new.get(k) != old.get(k)]
    for k in bad:
        print(f"DIFFERS {k}: rebuilt {new.get(k)!r} status file {old.get(k)!r}")
    if not bad:
        print(f"same identity: {module_json} and {status_json}")
    return 1 if bad else 0


def json_equal(a, b):
    same = load(a) == load(b)
    print(("same data: " if same else "DIFFERENT data: ") + f"{a} {b}")
    return 0 if same else 1


def notes(status_json):
    text = load(status_json)["entries"][0]["measure"]["notes"]
    marker = " truth rows without a call"
    if marker not in text:
        print(f"no marker {marker!r} in the notes of {status_json}", file=sys.stderr)
        return 1
    print(text.split(marker)[0])
    return 0


def measure(status_json, sens, spec):
    m = load(status_json)["entries"][0]["measure"]
    got_sens, got_spec = m["sensitivity"], m["specificity"]
    print(
        f"sensitivity {got_sens['value']:.3f} [{got_sens['lo']:.3f}, {got_sens['hi']:.3f}] "
        f"specificity {got_spec['value']:.3f} [{got_spec['lo']:.3f}, {got_spec['hi']:.3f}] "
        f"n_pos {m['n_pos']} n_neg {m['n_neg']}"
    )
    same = f"{got_sens['value']:.3f}" == sens and f"{got_spec['value']:.3f}" == spec
    print(f"expected sensitivity {sens} specificity {spec}: " + ("MATCH" if same else "DIFFERS"))
    return 0 if same else 3


def r0_check(status_json):
    entries = load(status_json)["entries"]
    bad = 0
    for e in entries:
        m = e["measure"]
        ok = "leakage_exception" in m and "negative_basis" in m
        bad += not ok
        print(
            f"taxa {e['taxa']} status {e['status']} set {m['calibration_set']} "
            f"leakage_exception {'leakage_exception' in m} negative_basis {'negative_basis' in m}"
        )
    return 1 if bad or not entries else 0


def models_searched(provenance_json, family_table):
    """Compare the model list recorded by pfam_hmmsearch.sbatch with the rows of the family table.

    A record is a list of {"acc": ..., "name": ...}. Models in the record that are not in the table are allowed.
    """
    rec = load(provenance_json)
    if not isinstance(rec.get("models"), list):
        print(f"NO RECORDED MODEL LIST in {provenance_json}: keys are {sorted(rec)}")
        return 4
    recorded = {(m["acc"], m["name"]) for m in rec["models"]}
    with open(family_table, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    wanted = {(r["pfam_acc"], r["name"]) for r in rows}
    missing = sorted(wanted - recorded)
    for acc, name in missing:
        print(
            f"NOT SEARCHED {acc} {name}: a row of {family_table} is not in the model list of {provenance_json}"
        )
    if not missing:
        print(f"all {len(wanted)} models of the family table were searched: {provenance_json}")
    return 1 if missing else 0


def main(argv):
    cmd, args = (argv[0] if argv else ""), argv[1:]
    table = {
        "identity-equal": (identity_equal, 2),
        "json-equal": (json_equal, 2),
        "notes": (notes, 1),
        "measure": (measure, 3),
        "models-searched": (models_searched, 2),
        "r0-check": (r0_check, 1),
    }
    if cmd not in table or len(args) != table[cmd][1]:
        print(__doc__, file=sys.stderr)
        return 2
    return table[cmd][0](*args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
