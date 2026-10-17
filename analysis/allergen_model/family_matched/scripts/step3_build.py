"""Step 4 prep: fetch sequences, names, verify filters, write family_negatives.faa/tsv."""

import csv
import hashlib
import json
import os
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(__file__))
from uni import get_all

OUT = "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/family_matched"
st = json.load(open(OUT + "/scripts/step2_state.json"))
sel = st["selected"]
fam_names = {}
for fam in sorted({v["family"] for v in sel.values()}):
    for _ in range(4):
        try:
            with urllib.request.urlopen(
                "https://www.ebi.ac.uk/interpro/api/entry/pfam/" + fam, timeout=60
            ) as r:
                fam_names[fam] = json.loads(r.read())["metadata"]["name"]["name"]
            break
        except Exception:
            time.sleep(3)
    else:
        fam_names[fam] = ""
        print("name lookup failed", fam, file=sys.stderr)
accs = sorted(sel)
CAND = {
    r["accession"]
    for r in csv.DictReader(
        open(
            "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/allergen_model/candidates/candidates.tsv"
        ),
        delimiter="\t",
    )
}
seqs = {}
for i in range(0, len(accs), 100):
    chunk = accs[i : i + 100]
    url = "https://rest.uniprot.org/uniprotkb/accessions?" + urllib.parse.urlencode(
        {"accessions": ",".join(chunk), "format": "fasta"}
    )
    with urllib.request.urlopen(url, timeout=120) as r:
        txt = r.read().decode()
    cur = None
    for ln in txt.splitlines():
        if ln.startswith(">"):
            cur = ln[1:].split("|")[1]
            seqs[cur] = []
        elif cur:
            seqs[cur].append(ln.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
missing = [a for a in accs if a not in seqs]
print("fasta missing", len(missing), missing[:5], file=sys.stderr)
# verify filters: re-query keyword/allergen exclusion and reviewed status for the final set
bad = []
for i in range(0, len(accs), 40):
    chunk = accs[i : i + 40]
    q = (
        "("
        + " OR ".join(f"accession:{a}" for a in chunk)
        + ") AND (keyword:KW-0020 OR cc_allergen:allergic)"
    )
    rows, _ = get_all(q, "accession")
    bad += [r["Entry"] for r in rows]
print("verify failures (KW-0020 or allergic):", len(bad), file=sys.stderr)
rev = {}
for i in range(0, len(accs), 40):
    chunk = accs[i : i + 40]
    q = "(" + " OR ".join(f"accession:{a}" for a in chunk) + ")"
    rows, _ = get_all(q, "accession,reviewed")
    for r in rows:
        rev[r["Entry"]] = r["Reviewed"]
with (
    open(OUT + "/family_negatives.faa", "w") as f,
    open(OUT + "/family_negatives.tsv", "w", newline="") as t,
):
    w = csv.writer(t, delimiter="\t", lineterminator="\n")
    w.writerow(
        [
            "accession",
            "family",
            "family_name",
            "organism",
            "length",
            "sha256",
            "reviewed",
            "overlap_with_candidates",
        ]
    )
    for a in accs:
        s = seqs[a]
        r = sel[a]["row"]
        fam = sel[a]["family"]
        f.write(f">F_{a}\n")
        for j in range(0, len(s), 60):
            f.write(s[j : j + 60] + "\n")
        w.writerow(
            [
                a,
                fam,
                fam_names.get(fam, ""),
                r["Organism"],
                len(s),
                hashlib.sha256(s.encode()).hexdigest(),
                rev.get(a, ""),
                "true" if a in CAND else "false",
            ]
        )
json.dump(
    {
        "fam_names": fam_names,
        "missing": missing,
        "verify_bad": bad,
        "rev": rev,
        "seqlen": {a: len(seqs[a]) for a in accs},
    },
    open(OUT + "/scripts/step3_state.json", "w"),
)
print("written", len(accs))
