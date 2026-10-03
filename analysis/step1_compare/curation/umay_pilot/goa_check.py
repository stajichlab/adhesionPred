import collections
import csv
import sys

sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/adhesionPred/analysis/step1_compare")
import gaf
import go_obo
import labels

D = "/bigdata/stajichlab/jstajich/projects/adhesionPred/_workdir/step1_compare/downloads/"
onto = go_obo.parse_obo(D + "go-basic.obo")
fl = gaf.filter_gaf(D + "MYCMD-uniprot.gaf.gz", onto, "5270")
print("GOA cc rows", len(fl.cc_rows), "dropped", dict(fl.dropped), file=sys.stderr)
goa = collections.defaultdict(set)
for r in fl.cc_rows:
    goa[r.gene_id].add((r.term, r.evidence))
cur = collections.defaultdict(set)
exp = {}
obs = []
for r in csv.DictReader(open("draft_rows.tsv"), delimiter="\t"):
    if not onto.known(r["go_term"]) or r["go_term"] in onto.obsolete:
        obs.append(r["go_term"])
    cur[r["uniprot_accession"]].add((r["go_term"], r["evidence_code"]))
    exp.setdefault(r["uniprot_accession"], (r["symbol"], set()))[1].add(r["expected_label"])


def lab(rows, pol):
    anyt = set()
    lt = set()
    for t, e in rows:
        anyt |= onto.ancestors(t)
        if labels.policy_accepts(pol, e):
            lt |= onto.ancestors(t)
    return labels.classify(frozenset(lt), frozenset(anyt))


out = []
print(
    "gene\tsymbol\texpected\tGOA_non_iea\tGOA_no_homology\tGOA_experimental\tcurated_only_non_iea\tmerged_non_iea\tmerged_no_homology\tmerged_experimental\tconflict\tgoa_rows"
)
for a, (sym, ex) in exp.items():
    g = goa.get(a, set())
    c = cur[a]
    m = g | c
    L = [lab(g, p) for p in labels.POLICIES]
    C = lab(c, "non_iea")
    M = [lab(m, p) for p in labels.POLICIES]
    e = ",".join(sorted(ex))
    conflict = "yes" if any(x != e for x in M) else "no"
    print(
        "\t".join(
            [a, sym, e, *L, C, *M, conflict, ";".join(sorted(f"{t}:{ev}" for t, ev in g)) or "none"]
        )
    )
print("obsolete/unknown terms:", obs, file=sys.stderr)
