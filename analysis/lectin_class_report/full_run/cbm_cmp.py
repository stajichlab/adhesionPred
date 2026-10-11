# ruff: noqa: UP031
import csv

new = {}
for r in csv.DictReader(open("out/lectin_report_hits.tsv"), delimiter="\t"):
    new.setdefault(r["model_acc"], {})[(r["proteome"], r["protein"])] = float(r["seq_score"])
com = {}
for r in csv.DictReader(open("cbm_whole_proteome_hits.tsv"), delimiter="\t"):
    com.setdefault(r["model_acc"].split(".")[0], {})[(r["proteome"], r["protein"])] = float(
        r["seq_score"]
    )
print(
    "model\tnew_proteins\tnew_genera\tcommitted_proteins\tcommitted_genera\tsame_set\tmax_score_diff"
)
for acc in sorted(new):
    n = new[acc]
    c = com.get(acc)
    g = len({k[0] for k in n})
    if c is None:
        print("%s\t%d\t%d\tNA\tNA\tNA\tNA" % (acc, len(n), g))
        continue
    d = max((abs(n[k] - c[k]) for k in set(n) & set(c)), default=0)
    print(
        "%s\t%d\t%d\t%d\t%d\t%s\t%.1f"
        % (acc, len(n), g, len(c), len({k[0] for k in c}), set(n) == set(c), d)
    )
