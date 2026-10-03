import sys

n = sum(1 for _ in open("query_log.tsv"))
with open("query_log.tsv", "a") as f:
    f.write("\t".join([str(n), "2026-10-02"] + sys.argv[1:]) + "\n")
