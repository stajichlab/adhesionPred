def log(tool, q, n, status, note):
    k = sum(1 for _ in open("query_log.tsv"))
    with open("query_log.tsv", "a") as f:
        f.write(f"{k}\t2026-10-02\t{tool}\t{q}\t{n}\t{status}\t{note}\n")
